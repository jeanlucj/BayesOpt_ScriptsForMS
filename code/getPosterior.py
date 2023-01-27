import sys
# I am going to get the initialization number for an argument passed to the script
init_num = sys.argv[1]
# To know how the parameters work, know how many stages are in the breeding scheme
n_stages = int(sys.argv[2])
# Based on how many iterations do you want the posterior?
n_iter = int(sys.argv[3])
# Get max or make predictions for a number of points?
# If the latter, give the file name, "file_name.rds"
pred_points = sys.argv[4]

print(init_num)
print(n_stages)
print(n_iter)
print(pred_points)

# ### Get posterior means and variances from optimization iterations
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using {} device".format(device))
dtype = torch.double
torch.set_default_dtype(dtype)

# ### Setup the dependent R instance
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
numpy2ri.activate()

# Put the initialization number also in the R instance
ro.globalenv['init_num'] = init_num
ro.globalenv['n_iter'] = n_iter

# ### Call R to get the optimization iterations
ro.r.source('getOptimizationIterations.R')
train_x = torch.tensor(ro.globalenv['budgets'])
train_obj = torch.tensor(ro.globalenv['gains']).unsqueeze(-1)

# ### Call R to get the budget for which you want predictions
if pred_points != "0":
    ro.globalenv['predFile'] = pred_points
    ro.r.source('getPredPoints.R')
    pred_x = torch.tensor(ro.globalenv['predBudgets'])

# ### Model initialization
from botorch.models.gp_regression import SingleTaskGP
from botorch.models.transforms.outcome import Standardize
from gpytorch.mlls.exact_marginal_log_likelihood import ExactMarginalLogLikelihood
    
def initialize_model(train_x, train_obj):
    # define model for objective
    surrogate = SingleTaskGP(train_x, train_obj, outcome_transform=Standardize(m=1))
    mll = ExactMarginalLogLikelihood(surrogate.likelihood, surrogate)
    # fit the model
    fit_gpytorch_model(mll)
    return surrogate

# ### Get the model all parameterized in R -- will need this for inequality constraints
from botorch import fit_gpytorch_model
from botorch.acquisition.monte_carlo import qExpectedImprovement
from botorch.optim import optimize_acqf
from botorch.exceptions import BadInitialCandidatesWarning

import warnings

warnings.filterwarnings('ignore', category=BadInitialCandidatesWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

MC_SAMPLES = 256
NUM_RESTARTS = 30 # 10 * input_dim
RAW_SAMPLES = 600 # 200 * input_dim

# ### Initialize the scheme and run the burn-in cycles
ro.r.source('breedSchemeInitialize.R')
budget_constraints = torch.tensor(ro.globalenv['budget_constraints'])
# Set this up so that only n_stages parameters are being captured
# So, they have to sum to less than 1.0
mb = 1.0 - budget_constraints[0] - budget_constraints[1]
bounds = torch.tensor(
    [[budget_constraints[0], 0.0, 0.0], [mb] * n_stages],
    device=device, dtype=dtype)

lr = budget_constraints[4]
inequality_constraints = [
    (torch.tensor([i for i in range(n_stages)]), torch.tensor([-1.0] * (n_stages)), -(1.0 - budget_constraints[1])), 
    (torch.tensor([0, 1]), torch.tensor([1.0, -budget_constraints[2]]), 0.0),
    (torch.tensor([1, 2]), torch.tensor([1.0, -budget_constraints[3]]), 0.0),
    (torch.tensor([i for i in range(n_stages)]), torch.tensor([lr, lr, 1+lr]), lr),
    ]

# ### Get parameter values for highest posterior mean
from botorch.acquisition.analytic import PosteriorMean

surrogate = initialize_model(train_x, train_obj)
posterior_mean_aqcf = PosteriorMean(model=surrogate)

x_opt, _ = optimize_acqf(
        acq_function=posterior_mean_aqcf,
        bounds=bounds,
        inequality_constraints=inequality_constraints,
        q=1,
        num_restarts=NUM_RESTARTS,
        raw_samples=RAW_SAMPLES,  # used for intialization heuristic
        options={"batch_limit": 5, "maxiter": 200},
    )

# ### Note: you can also send in a matrix of parameter values in x_opt here
posterior = surrogate.posterior(x_opt)
expected_gain = posterior.mean
post_var_gain = posterior.variance
print(x_opt)
print(expected_gain)
print(post_var_gain)

# ### Send the optimum X back to R
ro.globalenv['bestBudget'] = x_opt.numpy()
ro.globalenv['maxPredGain'] = expected_gain.detach().numpy()
ro.globalenv['postVarAtMax'] = post_var_gain.detach().numpy()

# ### In case I have specific X for which I want the posterior mean and variance
if pred_points != "0":
    post_points = surrogate.posterior(pred_x)
    pred_gains = post_points.mean
    pred_vars = post_points.variance
    ro.globalenv['predGains'] = pred_gains.detach().numpy()
    ro.globalenv['predVars'] = pred_vars.detach().numpy()

# ### Load the packages needed in the R instance
ro.r.source('processPosteriors.R')

import time
time.sleep(1)
