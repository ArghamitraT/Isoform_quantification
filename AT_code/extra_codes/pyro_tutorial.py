import torch
import torch.optim as optim
import torch.nn.functional as F
from torch.distributions.dirichlet import Dirichlet

# Sample data: 2D tensor where each row is a sample from a Dirichlet distribution
# data = torch.tensor([
#     [0.2, 0.5, 0.3],
#     [0.1, 0.3, 0.6],
#     [0.3, 0.4, 0.3]
# ], dtype=torch.float32)

# # Initialize alpha parameters (must be positive)
# alpha = torch.tensor([1.0, 1.0, 1.0], dtype=torch.float32, requires_grad=True)

# # Optimizer
# optimizer = optim.Adam([alpha], lr=0.01)

# # Number of iterations
# num_iterations = 1000

# Gradient descent loop
# for iteration in range(num_iterations):
#     optimizer.zero_grad()
    
#     # Calculate log-likelihood
#     dirichlet = Dirichlet(alpha)
#     log_likelihood = dirichlet.log_prob(data).mean()
    
#     # Since we want to maximize the log-likelihood, we minimize the negative log-likelihood
#     loss = -log_likelihood
    
#     # Backpropagation
#     loss.backward()
    
#     # Gradient descent step
#     optimizer.step()
    
#     # Ensure alpha remains positive
#     with torch.no_grad():
#         alpha.clamp_(min=1e-3)
    
#     # Print the loss and alpha parameters every 100 iterations
#     #print(f'Iteration {iteration + 1}/{num_iterations}, Loss: {loss.item()}, Alpha: {alpha.data}')

#     if (iteration + 1) % 100 == 0:
#         print(f'Iteration {iteration + 1}/{num_iterations}, Loss: {loss.item()}, Alpha: {alpha.data}')

# # Final alpha parameters
# print(f'Estimated Alpha: {alpha.data}')

# Example data in dictionary format
all_theta = {
    'sample1': {'isoform1': 0.2, 'isoform2': 0.5, 'isoform3': 0.3},
    'sample2': {'isoform1': 0.1, 'isoform2': 0.3, 'isoform3': 0.6},
    'sample3': {'isoform1': 0.3, 'isoform2': 0.4, 'isoform3': 0.3}
}

# Initial alpha values in dictionary format
all_alpha = {'isoform1': 1.0, 'isoform2': 1.0, 'isoform3': 1.0}

# Convert the dictionaries to tensors
isoforms = list(all_alpha.keys())
data = torch.tensor([[sample[isoform] for isoform in isoforms] for sample in all_theta.values()], dtype=torch.float32)
alpha = torch.tensor([all_alpha[isoform] for isoform in isoforms], dtype=torch.float32, requires_grad=True)

# Optimizer
optimizer = optim.Adam([alpha], lr=0.01)

# Number of iterations
num_iterations = 1000

def log_likelihood(data, alpha):
    term1 = torch.lgamma(torch.sum(alpha))
    term2 = torch.sum(torch.lgamma(alpha))
    term3 = torch.sum((alpha - 1) * torch.log(data), dim=1)
    return term1 - term2 + term3.mean()

# Gradient descent loop
for iteration in range(num_iterations):
    optimizer.zero_grad()
    
    # Calculate log-likelihood
    ll = log_likelihood(data, alpha)
    
    # Since we want to maximize the log-likelihood, we minimize the negative log-likelihood
    loss = -ll
    
    # Backpropagation
    loss.backward()
    
    # Gradient descent step
    optimizer.step()
    
    # Ensure alpha remains positive
    with torch.no_grad():
        alpha.clamp_(min=1e-3)
    
    # Print the loss and alpha parameters every 100 iterations
    if (iteration + 1) % 100 == 0:
        print(f'Iteration {iteration + 1}/{num_iterations}, Loss: {loss.item()}, Alpha: {alpha.data}')

# Final alpha parameters
print(f'Estimated Alpha: {alpha.data}')



# import torch
# import pyro
# import pyro.distributions as dist
# import pyro.optim as optim
# from pyro.infer import SVI, Trace_ELBO
# from pyro.infer.autoguide import AutoDelta

# class DirichletModel:
#     def __init__(self, all_theta, all_alpha, expectation_log_theta, GD_lr):
#         self.all_theta = all_theta
#         self.GD_lr = GD_lr
#         self.initial_alpha = all_alpha

#         # Initialize alpha parameters as a list of tensors with positive constraints
#         self.alpha = pyro.param("alpha", torch.tensor(all_alpha, dtype=torch.float32), constraint=dist.constraints.positive)

#         self.optimizer = optim.Adam({"lr": self.GD_lr})

#     def model(self, theta):
#         alpha = self.alpha.expand(len(theta), len(self.initial_alpha))
#         with pyro.plate("data", len(theta)):
#             pyro.sample("obs", dist.Dirichlet(alpha), obs=theta)

#     def guide(self, theta):
#         alpha_q = pyro.param("alpha_q", torch.tensor(self.initial_alpha, dtype=torch.float32), constraint=dist.constraints.positive)
#         alpha = alpha_q.expand(len(theta), len(self.initial_alpha))
#         with pyro.plate("data", len(theta)):
#             pyro.sample("obs_guide", dist.Dirichlet(alpha).to_event(1))

#     def update_alpha(self, max_iterations=10, tolerance=1e-6):
#         """
#         Gradient Descent of alpha considering all samples with convergence criteria.
#         :param max_iterations: int, maximum number of iterations to perform
#         :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
#         """
#         svi = SVI(self.model, self.guide, self.optimizer, loss=Trace_ELBO())

#         theta = self._theta_to_tensor(self.all_theta)  # Prepare the tensor once

#         previous_loss = float('inf')  # Initialize previous loss as infinity for comparison
#         for iteration in range(max_iterations):
#             loss = svi.step(theta)

#             # Print gradients and parameters for debugging
#             for name, value in pyro.get_param_store().items():
#                 print(name, pyro.param(name).data.cpu().numpy())

#             # Check for convergence: if the change in loss is less than the tolerance, stop the loop
#             if abs(previous_loss - loss) < tolerance:
#                 print(f"GD Convergence reached after {iteration + 1} iterations.")
#                 break
#             previous_loss = loss  # Update previous loss to current loss

#             # Optionally, print the current loss every few iterations to monitor progress
#             print(f"GD_Iteration {iteration}")
#             print(f"GD_Current_Loss = {loss}")

#         # After updating self.alpha, convert to non-log space and non-tensor
#         alpha_nontensor = self.alpha.detach().numpy()

#         return alpha_nontensor

#     def _theta_to_tensor(self, all_theta):
#         """
#         Convert all_theta (a dictionary of Counters) to a tensor for Pyro processing.

#         Parameters:
#         all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

#         Returns:
#         torch.Tensor: Tensor of theta values.
#         """
#         # Get the isoform names from all_theta keys
#         isoforms = sorted(set(isoform for sample in all_theta.values() for isoform in sample))
#         theta_list = []

#         for sample_key, isoform_dict in all_theta.items():
#             sample_list = []
#             for isoform in isoforms:
#                 if isoform in isoform_dict:
#                     sample_list.append(torch.tensor(isoform_dict[isoform], dtype=torch.float32))
#                 else:
#                     sample_list.append(torch.tensor(0.0, dtype=torch.float32))  # Missing isoforms are set to 0
#             theta_list.append(torch.stack(sample_list))

#         return torch.stack(theta_list)

# # Example usage
# # all_theta = {
# #     'sample1': {'isoform1': 0.2, 'isoform2': 0.8, 'isoform3': 0.0},  # isoform3 is missing in sample1
# #     'sample2': {'isoform1': 0.1, 'isoform2': 0.4, 'isoform3': 0.5}
# # }
# # all_alpha = [50, 5, 15]  # List of initial alpha values
# # expectation_log_theta = None
# # GD_lr = 0.01

# # model = DirichletModel(all_theta, all_alpha, expectation_log_theta, GD_lr)
# # optimized_alpha = model.update_alpha(max_iterations=1000)
# # print("Optimized alpha parameters:", optimized_alpha)






# # import torch
# # import pyro
# # import pyro.distributions as dist
# # import pyro.optim as optim
# # from pyro.infer import SVI, Trace_ELBO
# # from pyro.infer.autoguide import AutoDelta

# # class DirichletModel:
# #     def __init__(self, all_theta, all_alpha, expectation_log_theta, GD_lr):
# #         self.all_theta = all_theta
# #         self.GD_lr = GD_lr
# #         self.initial_alpha = all_alpha

# #         # Initialize alpha parameters with Pyro param and positive constraint
# #         self.alpha = {
# #             k: pyro.param(f"alpha_{k}", torch.tensor(v, dtype=torch.float32), constraint=dist.constraints.positive)
# #             for k, v in all_alpha.items()
# #         }

# #         #self.optimizer = optim.Adam({"lr": self.GD_lr})
# #         self.optimizer = optim.Adam({"lr": self.GD_lr})

# #     def model(self):
# #         alpha = torch.stack([self.alpha[isoform] for isoform in sorted(self.initial_alpha.keys())])
# #         theta = self._theta_to_tensor(self.all_theta)  # Prepare the tensor once
# #         with pyro.plate("data", len(theta)):
# #             pyro.sample("obs", dist.Dirichlet(alpha), obs=theta)

# #     def guide(self):
# #         alpha = torch.stack([self.alpha[isoform] for isoform in sorted(self.initial_alpha.keys())])
# #         pyro.sample("obs_guide", dist.Dirichlet(alpha))
# #         # AutoDelta automatically creates the guide using delta distributions
# #         #guide = AutoDelta(self.model)
# #         #guide(theta)

# #     def update_alpha(self, max_iterations=10, tolerance=1e-6):
# #         """
# #         Gradient Descent of alpha considering all samples with convergence criteria.
# #         :param max_iterations: int, maximum number of iterations to perform
# #         :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
# #         """
# #         svi = SVI(self.model, self.guide, self.optimizer, loss=Trace_ELBO())
# #         theta = self._theta_to_tensor(self.all_theta)  # Prepare the tensor once

# #         previous_loss = float('inf')  # Initialize previous loss as infinity for comparison
# #         for iteration in range(max_iterations):
# #             loss = svi.step()

# #             # Print gradients and parameters for debugging
# #             for name, value in pyro.get_param_store().items():
# #                 print(name, pyro.param(name).data.cpu().numpy())

# #             # Check for convergence: if the change in loss is less than the tolerance, stop the loop
# #             if abs(previous_loss - loss) < tolerance:
# #                 print(f"GD Convergence reached after {iteration + 1} iterations.")
# #                 break
# #             previous_loss = loss  # Update previous loss to current loss

# #             # Optionally, print the current loss every few iterations to monitor progress
# #             print(f"GD_Iteration {iteration}")
# #             print(f"GD_Current_Loss = {loss}")

# #         # After updating self.alpha, convert to non-log space and non-tensor
# #         alpha_nontensor = {k: v.item() for k, v in self.alpha.items()}

# #         return alpha_nontensor

# #     def _theta_to_tensor(self, all_theta):
# #         """
# #         Convert all_theta (a dictionary of Counters) to a tensor for Pyro processing.

# #         Parameters:
# #         all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

# #         Returns:
# #         torch.Tensor: Tensor of theta values.
# #         """
# #         isoforms = sorted(self.initial_alpha.keys())
# #         theta_list = []

# #         for sample_key, isoform_dict in all_theta.items():
# #             sample_list = []
# #             for isoform in isoforms:
# #                 if isoform in isoform_dict:
# #                     sample_list.append(torch.tensor(isoform_dict[isoform], dtype=torch.float32))
# #                 else:
# #                     sample_list.append(torch.tensor(0.0, dtype=torch.float32))  # Missing isoforms are set to 0
# #             theta_list.append(torch.stack(sample_list))

# #         return torch.stack(theta_list)



# import torch
# import pyro
# import pyro.distributions as dist
# from pyro.infer import SVI, Trace_ELBO
# from pyro.optim import Adam


# pyro.enable_validation(True)
# pyro.clear_param_store()

# data = torch.cat((torch.zeros(9), torch.ones(7), torch.empty(4).fill_(2.)))

# def model(data):
#   alpha = torch.tensor(6.0)
#   beta = torch.tensor(10.0)
#   #pay_probs = pyro.sample('pay_probs', dist.Beta(alpha, beta).expand(3).independent(1))
#   pay_probs = pyro.sample('pay_probs', dist.Beta(alpha, beta).expand((3,)).independent(1))
#   normalized_pay_probs = pay_probs / torch.sum(pay_probs)

#   with pyro.iarange('data_loop', len(data)):
#     pyro.sample('obs', dist.Categorical(probs=normalized_pay_probs), obs=data)


# def guide(data):
#     alphas = pyro.param('alphas', torch.tensor(6.).expand(3), constraint=dist.constraints.positive)
#     betas = pyro.param('betas', torch.tensor(10.).expand(3), constraint=dist.constraints.positive) 

#     for name, value in pyro.get_param_store().items():
#         print(name+"_GUIDE", pyro.param(name).data.cpu().numpy())

#     pyro.sample('pay_probs', dist.Beta(alphas, betas).independent(1))


# def print_progress():
#     alphas = pyro.param("alphas")
#     betas = pyro.param("betas")
#     # Print gradients and parameters for debugging
#     for name, value in pyro.get_param_store().items():
#         print(name, pyro.param(name).data.cpu().numpy())


#     if torch.cuda.is_available():
#         alphas.cuda()
#         betas.cuda()

#     means = alphas / (alphas + betas)
#     normalized_means = means / torch.sum(means)
#     factors = betas / (alphas * (1.0 + alphas + betas))
#     stdevs = normalized_means * torch.sqrt(factors)

#     tiger_pays_string = "probability Tiger pays: {0:.3f} +/- {1:.2f}".format(normalized_means[0], stdevs[0])
#     jason_pays_string = "probability Jason pays: {0:.3f} +/- {1:.2f}".format(normalized_means[1], stdevs[1])
#     james_pays_string = "probability James pays: {0:.3f} +/- {1:.2f}".format(normalized_means[2], stdevs[2])
#     print("[", step, "|", tiger_pays_string, "|", jason_pays_string, "|", james_pays_string, "]")


# adam_params = {"lr": 0.0005}
# optimizer = Adam(adam_params)
# svi = SVI(model, guide, optimizer, loss=Trace_ELBO())

# n_steps = 2501
# for step in range(n_steps):
#   svi.step(data)
#   for name, value in pyro.get_param_store().items():
#         print(name, pyro.param(name).data.cpu().numpy())

#   if step % 100 == 0:
#     print_progress()




# # import logging
# # import os

# # import torch
# # import numpy as np
# # import pandas as pd
# # #import seaborn as sns
# # import matplotlib.pyplot as plt

# # import pyro
# # import pyro.distributions as dist
# # import pyro.distributions.constraints as constraints


# # DATA_URL = "https://d2hg8soec8ck9v.cloudfront.net/datasets/rugged_data.csv"
# # data = pd.read_csv(DATA_URL, encoding="ISO-8859-1")
# # df = data[["cont_africa", "rugged", "rgdppc_2000"]]

# # df = df[np.isfinite(df.rgdppc_2000)]
# # df["rgdppc_2000"] = np.log(df["rgdppc_2000"])

# # train = torch.tensor(df.values, dtype=torch.float)
# # is_cont_africa, ruggedness, log_gdp = train[:, 0], train[:, 1], train[:, 2]

# # def model(is_cont_africa, ruggedness, log_gdp=None):
# #     a = pyro.sample("a", dist.Normal(0., 10.))
# #     b_a = pyro.sample("bA", dist.Normal(0., 1.))
# #     b_r = pyro.sample("bR", dist.Normal(0., 1.))
# #     b_ar = pyro.sample("bAR", dist.Normal(0., 1.))
# #     sigma = pyro.sample("sigma", dist.Uniform(0., 10.))

# #     mean = a + b_a * is_cont_africa + b_r * ruggedness + b_ar * is_cont_africa * ruggedness

# #     with pyro.plate("data", len(ruggedness)):
# #         return pyro.sample("obs", dist.Normal(mean, sigma), obs=log_gdp)

# # print()

