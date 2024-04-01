import torch
from scipy.special import gammaln
from collections import Counter
import torch_optimizer as optim



class DirichletModel:
    def __init__(self, all_theta):
        self.all_theta = all_theta
        #self.alpha = None  # Will be initialized later based on the number of unique isoforms
        self.isoforms = sorted(list(set(
            isoform for theta in self.all_theta.values() for isoform in theta.keys()
        )))
        self.isoform_to_index = {isoform: i for i, isoform in enumerate(self.isoforms)}
        self.alpha = torch.ones(len(self.isoforms), requires_grad=True)
        self.optimizer = torch.optim.Adam([self.alpha], lr=0.01)  # lr is the learning rate

    def initialize_alpha(self):
        # # Find the total number of unique isoforms across all samples
        # all_isoforms = set()
        # for sample_theta in self.all_theta.values():
        #     all_isoforms.update(sample_theta.keys())
        # # Initialize alpha
        # self.alpha = torch.ones(len(all_isoforms), requires_grad=True)
        # return self.alpha
        # self.alpha = torch.ones(len(self.isoforms), requires_grad=True)
        return self.alpha.detach().numpy(), self.isoform_to_index

    # this function does not use torch optimizer
    # def update_alpha(self,  learning_rate=0.01):
    #     # Reset gradients (important for multiple calls)
    #     if self.alpha.grad is not None:
    #         self.alpha.grad.zero_()
    #
    #     # Compute log likelihood
    #     log_likelihood = self.compute_log_likelihood()
    #
    #     # Compute gradients
    #     log_likelihood.backward()
    #
    #     # Update alpha using gradient descent
    #     with torch.no_grad():  # Temporarily disable gradient tracking
    #         self.alpha -= learning_rate * self.alpha.grad

    def update_alpha(self):
        self.optimizer.zero_grad()  # Clear the gradients from the previous step
        log_likelihood = self.compute_log_likelihood()
        loss = -log_likelihood  # We want to maximize log likelihood, hence minimize the negative log likelihood
        loss.backward()  # Compute the gradient of the loss w.r.t. the parameters (alpha)
        self.optimizer.step()  # Perform a gradient descent step to update alpha

    def compute_log_likelihood(self):
        # Convert all_theta Counters into a tensor for PyTorch processing
        theta_tensor = self._counters_to_tensor(self.all_theta)

        # Compute the log likelihood
        # alpha = self.alpha.detach().numpy() # to perform gammaln function you need to detach
        # log_B_alpha = torch.sum(gammaln(alpha)) - gammaln(torch.sum(alpha))
        log_B_alpha = torch.sum(torch.lgamma(self.alpha)) - torch.lgamma(torch.sum(self.alpha))
        log_likelihood = -log_B_alpha
        for f in range(theta_tensor.size(0)):  # Loop over each sample
            for i, alpha_i in enumerate(self.alpha):  # Loop over each alpha_i
                log_likelihood += (alpha_i - 1) * torch.log(theta_tensor[f, i])
        return log_likelihood

    def _counters_to_tensor(self, all_theta):
        # First, identify all unique isoforms across all samples.
        isoforms = set()
        for counter in all_theta.values():
            isoforms.update(counter.keys())
        isoforms = sorted(list(isoforms))

        # Create a mapping of isoform to column index
        # isoform_to_index = {isoform: i for i, isoform in enumerate(isoforms)}

        # Prepare an empty list to store the tensor rows
        tensor_list = []

        # Iterate over each sample and convert its Counter to a row in the tensor
        for sample_key, counter in all_theta.items():

            # Initialize a row with zeros
            #row = torch.zeros(len(isoforms))
            row = torch.full((len(isoforms),), torch.finfo(torch.float).tiny)

            # Fill in the row with the theta values from the Counter
            for isoform, theta_value in counter.items():
                i = self.isoform_to_index[isoform]
                row[i] = theta_value
                # Replace 0 with a very small positive value
                #row[i] = theta_value if theta_value > 0 else torch.finfo(row.dtype).tiny

            # Add the row to our list
            tensor_list.append(row)

        # Stack all rows to create a 2D tensor
        theta_tensor = torch.stack(tensor_list)

        return theta_tensor


# This method should convert your Counter objects into a tensor that can be processed by PyTorch
# The implementation of this will depend on the specific structure of your Counter objects and how you want to handle them
# ...

# optimizer = DirichletModel(all_theta=self.all_theta)
# optimizer.initialize_alpha()
# log_likelihood = optimizer.compute_log_likelihood()
#
# import torch
# from torch.autograd import Variable
# from torch.optim import Adam
# from scipy.special import gammaln, psi
#
#
# def dirichlet_log_likelihood(alpha, theta):
#     """
#     Compute the log likelihood of the Dirichlet distribution.
#
#     Parameters:
#     alpha : torch.Tensor
#         The concentration parameters of the Dirichlet distribution.
#     theta : torch.Tensor
#         The observed data for which we want to maximize the likelihood.
#
#     Returns:
#     log_likelihood : torch.Tensor
#         The log likelihood of the Dirichlet distribution for the data.
#     """
#     # Calculate the log of the Beta function using the gammaln function
#     log_B_alpha = torch.sum(gammaln(alpha)) - gammaln(torch.sum(alpha))
#
#     # Compute the log likelihood for all F observations
#     log_likelihood = -log_B_alpha  # No multiplication by F here
#     for f in range(theta.size(0)):  # Loop over each observation
#         for i in range(len(alpha)):  # Loop over each alpha_i
#             log_likelihood += (alpha[i] - 1) * torch.log(theta[f, i])
#     return log_likelihood
#
#
# #     F = theta.size(0)
# #     I = alpha.size(0)
# #
# #     # Calculate the log of the Beta function using the gammaln function
# #     log_B_alpha = torch.sum(gammaln(alpha)) - gammaln(torch.sum(alpha))
# #
# #     # Compute the log likelihood
# #     log_likelihood = F * (torch.log(torch.tensor(1.0)) - log_B_alpha)
# #     for i in range(I):
# #         log_likelihood += (alpha[i] - 1) * torch.sum(torch.log(theta[:, i]))
# #
# #     return log_likelihood
# #
# #
# # def dirichlet_log_likelihood(alpha, theta):
# #     # Calculate the log of the Beta function using the gammaln function
# #     log_B_alpha = torch.sum(gammaln(alpha)) - gammaln(torch.sum(alpha))
# #
# #     # Compute the log likelihood for all F observations
# #     log_likelihood = -log_B_alpha  # No multiplication by F here
# #     for f in range(theta.size(0)):  # Loop over each observation
# #         for i in range(len(alpha)):  # Loop over each alpha_i
# #             log_likelihood += (alpha[i] - 1) * torch.log(theta[f, i])
# #     return log_likelihood
#
#
# # Initialize alpha as a Variable with requires_grad=True to learn it
# # You should initialize it with some positive values (e.g., ones, or random positive numbers)
# alpha = Variable(torch.ones(I), requires_grad=True)
#
# # Define optimizer
# optimizer = Adam([alpha], lr=0.01)
#
# # Assume theta is your observed data and it's a tensor of shape (F, I)
#
# n_iterations = 10000  # Define the number of iterations
# for iteration in range(n_iterations):
#     optimizer.zero_grad()  # Clear the gradients
#     loss = -dirichlet_log_likelihood(alpha, theta)  # Negative log likelihood
#     loss.backward()  # Compute the gradient
#     optimizer.step()  # Perform a gradient descent step
#
#     # Optionally print the loss
#     if iteration % 1000 == 0:
#         print(f"Iteration {iteration}: loss {loss.item()}")
#
# # Check for convergence and print the result
# print(f"Optimized alpha: {alpha.data}")
