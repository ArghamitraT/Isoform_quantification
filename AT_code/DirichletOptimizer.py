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

    def return_alpha(self):
        """
        returns initialized alpha and isoform index
        """
        return self.alpha.detach().numpy(), self.isoform_to_index

    """ this GD of alpha function does not use torch optimizer """
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
        """
        Gradient Descent of alpha
        """
        self.optimizer.zero_grad()  # Clear the gradients from the previous step
        log_likelihood = self.compute_log_likelihood()
        loss = -log_likelihood  # We want to maximize log likelihood, hence minimize the negative log likelihood
        loss.backward()  # Compute the gradient of the loss w.r.t. the parameters (alpha)
        self.optimizer.step()  # Perform a gradient descent step to update alpha

    def compute_log_likelihood(self):
        """
        Calculates the log PDF of dirichlet function
        """

        # Convert all_theta Counters into a tensor for PyTorch processing
        theta_tensor = self._counters_to_tensor(self.all_theta)

        """ EXPLANATION:
            * Eqn --> \hat\alpha = \arg \max_{\alpha} \sum_{f=1}^{F}log(\frac{1}{B(\alpha)}) 
            \sum_{i=1}^{I} {(\alpha_i - 1)} log(\theta_{fi})} """

        # Compute the log likelihood
        log_B_alpha = torch.sum(torch.lgamma(self.alpha)) - torch.lgamma(torch.sum(self.alpha))
        log_likelihood = -log_B_alpha # we need log(1/B(alpha))

        """DAVID: try to avoid the for loop"""
        for f in range(theta_tensor.size(0)):  # Loop over each sample
            for i, alpha_i in enumerate(self.alpha):  # Loop over each alpha_i
                temp_val = (alpha_i - 1) * torch.log(theta_tensor[f, i]+1e-10)
                if torch.isnan(temp_val).any():
                    print('nan')
                    log_likelihood +=(alpha_i - 1) * torch.log(torch.zeros(1).squeeze(0)+1e-10)
                else:
                    log_likelihood += temp_val
        return log_likelihood

    def _counters_to_tensor(self, all_theta):
        # First, identify all unique isoforms across all samples.
        isoforms = set()
        for counter in all_theta.values():
            isoforms.update(counter.keys())
        isoforms = sorted(list(isoforms))

        # Prepare an empty list to store the tensor rows
        tensor_list = []

        # Iterate over each sample and convert its Counter to a row in the tensor
        for sample_key, counter in all_theta.items():

            # Initialize a row with zeros
            row = torch.full((len(isoforms),), torch.finfo(torch.float).tiny)

            # Fill in the row with the theta values from the Counter
            for isoform, theta_value in counter.items():
                i = self.isoform_to_index[isoform]
                row[i] = theta_value

            # Add the row to our list
            tensor_list.append(row)

        # Stack all rows to create a 2D tensor
        theta_tensor = torch.stack(tensor_list)

        return theta_tensor

