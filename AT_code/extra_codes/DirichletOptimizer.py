import torch
import torch.optim as optim
import torch.nn.functional as F
import numpy as np


class DirichletModel:
    #def __init__(self, all_theta, all_alpha, expectation_log_theta, GD_lr, process='expectation_log_theta'):
    def __init__(self, all_alpha, GD_lr, process='expectation_log_theta'):
    
        
        ## onvert the dictionaries to tensors
        self.isoforms = list(all_alpha.keys())
        self.process = process
        self.alpha = torch.tensor(([all_alpha[isoform] for isoform in self.isoforms]), dtype=torch.float32)
        self.GD_lr = GD_lr
        self.log_alpha = torch.tensor(np.log([all_alpha[isoform] for isoform in self.isoforms]), dtype=torch.float32, requires_grad=True)
        
        # Optimizer
        self.optimizer = optim.Adam([self.log_alpha], lr=self.GD_lr)

    def update_alpha(self, expectation_log_theta, all_theta, max_iterations=10, tolerance=1e-6,):
        """
        Gradient Descent of alpha considering all samples with convergence criteria.
        :param max_iterations: int, maximum number of iterations to perform
        :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
        """
        # Create the tensor for expectation_log_theta
        epsilon = 1e-10
        if self.process == 'expectation_log_theta':
            self.data = torch.tensor(
                [
                    [sample.get(isoform, np.log(epsilon)) for isoform in self.isoforms] 
                    for sample in expectation_log_theta.values()
                ], 
                dtype=torch.float32)
        # Create the tensor for log_expectation_theta
        elif self.process == 'log_expectation_theta':
            self.data = torch.tensor(
                [
                    [sample.get(isoform, epsilon) for isoform in self.isoforms] 
                    for sample in all_theta.values()
                ], 
                dtype=torch.float32)
            
        num_iterations = max_iterations
        loss_history = []  # Initialize a list to store the loss at each iteration
    
        for iteration in range(num_iterations):
            self.optimizer.zero_grad()

            # Transform log_alpha back to alpha
            alpha = torch.exp(self.log_alpha)
            
            # Calculate log-likelihood
            ll = self.log_likelihood(self.data, alpha)
            
            # Since we want to maximize the log-likelihood, we minimize the negative log-likelihood
            loss = -ll
            
            # Backpropagation
            loss.backward()
            
            # Gradient descent step
            self.optimizer.step()
            # Store the current loss in the history list
            loss_history.append(loss.item())

            print(f"GD_Iteration {iteration}")
            print(f"GD_Current_Loss = {loss.item()}")

        # After updating self.log_alpha, convert to non-log space and non-tensor
        alpha_nontensor = {self.isoforms[i]: torch.exp(self.log_alpha[i]).item() for i in range(len(self.isoforms))}

        return alpha_nontensor, loss_history
    
    def log_likelihood(self, data, alpha):
        term1 = torch.lgamma(torch.sum(alpha))
        term2 = torch.sum(torch.lgamma(alpha))
        if self.process == 'expectation_log_theta':
            term3 = torch.sum((alpha - 1) * data, dim=1)
        elif self.process == 'log_expectation_theta':
            term3 = torch.sum((alpha - 1) * torch.log(data), dim=1)
        
        return term1 - term2 + term3.sum()

