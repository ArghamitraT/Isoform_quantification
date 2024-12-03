import torch
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
from torch.distributions import Dirichlet

class DirichletModel:
    def __init__(self, all_alpha, GD_lr,  built_in=0, process='expectation_log_theta'):
        self.built_in = built_in
        if built_in:
            self.process = process
            self.alpha = torch.tensor(all_alpha, dtype=torch.float32, requires_grad=True)
            self.GD_lr = GD_lr
            # Optimizer
            self.optimizer = optim.Adam([self.alpha], lr=self.GD_lr)

        else:
            self.process = process
            self.alpha = torch.tensor(all_alpha, dtype=torch.float32)
            self.GD_lr = GD_lr
            self.log_alpha = torch.tensor(np.log(all_alpha), dtype=torch.float32, requires_grad=True)

            # Optimizer
            self.optimizer = optim.Adam([self.log_alpha], lr=self.GD_lr)

    def update_alpha(self, all_unique_isoforms, theta_names, expectation_log_theta, all_theta, max_iterations=10, tolerance=1e-6,):
        """
        Gradient Descent of alpha considering all samples with convergence criteria.
        :param max_iterations: int, maximum number of iterations to perform
        :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
        """
    
        if self.process == 'expectation_log_theta':
            working_data = expectation_log_theta
        # Create the tensor for log_expectation_theta
        elif self.process == 'log_expectation_theta':
            working_data = all_theta
        
        epsilon = 1e-10  # Small value to use when an isoform is missing
        unique_isoforms = np.array(all_unique_isoforms)  # Array of all unique isoforms
        num_isoforms = len(unique_isoforms)
        num_samples = len(working_data)

        if self.process == 'expectation_log_theta':
            aligned_data = np.full((num_samples, num_isoforms), np.log(epsilon), dtype=np.float32)
        elif self.process == 'log_expectation_theta':
            aligned_data = np.full((num_samples, num_isoforms), (epsilon), dtype=np.float32)

        # Step 1: Sort `unique_isoforms` and create an index mapping using `np.searchsorted`
        sort_indices = np.argsort(unique_isoforms)
        sorted_unique_isoforms = unique_isoforms[sort_indices]

        # Step 2: Create a list of sample isoforms and their values for vectorized processing
        sample_keys = list(working_data.keys())
        all_sample_isoforms = [theta_names[sample_key] for sample_key in sample_keys]
        all_theta_values = [working_data[sample_key] for sample_key in sample_keys]

        # Step 3: Flatten the isoform names and values across all samples for vectorized processing
        flat_isoforms = np.concatenate(all_sample_isoforms)
        flat_values = np.concatenate(all_theta_values)

        # Step 4: Use `np.searchsorted` to find indices of `flat_isoforms` in `sorted_unique_isoforms`
        flat_indices = np.searchsorted(sorted_unique_isoforms, flat_isoforms)

        # Step 5: Generate row indices for assignment
        sample_indices = np.concatenate([[i] * len(isoforms) for i, isoforms in enumerate(all_sample_isoforms)])

        # Step 6: Perform the vectorized assignment using `sort_indices` to map back to the original order
        aligned_data[sample_indices, sort_indices[flat_indices]] = flat_values
        for i in range(len(aligned_data)):
            aligned_data[i] = aligned_data[i]/np.sum(aligned_data[i])

        # Convert the final array to a PyTorch tensor
        self.data = torch.tensor(aligned_data, dtype=torch.float32)
        print(f"sum_sample1 {torch.sum(self.data[0])}")
        print(f"sum_sample2 {torch.sum(self.data[1])}")
            
        num_iterations = max_iterations
        loss_history = []  # Initialize a list to store the loss at each iteration

        if self.built_in:
            for iteration in range(num_iterations):
                self.optimizer.zero_grad()
                dirichlet = Dirichlet(self.alpha)
                log_likelihood = torch.zeros(1)
                for i in range(self.data.shape[0]):
                    log_likelihood += dirichlet.log_prob(self.data[i])
                # x= dirichlet.log_prob(self.data).sum()
                loss = -log_likelihood  # Minimize negative log-likelihood
                # Backpropagation
                loss.backward()
                
                # Gradient descent step
                self.optimizer.step()
                # Store the current loss in the history list
                loss_history.append(loss.item())

                print(f"GD_Iteration {iteration}")
                print(f"GD_Current_Loss = {loss.item()}")
            
            alpha_nontensor = (self.alpha).detach().numpy()
        else:
            for iteration in range(num_iterations):
                self.optimizer.zero_grad()
                # Transform log_alpha back to alpha
                alpha = torch.exp(self.log_alpha)
                # Calculate log-likelihood
                ll = self.log_likelihood(self.data, alpha)
                # Since we want to maximize the log-likelihood, we minimize the negative log-likelihood
                loss = -ll
                loss.backward()
                self.optimizer.step()
                loss_history.append(loss.item())
                print(f"GD_Iteration {iteration}")
                print(f"GD_Current_Loss = {loss.item()}")

            # After updating self.log_alpha, convert to non-log space and non-tensor
            alpha_nontensor = torch.exp(self.log_alpha).detach().numpy()

        return alpha_nontensor, loss_history
    
    def log_likelihood(self, data, alpha):
        term1 = torch.lgamma(torch.sum(alpha))
        term2 = torch.sum(torch.lgamma(alpha))
        if self.process == 'expectation_log_theta':
            term3 = torch.sum((alpha - 1) * data, dim=1)
        elif self.process == 'log_expectation_theta':
            term3 = torch.sum((alpha - 1) * torch.log(data), dim=1)
        
        log_liklihood = term3.shape[0]*(term1 - term2) + torch.sum(term3)
        
        return log_liklihood
