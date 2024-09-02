import torch
import torch.optim as optim
import torch.nn.functional as F



class DirichletModel:
    def __init__(self, all_theta, all_alpha, expectation_log_theta, GD_lr):
        
        # Convert the dictionaries to tensors
        isoforms = list(all_alpha.keys())
        self.data = torch.tensor([[sample[isoform] for isoform in isoforms] for sample in all_theta.values()], dtype=torch.float32)
        self.alpha = torch.tensor([all_alpha[isoform] for isoform in isoforms], dtype=torch.float32, requires_grad=True)
        self.GD_lr = GD_lr

        # Optimizer
        self.optimizer = optim.Adam([self.alpha], lr=self.GD_lr)

    def update_alpha(self, max_iterations=10, tolerance=1e-6):
        """
        Gradient Descent of alpha considering all samples with convergence criteria.
        :param max_iterations: int, maximum number of iterations to perform
        :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
        """
        
        num_iterations = 10
        for iteration in range(num_iterations):
            self.optimizer.zero_grad()
            
            # Calculate log-likelihood
            ll = self.log_likelihood(self.data, self.alpha)
            
            # Since we want to maximize the log-likelihood, we minimize the negative log-likelihood
            loss = -ll
            
            # Backpropagation
            loss.backward()
            
            # Gradient descent step
            self.optimizer.step()
            
            # Ensure alpha remains positive
            with torch.no_grad():
                self.alpha.clamp_(min=1e-3)
            
            # Print the loss and alpha parameters every 100 iterations
            #if (iteration + 1) % 100 == 0:
            print(f'Iteration {iteration + 1}/{num_iterations}, Loss: {loss.item()}, Alpha: {self.alpha.data}')

        # Final alpha parameters
        print(f'Estimated Alpha: {self.alpha.data}')
        # After updating self.log_alpha, convert to non-log space and non-tensor
        alpha_nontensor = {k: torch.exp(v).item() for k, v in self.log_alpha.items()}

        return alpha_nontensor
    
    def log_likelihood(self, data, alpha):
        term1 = torch.lgamma(torch.sum(alpha))
        term2 = torch.sum(torch.lgamma(alpha))
        term3 = torch.sum((alpha - 1) * torch.log(data), dim=1)
        return term1 - term2 + term3.mean()


# Example usage
# all_theta = {
#     'sample1': {'isoform1': 0.2, 'isoform2': 0.5, 'isoform3': 0.3},  # Example with three isoforms
#     'sample2': {'isoform1': 0.1, 'isoform2': 0.3, 'isoform3': 0.6},
#     'sample3': {'isoform1': 0.3, 'isoform2': 0.4, 'isoform3': 0.3}
# }
# all_alpha = {'isoform1': 1, 'isoform2': 1, 'isoform3': 1}  # Dictionary of initial alpha values
# expectation_log_theta = None
# GD_lr = 0.01

# model = DirichletModel(all_theta, all_alpha, expectation_log_theta, GD_lr)
# optimized_alpha = model.update_alpha(max_iterations=1000)
# print("Optimized alpha parameters:", optimized_alpha)

    # def compute_log_likelihood_LogTheta(self):
    #     """
    #     Calculates the log PDF of Dirichlet function considering all samples.
    #     """
    #     log_likelihood = 0
    #     # Convert all_theta Counters into a tensor for PyTorch processing
    #     theta_tensor = self._counters_to_tensor(self.all_theta)
    #     all_alpha = {k: torch.exp(v) for k, v in self.log_alpha.items()}

    #     # Compute the log likelihood
    #     log_B_alpha = torch.sum(torch.lgamma(torch.tensor(list(all_alpha.values())))) - torch.lgamma(torch.sum(torch.tensor(list(all_alpha.values()))))
    #     log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

    #     for sample_key, isoform_dict in theta_tensor.items():  # Loop over each sample
    #         for isoform, alpha_i in all_alpha.items():  # Loop over each alpha_i
    #             # Handle case where isoform is not in the isoform_dict or theta_value is zero
    #             theta_value = isoform_dict.get(isoform, torch.tensor(1e-10, dtype=torch.float32))
    #             temp_val = (alpha_i - 1) * torch.log(theta_value + 1e-10)
    #             if torch.isnan(temp_val).any():
    #                 print('nan')
    #                 log_likelihood += (alpha_i - 1) * torch.log(torch.zeros(1).squeeze(0) + 1e-10)
    #             else:
    #                 log_likelihood += temp_val

    #     return log_likelihood



    # def compute_log_likelihood_ExpLogTheta(self):
    #     """
    #     Calculates the log PDF of Dirichlet function considering all samples.
    #     """
    #     log_likelihood = 0
    #     # Convert all_theta Counters into a tensor for PyTorch processing
    #     theta_tensor = self._counters_to_tensor(self.expectation_log_theta)
    #     all_alpha = {k: torch.exp(v) for k, v in self.log_alpha.items()}
        

    #     # Compute the log likelihood
    #     log_B_alpha = torch.sum(torch.lgamma(torch.tensor(list(all_alpha.values())))) - torch.lgamma(torch.sum(torch.tensor(list(all_alpha.values()))))
    #     log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

    #     for sample_key, isoform_dict in theta_tensor.items():  # Loop over each sample
    #         for isoform, alpha_i in all_alpha.items():  # Loop over each alpha_i
    #             # Handle case where isoform is not in the isoform_dict or theta_value is zero
    #             theta_value = isoform_dict.get(isoform, torch.tensor(-1e10, dtype=torch.float32))
    #             temp_val = (alpha_i - 1) * (theta_value + 1e-10)
    #             if torch.isnan(temp_val).any():
    #                 print('nan')
    #                 log_likelihood += (alpha_i - 1) * (torch.zeros(1).squeeze(0) + 1e-10)
    #             else:
    #                 log_likelihood += temp_val

    #     return log_likelihood

    
    
    # def _counters_to_tensor(self, all_theta):
    #     """
    #     Convert all_theta (a dictionary of Counters) to a nested dictionary with tensors.
        
    #     Parameters:
    #     all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

    #     Returns:
    #     dict: Nested dictionary where the outer keys are sample names and inner keys are isoform names with theta values as tensors.
    #     """
    #     # Use alpha to get all unique isoform names
    #     isoforms = sorted(list(self.initial_alpha.keys()))

    #     # Prepare the result dictionary
    #     result_dict = {}

    #     # Iterate over each sample and convert its Counter to a dictionary of isoform: theta_value pairs
    #     for sample_key, counter in all_theta.items():
    #         # Initialize a dictionary for the current sample
    #         sample_dict = {}

    #         # Fill in the dictionary with the theta values from the Counter
    #         for isoform in isoforms:
    #             theta_value = counter.get(isoform, torch.tensor(float(torch.finfo(torch.float).tiny)))
    #             sample_dict[isoform] = torch.tensor(theta_value, dtype=torch.float32)

    #         # Add the sample dictionary to the result dictionary
    #         result_dict[sample_key] = sample_dict

    #     return result_dict

