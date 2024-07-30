import torch
import torch.optim as optim


class DirichletModel:
    def __init__(self, all_theta, all_alpha, expectation_log_theta):
        self.all_theta = all_theta
        self.all_alpha = {k: torch.tensor(v, dtype=torch.float32, requires_grad=True) if not isinstance(v, torch.Tensor) else v for k, v in all_alpha.items()}
        self.expectation_log_theta = expectation_log_theta

        # Initialize a global optimizer for all alpha parameters
        self.optimizer = optim.Adam([param for param in self.all_alpha.values()], lr=0.01)  # lr is the learning rate

    def return_alpha(self):
        """
        Returns initialized alpha and isoform index for each sample.
        """
        return {sample_key: alpha.detach().numpy() for sample_key, alpha in
                self.all_alpha.items()}, self.all_isoform_indices

    def update_alpha(self, max_iterations=10, tolerance=1e-6):
        """
        Gradient Descent of alpha considering all samples with convergence criteria.
        :param max_iterations: int, maximum number of iterations to perform
        :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
        """
        previous_loss = float('inf')  # Initialize previous loss as infinity for comparison
        for iteration in range(max_iterations):
            self.optimizer.zero_grad()  # Clear the gradients from the previous step
            log_likelihood = self.compute_log_likelihood_ExpLogTheta() # (AT)
            #log_likelihood = self.compute_log_likelihood_LogTheta() # (AT)
            loss = -log_likelihood  # We want to maximize log likelihood, hence minimize the negative log likelihood
            loss.backward()  # Compute the gradient of the loss w.r.t. the parameters (alpha)
            self.optimizer.step()  # Perform a gradient descent step to update alpha
            
            # Check for convergence: if the change in loss is less than the tolerance, stop the loop
            if abs(previous_loss - loss.item()) < tolerance:
                print(f"GD Convergence reached after {iteration + 1} iterations.")
                break
            previous_loss = loss.item()  # Update previous loss to current loss
            
            # Optionally, print the current loss every few iterations to monitor progress
            #if iteration % 100 == 0:
            print(f"GD_Iteration {iteration}")
            print(f"GD_Current_Loss = {loss.item()}")
        
        alpha_nontensor = {k: v.item() for k, v in self.all_alpha.items()}

        return alpha_nontensor



    def compute_log_likelihood_LogTheta(self):
        """
        Calculates the log PDF of Dirichlet function considering all samples.
        """
        log_likelihood = 0
        # Convert all_theta Counters into a tensor for PyTorch processing
        theta_tensor = self._counters_to_tensor(self.all_theta)

        # Compute the log likelihood
        log_B_alpha = torch.sum(torch.lgamma(torch.tensor(list(self.all_alpha.values())))) - torch.lgamma(torch.sum(torch.tensor(list(self.all_alpha.values()))))
        log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

        for sample_key, isoform_dict in theta_tensor.items():  # Loop over each sample
            for isoform, alpha_i in self.all_alpha.items():  # Loop over each alpha_i
                # Handle case where isoform is not in the isoform_dict or theta_value is zero
                theta_value = isoform_dict.get(isoform, torch.tensor(1e-10, dtype=torch.float32))
                temp_val = (alpha_i - 1) * torch.log(theta_value + 1e-10)
                if torch.isnan(temp_val).any():
                    print('nan')
                    log_likelihood += (alpha_i - 1) * torch.log(torch.zeros(1).squeeze(0) + 1e-10)
                else:
                    log_likelihood += temp_val

        return log_likelihood
    


    def compute_log_likelihood_ExpLogTheta(self):
        """
        Calculates the log PDF of Dirichlet function considering all samples.
        """
        log_likelihood = 0
        # Convert all_theta Counters into a tensor for PyTorch processing
        theta_tensor = self._counters_to_tensor(self.expectation_log_theta)

        # Compute the log likelihood
        log_B_alpha = torch.sum(torch.lgamma(torch.tensor(list(self.all_alpha.values())))) - torch.lgamma(torch.sum(torch.tensor(list(self.all_alpha.values()))))
        log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

        for sample_key, isoform_dict in theta_tensor.items():  # Loop over each sample
            for isoform, alpha_i in self.all_alpha.items():  # Loop over each alpha_i
                # Handle case where isoform is not in the isoform_dict or theta_value is zero
                theta_value = isoform_dict.get(isoform, torch.tensor(1e-10, dtype=torch.float32))
                temp_val = (alpha_i - 1) * (theta_value + 1e-10)
                if torch.isnan(temp_val).any():
                    print('nan')
                    log_likelihood += (alpha_i - 1) * (torch.zeros(1).squeeze(0) + 1e-10)
                else:
                    log_likelihood += temp_val

        return log_likelihood

    
    
    def _counters_to_tensor(self, all_theta):
        """
        Convert all_theta (a dictionary of Counters) to a nested dictionary with tensors.
        
        Parameters:
        all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

        Returns:
        dict: Nested dictionary where the outer keys are sample names and inner keys are isoform names with theta values as tensors.
        """
        # Use alpha to get all unique isoform names
        isoforms = sorted(list(self.all_alpha.keys()))

        # Prepare the result dictionary
        result_dict = {}

        # Iterate over each sample and convert its Counter to a dictionary of isoform: theta_value pairs
        for sample_key, counter in all_theta.items():
            # Initialize a dictionary for the current sample
            sample_dict = {}

            # Fill in the dictionary with the theta values from the Counter
            for isoform in isoforms:
                theta_value = counter.get(isoform, torch.tensor(float(torch.finfo(torch.float).tiny)))
                sample_dict[isoform] = torch.tensor(theta_value, dtype=torch.float32)

            # Add the sample dictionary to the result dictionary
            result_dict[sample_key] = sample_dict

        return result_dict

