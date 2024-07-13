import torch
import torch.optim as optim


class DirichletModel:
    def __init__(self, all_theta, all_isoform_indices, all_alpha, exp_log_theta):
        self.all_theta = all_theta
        self.all_alpha = {k: torch.tensor(v, dtype=torch.float32, requires_grad=True) if isinstance(v, list) else v for
                          k, v in all_alpha.items()}
        self.exp_log_theta = exp_log_theta

        self.all_isoform_indices = all_isoform_indices

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
            log_likelihood = self.compute_log_likelihood_ExpLogTheta()
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

    def compute_log_likelihood(self):
        """
        Calculates the log PDF of Dirichlet function considering all samples.
        """
        log_likelihood = 0

        for sample_key in self.all_theta:
            # Convert all_theta Counters into a tensor for PyTorch processing
            theta_tensor = self._counters_to_tensor(self.all_theta[sample_key], self.all_isoform_indices[sample_key])

            # Get alpha for the specific sample
            alpha = self.all_alpha[sample_key]

            # Compute the log likelihood for the current sample
            log_B_alpha = torch.sum(torch.lgamma(alpha)) - torch.lgamma(torch.sum(alpha))
            sample_log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

            for isoform, index in self.all_isoform_indices[sample_key].items():
                alpha_i = alpha[index]
                theta_fi = theta_tensor[index]
                temp_val = (alpha_i - 1) * torch.log(theta_fi + 1e-10) # takes log[E(theta)]
                if torch.isnan(temp_val).any():
                    sample_log_likelihood += (alpha_i - 1) * torch.log(torch.zeros(1).squeeze(0) + 1e-10)
                else:
                    sample_log_likelihood += temp_val

            log_likelihood += sample_log_likelihood

        return log_likelihood

    def compute_log_likelihood_ExpLogTheta(self):
        """
        Calculates the log PDF of Dirichlet function considering all samples.
        """
        log_likelihood = 0

        for sample_key in self.all_theta:
            # Convert exp_log_thetaa Counters into a tensor for PyTorch processing
            theta_tensor = self._counters_to_tensor(self.exp_log_theta[sample_key], self.all_isoform_indices[sample_key])

            # Get alpha for the specific sample
            alpha = self.all_alpha[sample_key]

            # Compute the log likelihood for the current sample
            log_B_alpha = torch.sum(torch.lgamma(alpha)) - torch.lgamma(torch.sum(alpha))
            sample_log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

            for isoform, index in self.all_isoform_indices[sample_key].items():
                alpha_i = alpha[index]
                theta_fi = theta_tensor[index]
                temp_val = (alpha_i - 1) * (theta_fi + 1e-10) # takes E[log(theta)]
                if torch.isnan(temp_val).any():
                    sample_log_likelihood += (alpha_i - 1) * torch.log(torch.zeros(1).squeeze(0) + 1e-10)
                else:
                    sample_log_likelihood += temp_val

            log_likelihood += sample_log_likelihood

        return log_likelihood

    def _counters_to_tensor(self, sample_theta, sample_isoform_indices):
        # Initialize a row with tiny values
        row = torch.full((len(sample_isoform_indices),), torch.finfo(torch.float).tiny)

        # Fill in the row with the theta values from the Counter
        for isoform, index in sample_isoform_indices.items():
            if isoform in sample_theta:
                row[index] = sample_theta[isoform]

        return row
