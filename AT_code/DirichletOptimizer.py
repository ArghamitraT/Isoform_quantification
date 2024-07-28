import torch
import torch.optim as optim


class DirichletModel:
    def __init__(self, all_theta, all_alpha, exp_log_theta):
        self.all_theta = all_theta
        self.all_alpha = {k: torch.tensor(v, dtype=torch.float32, requires_grad=True) if not isinstance(v, torch.Tensor) else v for k, v in all_alpha.items()}
        self.exp_log_theta = exp_log_theta

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



    def compute_log_likelihood_ExpLogTheta(self):
        """
        Calculates the log PDF of Dirichlet function considering all samples.
        """
        log_likelihood = 0
        # Convert all_theta Counters into a tensor for PyTorch processing
        theta_tensor = self._counters_to_tensor(self.exp_log_theta)

        # Compute the log likelihood
        log_B_alpha = torch.sum(torch.lgamma(torch.tensor(list(self.all_alpha.values())))) - torch.lgamma(torch.sum(torch.tensor(list(self.all_alpha.values()))))
        log_likelihood = -log_B_alpha  # we need log(1/B(alpha))

        for f in range(theta_tensor.size(0)):  # Loop over each sample
            for isoform, alpha_i in self.all_alpha.items():  # Loop over each alpha_i
                isoform_index = list(self.all_alpha.keys()).index(isoform)
                temp_val = (alpha_i - 1) * (theta_tensor[f, isoform_index] + 1e-10)
                if torch.isnan(temp_val).any():
                    print('nan')
                    log_likelihood += (alpha_i - 1) * torch.log(torch.zeros(1).squeeze(0) + 1e-10)
                else:
                    log_likelihood += temp_val

        return log_likelihood

    def _counters_to_tensor(self, all_theta):
        """
        Convert all_theta (a dictionary of Counters) to a 2D tensor.
        
        Parameters:
        all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

        Returns:
        torch.Tensor: 2D tensor where each row corresponds to a sample and each column corresponds to an isoform.
        """
        # Use alpha to get all unique isoform names
        isoforms = sorted(list(self.all_alpha.keys()))

        # Prepare an empty list to store the tensor rows
        tensor_list = []

        # Iterate over each sample and convert its Counter to a row in the tensor
        for sample_key, counter in all_theta.items():
            # Initialize a row with zeros
            row = torch.full((len(isoforms),), torch.finfo(torch.float).tiny)

            # Fill in the row with the theta values from the Counter
            for isoform, theta_value in counter.items():
                if isoform in self.all_alpha:
                    i = isoforms.index(isoform)
                    row[i] = theta_value

            # Add the row to our list
            tensor_list.append(row)

        # Stack all rows to create a 2D tensor
        theta_tensor = torch.stack(tensor_list)

        return theta_tensor
