import torch
import pyro
import pyro.distributions as dist
import pyro.optim as optim
from pyro.infer import SVI, Trace_ELBO

class DirichletModel:
    def __init__(self, all_theta, all_alpha, expectation_log_theta, GD_lr):
        self.all_theta = all_theta
        self.GD_lr = GD_lr
        self.initial_alpha = all_alpha

        # Initialize alpha parameters as a dictionary of tensors with positive constraints
        # self.alpha = {
        #     isoform: pyro.param(f"alpha_{isoform}", torch.tensor(value, dtype=torch.float32), constraint=dist.constraints.positive)
        #     for isoform, value in all_alpha.items()
        # }
        self.alpha = {
            isoform: torch.tensor(value, dtype=torch.float32, requires_grad=True)
            for isoform, value in all_alpha.items()
        }

        self.optimizer = optim.Adam({"lr": self.GD_lr})

    def model(self, theta):
        sorted_isoforms = sorted(self.initial_alpha.keys())
        alpha = torch.stack([self.alpha[isoform] for isoform in sorted_isoforms])
        with pyro.plate("data", len(theta)):
            pyro.sample("obs", dist.Dirichlet(alpha), obs=theta)

            # pyro.sample("obs", dist.Dirichlet(alpha).expand_by([len(theta)]).to_event(1), obs=theta)

        

    def guide(self, theta):
        sorted_isoforms = sorted(self.initial_alpha.keys())
        alpha_q = torch.stack([
            pyro.param(f"alpha_q_{isoform}", torch.tensor(self.initial_alpha[isoform], dtype=torch.float32), constraint=dist.constraints.positive)
            for isoform in sorted_isoforms
        ])
        # for name, value in pyro.get_param_store().items():
        #     print(name+"_GUIDE", pyro.param(name).data.cpu().numpy())

        with pyro.plate("data", len(theta)):
            pyro.sample("obs_guide", dist.Dirichlet(alpha_q))

            # pyro.sample("obs_guide", dist.Dirichlet(alpha_q).expand_by([len(theta)]).to_event(1))

    def update_alpha(self, max_iterations=10, tolerance=1e-6):
        svi = SVI(self.model, self.guide, self.optimizer, loss=Trace_ELBO())
        theta = self._theta_to_tensor(self.all_theta)  # Prepare the tensor once

        previous_loss = float('inf')  # Initialize previous loss as infinity for comparison
        for iteration in range(max_iterations):
            loss = svi.step(theta)

            # Print gradients and parameters for debugging
            # for name, value in pyro.get_param_store().items():
            #     print(name, pyro.param(name).data.cpu().numpy())

            # Check for convergence: if the change in loss is less than the tolerance, stop the loop
            if abs(previous_loss - loss) < tolerance:
                print(f"GD Convergence reached after {iteration + 1} iterations.")
                break
            previous_loss = loss  # Update previous loss to current loss

            # Optionally, print the current loss every few iterations to monitor progress
            print(f"GD_Iteration {iteration}")
            print(f"GD_Current_Loss = {loss}")

        # After updating self.alpha, convert to non-log space and non-tensor
        # Convert alpha_q to a dictionary and return
        alpha_nontensor = {isoform: pyro.param(f"alpha_q_{isoform}").detach().item() for isoform in self.initial_alpha}

        return alpha_nontensor

    def _theta_to_tensor(self, all_theta):
        """
        Convert all_theta (a dictionary of Counters) to a tensor for Pyro processing.

        Parameters:
        all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

        Returns:
        torch.Tensor: Tensor of theta values.
        """
        isoforms = sorted(self.initial_alpha.keys())
        theta_list = []

        for sample_key, isoform_dict in all_theta.items():
            sample_list = []
            for isoform in isoforms:
                if isoform in isoform_dict:
                    sample_list.append(torch.tensor(isoform_dict[isoform], dtype=torch.float32))
                else:
                    sample_list.append(torch.tensor(0.0, dtype=torch.float32))  # Missing isoforms are set to 0
            
            # Normalize the sample list to ensure it sums to 1
            sample_tensor = torch.stack(sample_list)
            sample_tensor /= sample_tensor.sum()
            theta_list.append(sample_tensor)

        return torch.stack(theta_list)

# Example usage
# all_theta = {
#     'sample1': {'isoform1': 0.2, 'isoform2': 0.3, 'isoform3': 0.5},  # Example with three isoforms
#     'sample2': {'isoform1': 0.1, 'isoform2': 0.4, 'isoform3': 0.5}
# }
# all_alpha = {'isoform1': 0.5, 'isoform2': 0.5, 'isoform3': 0.5}  # Dictionary of initial alpha values
# expectation_log_theta = None
# GD_lr = 0.01

# model = DirichletModel(all_theta, all_alpha, expectation_log_theta, GD_lr)
# optimized_alpha = model.update_alpha(max_iterations=1000)
# print("Optimized alpha parameters:", optimized_alpha)





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





# # # import torch
# # # import pyro
# # # import pyro.distributions as dist
# # # import pyro.optim as optim
# # # from pyro.infer import SVI, Trace_ELBO
# # # from pyro.infer.autoguide import AutoDelta

# # # class DirichletModel:
# # #     def __init__(self, all_theta, all_alpha, expectation_log_theta, GD_lr):
# # #         self.all_theta = all_theta
# # #         self.GD_lr = GD_lr
# # #         self.initial_alpha = all_alpha

# # #         # Initialize alpha parameters with Pyro param and positive constraint
# # #         self.alpha = {
# # #             k: pyro.param(f"alpha_{k}", torch.tensor(v, dtype=torch.float32), constraint=dist.constraints.positive)
# # #             for k, v in all_alpha.items()
# # #         }

# # #         self.optimizer = optim.Adam({"lr": self.GD_lr})

# # #     def model(self, theta):
# # #         # alpha = torch.tensor([self.alpha[isoform].detach() for isoform in sorted(self.initial_alpha.keys())])
# # #         # i = 0
# # #         # with pyro.plate("data", len(self.all_theta)):
# # #         #     if i == 0:
# # #         #         theta = self._theta_to_tensor(self.all_theta)
# # #         #         pyro.sample("obs", dist.Dirichlet(alpha), obs=theta[0])
# # #         #     else:
# # #         #         theta = self._theta_to_tensor(self.all_theta)
# # #         #         pyro.sample("obs", dist.Dirichlet(alpha), obs=theta[1])
# # #         #     i=+1
# # #         alpha = torch.stack([self.alpha[isoform] for isoform in sorted(self.initial_alpha.keys())])
# # #         with pyro.plate("data", len(theta)):
# # #             pyro.sample("obs", dist.Dirichlet(alpha), obs=theta)


# # #     def guide(self, theta):
# # #         guide = AutoDelta(self.model)
# # #         guide(theta)


# # #     def update_alpha(self, max_iterations=10, tolerance=1e-6):
# # #         """
# # #         Gradient Descent of alpha considering all samples with convergence criteria.
# # #         :param max_iterations: int, maximum number of iterations to perform
# # #         :param tolerance: float, tolerance to determine convergence (stop if change in loss is below this threshold)
# # #         """
# # #         svi = SVI(self.model, self.guide, self.optimizer, loss=Trace_ELBO())

# # #         alpha = torch.tensor([self.alpha[isoform].detach() for isoform in sorted(self.initial_alpha.keys())])
# # #         theta = self._theta_to_tensor(self.all_theta)

# # #         previous_loss = float('inf')  # Initialize previous loss as infinity for comparison
# # #         for iteration in range(max_iterations):
# # #             loss = svi.step(theta)

# # #             # Print gradients and parameters for debugging
# # #             for name, value in pyro.get_param_store().items():
# # #                 print(name, pyro.param(name).data.cpu().numpy())
            
# # #             # Check for convergence: if the change in loss is less than the tolerance, stop the loop
# # #             # if abs(previous_loss - loss) < tolerance:
# # #             #     print(f"GD Convergence reached after {iteration + 1} iterations.")
# # #             #     break
# # #             previous_loss = loss  # Update previous loss to current loss
            
# # #             # Optionally, print the current loss every few iterations to monitor progress
# # #             print(f"GD_Iteration {iteration}")
# # #             print(f"GD_Current_Loss = {loss}")
        
# # #         # After updating self.alpha, convert to non-log space and non-tensor
# # #         alpha_nontensor = {k: v.item() for k, v in self.alpha.items()}

# # #         return alpha_nontensor

# # #     def _theta_to_tensor(self, all_theta):
# # #         """
# # #         Convert all_theta (a dictionary of Counters) to a tensor for Pyro processing.

# # #         Parameters:
# # #         all_theta (dict): Dictionary where keys are sample names and values are Counters of isoforms and their theta values.

# # #         Returns:
# # #         torch.Tensor: Tensor of theta values.
# # #         """
# # #         isoforms = sorted(self.initial_alpha.keys())
# # #         theta_list = []

# # #         for sample_key, isoform_dict in all_theta.items():
# # #             sample_list = []
# # #             for isoform in isoforms:
# # #                 if isoform in isoform_dict:
# # #                     sample_list.append(torch.tensor(isoform_dict[isoform], dtype=torch.float32))
# # #                 else:
# # #                     sample_list.append(torch.tensor(0.0, dtype=torch.float32))  # Missing isoforms are set to 0
# # #             theta_list.append(torch.stack(sample_list))

# # #         return torch.stack(theta_list)

# # # # Example usage
# # # # all_theta = {
# # # #     'sample1': {'isoform1': 0.2, 'isoform2': 0.3, 'isoform3': 0.5},
# # # #     'sample2': {'isoform1': 0.1, 'isoform2': 0.4, 'isoform3': 0.5}
# # # # }
# # # # all_alpha = {'isoform1': 0.5, 'isoform2': 0.5, 'isoform3': 0.5}
# # # # expectation_log_theta = {
# # # #     'sample1': {'isoform1': torch.log(torch.tensor(0.2)), 'isoform2': torch.log(torch.tensor(0.3)), 'isoform3': torch.log(torch.tensor(0.5))},
# # # #     'sample2': {'isoform1': torch.log(torch.tensor(0.1)), 'isoform2': torch.log(torch.tensor(0.4)), 'isoform3': torch.log(torch.tensor(0.5))}
# # # # }
# # # # GD_lr = 0.01

# # # # model = DirichletModel(all_theta, all_alpha, expectation_log_theta, GD_lr)
# # # # optimized_alpha = model.update_alpha(max_iterations=1000)
# # # # print("Optimized alpha parameters:", optimized_alpha)





# # # # import torch
# # # # import pyro
# # # # import pyro.distributions as dist
# # # # import pyro.optim as optim
# # # # from pyro.infer import SVI, Trace_ELBO

# # # # def model(data, alpha_prior):
# # # #     alpha = pyro.param("alpha", alpha_prior)
# # # #     with pyro.plate("data", len(data)):
# # # #         pyro.sample("obs", dist.Dirichlet(alpha), obs=data)

# # # # def guide(data, alpha_prior):
# # # #     alpha_q = pyro.param("alpha_q", torch.ones_like(alpha_prior), constraint=dist.constraints.positive)
# # # #     pyro.sample("alpha", dist.Delta(alpha_q))

# # # # # Example data: 5 samples from a Dirichlet distribution
# # # # data = torch.tensor([
# # # #     [0.2, 0.5, 0.3],
# # # #     [0.1, 0.7, 0.2],
# # # #     [0.3, 0.4, 0.3],
# # # #     [0.2, 0.6, 0.2],
# # # #     [0.25, 0.5, 0.25]
# # # # ], dtype=torch.float32)

# # # # # Prior for alpha
# # # # alpha_prior = torch.tensor([1.0, 1.0, 1.0], dtype=torch.float32)

# # # # # Define the optimizer
# # # # adam_params = {"lr": 0.01}
# # # # optimizer = optim.Adam(adam_params)

# # # # # Define the SVI object
# # # # svi = SVI(model, guide, optimizer, loss=Trace_ELBO())

# # # # # Perform optimization
# # # # num_iterations = 1000
# # # # for step in range(num_iterations):
# # # #     loss = svi.step(data, alpha_prior)
# # # #     if step % 100 == 0:
# # # #         print(f"Step {step} : loss = {loss}")

# # # # # Retrieve the optimized alpha
# # # # optimized_alpha = pyro.param("alpha_q").detach().numpy()
# # # # print(f"Optimized alpha: {optimized_alpha}")
