import torch
from torch.distributions import Dirichlet
import torch.optim as optim
import numpy as np

torch.manual_seed(42)
np.random.seed(42)

def dirichlet_log_likelihood(alpha, data):
    """
    Compute the log likelihood of the Dirichlet distribution manually.

    Parameters:
    alpha (torch.Tensor): The parameters of the Dirichlet distribution (K-dimensional).
    data (torch.Tensor): The data samples (N x K, where N is the number of samples, K is the dimension).

    Returns:
    log_likelihood (torch.Tensor): The log-likelihood of the data given the parameters alpha.
    """
    alpha_sum = torch.sum(alpha)
    log_likelihood = 2* (torch.lgamma(alpha_sum) - torch.sum(torch.lgamma(alpha)) ) # First part
    
    # Second part: sum over the samples
    for x in data:
        log_likelihood += torch.sum((alpha - 1) * torch.log(x))
    
    return log_likelihood


# Generate some synthetic data sampled from a Dirichlet distribution
true_alpha = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float32)  # True alpha
dirichlet_dist = Dirichlet(true_alpha)
data = dirichlet_dist.sample((2,))  # Generate 2 samples from the true distribution
print("Data samples:\n", data)

# Initialize alpha (parameters to be learned)
alpha = torch.tensor([2.0, 2.0, 2.0], requires_grad=True)  # Initial guess for alpha

# Define the optimizer (Gradient Descent)
optimizer = optim.Adam([alpha], lr=0.01)  # Using Adam optimizer for faster convergence

# Number of iterations for gradient descent
num_iterations = 1000

# Gradient descent loop
for iteration in range(num_iterations):
    optimizer.zero_grad()

    # Create a Dirichlet distribution with current alpha estimate
    dirichlet = Dirichlet(alpha)

    # Compute negative log-likelihood (we want to minimize this)
    # log_likelihood = dirichlet.log_prob(data).sum()
    log_likelihood = dirichlet_log_likelihood(alpha, data)

    # log_likelihood = dirichlet_log_likelihood(alpha, data).sum()
    loss = -log_likelihood  # Minimize negative log-likelihood

    # Backpropagation
    loss.backward()

    # Gradient descent step
    optimizer.step()

    # Ensure that alpha remains positive
    with torch.no_grad():
        alpha.clamp_(min=0.01)  # Avoid alpha becoming zero or negative
    # print(f"Iteration {iteration}, Loss: {loss.item()}, Alpha: {alpha.data}, Alpha_sum: {np.sum(alpha.data.numpy())}")

    # Print the loss every 100 iterations
    if iteration % 100 == 0:
        print(f"Iteration {iteration}, Loss: {loss.item()}, Alpha: {alpha.data}, Alpha_sum: {np.sum(alpha.data.numpy())}")

# Final estimated alpha
print("\nEstimated Alpha:\n", alpha.data)
