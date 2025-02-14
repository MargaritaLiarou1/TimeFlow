import numpy as np
import torch
import torch.nn as nn

np.random.seed(1)
torch.manual_seed(1)

class DensityEstimation(nn.Module):

    """
    DensityEstimation class implements the Real NVP transform and estimates the probability density
    of the cell states as a function of the CD markers.
    """

    def __init__(self, D, base_dist, num_coupling_layers, shift_networks, scale_networks):
        super(DensityEstimation, self).__init__()

        # Number of CD markers
        self.D = D
        # Base distribution for the Z variable in the manuscript
        self.base_dist = base_dist
        # Number of coupling layers
        self.num_coupling_layers = num_coupling_layers
        # List of neural networks for the shift operation  in each coupling layer
        self.shift_nn = torch.nn.ModuleList([shift_networks() for _ in range(num_coupling_layers)])
        # List of neural networks for the scale operation in each coupling layer
        self.scale_nn = torch.nn.ModuleList([scale_networks() for _ in range(num_coupling_layers)])

    def coupling(self, x, i, forward=True):

        """
        Coupling transformation for the given input.

        Args:
            x: input data for the coupling layer (starts with x)
            i: index of the current coupling layer
            forward: direction of the transformation (forward for inference/reverse for data generation, not used here)

        Returns:
            torch.Tensor: transformed data (z)
            torch.Tensor: output of the scale network, used for the change of variables formula
        """

        # x: input data for each coupling layer
        # x1: d-dimensional partition of x
        # x2: (d+1:D)-dimensional partition of x

        (x1, x2) = torch.chunk(x, 2, 1)

        # Affine transformation applied on the unchanged data partition x1
        # Forward direction for inference: x -> z  (z = g^-1(x))
        # z1 = g^-1(x1)

        if forward:
            z2 = (x2 - self.shift_nn[i](x1)) * torch.exp(-self.scale_nn[i](x1))
        else:
            # Reverse direction for data generation: z -> x  (x = g(z))
            # x2 = g(z2)
            z2 = torch.exp(self.scale_nn[i](x1)) * x2 + self.shift_nn[i](x1)

        # Concatenate the partitions (z1=x1)
        # Return z data and scale networks to use later for the change of variables formula
        return torch.cat((x1, z2), 1), self.scale_nn[i](x1)

    def permute_dimensions(self, x):

        """
        Permutation of dimensions of input data.

        Args:
            x: input data

        Returns:
            torch.Tensor: Data with reversed order of dimensions.
        """

        # Reverse the order of input dimensions
        return x.flip(1)

    def function_g(self, x):

        """
        Application of coupling transformations and permutations on the input data.

        Args:
            x: input data

        Returns:
            torch.Tensor: transformed data (z)
            torch.Tensor: log-determinant of the Jacobian of the transformation
        """

        # g: R^D -> R^D, z=g(x)
        # Log-determinant of the Jacobian initialized with zeros
        log_det_jac, z = x.new_zeros(x.shape[0]), x

        for i in range(self.num_coupling_layers):
            # Apply the transformation
            z, s = self.coupling(z, i, forward=True)
            # Reverse the order of input dimensions for the next coupling layer
            z = self.permute_dimensions(z)
            # Sum inputs of all scaling neural networks from all coupling layers
            log_det_jac = log_det_jac - s.sum(dim=1)

        return z, log_det_jac

    def forward(self, x, loss='avg'):

        """
        Forward pass through the density estimation model.

        Args:
            x: input data
            loss: average or sum of negative log-likelihood function (loss function)

        Returns:
            torch.Tensor: final loss value
        """

        # Transform input x using the function_g method, z=g(x)
        z, log_det_jac = self.function_g(x)

        # Sum or average of the negative log-likelihood function (loss function)
        if loss == 'sum':
            return -(self.base_dist.log_prob(z) + log_det_jac).sum()
        else:
            return -(self.base_dist.log_prob(z) + log_det_jac).mean()

    def log_probability_outputs(self, x):

        """
        Computes the log-probability of the input data.

        Args:
            x: input data

        Returns:
            torch.Tensor: Log-probability of the input data.
        """

        # Transform input x using the function_g method, z=g(x)
        z, log_det_jac = self.function_g(x)

        # Use change of variables formula to compute the log-probability
        log_prob = -(self.base_dist.log_prob(z) + log_det_jac)

        # Return the log-probability of the data
        return log_prob
