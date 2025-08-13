"""
Functions for generating noise for the Noisy EM algorithm.
"""
import numpy as np
from scipy.stats import truncnorm

def advanced_nem_condition(data, means, sigma):
    """
    Applies the Advanced NEM Condition to generate noise for each data point.

    This corresponds to the `advNemCond` function in the Mathematica code.
    For each data point, it calculates a valid interval for the noise
    based on the NEM condition and draws a noise value from a normal
    distribution truncated to that interval. The base normal distribution
    is N(0, sigma).

    The NEM condition ensures that the noise added to a data point does not
    change its classification with respect to the current component means.

    Args:
        data (np.ndarray): The input data, a 1D array.
        means (np.ndarray): The current means of the mixture components (e.g., [mu0, mu1]).
        sigma (float): The standard deviation of the noise distribution.

    Returns:
        np.ndarray: A 1D array of noise values, one for each data point.
    """
    noise = np.zeros_like(data, dtype=float)

    # Ensure means is a numpy array for vectorized operations
    means = np.asarray(means)

    for i, y in enumerate(data):
        # This logic determines the interval for the truncated normal distribution
        # based on the data point's position relative to the means.
        diffs = 2 * (means - y)
        max_val = np.max(diffs)
        min_val = np.min(diffs)

        lower_bound, upper_bound = 0.0, 0.0
        if min_val > 0:
            # y is to the left of both means
            lower_bound = 0
            upper_bound = min_val
        elif max_val < 0:
            # y is to the right of both means
            lower_bound = max_val
            upper_bound = 0

        # If the interval is not zero, draw noise from a truncated normal distribution
        if lower_bound != upper_bound and sigma > 0:
            # The truncnorm function in scipy takes the bounds (a, b) in terms of
            # standard deviations from the mean.
            a = (lower_bound - 0) / sigma
            b = (upper_bound - 0) / sigma
            noise[i] = truncnorm.rvs(a, b, loc=0, scale=sigma)

    return noise


def nd_advanced_nem_condition(data, means, sigma):
    """
    Applies the Advanced NEM Condition to each dimension of N-D data.

    This corresponds to `ndNemCond` in the Mathematica code, which applies
    the 1D NEM condition column-wise to the data.

    Args:
        data (np.ndarray): The input data, shape (n_samples, n_dimensions).
        means (np.ndarray): The current means, shape (n_components, n_dimensions).
        sigma (float): The standard deviation of the noise distribution.

    Returns:
        np.ndarray: The NEM-conditioned noise, shape (n_samples, n_dimensions).
    """
    n_samples, n_dims = data.shape
    full_noise = np.zeros_like(data, dtype=float)

    # Apply the 1D NEM condition to each dimension independently
    for d in range(n_dims):
        data_d = data[:, d]
        means_d = means[:, d]
        # Note: advanced_nem_condition returns noise, not noisy_data
        noise_d = advanced_nem_condition(data_d, means_d, sigma)
        full_noise[:, d] = noise_d

    return full_noise
