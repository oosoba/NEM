"""
Core statistical distributions and helper functions.
"""
import numpy as np
from scipy.stats import norm, gamma, cauchy, multivariate_normal

def multivariate_normal_pdf(x, mu, sigma):
    """
    Calculates the probability density function of a multivariate normal distribution.

    This is a wrapper around scipy.stats.multivariate_normal.pdf.

    Args:
        x: The point(s) at which to evaluate the PDF.
        mu: The mean vector of the distribution.
        sigma: The covariance matrix of the distribution.

    Returns:
        The value of the PDF at x.
    """
    return multivariate_normal.pdf(x, mean=mu, cov=sigma)


def cauchy_pdf(x, loc, scale):
    """
    Calculates the probability density function of a Cauchy distribution.
    """
    return cauchy.pdf(x, loc=loc, scale=scale)


def cmm_mixture_pdf(y, alpha, locs, scales):
    """
    Calculates the PDF of a 2-component Cauchy Mixture Model.
    """
    return alpha * cauchy_pdf(y, locs[0], scales[0]) + (1 - alpha) * cauchy_pdf(y, locs[1], scales[1])


def cmm_joint_probability(y, z, alpha, locs, scales):
    """
    Calculates the joint probability p(y, z) for a 2-component CMM.
    """
    p_z = np.array([alpha, 1 - alpha])
    if z < 0 or z > 1:
        raise ValueError("Component index z must be 0 or 1.")
    return p_z[z] * cauchy_pdf(y, locs[z], scales[z])


def cmm_posterior_probability(y, z, alpha, locs, scales):
    """
    Calculates the posterior probability p(z | y) for a 2-component CMM.
    """
    joint_prob = cmm_joint_probability(y, z, alpha, locs, scales)
    mixture_pdf = cmm_mixture_pdf(y, alpha, locs, scales)
    return np.divide(joint_prob, mixture_pdf, out=np.zeros_like(joint_prob), where=mixture_pdf!=0)


def normal_pdf(x, mu, sigma):
    """
    Calculates the probability density function of a normal distribution.

    This is a wrapper around scipy.stats.norm.pdf for consistency with the
    interfaces required by the EM algorithms.

    Args:
        x: The point(s) at which to evaluate the PDF.
        mu: The mean of the normal distribution.
        sigma: The standard deviation of the normal distribution.

    Returns:
        The value of the PDF at x.
    """
    return norm.pdf(x, loc=mu, scale=sigma)


def gmm_mixture_pdf(y, alpha, mus, sigmas):
    """
    Calculates the probability density function of a 2-component Gaussian Mixture Model.

    This corresponds to the `Fy` function in the Mathematica code.

    Args:
        y: The point(s) at which to evaluate the PDF.
        alpha: The mixing weight of the first component.
        mus: A list or array of the means of the two Gaussian components, e.g., [mu0, mu1].
        sigmas: A list or array of the standard deviations of the two components, e.g., [sigma0, sigma1].

    Returns:
        The value of the GMM PDF at y.
    """
    return alpha * normal_pdf(y, mus[0], sigmas[0]) + (1 - alpha) * normal_pdf(y, mus[1], sigmas[1])


def gmm_joint_probability(y, z, alpha, mus, sigmas):
    """
    Calculates the joint probability p(y, z) for a 2-component GMM.

    This corresponds to the `fyzi` function in the Mathematica code.
    p(y, z) = p(z) * p(y|z)

    Args:
        y: The point(s) at which to evaluate the joint probability.
        z: The component index (0 or 1).
        alpha: The mixing weight of the first component.
        mus: A list or array of the means, e.g., [mu0, mu1].
        sigmas: A list or array of the standard deviations, e.g., [sigma0, sigma1].

    Returns:
        The joint probability p(y, z).
    """
    p_z = np.array([alpha, 1 - alpha])
    if z < 0 or z > 1:
        raise ValueError("Component index z must be 0 or 1.")

    return p_z[z] * normal_pdf(y, mus[z], sigmas[z])


def gmm_posterior_probability(y, z, alpha, mus, sigmas):
    """
    Calculates the posterior probability p(z | y) for a 2-component GMM.

    This is also known as the responsibility of component z for data point y.
    This corresponds to the `pzjyi` function in the Mathematica code.
    p(z | y) = p(y, z) / p(y)

    Args:
        y: The point(s) at which to evaluate the posterior.
        z: The component index (0 or 1).
        alpha: The mixing weight of the first component.
        mus: A list or array of the means, e.g., [mu0, mu1].
        sigmas: A list or array of the standard deviations, e.g., [sigma0, sigma1].

    Returns:
        The posterior probability p(z | y).
    """
    joint_prob = gmm_joint_probability(y, z, alpha, mus, sigmas)
    mixture_pdf = gmm_mixture_pdf(y, alpha, mus, sigmas)

    # Avoid division by zero for points with zero probability
    # (e.g., far in the tails)
    return np.divide(joint_prob, mixture_pdf, out=np.zeros_like(joint_prob), where=mixture_pdf!=0)
