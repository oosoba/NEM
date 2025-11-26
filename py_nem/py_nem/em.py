"""
Implementations of the Expectation-Maximization (EM) algorithm and its variants.
"""
import numpy as np
from scipy import special
from scipy.optimize import minimize
from . import distributions
from . import noise


def _expected_value_censored_gamma(T, alpha, theta):
    """
    Calculates E[z | z >= T] for a Gamma(alpha, theta) distribution.
    """
    # Using the relationship between Mathematica's Gamma[a, z] and SciPy's gammaincc
    # Gamma[a, z] = Gamma(a) * gammaincc(a, z)
    z = T / theta
    numerator = theta * special.gamma(alpha + 1) * special.gammaincc(alpha + 1, z)
    denominator = special.gamma(alpha) * special.gammaincc(alpha, z)
    if denominator == 0:
        return T # Fallback to the threshold if probability is tiny
    return numerator / denominator


def gamma_em_censored(data, alpha, T, initial_theta, tol=1e-6, max_iter=100):
    """
    Performs EM for censored Gamma-distributed data.

    This corresponds to `GammaEM` in the Mathematica code. It estimates the
    scale parameter `theta` for a Gamma distribution where the shape `alpha`
    is known and data is censored at a threshold `T`.

    Args:
        data (np.ndarray): The input data, which may be censored.
        alpha (float): The known shape parameter of the Gamma distribution.
        T (float): The censoring threshold. Values >= T are censored.
        initial_theta (float): The initial guess for the scale parameter theta.
        tol (float): The tolerance for convergence.
        max_iter (int): The maximum number of iterations.

    Returns:
        list: A list of the estimated theta at each iteration.
    """
    data = np.asarray(data)
    theta = initial_theta
    theta_history = [theta]
    n_samples = len(data)

    for i in range(max_iter):
        # E-step: Calculate the expected value of each observation
        expected_zs = np.copy(data)
        censored_mask = (data >= T)

        if np.any(censored_mask):
            e_censored = _expected_value_censored_gamma(T, alpha, theta)
            expected_zs[censored_mask] = e_censored

        # M-step: Update theta
        theta_new = np.sum(expected_zs) / (n_samples * alpha)

        # Check for convergence
        if np.abs(theta_new - theta) < tol:
            theta_history.append(theta_new)
            break

        theta = theta_new
        theta_history.append(theta)

    return theta_history


def _cauchy_m_step_objective(params, data, weights):
    """Objective function for the M-step of Cauchy EM (to be minimized)."""
    loc, scale = params
    if scale <= 0:
        return np.inf # Penalize non-positive scale

    # Negative of the log-likelihood for the component
    log_likelihood = np.sum(weights * distributions.cauchy.logpdf(data, loc=loc, scale=scale))

    return -log_likelihood


def cauchy_em_1d(data, initial_locs, initial_scales, initial_alpha, tol=1e-6, max_iter=100):
    """
    Performs EM for a 1D Cauchy Mixture Model.

    This corresponds to `CmmEM`. It uses numerical optimization for the M-step.
    """
    data = np.asarray(data)
    locs = np.array(initial_locs)
    scales = np.array(initial_scales)
    alpha = initial_alpha

    params_history = [{'locs': locs.copy(), 'scales': scales.copy(), 'alpha': alpha}]

    for i in range(max_iter):
        # E-step: Calculate responsibilities
        resp_0 = distributions.cmm_posterior_probability(data, 0, alpha, locs, scales)
        resp_1 = distributions.cmm_posterior_probability(data, 1, alpha, locs, scales)
        responsibilities = np.vstack([resp_0, resp_1]).T

        # M-step: Update parameters
        # Update alpha
        sum_resp = np.sum(responsibilities, axis=0)
        alpha = sum_resp[0] / len(data)

        # Update locs and scales via numerical optimization
        for k in range(2):
            res = minimize(
                _cauchy_m_step_objective,
                x0=[locs[k], scales[k]],
                args=(data, responsibilities[:, k]),
                bounds=[(None, None), (1e-6, None)] # scale must be > 0
            )
            locs[k], scales[k] = res.x

        # Store history and check for convergence
        current_params = {'locs': locs.copy(), 'scales': scales.copy(), 'alpha': alpha}
        prev_p = np.concatenate([
            params_history[-1]['locs'],
            params_history[-1]['scales'],
            [params_history[-1]['alpha']]
        ])
        curr_p = np.concatenate([
            current_params['locs'],
            current_params['scales'],
            [current_params['alpha']]
        ])

        param_change = np.max(np.abs(prev_p - curr_p))
        params_history.append(current_params)

        if param_change < tol:
            break

    return params_history


def gamma_nem_censored(data, alpha, T, initial_theta, noise_sigma, cooling_schedule=2.0, tol=1e-6, max_iter=100):
    """
    Performs Noisy EM for censored Gamma-distributed data.

    This corresponds to `GammaNEM` in the Mathematica code.

    Args:
        data (np.ndarray): The input data.
        alpha (float): The known shape parameter.
        T (float): The censoring threshold.
        initial_theta (float): Initial guess for the scale parameter.
        noise_sigma (float): The base standard deviation for noise generation.
        cooling_schedule (float): The exponent for the cooling schedule.
        tol (float): Tolerance for convergence.
        max_iter (int): Maximum number of iterations.

    Returns:
        list: A list of the estimated theta at each iteration.
    """
    data = np.asarray(data)
    theta = initial_theta
    theta_history = [theta]
    n_samples = len(data)

    for i in range(max_iter):
        # 1. Add cooled IID noise
        iter_num = i + 1
        current_noise_sigma = noise_sigma / (iter_num ** cooling_schedule)
        noise = np.random.normal(0, current_noise_sigma, n_samples)
        noisy_data = data + noise
        noisy_data[noisy_data <= 0] = 1e-9 # Ensure positivity

        # E-step: Calculate expected values
        # The original code uses the clean data to determine censoring, but
        # the noisy data for the uncensored values.
        expected_zs = np.copy(noisy_data)
        censored_mask = (data >= T)

        if np.any(censored_mask):
            e_censored = _expected_value_censored_gamma(T, alpha, theta)
            expected_zs[censored_mask] = e_censored

        # M-step: Update theta
        theta_new = np.sum(expected_zs) / (n_samples * alpha)

        # Check for convergence
        if np.abs(theta_new - theta) < tol:
            theta_history.append(theta_new)
            break

        theta = theta_new
        theta_history.append(theta)

    return theta_history


def gmm_nem_1d(data, initial_mus, initial_sigmas, initial_alpha, noise_sigma, cooling_schedule=2.0, tol=1e-6, max_iter=100):
    """
    Performs the Noisy Expectation-Maximization algorithm for a 1D GMM.

    This corresponds to the `GmmNEM` function in the Mathematica code. It uses
    the NEM condition to add noise only during the sigma update step, which is
    a form of Expectation Conditional Maximization (ECM).

    Args:
        data (np.ndarray): The input data, a 1D array of floats.
        initial_mus (list): Initial values for the means of the GMM components.
        initial_sigmas (list): Initial values for the standard deviations.
        initial_alpha (float): Initial value for the mixing weight.
        noise_sigma (float): The base standard deviation for the noise generation.
        cooling_schedule (float): The exponent for the cooling schedule.
                                  Noise is cooled by `iter_num ** schedule`.
        tol (float): The tolerance for convergence.
        max_iter (int): The maximum number of iterations.

    Returns:
        list: A list of parameter dictionaries from each iteration.
    """
    data = np.asarray(data)

    mus = np.array(initial_mus)
    sigmas = np.array(initial_sigmas)
    alpha = initial_alpha

    params_history = [{'mus': mus.copy(), 'sigmas': sigmas.copy(), 'alpha': alpha}]

    for i in range(max_iter):
        # E-step: Calculate responsibilities using clean data
        resp_0 = distributions.gmm_posterior_probability(data, 0, alpha, mus, sigmas)
        resp_1 = distributions.gmm_posterior_probability(data, 1, alpha, mus, sigmas)
        responsibilities = np.vstack([resp_0, resp_1]).T
        sum_resp = np.sum(responsibilities, axis=0)

        # M-step (alpha, mus) using clean data
        alpha = sum_resp[0] / len(data)
        mus = np.sum(responsibilities * data[:, np.newaxis], axis=0) / sum_resp

        # M-step (sigmas) using noisy data
        # 1. Calculate cooled noise sigma
        # The (i+2) makes it start from 1 on the first iter and avoid divide by zero,
        # similar to how Mathematica's k++ works in the original code's loop.
        # Let's use i+1 to match the iteration number (1-based).
        iter_num = i + 1
        current_noise_sigma = noise_sigma / (iter_num ** cooling_schedule)

        # 2. Generate NEM-conditioned noise
        nem_noise = noise.advanced_nem_condition(data, mus, current_noise_sigma)
        noisy_data = data + nem_noise

        # 3. Update sigmas with noisy data
        diff = noisy_data[:, np.newaxis] - mus
        sigmas_sq = np.sum(responsibilities * diff**2, axis=0) / sum_resp
        sigmas = np.sqrt(sigmas_sq)

        # Store current parameters and check for convergence
        current_params = {'mus': mus.copy(), 'sigmas': sigmas.copy(), 'alpha': alpha}
        prev_p = np.concatenate([
            params_history[-1]['mus'],
            params_history[-1]['sigmas'],
            [params_history[-1]['alpha']]
        ])
        curr_p = np.concatenate([
            current_params['mus'],
            current_params['sigmas'],
            [current_params['alpha']]
        ])

        param_change = np.max(np.abs(prev_p - curr_p))
        params_history.append(current_params)

        if param_change < tol:
            break

    return params_history


def gmm_nem_nd(data, initial_mus, initial_sigmas, initial_weights, noise_sigma, cooling_schedule=2.0, tol=1e-6, max_iter=100):
    """
    Performs the Noisy EM algorithm for an N-dimensional GMM.

    This corresponds to the `GmmNDEMFLaw` function with a noise schedule.
    It uses the NEM condition to add noise only during the sigma update step.

    Args:
        data (np.ndarray): The input data, shape (n_samples, n_dimensions).
        initial_mus (np.ndarray): Initial means, shape (n_components, n_dimensions).
        initial_sigmas (np.ndarray): Initial covariance matrices.
        initial_weights (np.ndarray): Initial mixing weights.
        noise_sigma (float): The base standard deviation for the noise generation.
        cooling_schedule (float): The exponent for the cooling schedule.
        tol (float): The tolerance for convergence.
        max_iter (int): The maximum number of iterations.

    Returns:
        list: A list of parameter dictionaries from each iteration.
    """
    n_samples, n_dims = data.shape
    n_components = len(initial_mus)

    mus = np.copy(initial_mus)
    sigmas = np.copy(initial_sigmas)
    weights = np.copy(initial_weights)

    params_history = [{'mus': mus.copy(), 'sigmas': sigmas.copy(), 'weights': weights.copy()}]

    for i in range(max_iter):
        # E-step with clean data
        responsibilities = np.zeros((n_samples, n_components))
        for k in range(n_components):
            regularized_cov = sigmas[k] + np.eye(n_dims) * 1e-6
            responsibilities[:, k] = weights[k] * distributions.multivariate_normal_pdf(data, mu=mus[k], sigma=regularized_cov)

        resp_sum = np.sum(responsibilities, axis=1)[:, np.newaxis]
        resp_sum[resp_sum == 0] = 1e-9
        responsibilities /= resp_sum

        # M-step (weights, mus) with clean data
        sum_resp = np.sum(responsibilities, axis=0)
        weights = sum_resp / n_samples
        mus = responsibilities.T @ data / sum_resp[:, np.newaxis]

        # M-step (sigmas) with noisy data
        iter_num = i + 1
        current_noise_sigma = noise_sigma / (iter_num ** cooling_schedule)

        nem_noise = noise.nd_advanced_nem_condition(data, mus, current_noise_sigma)
        noisy_data = data + nem_noise

        for k in range(n_components):
            diff = noisy_data - mus[k]
            resp_col = responsibilities[:, k]
            sigmas[k] = (diff.T * resp_col) @ diff / sum_resp[k]

        # Store history and check for convergence
        current_params = {'mus': mus.copy(), 'sigmas': sigmas.copy(), 'weights': weights.copy()}

        prev_p_mus = params_history[-1]['mus']
        prev_p_sigmas = params_history[-1]['sigmas']
        prev_p_weights = params_history[-1]['weights']

        mus_change = np.max(np.abs(prev_p_mus - mus))
        sigmas_change = np.max(np.abs(prev_p_sigmas - sigmas))
        weights_change = np.max(np.abs(prev_p_weights - weights))
        param_change = max(mus_change, sigmas_change, weights_change)

        params_history.append(current_params)

        if param_change < tol:
            break

    return params_history


def gmm_em_nd(data, initial_mus, initial_sigmas, initial_weights, tol=1e-6, max_iter=100):
    """
    Performs the EM algorithm for an N-dimensional GMM with full covariance.

    This corresponds to the `GmmNDEMFLaw` function in the Mathematica code.

    Args:
        data (np.ndarray): The input data, shape (n_samples, n_dimensions).
        initial_mus (np.ndarray): Initial means, shape (n_components, n_dimensions).
        initial_sigmas (np.ndarray): Initial covariance matrices, shape (n_components, n_dimensions, n_dimensions).
        initial_weights (np.ndarray): Initial mixing weights, shape (n_components,).
        tol (float): The tolerance for convergence.
        max_iter (int): The maximum number of iterations.

    Returns:
        list: A list of parameter dictionaries from each iteration.
    """
    n_samples, n_dims = data.shape
    n_components = len(initial_mus)

    # Initialize parameters
    mus = np.copy(initial_mus)
    sigmas = np.copy(initial_sigmas)
    weights = np.copy(initial_weights)

    params_history = [{'mus': mus.copy(), 'sigmas': sigmas.copy(), 'weights': weights.copy()}]

    for i in range(max_iter):
        # E-step: Calculate responsibilities
        responsibilities = np.zeros((n_samples, n_components))
        for k in range(n_components):
            # Add a small regularization term to the diagonal of the covariance matrix
            # to prevent it from becoming singular.
            regularized_cov = sigmas[k] + np.eye(n_dims) * 1e-6
            responsibilities[:, k] = weights[k] * distributions.multivariate_normal_pdf(data, mu=mus[k], sigma=regularized_cov)

        # Normalize responsibilities
        resp_sum = np.sum(responsibilities, axis=1)[:, np.newaxis]
        # handle case where resp_sum is 0
        resp_sum[resp_sum == 0] = 1e-9
        responsibilities /= resp_sum

        # M-step: Update parameters
        sum_resp = np.sum(responsibilities, axis=0)

        # Update weights
        weights = sum_resp / n_samples

        # Update means
        mus = responsibilities.T @ data / sum_resp[:, np.newaxis]

        # Update covariance matrices
        for k in range(n_components):
            diff = data - mus[k]
            resp_col = responsibilities[:, k]
            sigmas[k] = (diff.T * resp_col) @ diff / sum_resp[k]

        # Store history and check for convergence
        current_params = {'mus': mus.copy(), 'sigmas': sigmas.copy(), 'weights': weights.copy()}

        # Check convergence on all parameters
        prev_p_mus = params_history[-1]['mus']
        prev_p_sigmas = params_history[-1]['sigmas']
        prev_p_weights = params_history[-1]['weights']

        mus_change = np.max(np.abs(prev_p_mus - mus))
        sigmas_change = np.max(np.abs(prev_p_sigmas - sigmas))
        weights_change = np.max(np.abs(prev_p_weights - weights))

        param_change = max(mus_change, sigmas_change, weights_change)

        params_history.append(current_params)

        if param_change < tol:
            break

    return params_history


def gmm_em_1d(data, initial_mus, initial_sigmas, initial_alpha, tol=1e-6, max_iter=100):
    """
    Performs the Expectation-Maximization algorithm for a 1D Gaussian Mixture Model.

    This corresponds to the `GmmEM` function in the Mathematica code.

    Args:
        data (np.ndarray): The input data, a 1D array of floats.
        initial_mus (list): Initial values for the means of the GMM components.
        initial_sigmas (list): Initial values for the standard deviations of the GMM components.
        initial_alpha (float): Initial value for the mixing weight of the first component.
        tol (float): The tolerance for convergence. The algorithm stops when the
                     maximum absolute change in parameters is less than this value.
        max_iter (int): The maximum number of iterations to perform.

    Returns:
        list: A list of parameter dictionaries from each iteration. Each dictionary
              contains 'mus', 'sigmas', and 'alpha'.
    """
    # Ensure data is a numpy array
    data = np.asarray(data)

    # Initialize parameters
    mus = np.array(initial_mus)
    sigmas = np.array(initial_sigmas)
    alpha = initial_alpha

    params_history = [{'mus': mus.copy(), 'sigmas': sigmas.copy(), 'alpha': alpha}]

    for i in range(max_iter):
        # E-step: Calculate responsibilities (posterior probabilities)
        resp_0 = distributions.gmm_posterior_probability(data, 0, alpha, mus, sigmas)
        resp_1 = distributions.gmm_posterior_probability(data, 1, alpha, mus, sigmas)

        responsibilities = np.vstack([resp_0, resp_1]).T

        # M-step: Update parameters

        # Update mixing weight (alpha)
        sum_resp = np.sum(responsibilities, axis=0)
        alpha = sum_resp[0] / len(data)

        # Update means
        mus = np.sum(responsibilities * data[:, np.newaxis], axis=0) / sum_resp

        # Update standard deviations
        diff = data[:, np.newaxis] - mus
        sigmas_sq = np.sum(responsibilities * diff**2, axis=0) / sum_resp
        sigmas = np.sqrt(sigmas_sq)

        # Store current parameters
        current_params = {'mus': mus.copy(), 'sigmas': sigmas.copy(), 'alpha': alpha}

        # Check for convergence
        prev_p = np.concatenate([
            params_history[-1]['mus'],
            params_history[-1]['sigmas'],
            [params_history[-1]['alpha']]
        ])
        curr_p = np.concatenate([
            current_params['mus'],
            current_params['sigmas'],
            [current_params['alpha']]
        ])

        param_change = np.max(np.abs(prev_p - curr_p))
        params_history.append(current_params)

        if param_change < tol:
            break

    return params_history
