import numpy as np
from py_nem import em

def test_gmm_em_1d_basic():
    """
    Tests the gmm_em_1d function with synthetic data.
    """
    # 1. Generate synthetic data from a known GMM
    np.random.seed(42)  # for reproducibility
    true_mus = [-2.0, 2.0]
    true_sigmas = [0.5, 0.5]
    true_alpha = 0.5
    n_samples = 2000

    # Generate data from each component
    n0 = int(n_samples * true_alpha)
    n1 = n_samples - n0
    data_0 = np.random.normal(true_mus[0], true_sigmas[0], n0)
    data_1 = np.random.normal(true_mus[1], true_sigmas[1], n1)
    data = np.concatenate([data_0, data_1])

    # 2. Define initial parameters (should be different from true parameters)
    initial_mus = [-1.0, 1.0]
    initial_sigmas = [1.0, 1.0]
    initial_alpha = 0.5

    # 3. Run the EM algorithm
    params_history = em.gmm_em_1d(
        data,
        initial_mus,
        initial_sigmas,
        initial_alpha,
        tol=1e-6,
        max_iter=100
    )

    # 4. Get the final estimated parameters
    final_params = params_history[-1]
    est_mus = final_params['mus']
    est_sigmas = final_params['sigmas']
    est_alpha = final_params['alpha']

    # 5. Assert that the estimated parameters are close to the true parameters
    # Sort by means to handle label switching
    sort_indices = np.argsort(est_mus)
    est_mus_sorted = est_mus[sort_indices]
    est_sigmas_sorted = est_sigmas[sort_indices]

    # If the first estimated mean is negative, we assume it corresponds to the first true component
    if est_mus_sorted[0] < 0:
        est_alpha_sorted = est_alpha if sort_indices[0] == 0 else 1 - est_alpha
    else: # This case should not happen with this seed and data
        est_alpha_sorted = est_alpha if sort_indices[0] == 1 else 1 - est_alpha


    assert np.allclose(est_mus_sorted, true_mus, atol=0.1)
    assert np.allclose(est_sigmas_sorted, true_sigmas, atol=0.1)
    assert np.allclose(est_alpha_sorted, true_alpha, atol=0.1)


def test_gmm_em_nd_basic():
    """
    Tests the gmm_em_nd function with synthetic 2D data.
    """
    # 1. Generate synthetic data
    np.random.seed(42)
    true_mus = np.array([[-2, -2], [2, 2]])
    true_sigmas = np.array([[[1, 0.5], [0.5, 1]], [[1, -0.5], [-0.5, 1]]])
    true_weights = np.array([0.5, 0.5])
    n_samples = 2000

    n0 = int(n_samples * true_weights[0])
    n1 = n_samples - n0
    data_0 = np.random.multivariate_normal(true_mus[0], true_sigmas[0], n0)
    data_1 = np.random.multivariate_normal(true_mus[1], true_sigmas[1], n1)
    data = np.concatenate([data_0, data_1])

    # 2. Define initial parameters
    initial_mus = np.array([[-1, -1], [1, 1]])
    initial_sigmas = np.array([np.eye(2), np.eye(2)])
    initial_weights = np.array([0.5, 0.5])

    # 3. Run the EM algorithm
    params_history = em.gmm_em_nd(
        data,
        initial_mus,
        initial_sigmas,
        initial_weights,
        tol=1e-4, # N-D is more sensitive, relax tolerance slightly
        max_iter=100
    )

    # 4. Get final estimated parameters
    final_params = params_history[-1]
    est_mus = final_params['mus']
    est_sigmas = final_params['sigmas']
    est_weights = final_params['weights']

    # 5. Assert that estimated parameters are close to true parameters
    # Sort by the first dimension of the means to handle label switching
    sort_indices = np.argsort(est_mus[:, 0])
    est_mus_sorted = est_mus[sort_indices]
    est_sigmas_sorted = est_sigmas[sort_indices]
    est_weights_sorted = est_weights[sort_indices]

    assert np.allclose(est_mus_sorted, true_mus, atol=0.2)
    assert np.allclose(est_sigmas_sorted, true_sigmas, atol=0.2)
    assert np.allclose(est_weights_sorted, true_weights, atol=0.1)


def test_gamma_em_censored():
    """
    Tests the EM algorithm for censored Gamma-distributed data.
    """
    # 1. Generate synthetic data from a Gamma distribution
    np.random.seed(42)
    true_alpha = 2.0
    true_theta = 1.5
    n_samples = 5000

    # In scipy, scale parameter is theta
    data = np.random.gamma(shape=true_alpha, scale=true_theta, size=n_samples)

    # 2. Censor the data
    T = np.percentile(data, 80) # Censor the top 20%
    data[data > T] = T

    # 3. Run the Gamma EM algorithm
    initial_theta = 1.0
    theta_history = em.gamma_em_censored(
        data,
        alpha=true_alpha,
        T=T,
        initial_theta=initial_theta,
        max_iter=100
    )

    est_theta = theta_history[-1]

    # 4. Assert that the estimated theta is close to the true theta
    assert np.allclose(est_theta, true_theta, atol=0.1)


def test_gamma_nem_censored():
    """
    Tests the Noisy EM algorithm for censored Gamma-distributed data.
    """
    # 1. Generate synthetic data
    np.random.seed(42)
    true_alpha = 2.0
    true_theta = 1.5
    n_samples = 5000

    data = np.random.gamma(shape=true_alpha, scale=true_theta, size=n_samples)

    # 2. Censor the data
    T = np.percentile(data, 80)
    data[data > T] = T

    # 3. Run the Gamma NEM algorithm
    initial_theta = 1.0
    theta_history = em.gamma_nem_censored(
        data,
        alpha=true_alpha,
        T=T,
        initial_theta=initial_theta,
        noise_sigma=0.1,
        max_iter=100
    )

    est_theta = theta_history[-1]

    # 4. Assert that the estimated theta is close to the true theta
    assert np.allclose(est_theta, true_theta, atol=0.1)


def test_cauchy_em_1d():
    """
    Tests the EM algorithm for a 1D Cauchy Mixture Model.
    """
    # 1. Generate synthetic data
    np.random.seed(42)
    true_locs = [-3.0, 3.0]
    true_scales = [0.5, 0.5]
    true_alpha = 0.5
    n_samples = 3000

    n0 = int(n_samples * true_alpha)
    n1 = n_samples - n0
    data_0 = np.random.standard_cauchy(n0) * true_scales[0] + true_locs[0]
    data_1 = np.random.standard_cauchy(n1) * true_scales[1] + true_locs[1]
    data = np.concatenate([data_0, data_1])

    # 2. Define initial parameters
    initial_locs = [-1.0, 1.0]
    initial_scales = [1.0, 1.0]
    initial_alpha = 0.5

    # 3. Run the Cauchy EM algorithm
    params_history = em.cauchy_em_1d(
        data,
        initial_locs,
        initial_scales,
        initial_alpha,
        tol=1e-4,
        max_iter=100
    )

    # 4. Get final estimated parameters
    final_params = params_history[-1]
    est_locs = final_params['locs']
    est_scales = final_params['scales']
    est_alpha = final_params['alpha']

    # 5. Assert that estimated parameters are close to true parameters
    sort_indices = np.argsort(est_locs)
    est_locs_sorted = est_locs[sort_indices]
    est_scales_sorted = est_scales[sort_indices]

    if est_locs_sorted[0] < 0:
        est_alpha_sorted = est_alpha if sort_indices[0] == 0 else 1 - est_alpha
    else:
        est_alpha_sorted = est_alpha if sort_indices[0] == 1 else 1 - est_alpha

    assert np.allclose(est_locs_sorted, true_locs, atol=0.3)
    assert np.allclose(est_scales_sorted, true_scales, atol=0.3)
    assert np.allclose(est_alpha_sorted, true_alpha, atol=0.2)


def test_gmm_nem_nd_basic():
    """
    Tests the gmm_nem_nd function with synthetic 2D data.
    """
    # 1. Generate synthetic data
    np.random.seed(42)
    true_mus = np.array([[-2, -2], [2, 2]])
    true_sigmas = np.array([[[1, 0.5], [0.5, 1]], [[1, -0.5], [-0.5, 1]]])
    true_weights = np.array([0.5, 0.5])
    n_samples = 2000

    n0 = int(n_samples * true_weights[0])
    n1 = n_samples - n0
    data_0 = np.random.multivariate_normal(true_mus[0], true_sigmas[0], n0)
    data_1 = np.random.multivariate_normal(true_mus[1], true_sigmas[1], n1)
    data = np.concatenate([data_0, data_1])

    # 2. Define initial parameters
    initial_mus = np.array([[-1, -1], [1, 1]])
    initial_sigmas = np.array([np.eye(2), np.eye(2)])
    initial_weights = np.array([0.5, 0.5])

    # 3. Run the NEM algorithm
    params_history = em.gmm_nem_nd(
        data,
        initial_mus,
        initial_sigmas,
        initial_weights,
        noise_sigma=0.1,
        tol=1e-4,
        max_iter=100
    )

    # 4. Get final estimated parameters
    final_params = params_history[-1]
    est_mus = final_params['mus']
    est_sigmas = final_params['sigmas']
    est_weights = final_params['weights']

    # 5. Assert that estimated parameters are close to true parameters
    sort_indices = np.argsort(est_mus[:, 0])
    est_mus_sorted = est_mus[sort_indices]
    est_sigmas_sorted = est_sigmas[sort_indices]
    est_weights_sorted = est_weights[sort_indices]

    assert np.allclose(est_mus_sorted, true_mus, atol=0.2)
    assert np.allclose(est_sigmas_sorted, true_sigmas, atol=0.2)
    assert np.allclose(est_weights_sorted, true_weights, atol=0.1)


def test_gmm_nem_1d_basic():
    """
    Tests the gmm_nem_1d function with synthetic data.
    """
    # 1. Generate synthetic data from a known GMM
    np.random.seed(42)  # for reproducibility
    true_mus = [-2.0, 2.0]
    true_sigmas = [0.5, 0.5]
    true_alpha = 0.5
    n_samples = 2000

    n0 = int(n_samples * true_alpha)
    n1 = n_samples - n0
    data_0 = np.random.normal(true_mus[0], true_sigmas[0], n0)
    data_1 = np.random.normal(true_mus[1], true_sigmas[1], n1)
    data = np.concatenate([data_0, data_1])

    # 2. Define initial parameters
    initial_mus = [-1.0, 1.0]
    initial_sigmas = [1.0, 1.0]
    initial_alpha = 0.5

    # 3. Run the NEM algorithm
    params_history = em.gmm_nem_1d(
        data,
        initial_mus,
        initial_sigmas,
        initial_alpha,
        noise_sigma=0.1,
        cooling_schedule=2.0,
        tol=1e-6,
        max_iter=100
    )

    # 4. Get the final estimated parameters
    final_params = params_history[-1]
    est_mus = final_params['mus']
    est_sigmas = final_params['sigmas']
    est_alpha = final_params['alpha']

    # 5. Assert that the estimated parameters are close to the true parameters
    sort_indices = np.argsort(est_mus)
    est_mus_sorted = est_mus[sort_indices]
    est_sigmas_sorted = est_sigmas[sort_indices]

    if est_mus_sorted[0] < 0:
        est_alpha_sorted = est_alpha if sort_indices[0] == 0 else 1 - est_alpha
    else:
        est_alpha_sorted = est_alpha if sort_indices[0] == 1 else 1 - est_alpha

    assert np.allclose(est_mus_sorted, true_mus, atol=0.1)
    assert np.allclose(est_sigmas_sorted, true_sigmas, atol=0.1)
    assert np.allclose(est_alpha_sorted, true_alpha, atol=0.1)
