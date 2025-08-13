"""
Implementations of clustering algorithms from the original Mathematica library.
"""
import numpy as np

def kmeans(data, n_clusters, max_iter=100, initial_centroids=None):
    """
    Performs K-Means clustering on the input data.

    This implementation also returns the history of the centroids at each
    iteration, similar to the `KMeansEvolution` function in the
    Mathematica codebase.

    Args:
        data (np.ndarray): The input data, shape (n_samples, n_dimensions).
        n_clusters (int): The number of clusters (K).
        max_iter (int): The maximum number of iterations to perform.
        initial_centroids (np.ndarray, optional): The initial centroids, shape
            (n_clusters, n_dimensions). If None, centroids are initialized by
            randomly sampling from the data.

    Returns:
        tuple: A tuple containing:
            - labels (np.ndarray): The cluster label for each data point.
            - centroid_history (list): A list of the centroids at each iteration.
    """
    data = np.asarray(data)
    n_samples, n_dims = data.shape

    # Initialize centroids
    if initial_centroids is not None:
        centroids = np.copy(initial_centroids)
    else:
        # Randomly choose K data points as initial centroids
        random_indices = np.random.choice(n_samples, n_clusters, replace=False)
        centroids = data[random_indices]

    centroid_history = [centroids.copy()]

    for i in range(max_iter):
        # Assignment step: assign each data point to the closest centroid
        distances = np.zeros((n_samples, n_clusters))
        for k in range(n_clusters):
            # Euclidean distance squared (for efficiency, as sqrt is monotonic)
            distances[:, k] = np.sum((data - centroids[k])**2, axis=1)

        labels = np.argmin(distances, axis=1)

        # Update step: recalculate centroids
        new_centroids = np.zeros((n_clusters, n_dims))
        for k in range(n_clusters):
            points_in_cluster = data[labels == k]
            if len(points_in_cluster) > 0:
                new_centroids[k] = np.mean(points_in_cluster, axis=0)
            else:
                # Re-initialize centroid if a cluster becomes empty
                new_centroids[k] = data[np.random.choice(n_samples)]

        centroid_history.append(new_centroids.copy())

        # Check for convergence
        if np.allclose(centroids, new_centroids):
            break

        centroids = new_centroids

    # Final assignment with the converged centroids
    distances = np.zeros((n_samples, n_clusters))
    for k in range(n_clusters):
        distances[:, k] = np.sum((data - centroids[k])**2, axis=1)
    labels = np.argmin(distances, axis=1)

    return labels, centroid_history


def kmeans_nem(data, n_clusters, noise_sigma, cooling_schedule=2.0, max_iter=100, initial_centroids=None):
    """
    Performs K-Means clustering with NEM-conditioned noise.

    This corresponds to `KMeansNEMNoisyEvolution` in the Mathematica code.
    At each iteration, NEM-conditioned noise is added to the data before the
    assignment step.

    Args:
        data (np.ndarray): The input data.
        n_clusters (int): The number of clusters (K).
        noise_sigma (float): The base standard deviation for the noise generation.
        cooling_schedule (float): The exponent for the cooling schedule.
        max_iter (int): The maximum number of iterations.
        initial_centroids (np.ndarray, optional): The initial centroids.

    Returns:
        tuple: A tuple containing:
            - labels (np.ndarray): The final cluster labels.
            - centroid_history (list): The history of centroids.
    """
    data = np.asarray(data)
    n_samples, n_dims = data.shape

    # Need to import noise functions
    from . import noise

    # Initialize centroids
    if initial_centroids is not None:
        centroids = np.copy(initial_centroids)
    else:
        random_indices = np.random.choice(n_samples, n_clusters, replace=False)
        centroids = data[random_indices]

    centroid_history = [centroids.copy()]

    for i in range(max_iter):
        # Add NEM-conditioned noise
        iter_num = i + 1
        current_noise_sigma = noise_sigma / (iter_num ** cooling_schedule)

        nem_noise = noise.nd_advanced_nem_condition(data, centroids, current_noise_sigma)
        noisy_data = data + nem_noise

        # Assignment step with noisy data
        distances = np.zeros((n_samples, n_clusters))
        for k in range(n_clusters):
            distances[:, k] = np.sum((noisy_data - centroids[k])**2, axis=1)

        labels = np.argmin(distances, axis=1)

        # Update step with clean data
        new_centroids = np.zeros((n_clusters, n_dims))
        for k in range(n_clusters):
            points_in_cluster = data[labels == k]
            if len(points_in_cluster) > 0:
                new_centroids[k] = np.mean(points_in_cluster, axis=0)
            else:
                new_centroids[k] = data[np.random.choice(n_samples)]

        centroid_history.append(new_centroids.copy())

        # Check for convergence
        if np.allclose(centroids, new_centroids):
            break

        centroids = new_centroids

    # Final assignment with the converged centroids
    distances = np.zeros((n_samples, n_clusters))
    for k in range(n_clusters):
        distances[:, k] = np.sum((data - centroids[k])**2, axis=1)
    labels = np.argmin(distances, axis=1)

    return labels, centroid_history


def kmeans_noisy(data, n_clusters, noise_sigma, cooling_schedule=2.0, max_iter=100, initial_centroids=None):
    """
    Performs K-Means clustering with noise added at each iteration.

    This corresponds to `KMeansNoisyEvolution` in the Mathematica code.
    At each iteration, IID Gaussian noise is added to the data before the
    assignment step. The centroids are then updated based on the original
    clean data.

    Args:
        data (np.ndarray): The input data, shape (n_samples, n_dimensions).
        n_clusters (int): The number of clusters (K).
        noise_sigma (float): The base standard deviation for the noise generation.
        cooling_schedule (float): The exponent for the cooling schedule.
        max_iter (int): The maximum number of iterations.
        initial_centroids (np.ndarray, optional): The initial centroids.

    Returns:
        tuple: A tuple containing:
            - labels (np.ndarray): The final cluster labels.
            - centroid_history (list): The history of centroids.
    """
    data = np.asarray(data)
    n_samples, n_dims = data.shape

    # Initialize centroids
    if initial_centroids is not None:
        centroids = np.copy(initial_centroids)
    else:
        random_indices = np.random.choice(n_samples, n_clusters, replace=False)
        centroids = data[random_indices]

    centroid_history = [centroids.copy()]

    for i in range(max_iter):
        # Add noise to data for the assignment step
        iter_num = i + 1
        current_noise_sigma = noise_sigma / (iter_num ** cooling_schedule)
        noise = np.random.normal(0, current_noise_sigma, data.shape)
        noisy_data = data + noise

        # Assignment step with noisy data
        distances = np.zeros((n_samples, n_clusters))
        for k in range(n_clusters):
            distances[:, k] = np.sum((noisy_data - centroids[k])**2, axis=1)

        labels = np.argmin(distances, axis=1)

        # Update step with clean data
        new_centroids = np.zeros((n_clusters, n_dims))
        for k in range(n_clusters):
            points_in_cluster = data[labels == k]
            if len(points_in_cluster) > 0:
                new_centroids[k] = np.mean(points_in_cluster, axis=0)
            else:
                new_centroids[k] = data[np.random.choice(n_samples)]

        centroid_history.append(new_centroids.copy())

        # Check for convergence
        if np.allclose(centroids, new_centroids):
            break

        centroids = new_centroids

    # Final assignment with the converged centroids
    distances = np.zeros((n_samples, n_clusters))
    for k in range(n_clusters):
        distances[:, k] = np.sum((data - centroids[k])**2, axis=1)
    labels = np.argmin(distances, axis=1)

    return labels, centroid_history
