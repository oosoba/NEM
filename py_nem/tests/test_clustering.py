import numpy as np
from py_nem import clustering

def test_kmeans_basic():
    """
    Tests the basic functionality of the kmeans algorithm.
    """
    # 1. Generate synthetic data with 3 distinct clusters
    np.random.seed(42)
    means = np.array([[0, 0], [5, 5], [-5, 5]])
    n_samples_per_cluster = 100
    n_clusters = len(means)

    data = np.vstack([
        np.random.randn(n_samples_per_cluster, 2) + mu for mu in means
    ])

    # 2. Run the kmeans algorithm
    labels, centroid_history = clustering.kmeans(data, n_clusters, max_iter=100)
    final_centroids = centroid_history[-1]

    # 3. Check the results
    # The exact labels are arbitrary, but the final centroids should be
    # close to the means of the original clusters.

    # Sort the estimated and true centroids to handle label switching
    sort_indices_est = np.argsort(final_centroids[:, 0])
    est_centroids_sorted = final_centroids[sort_indices_est]

    sort_indices_true = np.argsort(means[:, 0])
    true_means_sorted = means[sort_indices_true]

    assert est_centroids_sorted.shape == (n_clusters, 2)
    assert np.allclose(est_centroids_sorted, true_means_sorted, atol=0.5)


def test_kmeans_nem_basic():
    """
    Tests the functionality of the kmeans_nem algorithm.
    """
    # 1. Generate synthetic data with 3 distinct clusters
    np.random.seed(42)
    means = np.array([[0, 0], [5, 5], [-5, 5]])
    n_samples_per_cluster = 100
    n_clusters = len(means)

    data = np.vstack([
        np.random.randn(n_samples_per_cluster, 2) + mu for mu in means
    ])

    # 2. Run the NEM K-Means algorithm
    labels, centroid_history = clustering.kmeans_nem(
        data,
        n_clusters,
        noise_sigma=0.5,
        max_iter=100
    )
    final_centroids = centroid_history[-1]

    # 3. Check the results
    sort_indices_est = np.argsort(final_centroids[:, 0])
    est_centroids_sorted = final_centroids[sort_indices_est]

    sort_indices_true = np.argsort(means[:, 0])
    true_means_sorted = means[sort_indices_true]

    assert est_centroids_sorted.shape == (n_clusters, 2)
    assert np.allclose(est_centroids_sorted, true_means_sorted, atol=0.5)

    # Check that the labels are assigned correctly (indirectly) by seeing
    # if the number of points in each cluster is correct.
    unique, counts = np.unique(labels, return_counts=True)
    assert len(unique) == n_clusters
    # The counts should be roughly equal
    assert np.all(counts > n_samples_per_cluster * 0.8)
    assert np.all(counts < n_samples_per_cluster * 1.2)


def test_kmeans_noisy_basic():
    """
    Tests the functionality of the kmeans_noisy algorithm.
    """
    # 1. Generate synthetic data with 3 distinct clusters
    np.random.seed(42)
    means = np.array([[0, 0], [5, 5], [-5, 5]])
    n_samples_per_cluster = 100
    n_clusters = len(means)

    data = np.vstack([
        np.random.randn(n_samples_per_cluster, 2) + mu for mu in means
    ])

    # 2. Run the noisy kmeans algorithm
    labels, centroid_history = clustering.kmeans_noisy(
        data,
        n_clusters,
        noise_sigma=0.5,
        max_iter=100
    )
    final_centroids = centroid_history[-1]

    # 3. Check the results
    sort_indices_est = np.argsort(final_centroids[:, 0])
    est_centroids_sorted = final_centroids[sort_indices_est]

    sort_indices_true = np.argsort(means[:, 0])
    true_means_sorted = means[sort_indices_true]

    assert est_centroids_sorted.shape == (n_clusters, 2)
    assert np.allclose(est_centroids_sorted, true_means_sorted, atol=0.5)
