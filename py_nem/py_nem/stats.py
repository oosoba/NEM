"""
Implementations of statistical tools like bootstrapping.
"""
import numpy as np
import matplotlib.pyplot as plt

def bootstrap_ci(samples, n_bootstraps=200, sample_size=53, ci=95, outlier_trim=2):
    """
    Creates confidence intervals for each aggregated data point using bootstrapping.

    This corresponds to the `BootstrapCI` function in the Mathematica code.

    Args:
        samples (np.ndarray): A 2D array of samples. Each row corresponds to a
                              different parameter value or condition.
        n_bootstraps (int): The number of bootstrap samples to generate.
        sample_size (int): The size of each bootstrap sample.
        ci (int): The desired confidence interval (e.g., 95 for 95%).
        outlier_trim (int): The number of samples to trim from each end after
                            sorting to handle extreme outliers.

    Returns:
        np.ndarray: A 3D array of shape (n_params, n_bootstraps, 2) containing
                    the lower and upper bounds of the confidence interval for each
                    parameter value.
    """
    samples = np.asarray(samples)
    n_params = samples.shape[0]

    # Sort and trim outliers
    if outlier_trim > 0:
        trimmed_samples = np.sort(samples, axis=1)[:, outlier_trim:-outlier_trim]
    else:
        trimmed_samples = samples

    # Perform bootstrapping
    bootstrap_means = np.zeros((n_params, n_bootstraps))
    for k in range(n_params):
        for i in range(n_bootstraps):
            bootstrap_sample = np.random.choice(
                trimmed_samples[k],
                size=sample_size,
                replace=True
            )
            bootstrap_means[k, i] = np.mean(bootstrap_sample)

    # Calculate confidence intervals from bootstrap means
    lower_percentile = (100 - ci) / 2
    upper_percentile = 100 - lower_percentile

    ci_lower = np.percentile(bootstrap_means, lower_percentile, axis=1)
    ci_upper = np.percentile(bootstrap_means, upper_percentile, axis=1)

    # The original code returns the sorted bootstrap means within the CI.
    # Returning the lower and upper bounds is more standard in Python.
    # Let's return the full data for plotting, similar to the original.
    # bdata = Transpose[Sort/@bdata];
    # bdata = bdata[[Ceiling[0.025*bn[[1]]];;Floor[0.975*bn[[1]]],;;]]

    sorted_bootstrap_means = np.sort(bootstrap_means, axis=1)
    lower_index = int(np.ceil(n_bootstraps * (lower_percentile / 100)))
    upper_index = int(np.floor(n_bootstraps * (upper_percentile / 100)))

    # Return the bootstrap distributions, trimmed to the CI
    return sorted_bootstrap_means[:, lower_index:upper_index]


def plot_confidence_bands(bdata, index, output_path=None, primary_color='blue', band_color='green', trend_color='red'):
    """
    Plots the bootstrap confidence bands, similar to `BootstrapCIBands`.

    Args:
        bdata (np.ndarray): The bootstrap data, from bootstrap_ci.
        index (np.ndarray): The x-axis values.
        output_path (str, optional): Path to save the plot image. If None,
                                     the plot is shown interactively.
        primary_color (str): Color for the main trend line.
        band_color (str): Color for the scatter plot of bootstrap means.
        trend_color (str): Color for the noise-free threshold line.
    """
    bdata = np.asarray(bdata)
    index = np.asarray(index)

    # Calculate the trend line (mean of bootstrap distributions)
    trend = np.mean(bdata, axis=1)

    plt.figure(figsize=(10, 6))

    # Plot the individual bootstrap means as a scatter plot
    for k in range(bdata.shape[1]):
        plt.scatter(index, bdata[:, k], color=band_color, alpha=0.1, s=5)

    # Plot the main trend line
    plt.plot(index, trend, color=primary_color, linewidth=2, label='Mean Trend')

    # Plot the noise-free threshold line
    if len(trend) > 0:
        plt.axhline(y=trend[0], color=trend_color, linestyle='--', label='Noise-free Threshold')

    plt.xlabel("Parameter Index")
    plt.ylabel("Value")
    plt.title("Bootstrap Confidence Bands")
    plt.legend()
    plt.grid(True)

    if output_path:
        plt.savefig(output_path)
        plt.close() # Close the plot to free memory
    else:
        plt.show()
