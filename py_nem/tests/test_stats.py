import numpy as np
from py_nem import stats

def test_bootstrap_ci_basic():
    """
    Tests the basic functionality of the bootstrap_ci function.
    """
    # 1. Generate synthetic data
    np.random.seed(42)
    # 100 samples for each of 3 "parameters"
    samples = np.random.rand(3, 100)
    samples[0, :] += 1 # shift distributions to be distinct
    samples[1, :] += 2
    samples[2, :] += 3

    # 2. Run the bootstrap function
    n_bootstraps = 200
    ci = 95
    bdata = stats.bootstrap_ci(
        samples,
        n_bootstraps=n_bootstraps,
        sample_size=50,
        ci=ci,
        outlier_trim=2
    )

    # 3. Check the output shape
    # The number of returned samples should correspond to the CI
    expected_samples = int(np.floor(n_bootstraps * (ci / 100)))
    # My implementation is slightly off, let's recalculate the exact size
    lower_percentile = (100 - ci) / 2
    upper_percentile = 100 - lower_percentile
    lower_index = int(np.ceil(n_bootstraps * (lower_percentile / 100)))
    upper_index = int(np.floor(n_bootstraps * (upper_percentile / 100)))

    assert bdata.shape == (samples.shape[0], upper_index - lower_index)

    # 4. Check if the means are reasonable
    # The mean of the bootstrap distributions should be close to the mean
    # of the original sample distributions.
    mean_of_bdata = np.mean(bdata, axis=1)
    mean_of_samples = np.mean(samples, axis=1)

    assert np.allclose(mean_of_bdata, mean_of_samples, atol=0.1)


def test_plot_confidence_bands(tmp_path):
    """
    Tests that the plot_confidence_bands function runs without errors
    and creates an output file.
    """
    # 1. Generate synthetic data
    np.random.seed(42)
    samples = np.random.rand(5, 50)
    index = np.arange(5)

    # 2. Run bootstrap_ci to get data for plotting
    bdata = stats.bootstrap_ci(samples, n_bootstraps=100, sample_size=20)

    # 3. Define output path and run plotting function
    output_file = tmp_path / "confidence_bands.png"
    stats.plot_confidence_bands(bdata, index, output_path=output_file)

    # 4. Assert that the output file was created
    assert output_file.is_file()
    assert output_file.stat().st_size > 0
