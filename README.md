# PSD Calculation in Rust

This program essentially copies the `pwelch` function from MATLAB. It will:

1. Remove mean from signal.
2. Apply a [Hanning window](https://en.wikipedia.org/wiki/Window_function#Hann_and_Hamming_windows).
3. Compute a [Cooley-Tukey FFT](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm) on an n-length sample, where `n` is the power of 2 that results in the closest possible window to 1 second.
4. Compute the PSD up to the [Nyquist frequency](https://en.wikipedia.org/wiki/Nyquist_frequency) of the signal.
5. Convert to [1/6-octave band](https://en.wikipedia.org/wiki/Octave_band) (1/3 and 1/6 OB are by convention in vibration analysis).
6. Compute a final envelope for the entire signal.

My sample time series were free and came from [Mide](https://blog.mide.com/vibration-analysis-fft-psd-and-spectrogram), whose blog also provides a great into to vibration analysis.

My PSD and scale factor calculation came from equation 23 of [this paper](https://www.researchgate.net/publication/267956210_Spectrum_and_spectral_density_estimation_by_the_Discrete_Fourier_transform_DFT_including_a_comprehensive_list_of_window_functions_and_some_new_flat-top_windows).

## Instructions

Modify the file path in your `main.rs` to point to a CSV file with time in the first column and signal in the second. This program makes the following assumptions:

- Time is in seconds and **all points are evenly spaced**.
- Signal is stochastic and [stationary](https://en.wikipedia.org/wiki/Stationary_process).
- Data has been cleaned (viz, all time points have a corresponding signal point that can be read as an `f64`).
- Signal is in unit acceleration (g), although this is mostly related to plots.

```
cargo build --release
cargo run --release
```

The `--release` build will run significantly faster than the debug build.

## To Do

- Make this a Parity Substrate off-chain worker.
- Change this to a library that accepts a `Vec<_>` and returns the envelope.
- Add traits to allow more input types.
- Add more error handling (e.g. if signal contains `NaN` values).
- Add more features (e.g. input different units, plot options, window options).
- Implement FFT algorithms that accept non-power-of-two signal lengths.
- Performance comparisons against Python, MATLAB, etc.

## Additional inspiration

- [Tom Irvine's website](http://vibrationdata.com/random.htm) and his [C++ PSD Calculation](http://www.vibrationdata.com/software_alt/poweri_lite.cpp).
- [Scientific Computing in Rust](https://www.lpalmieri.com/posts/2019-02-23-scientific-computing-a-rust-adventure-part-0-vectors/)

## Friendly Warning

This has not been peer-reviewed at all. I have compared FFT results to other libraries and compared time-based and PSD-based RMS results, and everything looks reasonable. However, this was much more about learning more Rust and eventually integrating into Substrate as a _demo_.
