use std::fs;
use std::error::Error;
use std::f64::consts::PI;
use num_complex::Complex;
use gnuplot::{Figure, Color, AxesCommon};

fn main() {

	// Read in the data from a CSV.
	let data_file: &'static str = "./aircraft.csv";
	let (time, signal) = input_from_file(data_file).unwrap();

	// Sample Frequency (Hz)
	let f_s: f64 = 1. / (time[1] - time[0]);

	// Total points
	let ns: usize = signal.len();

	// Number of data points in a slice. Equals Fs * seconds/slice. Must be power of 2.
	let n_segment: usize = 2_usize.pow(f_s.log2() as u32);

	// Bandwidth
	let df: f64 = f_s / n_segment as f64;

	// Frequency resolution for FFT
	let f_fft: Vec<f64> = (0..n_segment).map(|ii| ii as f64 * df).collect();

	// Amount slices should overlap.
	let n_overlap: f64 = 0.5;

	// Set up indexing variables
	let idx_increment: usize = ((n_segment as f64) * n_overlap) as usize;
	let total_windows: usize = time.len() / n_segment;
	let mut start_idxs: Vec<usize> = (0..total_windows).map(|ii| ii * idx_increment).collect();
	start_idxs.push(ns - n_segment); // Add last segment

	// The window will reduce the amplitude of the FFT. We will use this scale factor on the results.
	let mut hann_scale: Vec<f64> = (0..n_segment).map(|_| 1.0_f64).collect();
	hanning_window(&mut hann_scale);
	let scale_factor: f64 = hann_scale.iter().map(|ii| ii * ii).sum();

	// Calculate output frequency vector
	let f_psd: Vec<f64> = (0..f_fft.len()/2).map(|ii| ii as f64 * f_s / f_fft.len() as f64).collect();
	let center_frequencies: Vec<f64> = octave_bands(&f_psd, 6);
	let mut psd_max: Vec<f64> = (0..center_frequencies.len()).map(|_| 0.0 as f64).collect();

	println!("\nSignal of {:.*} seconds.", 2, time[ns - 1] - time[0]);
	println!("Sample rate: {:.*} Hz", 0, f_s);
	println!("Frequency resolution: {:.*} Hz", 3, df);
	println!("Analyzing {} segments of {:.*} seconds each.\n", start_idxs.len(), 3, n_segment as f64 / f_s);

	// Process Data
	let mut ii: usize = 0;
	for _ in 0..start_idxs.len() {

		let idx_s = start_idxs[ii];
		let idx_e: usize = idx_s + n_segment;

		let mut analysis_segment: Vec<f64> = signal[idx_s..idx_e].to_vec();

		// Perform mean removal
		mean_removal(&mut analysis_segment);

		// Apply Hanning Window
		hanning_window(&mut analysis_segment);

		// Compute FFT
		let input: Vec<Complex<f64>> = (0..n_segment).map(|ii| Complex::new(analysis_segment[ii], 0.0_f64)).collect();
		let spectrum = fft(input);

		// Compute PSD
		let psd = psd(&spectrum, &f_s, &scale_factor);

		// Reband PSD
		let one_sixth_psd: Vec<f64> = reband(&psd, &f_psd, &df, &center_frequencies);

		// Max PSD
		assert!(psd_max.len() == one_sixth_psd.len());
		for jj in 0..psd_max.len() {
			if one_sixth_psd[jj] > psd_max[jj] {
				psd_max[jj] = one_sixth_psd[jj];
			}
		}

		// Sanity check, comment out if you don't want plots
		if ii == start_idxs.len()/2 {
			plot_time(&time[idx_s..idx_e].to_vec(), &analysis_segment.clone());
			let yfft: Vec<f64> = (0..spectrum.len()).map(|ii| spectrum[ii].norm_sqr().sqrt()).collect();
			plot_freq(&f_fft, &yfft, "FFT");
			plot_freq(&f_psd, &psd, "PSD (g^2/Hz)");
			plot_freq(&center_frequencies, &one_sixth_psd, "PSD in 1/6 Octave Band (g^2/Hz)");
			let rms_time: f64 = rms_from_time(&signal[idx_s..idx_e].to_vec());
			let rms_psd: f64 = rms_from_psd(&psd, &df);
			println!("RMS from time: {:.*} g", 3, rms_time);
			println!("RMS from PSD:  {:.*} g", 3, rms_psd);
		}

		ii += 1;
	}
	println!("Processed {} FFTs",ii);

	plot_freq(&center_frequencies, &psd_max, "Max PSD in 1/6 Octave Band (g^2/Hz)");
}

// Read in data from a CSV file.
fn input_from_file(fname: &'static str) -> Result<(Vec<f64>, Vec<f64>), Box<Error>> {

	let file = fs::File::open(fname)?;

	let mut rdr = csv::Reader::from_reader(file);
	// Stupid workaround for allocating vector length. Consumes rdr so need to read again.
	let n = rdr.records().count();

	let mut time: Vec<f64> = (0..n).map(|_| 0 as f64).collect();
	let mut signal: Vec<f64> = (0..n).map(|_| 0 as f64).collect();

	let file = fs::File::open(fname)?;
	let mut rdr = csv::Reader::from_reader(file);

	let mut ii: usize = 0;
	for result in rdr.deserialize() {
		let record: (f64, f64) = result?;
		time[ii] = record.0;
		signal[ii] = record.1;
		ii += 1
	}

	Ok((time, signal))
}

// Remove the mean from the signal.
fn mean_removal(v: &mut Vec<f64>) {
	let l = v.len();
	let s: f64 = v.iter().sum();
	let mean: f64 = s / (l as f64);

	for ii in 0..l {
		v[ii] -= mean;
	}
}

// Apply a hanning window to a signal.
fn hanning_window(v: &mut Vec<f64>) {
	let l = v.len();

	for ii in 0..l {
		v[ii] *= (PI * (ii as f64) / (l as f64)).sin();
	}
}

// Compute an FFT using the Cooley-Tukey algorithm. The length of the input
// must be a power of 2.
fn fft(mut v: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
	let m = v.len() / 2;

	if v.len() == 1 {
		return v;
	}

	let mut evens = fft(v.iter().step_by(2).cloned().collect());
	let mut odds = fft(v.iter().skip(1).step_by(2).cloned().collect());
	evens.append(&mut odds);
	v = evens;

	for ii in 0..m {
		let n = v.len();
		let theta: f64 = -2.0_f64 * PI * (ii as f64) / (n as f64);
		let twid: Complex<f64> = Complex::from_polar(&1.0_f64, &theta);
		let t = v[ii];
		v[ii] = t + twid * v[ii + m];
		v[ii + m] = t - twid * v[ii + m];
	}
	return v;
}

// Compute the Power Spectral Density based on an FFT as input. Will scale
// according to the sample frequency and window bias.
fn psd(ft: &Vec<Complex<f64>>, f_s: &f64, wsf: &f64) -> Vec<f64> {
	let n = ft.len();
	// Complex conjugate squared, up to Nyquist freq (n/2), scaled by hanning window
	// factor and sample frequency. Many references will note to scale by signal length
	// and sample frequency, but the hanning factor *includes* the signal length and the
	// amplitude adjustment necessary to account for the window.
	let mut psd: Vec<f64> = (0..n/2).map(|ii| ft[ii].norm_sqr() / (wsf * f_s)).collect();

	// Do not include end points because of symmetry.
	for ii in 1..psd.len()-1 {
		// Normalize by 4/2 to account for Nyquist cutoff and peak to rms.
		psd[ii] *= 2.0_f64;
	}

	return psd;
}

// Compute new octave bands.
//
// f0: Input frequency vector.
// ob: Octave band (3 = 1/3 octave band).
fn octave_bands(f0: &Vec<f64>, ob: usize) -> Vec<f64> {

	let f_min: f64 = f0[0].max(15.625); // This will make "1000 Hz" one of the points, which is by convention.
	let f_max: f64 = f0[f0.len() - 1];
	let n1: f64 = ob as f64 * (f_max / f_min).log2();
	let n2: usize = n1 as usize;

	let center_freq: Vec<f64> = (0..n2).map(|ii| f_min * 2.0_f64.powf(ii as f64 / ob as f64)).collect();

	assert!(center_freq[0] >= f0[0]);
	assert!(center_freq[center_freq.len() - 1] <= f0[f0.len() - 1]);

	return center_freq;
}

// Reband the PSD to a new frequency vector. The new PSD should retain the same
// overall G_rms, given by the area under the curve.
//
// psd: Input spectrum.
// f0:  Frequency vector for `psd`. Must be the same length as `psd`.
// df:  Bandwidth of `f0`. Must be constant.
// cf:  Center frequencies of output.
fn reband(psd: &Vec<f64>, f0: &Vec<f64>, df: &f64, cf: &Vec<f64>) -> Vec<f64> {
	// println!("f: {}, psd: {}", f0.len(), psd.len());
	assert!(psd.len() == f0.len());
	let n: usize = cf.len();

	let mut output: Vec<f64> = (0..n).map(|_| 0.0_f64).collect();

	for ii in 0..n {
		let mut lb: f64; // lower bound
		let mut ub: f64; // upper bound

		// Find the lower and upper frequency bounds for this center frequency.
		if ii == 0 {
			lb = cf[0];
			ub = 10.0_f64.powf((cf[ii].log10() + cf[ii+1].log10()) / 2.0_f64);
		} else if ii == n - 1 {
			lb = 10.0_f64.powf((cf[ii-1].log10() + cf[ii].log10()) / 2.0_f64);
			ub = cf[n - 1];
		} else {
			lb = 10.0_f64.powf((cf[ii-1].log10() + cf[ii].log10()) / 2.0_f64);
			ub = 10.0_f64.powf((cf[ii].log10() + cf[ii+1].log10()) / 2.0_f64);
		}

		// Collect the points from `f0` that will correspond to this band
		let idx_lb: usize = f0.iter().position(|&a| a > lb).unwrap(); // index of lb in f0
		let idx_ub: usize = f0.iter().position(|&a| a > ub).unwrap(); // index of ub in f0

		let psd_segment: Vec<f64> = psd[idx_lb..idx_ub].to_vec();
		let rms = rms_from_psd(&psd_segment, df);

		output[ii] = rms / (ub - lb);
	}

	return output;
}

// Compute the overall RMS of a spectrum.
fn rms_from_psd(p: &Vec<f64>, df: &f64) -> f64 {
	let psd_df: Vec<f64> = (0..p.len()).map(|ii| p[ii] * df).collect();
	let psd_sum: f64 = psd_df.iter().sum();
	let rms: f64 = psd_sum.sqrt();
	return rms;
}

// Compute the overall RMS of a time series.
fn rms_from_time(sig: &Vec<f64>) -> f64 {
	let n = sig.len();
	let mut rms: f64 = 0.0_f64;
	for ii in 0..n {
		rms += sig[ii] * sig[ii] / n as f64;
	}
	rms = rms.sqrt();
	return rms;
}

// Plot a signal with time on the x-axis.
fn plot_time(x: &Vec<f64>, y: &Vec<f64>) {
	let mut fg = Figure::new();
	fg.axes2d().lines(x, y, &[Color("red")])
		.set_x_label("Time (s)", &[])
		.set_y_label("Signal (g)", &[]);
	fg.show();
}

// Plot a spectrum in the frequency domain. Log-log axes.
fn plot_freq(f: &Vec<f64>, y: &Vec<f64>, label: &'static str) {
	let mut fg = Figure::new();
	fg.axes2d().lines(f, y, &[Color("red")])
		.set_x_label("Frequency (Hz)", &[])
		.set_y_label(label, &[])
		.set_x_log(Some(10.0_f64))
		.set_y_log(Some(10.0_f64));
	fg.show();
}
