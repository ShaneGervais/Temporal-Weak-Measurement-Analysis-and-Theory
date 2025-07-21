
# weak_value_analysis.py
# Final: Error bars for real & imaginary weak values and propagated errors for a,b

import os
import glob
import re
import numpy as np
import pandas as pd
from scipy.fft import fft, fftshift
from scipy.optimize import curve_fit
from scipy.signal import periodogram, welch, get_window
import matplotlib.pyplot as plt
from pathlib import Path
plt.rcParams['text.usetex'] = True
plt.rcParams['toolbar'] = 'toolbar2'
# --- SCRIPT DIRECTORY SETUP ---
script_dir = Path(__file__).resolve().parent
base_dir = script_dir.parent

# --- CONFIGURATION ---
calibration_folder = base_dir / 'calibration_0907'
measurement_folder = script_dir / 'measurement4_2706'
sampling_interval_ps = 4  # ps

time_step = sampling_interval_ps * 1e-12  # seconds

# Debug prints
print(f"Script directory: {script_dir}")
print(f"Calibration folder: {calibration_folder}")
print(f"Measurement folder: {measurement_folder}")

# --- HELPER FUNCTIONS ---

def my_fftfreq(n, d):
    """Generate frequency vector like MATLAB's fftfreq"""
    val = 1.0 / (n * d)
    if n % 2 == 0:
        k = np.arange(-n//2, n//2)
    else:
        k = np.arange(-(n-1)//2, (n-1)//2 + 1)
    return k * val


def extract_arrival_time(file_list, time_step):
    """
    Compute time-of-arrival via high-res gradient of poly4 fit
    """
    traces = []
    for fpath in file_list:
        data = np.loadtxt(fpath, delimiter=',', skiprows=2)
        traces.append(data[:, 1])
    avg = np.mean(traces, axis=0)
    t = np.arange(len(avg)) * time_step * 1e12  # in ps

    # 40% threshold window
    thr = 0.4 * np.max(avg)
    idxs = np.where(avg >= thr)[0]
    if idxs.size == 0:
        raise ValueError(f"No data above 40% threshold: {file_list}")
    t0, t1 = t[idxs[0]], t[idxs[-1]]
    mask = (t >= t0) & (t <= t1)
    t_fit, y_fit = t[mask], avg[mask]

    # 4th-order polynomial fit
    p = np.polyfit(t_fit, y_fit, 4)
    t_hr = np.linspace(t0, t1, 10000)
    y_hr = np.polyval(p, t_hr)

    # gradient and peak location
    grad_hr = np.gradient(y_hr, t_hr)
    idx_max = np.argmax(np.abs(grad_hr))
    return t_hr[idx_max]

def pspectrum(x, fs=1.0, mode='psd', window='hann', nperseg=None, noverlap=None, nfft=None):
    """
    Python equivalent of MATLAB’s pspectrum, returning frequencies and power.

    Parameters
    ----------
    x : 1D array
        Time-series signal.
    fs : float, optional
        Sampling frequency (Hz). Default is 1.
    mode : {'psd','spectrum'}, optional
        'psd'       → power spectral density (units²/Hz), like pspectrum default.
        'spectrum'  → power spectrum (units²), i.e. PSD * df.
    window : str or tuple or array_like, optional
        Window specification passed to scipy.signal.get_window.
    nperseg : int, optional
        Length of each segment for Welch’s method. Default = len(x).
    noverlap : int, optional
        Number of points to overlap between segments. Default = nperseg//2.
    nfft : int, optional
        Number of FFT points. Default = nperseg.

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxx : ndarray
        Power spectral density or spectrum of x.
    """
    # Choose method: Welch for PSD or spectrum
    # If user wants full-periodogram, set nperseg=len(x), noverlap=0, mode='spectrum'
    if nperseg is None:
        nperseg = len(x)
    if noverlap is None:
        noverlap = nperseg // 2
    if nfft is None:
        nfft = nperseg

    win = get_window(window, nperseg, fftbins=True)

    if mode == 'psd':
        # Welch's method for PSD
        f, Pxx = welch(
            x, fs=fs, window=win, nperseg=nperseg,
            noverlap=noverlap, nfft=nfft, detrend=False,
            return_onesided=True, scaling='density'
        )
    elif mode == 'spectrum':
        # Welch but scaled to total power per bin (power spectrum)
        f, Pxx = welch(
            x, fs=fs, window=win, nperseg=nperseg,
            noverlap=noverlap, nfft=nfft, detrend=False,
            return_onesided=True, scaling='spectrum'
        )
    else:
        raise ValueError("mode must be 'psd' or 'spectrum'")

    return f, Pxx

# --- CALIBRATION TIMES ---
if not calibration_folder.exists():
    raise FileNotFoundError(f"Calibration folder not found: {calibration_folder}")

H_files = glob.glob(str(calibration_folder / '*_3_*.csv'))
V_files = glob.glob(str(calibration_folder / '*_48_*.csv'))
print(f"Found H files: {H_files}")
print(f"Found V files: {V_files}")
if not H_files or not V_files:
    raise FileNotFoundError("H or V calibration files not found.")

t_H = extract_arrival_time(H_files, time_step)
t_V = extract_arrival_time(V_files, time_step)
print(f"Calibration: t_H={t_H:.3f} ps, t_V={t_V:.3f} ps")

# --- LOAD MEASUREMENT FILES ---
files = glob.glob(str(measurement_folder / '*.csv'))
angle_map = {}
for fpath in files:
    m = re.search(r'_(\d+)_deg_', os.path.basename(fpath))
    if not m:
        continue
    angle = int(m.group(1))
    angle_map.setdefault(angle, []).append(fpath)
if not angle_map:
    raise FileNotFoundError(f"No measurement files in {measurement_folder}")

# --- PREP ANGLES ---
angles = sorted(angle_map)
scaled_states = np.array([2 * a for a in angles])
ref_idx = angles.index(48)

# --- REAL PART ANALYSIS (with error bars) ---
arrival_means = []
arrival_stds = []
for angle in angles:
    paths = angle_map[angle]
    arrivals = [extract_arrival_time([p], time_step) for p in paths]
    arrival_means.append(np.mean(arrivals))
    arrival_stds.append(np.std(arrivals))
arrival_means = np.array(arrival_means)
arrival_stds = np.array(arrival_stds)
print(t_H, t_V)
print(arrival_means[0], arrival_means[19])
#path 4
t_H = t_H + 167
t_V = t_V + 167
print(t_H-t_V)
#path 5
#NOTE: THE DATA OSCILLATES FOR 20PS
#t_H = t_H - 180   # Adjust H calibration time
#t_V = t_V - 560  # Adjust V calibration time
#path 3
#t_H = t_H + 32
#t_V = t_V + 20
print(t_H, t_V)
print(t_H - t_V)


# Normalize real weak value between V and H
norm_real = (arrival_means - t_V) / (t_H - t_V)
norm_real = norm_real#/np.max(np.abs(norm_real))
norm_real_err = arrival_stds / abs(t_H - t_V)

# --- IMAGINARY PART ANALYSIS (with error bars) ---
plt.figure()
plt.title("Chirp Frequency Analysis")
plt.xlabel("Angle (degrees)")
plt.ylabel("Chirp Frequency (MHz)")
# Compute chirp frequencies for each angle
chirp_means = []
chirp_stds = []
for angle in angles:
    freqs = []
    for fpath in angle_map[angle]:
        data = np.loadtxt(fpath, delimiter=',', skiprows=2)
        trace = data[:,1]
        N = len(trace)
        S = fft(trace)
        """f_psd, Pxx = pspectrum(
            trace,
            fs=1.0/(time_step),      # sampling frequency in Hz
            mode='spectrum',         # or 'psd' if you want density
            window='hann',
            nperseg=len(trace),
            noverlap=0,
            nfft=len(trace)
        )
        # Pxx is already the power spectrum, so no need for abs**2

        #S = pspectrum(trace)
        #P = np.abs(S)**2
        #plt.plot(np.abs(S), label=f"{angle}°")
        freq = np.sum(f_psd * Pxx) / np.sum(Pxx)*1e-6
        #f = fftshift(my_fftfreq(N, time_step))
        P_shift = fftshift(Pxx)
        #pos = f >= 0
        #freq = np.sum(f[pos] * P_shift[pos]) / np.sum(P_shift[pos]) * 1e-6
        freqs.append(freq)
        """
        P = np.abs(S)**2
        plt.plot(np.abs(S), label=f"{angle}°")
        f = fftshift(my_fftfreq(N, time_step))
        P_shift = fftshift(P)
        pos = f >= 0
        freq = np.sum(f[pos] * P_shift[pos]) / np.sum(P_shift[pos]) 
        freqs.append(freq)
    chirp_means.append(np.mean(freqs))
    chirp_stds.append(np.std(freqs))
chirp_means = np.array(chirp_means)
chirp_stds = np.array(chirp_stds)

# Normalize imaginary part
f_ref = chirp_means[ref_idx]
norm_imag = (chirp_means - f_ref) #/ 16.61e3
print(chirp_means)
print(norm_imag)
print(np.max(np.abs(chirp_means)))
norm_imag_err = chirp_stds/ np.abs(f_ref)

#theorie_real = 1*np.cos(np.radians(scaled_states - 11))**2
theorie_real = 1*np.cos(np.radians(scaled_states))**2 #+ 0.05
#theorie_real = 0*np.cos(np.radians(scaled_states))**2 + 0.5
theorie_imag = 0.5*np.sin(np.radians(4*scaled_states - 11)) + 0.5

theorie_a_re = np.sqrt(theorie_real)
theorie_b_re = np.sqrt(1-np.abs(theorie_a_re)**2)
theorie_a_im = np.sqrt(theorie_imag)
theorie_b_im = 1 - theorie_a_im


re_a = (np.sqrt(norm_real)) # scale to ±1
re_b = np.sqrt(1-(re_a))#(1 - norm_real)

#re_a = re_a / np.max(re_a)
#re_b = re_b / np.max(re_b)
# propagate error: sigma_re = sigma_norm/2
#re_err_a = np.abs(re_a)*(1/2)*norm_real_err/np.abs(norm_real)
#re_err_b = np.abs(re_b)*(1/2)*re_err_a/np.abs(re_a)
re_err_a = norm_real_err / (2 * re_a)
re_err_b = norm_real_err / (2 * re_b)

# scale imag to ±1
max_abs = np.max(np.abs(norm_imag))
norm_imag_scaled = norm_imag / max_abs
norm_imag_scaled_err = norm_imag_err / max_abs
im_a = np.sqrt(norm_imag_scaled)
im_b = np.sqrt(1 - np.abs(im_a)**2)
# propagate error: sigma_im = sigma_scaled/2
im_err = norm_imag_scaled_err / 2

# --- PLOTTING ---
# Raw Real Weak Value
plt.figure()
plt.plot(scaled_states, theorie_real, 'k--', label=r'Théorie')
plt.errorbar(scaled_states, norm_real,
             yerr=norm_real_err, xerr=0.3,
             fmt='bs', ecolor='b', capsize=6,
             label=r'Données mesurées')
plt.xlabel(r"\textbf{État d'entrée } $|\psi(\theta)\rangle$ (degré)", fontsize=14)
plt.ylabel(r"\textbf{Valeur faible } $\mathcal{R}(\langle \hat{\pi} \rangle_W)$ (u.a.)", fontsize=14)
plt.grid(False)
plt.legend()
#plt.ylim(0, 1)
plt.savefig('real_weak_value_path_4.png', dpi=300) 

# Raw Imaginary Weak Value
plt.figure()
#plt.plot(scaled_states, theorie_imag, 'k--', label=r'Théorie')
plt.errorbar(scaled_states, norm_imag_scaled,
             yerr=norm_imag_scaled_err, xerr=0.3,
             fmt='rs', ecolor='r', capsize=6,
             label=r'Données mesurées')
plt.xlabel(r"\textbf{État d'entrée } $|\psi(\theta)\rangle$ (degré)", fontsize=14)
plt.ylabel(r"\textbf{Valeur faible } $\mathcal{I}(\langle \hat{\pi} \rangle_W)$ (u.a.)", fontsize=14)
plt.grid(False)
plt.legend()
#plt.ylim(0, 1)
plt.savefig('imag_weak_value_path_4.png', dpi=300) 

# Re(a) & Re(b)
plt.figure()
plt.plot(scaled_states, theorie_a_re, 'k--', label=r'Théorie $|a|$')
plt.plot(scaled_states, theorie_b_re, 'k--', label=r'Théorie $|b|$')
plt.errorbar(scaled_states, re_a,
             yerr=re_err_a, xerr=0.3,
             fmt='o', color='b', ecolor='b', capsize=6,
             label=r'$|a|$')
plt.errorbar(scaled_states, re_b,
             yerr=re_err_b, xerr=0.3,
             fmt='o', color='r', ecolor='r', capsize=6,
             label=r'$|b|$')
plt.xlabel(r"\textbf{État d'entrée } $|\psi(\theta)\rangle$ (degré)", fontsize=14)
plt.ylabel(r"\textbf{Amplitude de probabilité } $|a|, |b|$ (u.a.)", fontsize=14)
plt.grid(False)
plt.legend()
#plt.ylim(0, 1)
plt.savefig('real_probability_path_4.png', dpi=300)

# Im(a) & Im(b)
plt.figure()
#plt.plot(scaled_states, theorie_a_im, 'k--', label=r'Théorie $\mathcal{I}(a)$')
#plt.plot(scaled_states, theorie_b_im, 'k--', label=r'Théorie $\mathcal{I}(b)$')
plt.errorbar(scaled_states, im_a,
             yerr=im_err, xerr=0.3,
             fmt='o', color='b', ecolor='b', capsize=6,
             label=r'$|a|$')
plt.errorbar(scaled_states, im_b,
             yerr=im_err, xerr=0.3,
             fmt='o', color='r', ecolor='r', capsize=6,
             label=r'$|b|$')
plt.xlabel(r"\textbf{État d'entrée } $|\psi(\theta)\rangle$ (degré)", fontsize=14)
plt.ylabel(r"\textbf{Amplitude de probabilité } $|a|, |b|$ (u.a.)", fontsize=14)
plt.grid(False)
plt.legend()
#plt.ylim(0,1)
plt.savefig('imaginary_probability_path_4.png', dpi=300)

# Table output
results = pd.DataFrame({
    'Angle_deg': angles,
    'NormReal': norm_real,
    'NormRealErr': norm_real_err,
    'NormImag': norm_imag_scaled,
    'NormImagErr': norm_imag_scaled_err,
    'Re(a)': re_a,
    'Re(b)': re_b,
    'ReErr(a)': re_err_a,
    'ReErr(b)': re_err_b,
    'Im(a)': im_a,
    'Im(b)': im_b,
    'ImErr': im_err
})
print(results)
# Save results to CSV
results.to_csv('weak_values_results_path_4_real.csv', index=False)

plt.show()
