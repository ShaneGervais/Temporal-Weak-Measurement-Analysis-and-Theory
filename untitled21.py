import zipfile
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter


# Define extraction paths
data_folder = "imag_0302_DRAL"
data_files = os.listdir(data_folder)

# Function to load data and compute FFT with zero padding and windowing
def compute_avg_power_spectrum(file_paths, zero_padding_factor=4):
    if not file_paths:
        return None, None  # Handle empty file lists
    
    spectra = []
    freqs = None

    for file in file_paths:
        data = pd.read_csv(file, skiprows=2, header=None, names=["Time", "Voltage"], delimiter=",")
        time = data["Time"].values
        voltage = data["Voltage"].values
        print(np.mean(voltage))
        voltage = voltage - np.mean(voltage)
        
        # Apply a Hanning window to reduce spectral leakage
        window = np.hanning(len(voltage))
        voltage_windowed = voltage * window
        #voltage_windowed = voltage_windowed - np.mean(voltage_windowed)
        
        # Zero padding to improve frequency resolution
        n_fft = len(voltage) * zero_padding_factor
        freq = np.fft.fftfreq(n_fft, d=(time[1] - time[0]))
        spectrum = smooth_spectrum(np.abs(np.fft.fft(voltage_windowed, n=n_fft))**2)
        spectrum = (np.log(spectrum))
        spectra.append(spectrum)
        
        if freqs is None:
            freqs = freq[:len(spectrum)//2]  # Ensuring frequency array matches spectrum size

    avg_spectrum = np.mean(spectra, axis=0)[:len(freqs)]
    return freqs, avg_spectrum

def smooth_spectrum(spectrum, window_size=51, poly_order=3):
    return savgol_filter(spectrum, window_size, poly_order) if spectrum is not None else None

# Function to filter files by angle
def filter_files_by_angle(files, angle, background=False):
    if background:
        pattern = f"mesure_faible_background_5cm_{angle}_deg_"
    else:
        pattern = f"mesure_faible_5cm_{angle}_deg_"
    return [os.path.join(data_folder, file) for file in files if pattern in file]

# Define polarization states
angles = {"D": 3, "A": 41, "R": 28, "L": 67}
spectra = {}
bg_spectra = {}
diff_spectra = {}
peak_data = {}

# Compute power spectra for all polarization states and their backgrounds
for pol, angle in angles.items():
    pol_files = filter_files_by_angle(data_files, angle, background=False)
    bg_files = filter_files_by_angle(data_files, angle, background=True)
    
    freqs, avg_spectrum = compute_avg_power_spectrum(pol_files)
    bg_freqs, avg_bg_spectrum = compute_avg_power_spectrum(bg_files)
    
    spectra[pol] = (freqs, avg_spectrum)
    bg_spectra[pol] = (bg_freqs, avg_bg_spectrum)
    
    if avg_spectrum is not None and avg_bg_spectrum is not None:
        diff_spectra[pol] = smooth_spectrum(np.abs(avg_bg_spectrum - avg_spectrum))
    
    # Identify peaks in the power spectrum
    if diff_spectra[pol] is not None and len(diff_spectra[pol]) > 0:
        peaks, _ = find_peaks(diff_spectra[pol], height=0.001)
        peak_data[pol] = freqs[peaks] if len(peaks) > 0 else []

# Compute frequency differences compared to D
freq_diff = {}
for pol in ["A", "R", "L"]:
    if pol in diff_spectra and "D" in diff_spectra:
        freq_diff[pol] = np.abs(smooth_spectrum((diff_spectra[pol] - diff_spectra["D"])))

# Plot power spectra for all polarization states and their backgrounds
plt.figure(figsize=(10, 6))
for pol in angles.keys():
    freqs, avg_spectrum = spectra[pol]
    bg_freqs, avg_bg_spectrum = bg_spectra[pol]
    if avg_spectrum is not None and avg_bg_spectrum is not None:
        plt.plot(freqs, avg_spectrum, label=f"{pol} (With Interference)", alpha=0.7)
        plt.plot(bg_freqs, avg_bg_spectrum, linestyle="--", label=f"{pol} (Background)", alpha=0.7)

plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectrum (a.u.)")
plt.title("Power Spectra for DRAL States and Their Backgrounds")
plt.legend()
plt.grid()
plt.show()

# Plot Difference Spectra for DRAL states (Interference - Background)
plt.figure(figsize=(10, 6))
for pol in angles.keys():
    if pol in diff_spectra:
        plt.plot(freqs, diff_spectra[pol], label=f"{pol} - Background", alpha=0.7)

plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectrum Difference (a.u.)")
plt.title("Difference in Power Spectra for DRAL (With vs. Without Background)")
plt.legend()
plt.grid()
plt.show()

# Plot frequency differences compared to D
plt.figure(figsize=(10, 6))
for pol in ["A", "R", "L"]:
    if pol in freq_diff:
        plt.plot(freqs, freq_diff[pol], label=f"{pol} - D", alpha=0.7)

plt.xlabel("Frequency (Hz)")
plt.ylabel("Absolute Frequency Difference (a.u.)")
plt.title("Frequency Difference Compared to D")
plt.legend()
plt.grid()
plt.show()

# Compute statistical comparison of peak frequencies
linear_peaks = np.concatenate([peak_data["D"], peak_data["A"]]) if "D" in peak_data and "A" in peak_data else []
circular_peaks = np.concatenate([peak_data["R"], peak_data["L"]]) if "R" in peak_data and "L" in peak_data else []

linear_mean = np.mean(linear_peaks) if len(linear_peaks) > 0 else None
linear_std = np.std(linear_peaks) if len(linear_peaks) > 0 else None
circular_mean = np.mean(circular_peaks) if len(circular_peaks) > 0 else None
circular_std = np.std(circular_peaks) if len(circular_peaks) > 0 else None

# Create DataFrame for statistical results
comparison_df = pd.DataFrame({
    "Polarization Type": ["Linear (D, A)", "Circular (R, L)"],
    "Mean Peak Frequency (Hz)": [linear_mean, circular_mean],
    "Std Dev (Hz)": [linear_std, circular_std]
})

# Display results
print(comparison_df)
