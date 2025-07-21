#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 01:35:43 2025

@author: shanegervais
"""

import zipfile
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


# Function to load data and compute FFT with zero padding and windowing
def compute_avg_power_spectrum(file_paths, zero_padding_factor=4):
    spectra = []
    freqs = None

    for file in file_paths:
        data = pd.read_csv(file, skiprows=2, header=None, names=["Time", "Voltage"], delimiter=",")
        time = data["Time"].values
        voltage = data["Voltage"].values
        
        # Apply a Hanning window to reduce spectral leakage
        window = np.hanning(len(voltage))
        voltage_windowed = voltage * window
        
        # Zero padding to improve frequency resolution
        n_fft = len(voltage) * zero_padding_factor
        freq = np.fft.fftfreq(n_fft, d=(time[1] - time[0]))
        spectrum = np.abs(np.fft.fft(voltage_windowed, n=n_fft))**2

        spectra.append(spectrum)
        
        if freqs is None:
            freqs = freq

    avg_spectrum = np.mean(spectra, axis=0)
    return freqs, avg_spectrum, time, voltage

# Function to find peak differences
def find_peak_differences(freqs, spectrum_diff):
    peaks, _ = find_peaks(np.abs(spectrum_diff), height=0.001)
    return freqs[peaks], spectrum_diff[peaks]

# Function to pad arrays to the same length
def pad_with_nan(arr, length):
    return np.pad(arr, (0, length - len(arr)), constant_values=np.nan)

def filter_files_by_angle(files, angles):
    selected_files = []
    for file in files:
        for angle in angles:
            if f"_{angle}_" in file:
                selected_files.append(os.path.join(data_folder, file))
    return selected_files

# Define the paths for extracted files
data_folder = "imag_0127_copy/"
data_files = os.listdir(data_folder)

# Define polarization states
D_files = filter_files_by_angle(data_files, [3, 93])
A_files = filter_files_by_angle(data_files, [48])
R_files = filter_files_by_angle(data_files, [25])
L_files = filter_files_by_angle(data_files, [70])

# Compute power spectra with improved FFT
freqs_D, avg_spectrum_D, time_D, voltage_D = compute_avg_power_spectrum(D_files)
freqs_A, avg_spectrum_A, _, _ = compute_avg_power_spectrum(A_files)
freqs_R, avg_spectrum_R, _, _ = compute_avg_power_spectrum(R_files)
freqs_L, avg_spectrum_L, _, _ = compute_avg_power_spectrum(L_files)

# Plot Raw Signal Before FFT (for D polarization as an example)
plt.figure(figsize=(10, 5))
plt.plot(time_D, voltage_D, label="Raw Signal (D Polarization)", alpha=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Voltage (a.u.)")
plt.title("Raw Signal Before FFT")
plt.legend()
plt.grid()
plt.show()

# Plot Power Spectra for DRAL
plt.figure(figsize=(10, 6))
plt.plot(freqs_D[:len(freqs_D)//2], avg_spectrum_D[:len(freqs_D)//2], label="D (Diagonal)", alpha=0.7)
plt.plot(freqs_A[:len(freqs_A)//2], avg_spectrum_A[:len(freqs_A)//2], label="A (Anti-Diagonal)", alpha=0.7)
plt.plot(freqs_R[:len(freqs_R)//2], avg_spectrum_R[:len(freqs_R)//2], label="R (Right Circular, 25°)", alpha=0.7)
plt.plot(freqs_L[:len(freqs_L)//2], avg_spectrum_L[:len(freqs_L)//2], label="L (Left Circular, 70°)", alpha=0.7)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power Spectrum (a.u.)")
plt.title("Comparison of Power Spectra for DRAL Polarization States")
plt.legend()
plt.grid()
plt.show()

# Compute differences
spectrum_diff_D93_D3 = avg_spectrum_D[:len(freqs_D)//2] - avg_spectrum_D[:len(freqs_D)//2]
spectrum_diff_A_D3 = avg_spectrum_A[:len(freqs_A)//2] - avg_spectrum_D[:len(freqs_D)//2]
spectrum_diff_R_D3 = avg_spectrum_R[:len(freqs_R)//2] - avg_spectrum_D[:len(freqs_D)//2]
spectrum_diff_L_D3 = avg_spectrum_L[:len(freqs_L)//2] - avg_spectrum_D[:len(freqs_D)//2]

# Plot Power Spectrum Differences
for diff, label, color in zip([spectrum_diff_D93_D3, spectrum_diff_A_D3, spectrum_diff_R_D3, spectrum_diff_L_D3],
                              ["D(93°) - D(3°)", "A - D(3°)", "R - D(3°)", "L - D(3°)"],
                              ["purple", "blue", "green", "red"]):
    plt.figure(figsize=(10, 5))
    plt.plot(freqs_D[:len(freqs_D)//2], diff, label=label, color=color, alpha=0.7)
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power Spectrum Difference (a.u.)")
    plt.title(f"Difference in Power Spectrum ({label})")
    plt.legend()
    plt.grid()
    plt.show()

# Find Peaks in Difference Spectra
peak_freqs_D93_D3, peak_values_D93_D3 = find_peak_differences(freqs_D[:len(freqs_D)//2], spectrum_diff_D93_D3)
peak_freqs_A_D3, peak_values_A_D3 = find_peak_differences(freqs_A[:len(freqs_A)//2], spectrum_diff_A_D3)
peak_freqs_R_D3, peak_values_R_D3 = find_peak_differences(freqs_R[:len(freqs_R)//2], spectrum_diff_R_D3)
peak_freqs_L_D3, peak_values_L_D3 = find_peak_differences(freqs_L[:len(freqs_L)//2], spectrum_diff_L_D3)

# Adjust lists to have the same length
max_length = max(len(peak_freqs_D93_D3), len(peak_freqs_A_D3), len(peak_freqs_R_D3), len(peak_freqs_L_D3))
peak_freqs_D93_D3 = pad_with_nan(peak_freqs_D93_D3, max_length)
peak_values_D93_D3 = pad_with_nan(peak_values_D93_D3, max_length)
peak_freqs_A_D3 = pad_with_nan(peak_freqs_A_D3, max_length)
peak_values_A_D3 = pad_with_nan(peak_values_A_D3, max_length)
peak_freqs_R_D3 = pad_with_nan(peak_freqs_R_D3, max_length)
peak_values_R_D3 = pad_with_nan(peak_values_R_D3, max_length)
peak_freqs_L_D3 = pad_with_nan(peak_freqs_L_D3, max_length)
peak_values_L_D3 = pad_with_nan(peak_values_L_D3, max_length)