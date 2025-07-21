import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
plt.rcParams['text.usetex'] = True
plt.rcParams['toolbar'] = 'toolbar2'
# --- CONFIGURATION ---
# Path to folder containing your *_int_signal.csv files
data_dir = Path('./imag_plusieur')  
sampling_interval_ps = 4      # in picoseconds
dt = sampling_interval_ps * 1e-12  # in seconds per sample

# --- COLLECT SIGNAL FILES ---
signal_files = sorted(data_dir.glob('*_int_fft.csv'))
if not signal_files:
    raise FileNotFoundError(f"No *_int_signal.csv files found in {data_dir}")

# --- HELPER FUNCTION ---

def compute_centroid(trace, dt):
    
    #Compute the power-spectrum centroid of a real-valued time trace.
    #Returns centroid in MHz.
    
    N = len(trace)*100
    # FFT and power
    S = np.fft.fft(trace, n=N)
    P = np.abs(S)**2
    # Frequency axis in Hz
    freqs = np.fft.fftshift(np.fft.fftfreq(N, dt))
    plt.plot(freqs, P, label='Power Spectrum', marker='o', markersize=2, linestyle='-', linewidth=0.5)
    plt.xlabel(r'\textbf{Frequency} (Hz)')
    plt.ylabel('Power')
    plt.title('Power Spectrum of Time Trace')

    P_shift = np.fft.fftshift(P)
    # Keep only positive freqs
    mask = freqs >= 0
    freqs_pos = freqs[mask]
    
    P_pos = P_shift[mask]

    peak_idxs, _ = find_peaks(P_pos, height=1e7)

    # get their frequencies and powers
    peak_freqs = freqs_pos[peak_idxs]
    print(f"Peak Frequencies (Hz): {peak_freqs}")
    peak_powers = P_pos[peak_idxs]

    # compute centroid over only the peaks
    centroid_peaks_hz = np.sum(peak_freqs * peak_powers) / np.sum(peak_powers)
    centroid_peaks_mhz = centroid_peaks_hz * 1e-6
    # Weighted average (centroid)
    #centroid_hz = np.sum(freqs_pos * P_pos) / np.sum(P_pos)
    #return centroid_hz * 1e-6  # convert to MHz
    return centroid_peaks_mhz



# --- PROCESS ALL STATES ---
states = []
centroids = []

plt.figure(figsize=(10, 6))
for fpath in signal_files:
    # Derive label (e.g. 'A' from 'A_int_signal.csv')
    state = fpath.stem.split('_')[0].upper()
    states.append(state)

    # Load the time-domain trace
    data = np.loadtxt(fpath, delimiter=',', skiprows=0)
    trace = data[:,1] if data.ndim > 1 else data

    # Compute and store centroid
    centroids.append(compute_centroid(trace, dt))
#plt.show()
# Build DataFrame
df = pd.DataFrame({
    'State': states,
    'Centroid_MHz': centroids
})

# Normalize labels and ensure V present
df['State'] = df['State'].str.upper()
if 'V' not in df['State'].values:
    raise ValueError("Reference state 'V' not found — check your filenames.")

# Compute shift relative to V
v_centroid = df.loc[df['State']=='V', 'Centroid_MHz'].iat[0]
df['Shift_MHz'] = df['Centroid_MHz'] - v_centroid

# --- PLOT BAR CHART ---
plt.figure(figsize=(8,4))
plt.bar(df['State'], df['Shift_MHz'], color='skyblue')
plt.xlabel('Polarization State', fontsize=12)
plt.ylabel('Centroid Shift (MHz)', fontsize=12)
plt.title('Frequency Centroid Shift Relative to V State', fontsize=14)
plt.grid(True, axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# --- OPTIONAL: display the table ---
print(df.to_string(index=False))
