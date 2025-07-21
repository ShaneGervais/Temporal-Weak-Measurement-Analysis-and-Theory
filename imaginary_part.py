import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline

plt.rcParams['text.usetex'] = True
plt.rcParams['toolbar'] = 'toolbar2'

# ——— Parameters ———
polarization_states = ['h','v','r','a','l','d']
data_folder         = Path('./imag_plusieur')  # Folder with *_int_fft.csv files
threshold           = 100_000_000               # Min PSD to consider as a peak
search_half_width   = 8                        # Points around each peak for local fit
tolerance_hz        = 20e6                     # 20 MHz clustering tolerance
dt                  = 4e-12                   # 4 ps sampling interval
fsamp               = 1/dt                     # Sample rate

# Storage for fitted peaks (Hz)
peak_map    = {pol: [] for pol in polarization_states}
v_ref_freq  = None   # will hold the first 'V' peak we detect

# ——— Loop over each polarization file ———
plt.figure(figsize=(10,6))
colors = plt.cm.tab10.colors

for idx, pol in enumerate(polarization_states):
    fname = data_folder / f"{pol}_int_fft.csv"
    if not fname.exists():
        print(f"  ⚠️  Missing {fname}")
        continue

    # 1) Load the time‑trace and compute PSD via zero-padded FFT
    raw_y = np.loadtxt(fname).flatten()
    N     = len(raw_y) * 100
    S     = np.fft.fft(raw_y, n=N)
    y     = np.abs(S)**2

    # 2) Build frequency axis
    f_center = 703.12e6
    f_span   = 1.4062e9
    freq_hz  = np.linspace(f_center - f_span/2,
                           f_center + f_span/2,
                           y.size)

    # 3) Plot raw PSD
    plt.scatter(freq_hz/1e6, y, s=8, color=colors[idx], label=pol.upper())

    # 4) Find coarse peaks
    locs, props = find_peaks(y, height=threshold)
    for loc in locs:
        # Local window around the peak
        i1    = max(0, loc - search_half_width)
        i2    = min(len(y), loc + search_half_width + 1)
        x_fit = freq_hz[i1:i2]
        y_fit = y[i1:i2]
        if x_fit.size < 5:
            continue

        # 5) Lightly smooth & fit a smoothing spline
        y_sm  = pd.Series(y_fit).rolling(5, center=True, min_periods=1).mean().values
        w     = np.sqrt(y_sm) + 1e-6
        spline = UnivariateSpline(x_fit, y_sm, w=w, s=0.9 * x_fit.size)
        xg     = np.linspace(x_fit.min(), x_fit.max(), 2000)
        yg     = spline(xg)

        # 6) Offset spline to pass through raw peak
        x_data_peak = freq_hz[loc]
        y_data_peak = y[loc]
        yg += (y_data_peak - spline(x_data_peak))

        # 7) Extract "true" peak from fitted curve
        imax      = np.argmax(yg)
        true_freq = xg[imax]
        true_amp  = yg[imax]

        # store peak
        peak_map[pol].append(true_freq)

        # if this is the first V‑peak, record it
        if pol == 'v' and v_ref_freq is None:
            v_ref_freq = true_freq

        # compute and print delay if we have V reference
        if v_ref_freq is not None:
            delay_khz = (true_freq - v_ref_freq) 
            print(
                f"Pol {pol.upper():>1} → "
                f"{true_freq/1e6:7.3f} MHz  "
                f"(amp {true_amp:.1f}, delay {delay_khz:.3f} kHz)"
            )
        else:
            print(
                f"Pol {pol.upper():>1} → "
                f"{true_freq/1e6:7.3f} MHz  "
                f"(amp {true_amp:.1f})"
            )

        # plot fitted curve & marker
        #plt.plot(xg, yg, '-', color=colors[idx], lw=1.5)
        plt.axvline(true_freq/1e6, color=colors[idx], ls='--')

# Finalize plot
plt.xlabel(r'\textbf{Fréquence} (MHz)', fontsize=14)
plt.ylabel(r'\textbf{Amplitude} (u.a.)', fontsize=14)
plt.legend(title='Polarisations')
plt.grid(False)
plt.xlim(0, 710)
plt.ylim(-0.1e9, 1.75e9)
plt.tight_layout()
plt.savefig('PSD.png', dpi=300)
plt.show()


# ——— Cluster all peaks into “zones” ———
all_peaks = []
for pol_idx, pol in enumerate(polarization_states):
    for f in peak_map[pol]:
        all_peaks.append((f, pol_idx))
all_peaks.sort(key=lambda x: x[0])

zones = []
current_zone = [all_peaks[0]]
for f, pi in all_peaks[1:]:
    if abs(f - current_zone[-1][0]) <= tolerance_hz:
        current_zone.append((f, pi))
    else:
        zones.append(current_zone)
        current_zone = [(f, pi)]
zones.append(current_zone)


# ——— Build results table ———
rows = []
for z_idx, zone in enumerate(zones, start=1):
    # reference freq for V in this zone
    ref_freq = next((f for f, pi in zone if polarization_states[pi]=='v'), np.nan)
    for pol_idx, pol in enumerate(polarization_states):
        freq_pos = next((f for f, pi in zone if pi==pol_idx), np.nan)
        delay_hz = (freq_pos - ref_freq) if not np.isnan(ref_freq) and not np.isnan(freq_pos) else np.nan
        rows.append({
            'Zone':          z_idx,
            'Polarisation':  pol.upper(),
            'Frequence_MHz': freq_pos / 1e6,
            'Delai_kHz':     delay_hz / 1e3
        })

df = pd.DataFrame(rows, columns=['Zone','Polarisation','Frequence_MHz','Delai_kHz'])
print(df)
df.to_csv('delais_par_rapport_V.csv', index=False)
