import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sgn
import scipy
import pandas as pd

# Load data from file
with open("visibility_data.txt") as f:
    data = f.read()

# Split data into lines
data = data.split('\n')

# Extract columns from data
x = np.double([row.split()[0] for row in data])
y = np.double([row.split()[1] for row in data])
z = np.double([row.split()[2] for row in data])

# Constants
offset = 11.8e-3

# Convert x values from cm to inches (2.54 cm = 1 inch)
x_cm = (x-0) * 2.54 *2
# Voltages
vmin = y - offset
vmax = z - offset

# Calculate visibility
vis = (vmax - vmin) / (vmin + vmax)
# Create figure and subplots with shared x-axis

# Analysis
vis_part = vis[:330]
peaks, _ = sgn.find_peaks(vis_part,distance=5)
peaks = list(peaks)  # Convert peaks to a list (if it's not already)

# Append the extra element to peaks
peaks.append(0)
peaks.append(209)
#peaks.append(244)
peaks.append(248)
peaks.append(288)
peaks.append(364)
peaks.append(302)

main_peaks, _ = sgn.find_peaks(vis_part,distance=30,prominence=0.1)
main_peaks = list(main_peaks)
main_peaks.append(0)
main_peaks.append(323)
main_peaks.append(364)

# Fit
def monoExp(x, m, t, b):
    return m * np.exp(-t * x**2) +b

main_peak_x=x_cm[main_peaks]
main_peak_y=vis[main_peaks]
params, cv = scipy.optimize.curve_fit(monoExp, main_peak_x, main_peak_y)
m, t, b = params
print(params)
sample_x=np.linspace(0,6,200)
fit_val=monoExp(sample_x,m,t,b)

# Locations peaks
differences = []
all_peak_pos = sorted(x_cm[peaks])
for i in range(0, len(all_peak_pos)-2):
    diff = all_peak_pos[i+1] - all_peak_pos[i]
    #print(diff)
    differences.append(diff)
differences2=[]
for i in range(0, len(all_peak_pos)-4):
    diff2 = differences[i+1] - differences[i]
    print(diff2)
    differences2.append(diff2)
print(differences2)
"""
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9))  # 1 row, 2 columns

ax1.scatter(all_peak_pos[0:-2],differences)
ax1.hlines(np.mean(differences),0,6,label='Mean distance = {}cm'.format(np.round(np.mean(differences),4)))
ax1.set_xlabel('Path difference/cm')
ax1.set_ylabel('Distance between peaks/cm')
ax1.set_title('Locations of Peaks')
print(np.mean(differences))
ax1.legend()
#plt.show()

# Plot 2: Visibility
ax2.plot(x_cm, vis,color='b')
ax2.scatter(x_cm, vis,s=10,color='b')

ax2.scatter(x_cm[peaks],vis[peaks],s=20,color='orange')

ax2.scatter(main_peak_x,main_peak_y,color='red', marker='*', s=25)
ax2.plot(sample_x, fit_val, color='red', label=r'Fit: ${:.4f} \cdot e^{{-{:.4f}/ cm \cdot x}} $'.format(m, t))
ax2.set_ylabel('Visibility')
ax2.set_xlabel('Path difference/cm')
ax2.set_title('Visibility for Different Lens Positions')
ax2.legend()
plt.tight_layout()
plt.show()
"""


plt.scatter(x_cm, vis,s=10,color='b')
plt.scatter(x_cm[peaks],vis[peaks],s=20,color='orange')
plt.scatter(main_peak_x,main_peak_y,color='red', marker='*', s=25)
plt.xlabel(r"Distance (cm)")
plt.ylabel(r"Visibilité (u.a.)")
plt.grid(True)
plt.show()

results = pd.DataFrame(
    {
        'x_cm': x_cm,
        'vmin': vmin,
        'vmax': vmax,
        'visibility': vis
    }
)

results.to_csv('visibility_results.csv', index=False)

print("Results saved to visibility_results.csv")