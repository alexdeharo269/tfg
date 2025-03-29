import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import mode
import csv
import os
from macro_analysis import SimulationRunner


interactive=False
sim_folder = "ising_data/simulations"
csv_file = "./ising_data/simulation_index.csv"
runner = SimulationRunner(csv_file, sim_folder)# --- Step 1: Process the CSV using the csv module ---
def process_row(row):
    """
    Process a row from the CSV.
    Assumes:
      Fixed columns (indices 0-9):
         0: sim_id, 1: temperature, 2: avg_mag, 3: n, 4: R, 5: Rg,
         6: init_mag, 7: meanDomainLength, 8: frozentime, 9: num_peaks
      Then, if num_peaks > 0, extra columns:
         For i in 0...num_peaks-1: peak_heights at columns 10+i,
         and for i in 0...num_peaks-1: peak_widths at columns 10+num_peaks+i.
    """
    fixed = row[:10]
    try:
        num_peaks = int(fixed[9])
    except Exception:
        num_peaks = 0
    peak_heights = []
    peak_widths = []
    if num_peaks > 0:
        for i in range(num_peaks):
            try:
                peak_heights.append(float(row[10 + i]))
            except Exception:
                peak_heights.append(np.nan)
        for i in range(num_peaks):
            try:
                peak_widths.append(float(row[10 + num_peaks + i]))
            except Exception:
                peak_widths.append(np.nan)
    return {
        'sim_id': fixed[0],
        'temperature': float(fixed[1]),
        'avg_mag': float(fixed[2]),
        'n': int(fixed[3]),
        'R': int(fixed[4]),
        'Rg': int(fixed[5]),
        'init_mag': float(fixed[6]),
        'meanDomainLength': float(fixed[7]),
        'frozentime': int(fixed[8]),
        'num_peaks': num_peaks,
        'peak_heights': peak_heights,
        'peak_widths': peak_widths
    }


data_list = []
with open('./ising_data/simulation_index.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        processed = process_row(row)
        data_list.append(processed)

# Build a DataFrame from the processed rows.
df_processed = pd.DataFrame(data_list)
print(df_processed.head())

# --- Step 2: Build the phase diagram data ---
# Use:
#   x-axis: avg_mag (column 'avg_mag')
#   y-axis: R (column 'R')
#   Color: number of peaks (column 'num_peaks')
x = df_processed['avg_mag'].values
y = df_processed['R'].values
z = df_processed['num_peaks'].values.astype(int)

# Modified binning approach - create bins centered at each unique point
unique_x = np.sort(np.unique(x))
unique_y = np.sort(np.unique(y))

# Create bin edges centered at each point
x_edges = np.zeros(len(unique_x) + 1)
y_edges = np.zeros(len(unique_y) + 1)

# Calculate bin edges between points
for i in range(len(unique_x) - 1):
    x_edges[i+1] = (unique_x[i] + unique_x[i+1]) / 2
# Add outer edges
x_edges[0] = unique_x[0] - (unique_x[1] - unique_x[0]) / 2
x_edges[-1] = unique_x[-1] + (unique_x[-1] - unique_x[-2]) / 2

# Do the same for y
for i in range(len(unique_y) - 1):
    y_edges[i+1] = (unique_y[i] + unique_y[i+1]) / 2
y_edges[0] = unique_y[0] - (unique_y[1] - unique_y[0]) / 2
y_edges[-1] = unique_y[-1] + (unique_y[-1] - unique_y[-2]) / 2

# Better approach to get the most common number of peaks in each bin
z_mode = np.full((len(unique_y), len(unique_x)), -1)  # -1 will indicate no data


# Store simulation IDs for each bin to enable lookups on click
bin_sim_ids = [[[] for _ in range(len(unique_x))] for _ in range(len(unique_y))]

for i in range(len(unique_x)):
    for j in range(len(unique_y)):
        in_bin = ((x >= x_edges[i]) & (x < x_edges[i+1]) &
                  (y >= y_edges[j]) & (y < y_edges[j+1]))
        if np.any(in_bin):
            z_vals = z[in_bin]
            # Count occurrences
            unique_vals, counts = np.unique(z_vals, return_counts=True)
            if len(unique_vals) > 0:
                # Get the most frequent value
                z_mode[j, i] = unique_vals[np.argmax(counts)]
                
                
                # Store the simulation IDs in this bin
                indices = np.where(in_bin)[0]
                for idx in indices:
                    bin_sim_ids[j][i].append(df_processed.iloc[idx]['sim_id'])
    
# Determine maximum number of peaks in the data
max_peaks = int(np.max(z))
print(f"Maximum number of peaks found: {max_peaks}")

# Create a colormap that scales from 0 to max_peaks
# Use a color sequence from blue to red
colors = plt.cm.viridis(np.linspace(0, 1, max_peaks + 2))  # +2 for -1 (no data) and 0
colors_list = ['lightgray'] + [colors[i] for i in range(1, max_peaks + 2)]

# Create boundaries for the colormap
boundaries = np.arange(-1.5, max_peaks + 1.5, 1)
norm = mcolors.BoundaryNorm(boundaries, ncolors=len(colors_list))
discrete_cmap = mcolors.ListedColormap(colors_list)

# For colorbar ticks, limit to at most 5 ticks
if max_peaks <= 4:
    ticks = np.arange(0, max_peaks + 1)
else:
    ticks = np.linspace(0, max_peaks, 5, dtype=int)

# --- Step 5: Plot the phase diagram ---
plt.figure(figsize=(8, 6))
pc = plt.pcolormesh(x_edges, y_edges, z_mode, cmap=discrete_cmap, norm=norm, shading='auto')
plt.xlabel('Average Magnetization', fontsize=12)
plt.ylabel('R', fontsize=12)
plt.title('Phase Diagram: Number of Peaks', fontsize=14)
cbar = plt.colorbar(pc, ticks=[0, 1, 2, 3])
cbar.set_label('Number of Peaks', fontsize=12)
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(8, 6))


# Callback function for mouse clicks
def on_click(event):
    if event.inaxes != ax:
        return
        
    # Get the clicked coordinates
    x_click, y_click = event.xdata, event.ydata
        
    # Find which bin was clicked
    x_bin = np.digitize(x_click, x_edges) - 1
    y_bin = np.digitize(y_click, y_edges) - 1
        
    # Check if valid bin
    if 0 <= x_bin < len(unique_x) and 0 <= y_bin < len(unique_y):
        sim_ids = bin_sim_ids[y_bin][x_bin]
        if sim_ids:
            # If multiple simulations in the bin, pick the first one
            sim_id = sim_ids[0]
            print(f"Clicked bin contains simulation(s): {sim_ids}")
            print(f"Displaying peaks for simulation: {sim_id}")
                
            # Get the row corresponding to this simulation ID
            row = df_processed[df_processed['sim_id'] == sim_id].iloc[0]
            num_peaks = int(row['num_peaks'])
            peak_heights = row['peak_heights']
            peak_widths = row['peak_widths']
                
            # Use the actual SimulationRunner to plot
            sim_params = [sim_id]  # Build the parameters list as needed
            input_file = os.path.join(sim_folder, f"{sim_id}.dat")
                    
                # Run the necessary analysis and plot
                # Note: This part would need to be customized based on your specific needs
            temp_data = np.random.rand(100)  # Placeholder
            peaks = [int(i*10) for i in range(num_peaks)]  # Placeholder
            runner.plot_warning_details(temp_data, peaks, peak_heights, peak_widths, sim_id)
                
        else:
            print(f"No simulations in the selected bin at x={x_click:.2f}, y={y_click:.2f}")
    
    # Connect the click event
    fig.canvas.mpl_connect('button_press_event', on_click)
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax
