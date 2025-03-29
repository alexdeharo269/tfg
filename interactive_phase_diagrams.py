import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import csv
import os
import random
import tempfile
import subprocess
from macro_analysis import SimulationRunner
import time

class IsingVisualizer:
    def __init__(self, csv_file, sim_folder, incoherence_folder):
        """Initialize the visualizer with data sources."""
        self.csv_file = csv_file
        self.sim_folder = sim_folder
        self.incoherence_folder = incoherence_folder
        self.runner = SimulationRunner(csv_file, sim_folder, incoherence_folder)
        self.df_processed = None
        self.unique_x = None
        self.unique_y = None
        self.x_edges = None
        self.y_edges = None
        self.z_values = None
        self.bin_sim_ids = None
        self.max_z = None
        
        # Available parameters for plotting
        self.input_params = ['temperature', 'Rg', 'pp', 'n', 'avg_mag', 'R']
        self.output_params = ['num_peaks', 'meanDomainLength', 'frozentime', 'init_mag']
        # Additional derived parameters that can be calculated from peak data
        self.derived_params = ['mean_peak_height', 'max_peak_height', 'mean_peak_width', 'max_peak_width']
        
    def process_csv(self):
        """Process the CSV file and create the dataframe."""
        data_list = []
        with open(self.csv_file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                processed = self._process_row(row)
                data_list.append(processed)

        # Build a DataFrame from the processed rows
        self.df_processed = pd.DataFrame(data_list)
        
        # Add derived parameters
        self._add_derived_parameters()
        
        return self.df_processed
    
    def _process_row(self, row):
        """
        Process a row from the CSV.
        Assumes:
          Fixed columns (indices 0-11):
             0: sim_id, 1: pp, 2: temperature, 3: avg_mag, 4: n, 5: R, 6: Rg,
             7: init_mag, 8: meanDomainLength, 9: frozentime, 10: seed, 11: num_peaks
        """
        fixed = row[:12]
        try:
            num_peaks = int(fixed[11])
            seed = int(fixed[10])
        except Exception:
            num_peaks = 0
            seed = 0
            
        peak_heights = []
        peak_widths = []
        peak_positions = []
        
        if num_peaks > 0:
            for i in range(num_peaks):
                try:
                    peak_heights.append(float(row[12 + i]))
                except Exception:
                    peak_heights.append(np.nan)
            for i in range(num_peaks):
                try:
                    peak_widths.append(float(row[12 + num_peaks + i]))
                except Exception:
                    peak_widths.append(np.nan)
            for i in range(num_peaks):
                try:
                    peak_positions.append(int(row[12 + 2*num_peaks + i]))
                except:
                    peak_positions.append(np.nan)
                    
        return {
            'sim_id': fixed[0],
            'pp': float(fixed[1]),
            'temperature': float(fixed[2]),
            'avg_mag': float(fixed[3]),
            'n': int(fixed[4]),
            'R': int(fixed[5]),
            'Rg': int(fixed[6]),
            'init_mag': float(fixed[7]),
            'meanDomainLength': float(fixed[8]),
            'frozentime': int(fixed[9]),
            'seed': seed,
            'num_peaks': num_peaks,
            'peak_heights': peak_heights,
            'peak_widths': peak_widths,
            'peak_positions': peak_positions
        }
    
    def _add_derived_parameters(self):
        """Add derived parameters to the dataframe based on peak data."""
        # Calculate mean peak height for each simulation
        self.df_processed['mean_peak_height'] = self.df_processed.apply(
            lambda row: np.mean(row['peak_heights']) if row['peak_heights'] and row['num_peaks'] > 0 else np.nan, axis=1)
        
        # Calculate maximum peak height
        self.df_processed['max_peak_height'] = self.df_processed.apply(
            lambda row: np.max(row['peak_heights']) if row['peak_heights'] and row['num_peaks'] > 0 else np.nan, axis=1)
        
        # Calculate mean peak width
        self.df_processed['mean_peak_width'] = self.df_processed.apply(
            lambda row: np.mean(row['peak_widths']) if row['peak_widths'] and row['num_peaks'] > 0 else np.nan, axis=1)
        
        # Calculate max peak width
        self.df_processed['max_peak_width'] = self.df_processed.apply(
            lambda row: np.max(row['peak_widths']) if row['peak_widths'] and row['num_peaks'] > 0 else np.nan, axis=1)
        
    def filter_data(self, fixed_params=None):
        """
        Filter the data based on fixed parameter values.
        
        Args:
            fixed_params (dict): Dictionary of parameters to fix and their values
                                 Example: {'temperature': 2.0, 'n': 100}
        
        Returns:
            pd.DataFrame: Filtered dataframe
        """
        if self.df_processed is None:
            self.process_csv()
            
        if fixed_params is None or len(fixed_params) == 0:
            return self.df_processed
        
        # Start with full dataset
        filtered_df = self.df_processed
        
        # Apply filters for each fixed parameter
        for param, value in fixed_params.items():
            if param in filtered_df.columns:
                # For numeric parameters, allow for small floating-point differences
                if isinstance(value, (int, float)):
                    filtered_df = filtered_df[np.isclose(filtered_df[param], value)]
                else:
                    filtered_df = filtered_df[filtered_df[param] == value]
        
        return filtered_df
    
    def build_phase_diagram_data(self, x_param, y_param, z_param, fixed_params=None):
        """
        Build data for phase diagram with flexible parameter selection.
        
        Args:
            x_param (str): Parameter to use for x-axis
            y_param (str): Parameter to use for y-axis
            z_param (str): Parameter to use for z-axis (color)
            fixed_params (dict): Parameters to keep fixed
        """
        if self.df_processed is None:
            self.process_csv()
            
        # Filter data based on fixed parameters
        filtered_df = self.filter_data(fixed_params)
        
        if len(filtered_df) == 0:
            raise ValueError("No data points match the fixed parameter criteria")
            
        # Extract data for plotting
        x = filtered_df[x_param].values
        y = filtered_df[y_param].values
        z = filtered_df[z_param].values
        
        # Handle NaN values in z
        mask = ~np.isnan(z)
        x, y, z = x[mask], y[mask], z[mask]
        
        # Get unique values for binning
        self.unique_x = np.sort(np.unique(x))
        self.unique_y = np.sort(np.unique(y))
        
        # Create bin edges
        self.x_edges = np.zeros(len(self.unique_x) + 1)
        self.y_edges = np.zeros(len(self.unique_y) + 1)
        
        # Calculate bin edges between points
        for i in range(len(self.unique_x) - 1):
            self.x_edges[i+1] = (self.unique_x[i] + self.unique_x[i+1]) / 2
        # Add outer edges
        if len(self.unique_x) > 1:
            self.x_edges[0] = self.unique_x[0] - (self.unique_x[1] - self.unique_x[0]) / 2
            self.x_edges[-1] = self.unique_x[-1] + (self.unique_x[-1] - self.unique_x[-2]) / 2
        else:
            # Handle case with only one unique x value
            self.x_edges[0] = self.unique_x[0] - 0.5
            self.x_edges[-1] = self.unique_x[0] + 0.5
        
        # Do the same for y
        for i in range(len(self.unique_y) - 1):
            self.y_edges[i+1] = (self.unique_y[i] + self.unique_y[i+1]) / 2
        if len(self.unique_y) > 1:
            self.y_edges[0] = self.unique_y[0] - (self.unique_y[1] - self.unique_y[0]) / 2
            self.y_edges[-1] = self.unique_y[-1] + (self.unique_y[-1] - self.unique_y[-2]) / 2
        else:
            # Handle case with only one unique y value
            self.y_edges[0] = self.unique_y[0] - 0.5
            self.y_edges[-1] = self.unique_y[0] + 0.5
        
        # Calculate mean of z values in each bin
        self.z_values = np.full((len(self.unique_y), len(self.unique_x)), np.nan)  # NaN will indicate no data
        
        # Store simulation IDs for each bin
        self.bin_sim_ids = [[[] for _ in range(len(self.unique_x))] for _ in range(len(self.unique_y))]
        
        print(f"Processing {z_param} values for phase diagram...")
        
        # Map filtered_df indices to original dataframe for retrieval
        index_map = filtered_df.index.tolist()
        
        # In the build_phase_diagram_data method when calculating mean z values:
        for i in range(len(self.unique_x)):
            for j in range(len(self.unique_y)):
                in_bin = ((x >= self.x_edges[i]) & (x < self.x_edges[i+1]) &
                         (y >= self.y_edges[j]) & (y < self.y_edges[j+1]))
                if np.any(in_bin):
                    z_vals = z[in_bin]
                    # Only consider non-NaN values
                    valid_z = z_vals[~np.isnan(z_vals)]
                    if len(valid_z) > 0:
                        # Calculate mean of valid z values
                        self.z_values[j, i] = np.mean(valid_z)
                        
                        # Store the simulation IDs in this bin
                        bin_indices = np.where(in_bin)[0]
                        for idx in bin_indices:
                            original_idx = index_map[idx]
                            sim_id = filtered_df.iloc[idx]['sim_id']
                            self.bin_sim_ids[j][i].append(sim_id)
        
        # Determine maximum z value (for colormap)
        self.max_z = np.nanmax(self.z_values)
        print(f"Maximum {z_param} value: {self.max_z}")
        
        return {
            'x_param': x_param,
            'y_param': y_param,
            'z_param': z_param,
            'fixed_params': fixed_params
        }
    
    def create_colormap(self, z_param):
        """Create a suitable colormap for the visualization."""
        if self.z_values is None:
            raise ValueError("Phase diagram data must be built first")
            
        # Choose colormap based on parameter type
        if z_param in ['num_peaks', 'frozentime']:
            # For count/discrete data, use a sequential colormap
            cmap = plt.cm.viridis
        elif 'height' in z_param:
            # For height-related params, use a "hot" colormap
            cmap = plt.cm.hot
        elif 'width' in z_param:
            # For width-related params, use a different colormap
            cmap = plt.cm.cool
        else:
            # Default
            cmap = plt.cm.viridis
        
        # Define normalization for the colorbar
        vmin = np.nanmin(self.z_values)
        vmax = self.max_z
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        
        # For colorbar ticks, limit to 5 or 6 ticks
        if z_param == 'num_peaks' and vmax <= 4:
            ticks = np.arange(0, vmax + 1)
        else:
            ticks = np.linspace(vmin, vmax, 6)
            
        return cmap, norm, ticks
    
    def plot_phase_diagram(self, x_param, y_param, z_param, fixed_params=None):
        """
        Plot a static phase diagram with flexible parameter selection.
        
        Args:
            x_param (str): Parameter to use for x-axis
            y_param (str): Parameter to use for y-axis
            z_param (str): Parameter to use for z-axis (color)
            fixed_params (dict): Parameters to keep fixed (e.g., {'temperature': 2.0})
        """
        # Build data for the phase diagram
        plot_info = self.build_phase_diagram_data(x_param, y_param, z_param, fixed_params)
        
        # Create colormap
        cmap, norm, ticks = self.create_colormap(z_param)
        
        plt.figure(figsize=(10, 8))
        # Create a masked array to hide NaN values
        masked_z = np.ma.masked_invalid(self.z_values)
        pc = plt.pcolormesh(self.x_edges, self.y_edges, masked_z, cmap=cmap, norm=norm, shading='auto')
        
        # Add parameter information to labels
        plt.xlabel(f'{x_param}', fontsize=12)
        plt.ylabel(f'{y_param}', fontsize=12)
        
        # Create title that includes fixed parameters
        if fixed_params:
            fixed_str = ', '.join([f'{k}={v}' for k, v in fixed_params.items()])
            title = f'Phase Diagram: {z_param} (Fixed: {fixed_str})'
        else:
            title = f'Phase Diagram: {z_param}'
        
        plt.title(title, fontsize=14)
        cbar = plt.colorbar(pc, ticks=ticks)
        cbar.set_label(f'{z_param}', fontsize=12)
        plt.tight_layout()
        plt.show()
    
    def create_interactive_plot(self, x_param, y_param, z_param, fixed_params=None):
        """
        Create an interactive phase diagram with flexible parameter selection.
        
        Args:
            x_param (str): Parameter to use for x-axis
            y_param (str): Parameter to use for y-axis
            z_param (str): Parameter to use for z-axis (color)
            fixed_params (dict): Parameters to keep fixed (e.g., {'temperature': 2.0})
        """
        # Build data for the phase diagram
        plot_info = self.build_phase_diagram_data(x_param, y_param, z_param, fixed_params)
        
        # Create colormap
        cmap, norm, ticks = self.create_colormap(z_param)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        # Create a masked array to hide NaN values
        masked_z = np.ma.masked_invalid(self.z_values)
        pc = ax.pcolormesh(self.x_edges, self.y_edges, masked_z, cmap=cmap, norm=norm, shading='auto')
        
        # Add parameter information to labels
        ax.set_xlabel(f'{x_param}', fontsize=12)
        ax.set_ylabel(f'{y_param}', fontsize=12)
        
        # Create title that includes fixed parameters
        if fixed_params:
            fixed_str = ', '.join([f'{k}={v}' for k, v in fixed_params.items()])
            title = f'Phase Diagram: {z_param} (Fixed: {fixed_str})'
        else:
            title = f'Phase Diagram: {z_param}'
        
        ax.set_title(f'{title}\nClick to view random simulation', fontsize=14)
        cbar = fig.colorbar(pc, ax=ax, ticks=ticks)
        cbar.set_label(f'{z_param}', fontsize=12)
        
        # Dictionary to track which simulations have been shown for each bin
        shown_simulations = {}
        
        # Define callback function for mouse clicks
        def on_click(event):
            if event.inaxes != ax:
                return
                
            # Get the clicked coordinates
            x_click, y_click = event.xdata, event.ydata
                
            # Find which bin was clicked
            x_bin = np.digitize(x_click, self.x_edges) - 1
            y_bin = np.digitize(y_click, self.y_edges) - 1
                
            # Check if valid bin
            if 0 <= x_bin < len(self.unique_x) and 0 <= y_bin < len(self.unique_y):
                bin_key = f"{x_bin}_{y_bin}"
                sim_ids = self.bin_sim_ids[y_bin][x_bin]
                
                if sim_ids:
                    # Initialize tracking for this bin if it's the first click
                    if bin_key not in shown_simulations:
                        shown_simulations[bin_key] = {
                            'all_ids': sim_ids.copy(),
                            'remaining_ids': sim_ids.copy(),
                            'last_id': None
                        }
                    
                    bin_data = shown_simulations[bin_key]
                    
                    # If we've shown all simulations, reset the list
                    if not bin_data['remaining_ids']:
                        print("All simulations in this bin have been viewed. Resetting.")
                        bin_data['remaining_ids'] = bin_data['all_ids'].copy()
                        # Remove the last shown ID to avoid showing it again immediately
                        if bin_data['last_id'] in bin_data['remaining_ids'] and len(bin_data['remaining_ids']) > 1:
                            bin_data['remaining_ids'].remove(bin_data['last_id'])
                    
                    # Select a random simulation from remaining options
                    sim_id = random.choice(bin_data['remaining_ids'])
                    bin_data['remaining_ids'].remove(sim_id)
                    bin_data['last_id'] = sim_id
                    
                    print(f"Clicked bin at {x_param}={self.unique_x[x_bin]:.4f}, {y_param}={self.unique_y[y_bin]:.4f}")
                    print(f"{z_param} in this bin: {self.z_values[y_bin, x_bin]:.4f}")
                    print(f"Showing simulation {len(bin_data['all_ids']) - len(bin_data['remaining_ids'])} of {len(bin_data['all_ids'])}")
                    print(f"Displaying peaks for simulation: {sim_id}")
                        
                    # Get the row corresponding to this simulation ID
                    row = self.df_processed[self.df_processed['sim_id'] == sim_id].iloc[0]
                    num_peaks = int(row['num_peaks'])
                    peak_heights = row['peak_heights']
                    peak_widths = row['peak_widths']
                    peak_positions = row['peak_positions']
                    seed = row['seed']
                    
                    print(f"Using seed: {seed}, Number of peaks: {num_peaks}")
                    
                    # Load the simulation data
                    input_file = os.path.join(self.sim_folder, f"{sim_id}.dat")
                    incoherence_file = os.path.join(self.incoherence_folder, f"incoherence{sim_id}.dat")    
                    # Load incoherence data
                    incoherence_data = np.loadtxt(incoherence_file, usecols=1)
                    
                    # Call the runner's plot_warning_details with the right parameters
                    self.runner.plot_warning_details(sim_id, input_file, incoherence_data, peak_positions, peak_heights, peak_widths)
                    
                else:
                    print(f"No simulations in the selected bin at {x_param}={x_click:.4f}, {y_param}={y_click:.4f}")
        
        # Connect the click event
        fig.canvas.mpl_connect('button_press_event', on_click)
        
        plt.tight_layout()
        plt.show()
        return fig, ax
    
    def get_available_parameters(self):
        """Return lists of available parameters for plotting."""
        return {
            'input_params': self.input_params,
            'output_params': self.output_params,
            'derived_params': self.derived_params
        }


def main():
    # Define paths
    sim_folder = "ising_data/simulations"
    csv_file = "./ising_data/simulation_index.csv" 
    incoherence_folder = "ising_data/incoherence"
    
    # Create visualizer
    visualizer = IsingVisualizer(csv_file, sim_folder, incoherence_folder)
    
    # Process data
    visualizer.process_csv()
    
    # USER-CONFIGURABLE PARAMETERS:
    # Choose parameters for x, y, and z axes
    x_param = 'avg_mag'  # Parameter for x-axis
    y_param = 'R'        # Parameter for y-axis
    z_param = 'num_peaks'  # Parameter for z-axis (color)
    
    # Optional: Specify fixed parameters (all other input parameters must be constant)
    # For example, to filter for simulations with temperature=2.0, Rg=5, etc.
    fixed_params = {
        'temperature': 0.00,  # Example value
        'Rg': 1,             # Example value
        'pp': 1.00,           # Example value
        'n': 256             # Example value
    }
    
    # Alternative z-parameters that can be used:
    # - 'num_peaks': Number of peaks (default)
    # - 'meanDomainLength': Mean domain length
    # - 'frozentime': Time to freeze
    # - 'mean_peak_height': Mean height of peaks
    # - 'max_peak_height': Maximum peak height
    # - 'mean_peak_width': Mean width of peaks
    # - 'max_peak_width': Maximum peak width
    
    start_time = time.time()
    
    # Print available parameters
    params = visualizer.get_available_parameters()
    print("Available input parameters:", params['input_params'])
    print("Available output parameters:", params['output_params'])
    print("Available derived parameters:", params['derived_params'])
    
    # Create and display the interactive plot
    visualizer.create_interactive_plot(x_param, y_param, z_param, fixed_params)
    
    total_time = time.time() - start_time
    print(f"Total execution time: {total_time:.2f} seconds")

if __name__ == "__main__":
    # Run the main function
    main()