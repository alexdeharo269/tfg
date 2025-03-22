import os
import csv
import subprocess
import tempfile
import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import multiprocessing
import traceback
import functools
import warnings
#Usar initialize ring balance para reac dif. 
#para kawa T0 puede que vaya bien hacer solo number_changes y luego si no funciona bien probamos con la incoherence.
#Hay parametros que tocar en la función find_peacks. Otro requisito puede ser imponer un minimo de altura 
#La altura media de los picos está en grouped_heights, pero puede ser mejor ir a por la total. 
#También podemos ir a por si cambia el número de dominios dentro de una simulación, aunque habria que quitar tiempos iniciales. 

class SimulationRunner:
    def __init__(self, csv_path, sim_folder):
        """
        Inicializa el runner con la ruta del CSV y el folder de simulaciones.
        """
        self.csv_path = csv_path
        self.sim_folder = sim_folder
        self._compile_executables()

    def _compile_executables(self):
        # Compilar test_process
        if not os.path.exists("./test_process"):
            subprocess.run(["g++", "test_process_cpp.cpp", "-o", "test_process"], check=True)
        # Compilar changes_incoherence
        if not os.path.exists("./changes_incoherence"):
            subprocess.run(["g++", "chnages_incoherence.cpp", "-o", "changes_incoherence"], check=True)
    
    @staticmethod
    def run_simulation_parallel(runner,row):
        """Wrapper function to run a single simulation in parallel, handling errors."""
        sim_id = row[0]  # Assuming the first column is the simulation key
        try:
            sim_output = runner.run_simulation(row)
            row.extend(map(str, sim_output))
            return row
        except Exception as e:
            error_message = f"\nError in simulation {sim_id}:\n{traceback.format_exc()}"
            print(error_message, flush=True)  # Ensure error prints immediately

            return row  # Return original row (unchanged) to maintain CSV structure
    
    def chimeras(self, sim_params, temp_incoherence):
        # Load only the second column (number of changes)
        data = np.loadtxt(temp_incoherence, usecols=1)
        length = len(data)  # Total number of points

        # Detect peaks with constraints
        height_peaks=max(int(round(0.01*float(sim_params[8]))),10)
        peaks, properties = find_peaks(data, distance=length*0.1, prominence=max(data)*0.2,height=height_peaks)
        peak_heights = data[peaks]

        # Group nearby peaks into single maxima
        grouped_peaks = []
        grouped_heights = []
        if len(peaks) > 0:
            grouped_peaks.append(peaks[0])
            grouped_heights.append(peak_heights[0])
            for i in range(1, len(peaks)):
                if peaks[i] - grouped_peaks[-1] > 20:  # Avoid clustering too close
                    grouped_peaks.append(peaks[i])
                    grouped_heights.append(peak_heights[i])

        # --- Handling Periodic Boundary Conditions ---
        if len(grouped_peaks) > 1:
            first_peak = grouped_peaks[0]
            last_peak = grouped_peaks[-1]
            # Check if first and last peaks are close under PBC
            if (first_peak + (length - last_peak)) < 20:  
                # Merge the two: take an approximate middle position and average height
                avg_position = (first_peak + last_peak) // 2
                avg_height = (grouped_heights[0] + grouped_heights[-1]) / 2
                grouped_peaks = grouped_peaks[1:-1]  # Remove first and last
                grouped_peaks.append(avg_position)
                grouped_heights = grouped_heights[1:-1]  # Remove first and last
                grouped_heights.append(avg_height)

        # Widths at half-height for the grouped peaks
        results_half = peak_widths(data, grouped_peaks, rel_height=0.5)
        widths = results_half[0]
        min_width=10
        valid_peak=widths>=min_width
        # Filter grouped_peaks, grouped_heights, and widths
        grouped_peaks = np.array(grouped_peaks)[valid_peak].tolist()
        grouped_heights = np.array(grouped_heights)[valid_peak].tolist()
        widths = widths[valid_peak]
        
        num_big_maxima = len(grouped_peaks)
        analysis_results = [num_big_maxima] + grouped_heights + widths.tolist()
        return analysis_results, data, grouped_peaks, grouped_heights, widths

    def plot_peaks(self, data, peaks, heights, widths, sim_id=""):
        plt.figure(figsize=(10, 6))
        plt.plot(data, label="Activity (changes)")
        plt.plot(peaks, data[peaks], "ro", label="Detected Peaks")
        
        # Annotate each peak with its height and width
        for i, peak in enumerate(peaks):
            plt.annotate(f"H: {heights[i]:.2f}\nW: {widths[i]:.1f}",
                         xy=(peak, data[peak]), xytext=(5, 5),
                         textcoords="offset points", fontsize=8, color="darkblue")
        
        plt.xlabel("Position")
        plt.ylabel("Activity (number of changes)")
        plt.title(f"Simulation {sim_id} - Detected Peaks")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_warning_details(self, sim_id, input_file, data, peaks, heights, widths):
        # Attempt to load the simulation matrix from the input file.
        try:
           sim_matrix = np.loadtxt(input_file)
        except Exception as e:
           sim_matrix = None
           print(f"Could not load simulation matrix from {input_file}: {e}")

        # Create a figure with two subplots.
        fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True ,figsize=(10, 12))
    
        # Top subplot: Simulation matrix
        if sim_matrix is not None:
            im = ax1.imshow(sim_matrix, aspect='auto', interpolation='nearest', cmap=None,origin='lower')
            ax1.set_title(f"Simulation {sim_id} Matrix")
            fig.colorbar(im, ax=ax1)
        else:
            ax1.text(0.5, 0.5, "Simulation matrix not available",
                 ha='center', va='center', transform=ax1.transAxes)
            ax1.set_title(f"Simulation {sim_id} Matrix")
    
        ax1.set_xlim(0, len(data))
        ax2.set_xlim(0, len(data))
        # Bottom subplot: Peaks analysis
        ax2.plot(data, label="Activity (changes)")
        ax2.plot(peaks, data[peaks], "ro", label="Detected Peaks")
        for i, peak in enumerate(peaks):
            ax2.annotate(f"H: {heights[i]:.2f}\nW: {widths[i]:.1f}",
                     xy=(peak, data[peak]), xytext=(5, 5),
                     textcoords="offset points", fontsize=8, color="darkblue")
        ax2.set_xlabel("Position")
        ax2.set_ylabel("Activity (number of changes)")
        ax2.set_title("Peaks Analysis")
        ax2.legend()
    
        plt.tight_layout()
        plt.show()  # Blocks execution until the window is closed


    def run_simulation(self, sim_params):
        # La clave de la simulación es el primer elemento de la fila CSV
        sim_id = sim_params[0]
        # Construir el nombre del fichero de entrada basado en la clave
        input_file = os.path.join("./", self.sim_folder, f"{sim_id}.dat")
        
        # Ejecutar el primer comando (test_process) usando un fichero temporal para centroids
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            temp_centroids = tmp.name
        with tempfile.NamedTemporaryFile(delete=False) as tmp_2:
            temp_incoherence = tmp_2.name
        
        subprocess.run(
            ["./test_process", input_file, "-", "-", temp_centroids, "false"],
            check=True
        )
        # Ejecutar el segundo comando (changes_incoherence)
        subprocess.run(
            ["./changes_incoherence", input_file, temp_centroids, temp_incoherence, "-", "false"],
            check=True
        )
        os.unlink(temp_centroids)

        # Run the chimeras analysis, which now returns extra debugging information.
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim_output, data, big_peaks, big_heights, widths = self.chimeras(sim_params, temp_incoherence)
            if w:
                for warn in w:
                    print(f"Warning in simulation {sim_id}: {warn.message}")
                    # Freeze execution and display a double plot:
                    
                self.plot_warning_details(sim_id, input_file, data, big_peaks, big_heights, widths)


        os.unlink(temp_incoherence)
        # Return a tuple with both the simulation output and debug info
        return sim_output

    def run_all_simulations(self, debug_mode=False, num_workers=None):

        warnings.filterwarnings("ignore")

        updated_rows = []
        start_time = time.time()  # Start timing the entire process

        # Read all rows first
        with open(self.csv_path, "r", newline="") as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)

        total_simulations = len(rows)

        
        # Set up multiprocessing pool
        num_workers = num_workers or multiprocessing.cpu_count()
        with multiprocessing.Pool(processes=num_workers) as pool:
                # Use tqdm to track progress while handling errors
            run_sim_func = functools.partial(SimulationRunner.run_simulation_parallel, self)
            results = list(tqdm(pool.imap(
                            run_sim_func, rows), 
                            total=total_simulations, 
                            desc="Running Simulations\n",
                            ncols=80))

            updated_rows.extend(results)

        # Write results back to CSV after all simulations are done
        with open(self.csv_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(updated_rows)

        # Print final statistics
        total_time = time.time() - start_time
        print(f"\nTotal simulations: {total_simulations}. Total time: {total_time:.2f} s.")
        


# ----------------------------
# Ejemplo de uso:
# ----------------------------
if __name__ == "__main__":
    sim_folder = "ising_data/simulations"
    csv_file = "./ising_data/simulation_index.csv"
    debug_mode = False  # Cambia a True para ver la gráfica con datos y picos detectados
    
    runner = SimulationRunner(csv_file, sim_folder)
    runner.run_all_simulations(debug_mode=debug_mode)
