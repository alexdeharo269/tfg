import subprocess
import os

#Quiero que en el fichero del analisis salga N, Dinámica, Temperatura, Tiempo hasta estacionario en un fichero aparte resumen.txt


# config.py
range=40
# Path to the first input file
input_file =f"./ising_data/kawa_T0/raw/ising_dataR{range}.dat"

# Folder for intermediate and final output files
output_folder = "./ising_data/kawa_T0"

# Define output files
reduced = f"{output_folder}/reduced/reduced{range}.dat"
domains = f"{output_folder}/domains/domains{range}.dat"
centroids= f"{output_folder}/centroids/centroids{range}.dat"
number_changes= f"{output_folder}/changes/changes{range}.dat"


cpp_source = "compute.cpp"
cpp_executable = "compute"

if not os.path.exists(cpp_executable):  # Compile only if necessary
    print(f"Compiling {cpp_source}...")
    subprocess.run(["g++", cpp_source, "-o", cpp_executable], check=True)




# Run script2 (C++)
subprocess.run([f"./{cpp_executable}", input_file, reduced], check=True)

# Run script1 (Python)
subprocess.run(["python", "test_process.py", reduced, domains, centroids], check=True)

# Run script3 (Python)
subprocess.run(["python", "displaced_number_changes.py",input_file,centroids, number_changes], check=True)

print("Analysis complete.")
