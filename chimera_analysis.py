import subprocess
import os

#Quiero que en el fichero del analisis salga N, Dinámica, Temperatura, Tiempo hasta estacionario en un fichero aparte resumen.txt


# config.py
range=51 
# Path to the first input file
input_file =f"./ising_data/kawa_T0/Avg_mag=0.75/raw/ising_dataR{range}T0.0.dat"

# Folder for intermediate and final output files
output_folder = "./ising_data/kawa_T0/Avg_mag=0.75"

# Define output files
reduced = f"{output_folder}/reduced/reduced{range}.dat"
domains = f"{output_folder}/domains/domains{range}.dat"
centroids= f"{output_folder}/centroids/centroids{range}.dat"
number_changes= f"{output_folder}/changes/changes{range}.dat"
incoherence= f"{output_folder}/incoherence/incoherence{range}.dat"
incoherence_matrix= f"{output_folder}/incoherence_matrix/incoherence_matrix{range}.dat"

cpp_file_1 = "compute.cpp"
cpp_executable_1 = "compute"
cpp_file_2 = "test_process_cpp.cpp"
cpp_executable_2 = "test_process"
cpp_file_3 = "displaced_number_changes_cpp.cpp"
cpp_executable_3 = "displaced_number_changes"
cpp_file_4 = "clusters_centrois_incoherence.cpp"
cpp_executable_4 = "clusters_centroids_incoherence"
cpp_file_5="chnages_incoherence.cpp"
cpp_executable_5="changes_incoherence"


if not os.path.exists(cpp_executable_1):  # Compile only if necessary 
    print(f"Compiling {cpp_file_1}...")
    subprocess.run(["g++", cpp_file_1, "-o", cpp_executable_1], check=True)
if not os.path.exists(cpp_executable_2):  # Compile only if necessary
    print(f"Compiling {cpp_file_2}...")
    subprocess.run(["g++", cpp_file_2, "-o", cpp_executable_2], check=True)
if  not os.path.exists(cpp_executable_5):  # Compile only if necessary
    print(f"Compiling {cpp_file_5}...")
    subprocess.run(["g++", cpp_file_5, "-o", cpp_executable_5], check=True)


subprocess.run([f"./{cpp_executable_1}", input_file, reduced], check=True)

subprocess.run([f"./{cpp_executable_2}", reduced, domains, centroids], check=True)

subprocess.run([f"./{cpp_executable_5}", input_file,centroids, incoherence,incoherence_matrix], check=True)

print("Analysis complete.")
