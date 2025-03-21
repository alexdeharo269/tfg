import subprocess
import os
import tempfile

#Quiero que en el fichero del analisis salga N, Dinámica, Temperatura, Tiempo hasta estacionario en un fichero aparte resumen.txt
#Sistematizar analisis, incoherence_matrix tiene q salir perfe. 
#poner todo el output en el fichero del metadata usando el flag de la sim_id.
#comprobar bien bien aquí porque para sistematizarlo no hacen falta todos los ficheros de output, solo el parametro para el diagrama de fases.
#podria intentar hacer un diagrama de fases preliminar pero para que fuera bueno tendria que correrlo con varias semillas. 


# config.py
sim_id=0
debug=False
# Path to the first input file
input_file =f"./ising_data/simulations/{sim_id}.dat"

# Folder for intermediate and final output files
output_folder = "./ising_data/simulations"

# Define output files
reduced = f"{output_folder}/reduced/reduced{sim_id}.dat"
domains = f"{output_folder}/domains/domains{sim_id}.dat"
centroids= f"{output_folder}/centroids/centroids{sim_id}.dat"
number_changes= f"{output_folder}/changes/changes{sim_id}.dat"
incoherence= f"{output_folder}/incoherence/incoherence{sim_id}.dat"
incoherence_matrix= f"{output_folder}/incoherence_matrix/incoherence_matrix{sim_id}.dat"

# Create output directories for final outputs only
for directory in [os.path.dirname(p) for p in [incoherence]]:
    os.makedirs(directory, exist_ok=True)

# Create directories for intermediate files if we're not using pipes
if  debug:
    for directory in [os.path.dirname(p) for p in [reduced, domains, centroids, incoherence_matrix]]:
        os.makedirs(directory, exist_ok=True)

cpp_file_1 = "compute.cpp"   #Compila la matriz reducida, de momento innecesario porque lo meto en test_process_cpp
cpp_executable_1 = "compute"
cpp_file_2 = "test_process_cpp.cpp"
cpp_executable_2 = "test_process"
cpp_file_3 = "displaced_number_changes_cpp.cpp"
cpp_executable_3 = "displaced_number_changes" #Se centra en el cambio de valor en de las posiciones ajustadas a los centroides.
#esta bien pero tiene en cuenta las fronteras, por lo que tiene máximos ahí
cpp_file_4 = "clusters_centrois_incoherence.cpp"
cpp_executable_4 = "clusters_centroids_incoherence"
cpp_file_5="chnages_incoherence.cpp" #Se centra en el cambio respecto al modo de un cluster, por lo que a priori no tiene en cuenta las fronteras. 
cpp_executable_5="changes_incoherence"


if not os.path.exists(cpp_executable_2):  # Compile only if necessary
    subprocess.run(["g++", cpp_file_2, "-o", cpp_executable_2], check=True)
if  not os.path.exists(cpp_executable_5):  # Compile only if necessary
    subprocess.run(["g++", cpp_file_5, "-o", cpp_executable_5], check=True)

if not debug:
    debug_arg = "false"
    # Create a named temporary file for centroids
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        temp_centroids = tmp.name
    subprocess.run(
        [f"./{cpp_executable_2}", input_file, reduced, domains, temp_centroids, debug_arg],
        check=True
    )
    subprocess.run(
        [f"./{cpp_executable_5}", input_file, temp_centroids, incoherence, incoherence_matrix, debug_arg],
        check=True
    )
    os.unlink(temp_centroids)
else:
    print("Using file-based approach for debugging/testing...")
    # Traditional approach with intermediate files
    #subprocess.run([f"./{cpp_executable_1}", input_file, reduced], check=True)
    #subprocess.run([f"./{cpp_executable_2}", reduced, domains, centroids], check=True)
    debug_arg = "true" 

    subprocess.run([f"./{cpp_executable_2}", input_file, reduced,domains,centroids,debug_arg], check=True) #Poner que si debug no imprima
    subprocess.run([f"./{cpp_executable_5}", input_file, centroids, incoherence, incoherence_matrix,debug_arg], check=True)

    print("Analysis complete.")