import numpy as np
from dbscan1d.core import DBSCAN1D
'''Falta meter condiciones de contorno periodicas'''
def read_reduced_matrix(filename):
    """
    Read a reduced matrix from a file.
    Assumes the file contains whitespace‐separated floats.
    """
    return np.loadtxt(filename)

def process_row(row_data, eps, min_samples, threshold):
    """
    Process a single row (1D array) with DBSCAN1D.
    - Clusters are determined by DBSCAN1D.
    - For each cluster, the centroid is computed.
      If the centroid >= threshold, we assign the entire cluster domain 1,
      otherwise 0.
    - Noise points (label == -1) are assigned based on their own value.
    Returns a 1D binary array with the domain (0 or 1) for each point.
    """
    # Run DBSCAN1D on the row
    dbs = DBSCAN1D(eps=eps, min_samples=min_samples)
    labels = dbs.fit_predict(row_data)
    
    # For each detected cluster (ignoring noise, label == -1),
    # decide its domain based on its centroid.
    unique_labels = np.unique(labels)
    cluster_domain = {}
    for label in unique_labels:
        if label == -1:
            continue  # noise will be handled separately
        indices = np.where(labels == label)[0]
        if indices.size > 0:
            centroid = np.mean(row_data[indices])
            # Domain is 1 if centroid >= threshold; else 0.
            cluster_domain[label] = 1 if centroid >= threshold else 0

    # Build the binary output for the row.
    binary_row = np.empty_like(row_data, dtype=int)
    for idx, label in enumerate(labels):
        if label == -1:
            # For noise, assign based on the actual value.
            binary_row[idx] = 1 if row_data[idx] >= threshold else 0
        else:
            binary_row[idx] = cluster_domain[label]
    return binary_row

def process_reduced_matrix(reduced_matrix, eps, min_samples, threshold):
    """
    Process every row of the reduced matrix.
    Returns a binary matrix (same shape as reduced_matrix) where each element is 0 or 1.
    """
    num_rows, num_cols = reduced_matrix.shape
    binary_matrix = np.zeros((num_rows, num_cols), dtype=int)
    
    for i in range(num_rows):
        row_data = reduced_matrix[i, :]
        binary_row = process_row(row_data, eps=eps, min_samples=min_samples, threshold=threshold)
        binary_matrix[i, :] = binary_row
        
    return binary_matrix

def write_binary_matrix(binary_matrix, filename):
    """
    Write the binary matrix to a text file.
    """
    np.savetxt(filename, binary_matrix, fmt='%d')

if __name__ == "__main__":
    # File names (adjust the paths as needed)
    reduced_matrix_file = "./ising_data_R_low_temp/reduced_matrix.dat"
    output_binary_file = "./ising_data_R_low_temp/binary_domains_matrix.dat"
    
    # Read the reduced matrix from file.
    reduced_matrix = read_reduced_matrix(reduced_matrix_file)
    print(f"Read reduced matrix of shape: {reduced_matrix.shape}")
    
    # Process each row with DBSCAN1D.
    # Adjust eps and min_samples if needed so that clusters separate regions near 0 and near 1.
    binary_matrix = process_reduced_matrix(reduced_matrix, eps=0.1, min_samples=50, threshold=0.5)
    
    # Write the binary domains matrix to file.
    write_binary_matrix(binary_matrix, output_binary_file)
    print(f"Binary domains matrix written to {output_binary_file}")
