import numpy as np
import math

###########################################
# A custom DBSCAN1D for indices with periodic boundary conditions.
###########################################
def cluster_positions_periodic(positions, n_total, eps, min_samples):
    """
    Cluster a sorted array of 1D positions (indices) on a circle of length n_total.
    Distance is defined as:
      d = min(|a - b|, n_total - |a - b|)
    Uses a DBSCAN–like expansion.
    
    Parameters:
      positions   : 1D numpy array of (sorted) indices
      n_total     : total number of points (period)
      eps         : maximum distance (in index units) for neighbors
      min_samples : minimum neighbors required for a core point
      
    Returns:
      labels : numpy array of the same length as positions, with cluster labels.
               Points not assigned are labeled -1.
    """
    L = len(positions)
    if L == 0:
        return np.array([], dtype=int)
    
    labels = np.full(L, -1, dtype=int)
    visited = np.zeros(L, dtype=bool)
    cluster_id = 0

    def pdist(a, b):
        d = abs(a - b)
        return min(d, n_total - d)
    
    for i in range(L):
        if visited[i]:
            continue
        neighbors = []
        for j in range(L):
            if pdist(positions[i], positions[j]) <= eps:
                neighbors.append(j)
        if len(neighbors) < min_samples:
            visited[i] = True
            continue
        labels[i] = cluster_id
        visited[i] = True
        seed_set = set(neighbors)
        seed_set.discard(i)
        while seed_set:
            j = seed_set.pop()
            if not visited[j]:
                visited[j] = True
                neighbors_j = []
                for k in range(L):
                    if pdist(positions[j], positions[k]) <= eps:
                        neighbors_j.append(k)
                if len(neighbors_j) >= min_samples:
                    seed_set.update(neighbors_j)
            if labels[j] == -1:
                labels[j] = cluster_id
        cluster_id += 1
    return labels

###########################################
# Process a single row with periodic DBSCAN.
###########################################
def process_row_dbscan_periodic(row, eps=5, min_samples=3, threshold=0.5):
    """
    Process a single row (1D numpy array) using DBSCAN with periodic boundary conditions.
    
    Steps:
      1. Binarize the row: values >= threshold → 1, else 0.
      2. For each binary group (0 and 1), extract indices and run DBSCAN on these positions.
      3. For each found cluster, compute the circular centroid.
      4. Only accept clusters with size >= (n/10); reject smaller clusters.
      5. Build a domain label array (same length as row) where accepted clusters are assigned a unique domain label and
         all noise or small clusters are set to 0.
      6. Return also a list of cluster info dictionaries.
    
    Returns:
      domain_labels : 1D integer array with accepted domain labels (0 means no domain)
      clusters_info : list of dictionaries for each accepted cluster with keys:
                      'domain_label', 'centroid', 'domain_value', 'indices', 'length'
    """
    n = len(row)
    binary = (row >= threshold).astype(int)
    domain_labels = np.full(n, -1, dtype=int)
    clusters_info = []
    current_domain_label = 1  # Use positive numbers for accepted clusters; 0 will mean no cluster.

    # Process each binary group separately (for 0 and for 1)
    for domain_value in [0, 1]:
        idx = np.where(binary == domain_value)[0]
        if len(idx) == 0:
            continue
        cluster_labels = cluster_positions_periodic(idx, n, eps, min_samples)
        for cl in np.unique(cluster_labels):
            if cl == -1:
                continue  # Skip noise
            cluster_idx = idx[cluster_labels == cl]
            # Enforce minimum domain size
            if len(cluster_idx) < n/10 or len(cluster_idx)>n/1.5:
                continue
            # Compute circular centroid
            angles = 2 * np.pi * cluster_idx / n
            sin_mean = np.mean(np.sin(angles))
            cos_mean = np.mean(np.cos(angles))
            mean_angle = math.atan2(sin_mean, cos_mean)
            if mean_angle < 0:
                mean_angle += 2 * np.pi
            centroid = (mean_angle / (2 * np.pi)) * n
            
            clusters_info.append({
                'domain_label': current_domain_label,
                'centroid': centroid,
                'domain_value': domain_value,
                'indices': cluster_idx,
                'length': len(cluster_idx)
            })
            domain_labels[cluster_idx] = current_domain_label
            current_domain_label += 1
            
    # Replace any remaining -1 (noise or rejected small domains) with 0.
    domain_labels[domain_labels < 0] = 0
    return domain_labels, clusters_info

###########################################
# Process an entire matrix row-by-row.
###########################################
def process_matrix_dbscan_periodic(reduced_matrix, eps=5, min_samples=3, threshold=0.5):
    """
    Process every row of the reduced matrix using the periodic DBSCAN routine.
    
    Returns:
      domain_label_matrix : integer matrix (same shape as reduced_matrix) with domain labels (0 for none).
      centroids_all       : list (one per row) of clusters_info lists.
    """
    num_rows, num_cols = reduced_matrix.shape
    domain_label_matrix = np.full((num_rows, num_cols), 0, dtype=int)
    centroids_all = []
    
    for i in range(num_rows):
        row = reduced_matrix[i, :]
        dlabels, clusters_info = process_row_dbscan_periodic(row, eps=eps, min_samples=min_samples, threshold=threshold)
        domain_label_matrix[i, :] = dlabels
        centroids_all.append(clusters_info)
    return domain_label_matrix, centroids_all

###########################################
# Build a centroid matrix.
###########################################
def build_centroid_matrix(centroids_all, num_cols,RV):
    """
    Build a matrix (one row per time step) that is zero except at the (rounded) centroid positions.
    At each centroid position, place the corresponding domain label.
    
    Parameters:
      centroids_all : list (one per row) of clusters_info lists.
      num_cols      : number of columns in the reduced matrix.
      
    Returns:
      A binary (sparse) integer matrix of shape (num_rows, num_cols).
    """
    num_rows = len(centroids_all)
    centroid_matrix = np.zeros((num_rows, num_cols), dtype=int)
    
    for i, clusters_info in enumerate(centroids_all):
        for info in clusters_info:
            # Round the centroid to the nearest column index and wrap around using modulo.
            col = int(round(info['centroid'])) % num_cols
            for j in range(-RV,RV):
                centroid_matrix[i, col+j] = info['domain_label']
            
    return centroid_matrix

###########################################
# File I/O functions.
###########################################
def read_reduced_matrix(filename):
    """
    Read the reduced matrix from a text file.
    Assumes whitespace-separated floats.
    """
    return np.loadtxt(filename)

def write_domain_label_matrix(matrix, filename):
    """
    Write the domain label matrix to a text file.
    """
    np.savetxt(filename, matrix, fmt='%d')

def write_centroid_matrix(matrix, filename):
    """
    Write the centroid matrix to a text file.
    """
    np.savetxt(filename, matrix, fmt='%d')

###########################################
# Main.
###########################################
if __name__ == "__main__":
    # --- File names (adjust paths as needed) ---
    reduced_matrix_file   = "./ising_data_R_low_temp/reduced_matrix.dat"
    domain_label_file     = "./ising_data_R_low_temp/domain_labels_dbscan.dat"
    centroid_matrix_file  = "./ising_data_R_low_temp/domain_centroids_dbscan.dat"
    RV = 5
    
    # --- Read the reduced matrix ---
    reduced_matrix = read_reduced_matrix(reduced_matrix_file)
    num_rows, num_cols = reduced_matrix.shape
    print(f"Read reduced matrix with shape {reduced_matrix.shape}")
    
    # --- Process the matrix row-by-row using periodic DBSCAN ---
    # Adjust eps (in index units), min_samples, and threshold as needed.
    domain_label_matrix, centroids_all = process_matrix_dbscan_periodic(
        reduced_matrix, eps=5, min_samples=10, threshold=0.5)
    
    # --- Build the centroid matrix ---
    centroid_matrix = build_centroid_matrix(centroids_all, num_cols,RV)
    
    # --- Write the output files ---
    write_domain_label_matrix(domain_label_matrix, domain_label_file)
    write_centroid_matrix(centroid_matrix, centroid_matrix_file)
    
    print("Processing complete.")

