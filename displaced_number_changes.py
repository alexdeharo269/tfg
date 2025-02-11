import numpy as np
import sys

def read_matrix_from_file(input_file, dtype=int):
    """
    Reads a whitespace-separated matrix from file.
    Each line is split into numbers of type 'dtype'.
    Returns a NumPy array.
    """
    with open(input_file, 'r') as file:
        matrix = [list(map(dtype, line.strip().split()))
                  for line in file if line.strip()]
    return np.array(matrix)

def replicate_rows(matrix, replication_factor):
    """
    Replicates each row of the given matrix 'replication_factor' times.
    (Vertical replication.)
    """
    return np.repeat(matrix, replication_factor, axis=0)

def find_nearest_centroid(row, j):
    """
    Given a 1D NumPy array 'row' (from the replicated centroid matrix)
    and a column index j, return the column index of the nearest nonzero entry.
    If no nonzero entry exists, returns None.
    """
    indices = np.nonzero(row)[0]
    if len(indices) == 0:
        return None
    # Find the index with the smallest absolute difference from j.
    diffs = np.abs(indices - j)
    nearest = indices[np.argmin(diffs)]
    return nearest

def count_displaced_column_changes(original_matrix, centroid_matrix):
    """
    For each column (with the original column index) of the full–resolution original matrix,
    count the number of rows (from row 1 onward) where the value changes
    when comparing the current row to the previous row—but with the previous row's
    column index displaced by the change in the nearest centroid position.
    
    In detail, for each row i (starting at 1) and each column j:
      - Find the nearest centroid in row i–1 and row i from the centroid matrix.
      - Compute delta = (current row's nearest centroid) – (previous row's nearest centroid).
        (If either is not found, set delta = 0.)
      - Compare original_matrix[i, j] with original_matrix[i–1, (j+delta) mod num_columns].
      - If they differ, count that as one change for column j.
    
    Returns:
      A NumPy array (of length num_columns) with the change count for each column.
    """
    num_rows, num_cols = original_matrix.shape
    change_counts = np.zeros(num_cols, dtype=int)
    
    for i in range(1, num_rows):
        for j in range(num_cols):
            prev_centroid = find_nearest_centroid(centroid_matrix[i-1], j)
            curr_centroid = find_nearest_centroid(centroid_matrix[i], j)
            if prev_centroid is None or curr_centroid is None:
                delta = 0
            else:
                delta = curr_centroid - prev_centroid
            shifted_index = (j + delta) % num_cols
            if original_matrix[i, j] != original_matrix[i-1, shifted_index]:
                change_counts[j] += 1
    return change_counts

def write_changes_to_file(output_file, change_counts):
    """
    Writes the change counts to a file.
    Each line will have the column index (original index) and the change count.
    """
    with open(output_file, 'w') as file:
        num_columns = len(change_counts)
        for j in range(num_columns):
            file.write(f"{j} {change_counts[j]}\n")

if __name__ == "__main__":
    # --- File names (adjust paths as needed) ---
    # The full–resolution binary matrix file (e.g., ising_dataR39.dat)
    original_file = sys.argv[1]
    # The domain centroids file (produced from the reduced matrix via DBSCAN)
    centroid_file = sys.argv[2]
    # The output file for the displaced change counts
    output_file = sys.argv[3]
    
    # --- Read matrices ---
    # Original full-resolution binary matrix (0s and 1s)
    original_matrix = read_matrix_from_file(original_file, dtype=int)
    # Centroid matrix from DBSCAN (each row has zeros except at centroid columns)
    centroid_matrix_reduced = read_matrix_from_file(centroid_file, dtype=int)
    
    # --- Replicate the reduced centroid matrix to full resolution ---
    # (Each row in the reduced centroid file corresponds to 50 original rows.)
    replication_factor = 50
    centroid_matrix_full = replicate_rows(centroid_matrix_reduced, replication_factor)
    
    if centroid_matrix_full.shape[0] != original_matrix.shape[0]:
        print("Warning: The replicated centroid matrix rows do not match the original matrix rows.")
    
    # --- Count displaced column changes ---
    change_counts = count_displaced_column_changes(original_matrix, centroid_matrix_full)
    
    # --- Write the change counts to file ---
    write_changes_to_file(output_file, change_counts)
    print("Processing complete. Change counts written to", output_file)
