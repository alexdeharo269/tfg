#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cmath>

using namespace std;

//---------------------------------------------------------------------
// Function: readMatrixFromFile
//
// Reads a whitespace-separated matrix from a file. Each line is split
// into numbers (of type int). Returns a 2D vector representing the matrix.
//---------------------------------------------------------------------
vector<vector<int>> readMatrixFromFile(const string& input_file) {
    ifstream file(input_file);
    if (!file) {
        cerr << "Error opening file: " << input_file << endl;
        exit(EXIT_FAILURE);
    }
    vector<vector<int>> matrix;
    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        vector<int> row;
        int value;
        while (iss >> value)
            row.push_back(value);
        if (!row.empty())
            matrix.push_back(row);
    }
    file.close();
    return matrix;
}

//---------------------------------------------------------------------
// Function: replicateRows
//
// Replicates each row of the given matrix 'replication_factor' times.
// (Vertical replication.)
//---------------------------------------------------------------------
vector<vector<int>> replicateRows(const vector<vector<int>>& matrix, int replication_factor) {
    vector<vector<int>> replicated;
    for (const auto& row : matrix) {
        for (int i = 0; i < replication_factor; i++) {
            replicated.push_back(row);
        }
    }
    return replicated;
}

//---------------------------------------------------------------------
// Function: findNearestCentroid
//
// Given a row (from the replicated centroid matrix) and a column index j,
// returns the column index of the nearest nonzero entry. If no nonzero
// entry exists, returns -1.
//---------------------------------------------------------------------
int findNearestCentroid(const vector<int>& row, int j) {
    vector<int> indices;
    for (size_t i = 0; i < row.size(); i++) {
        if (row[i] != 0)
            indices.push_back(static_cast<int>(i));
    }
    if (indices.empty())
        return -1;  // Equivalent to Python's None.
    
    int nearest = indices[0];
    int minDiff = abs(nearest - j);
    for (int idx : indices) {
        int diff = abs(idx - j);
        if (diff < minDiff) {
            minDiff = diff;
            nearest = idx;
        }
    }
    return nearest;
}

//---------------------------------------------------------------------
// Function: countDisplacedColumnChanges
//
// For each column of the full–resolution original matrix, counts the number
// of rows (from row 1 onward) where the value changes when comparing the
// current row to the previous row, with the previous row's column index shifted
// by the change in the nearest centroid position.
//---------------------------------------------------------------------
vector<int> countDisplacedColumnChanges(const vector<vector<int>>& original_matrix,
                                          const vector<vector<int>>& centroid_matrix) {
    int num_rows = static_cast<int>(original_matrix.size());
    if (num_rows == 0)
        return vector<int>{};
    int num_cols = static_cast<int>(original_matrix[0].size());
    
    vector<int> change_counts(num_cols, 0);
    
    for (int i = 1; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            int prev_centroid = findNearestCentroid(centroid_matrix[i - 1], j);
            int curr_centroid = findNearestCentroid(centroid_matrix[i], j);
            int delta = 0;
            if (prev_centroid == -1 || curr_centroid == -1)
                delta = 0;
            else
                delta = curr_centroid - prev_centroid;
            
            int shifted_index = (j + delta) % num_cols;
            if (shifted_index < 0)  // Ensure non-negative index.
                shifted_index += num_cols;
            
            if (original_matrix[i][j] != original_matrix[i - 1][shifted_index])
                change_counts[j]++;
        }
        if (i % 1000 == 0){
            cout << "Processed row " << i << " of " << num_rows << "\r";
        }
    }
    return change_counts;
}

//---------------------------------------------------------------------
// Function: writeChangesToFile
//
// Writes the change counts to a file. Each line contains the column index
// and the corresponding change count.
//---------------------------------------------------------------------
void writeChangesToFile(const string& output_file, const vector<int>& change_counts) {
    ofstream file(output_file);
    if (!file) {
        cerr << "Error opening output file: " << output_file << endl;
        exit(EXIT_FAILURE);
    }
    int num_columns = static_cast<int>(change_counts.size());
    for (int j = 0; j < num_columns; j++) {
        file << j << " " << change_counts[j] << "\n";
    }
    file.close();
}

//---------------------------------------------------------------------
// Main function
//
// Expects three command–line arguments:
//   1. The full–resolution original binary matrix file (e.g., ising_dataR39.dat)
//   2. The domain centroids file (from the reduced matrix via DBSCAN)
//   3. The output file for the displaced change counts
//---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <original_file> <centroid_file> <output_file>" << endl;
        return EXIT_FAILURE;
    }
    cout << "Processing..." << endl;
    
    string original_file = argv[1];
    string centroid_file = argv[2];
    string output_file = argv[3];
    
    // Read matrices.
    vector<vector<int>> original_matrix = readMatrixFromFile(original_file);
    vector<vector<int>> centroid_matrix_reduced = readMatrixFromFile(centroid_file);
    
    // Replicate the reduced centroid matrix to full resolution.
    // (Each row in the reduced centroid file corresponds to 50 original rows.)
    int replication_factor = 50;
    vector<vector<int>> centroid_matrix_full = replicateRows(centroid_matrix_reduced, replication_factor);
    
    if (centroid_matrix_full.size() != original_matrix.size()) {
        cout << "Warning: The replicated centroid matrix rows do not match the original matrix rows." << endl;
    }
    
    // Count displaced column changes.
    vector<int> change_counts = countDisplacedColumnChanges(original_matrix, centroid_matrix_full);
    
    // Write the change counts to file.
    writeChangesToFile(output_file, change_counts);
    
    cout << "Processing complete. Change counts written to " << output_file << endl;
    return EXIT_SUCCESS;
}
