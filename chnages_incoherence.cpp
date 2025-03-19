#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <utility>

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
            cout << "Time: " << i/num_rows << "\r";
        }
    }
    return change_counts;
}

//---------------------------------------------------------------------
// Function: countIncoherencePerColumn
//
// For each row of the original simulation matrix (binary data assumed),
// groups the sites by their domain label (from domain_label_matrix). For each domain,
// computes the mode (if the average of values in that domain is >= 0.5, mode=1, else 0).
// Then, for every column in that row that belongs to a domain, if the value deviates
// from the domain's mode, an incoherence count is added for that column.
// The function returns a vector (one entry per column) with the total incoherence count.
//---------------------------------------------------------------------
vector<int> countDisplacedIncoherencePerColumn(const vector<vector<int>> &original_matrix,
                                               const vector<vector<int>> &centroid_matrix,
                                               const vector<vector<int>> &domain_label_matrix)
{
    int num_rows = static_cast<int>(original_matrix.size());
    if (num_rows == 0)
        return vector<int>{};
    int num_cols = static_cast<int>(original_matrix[0].size());

    vector<int> incoherence_counts(num_cols, 0);

    // Start from row 1 (since we need a previous row to compute displacement)
    for (int i = 1; i < num_rows; i++)
    {
        // For each row, we create a mapping from "shifted column index" to a vector of values.
        unordered_map<int, vector<int>> shiftedDomainValues;
        for (int j = 0; j < num_cols; j++)
        {
            // Compute displacement delta as in countDisplacedColumnChanges:
            int prev_centroid = findNearestCentroid(centroid_matrix[i - 1], j);
            int curr_centroid = findNearestCentroid(centroid_matrix[i], j);
            int delta = (prev_centroid == -1 || curr_centroid == -1) ? 0 : (curr_centroid - prev_centroid);
            int shifted_index = (j + delta) % num_cols;
            if (shifted_index < 0)
                shifted_index += num_cols;

            // Use the domain label from the current row (or, if desired, from the previous row)
            int label = domain_label_matrix[i][j];
            // Group the current row's value into the bucket corresponding to the shifted index.
            // (This aligns the domain spatially with the previous row.)
            shiftedDomainValues[shifted_index].push_back(original_matrix[i][j]);
        }
        // Now, for each shifted group, compute the domain mode and count incoherent entries.
        for (const auto &kv : shiftedDomainValues)
        {
            int col = kv.first; // the (shifted) column index
            const vector<int> &values = kv.second;
            double sum = 0.0;
            for (int v : values)
                sum += v;
            double avg = sum / values.size();
            // Ensure that if the domain is completely zero, mode is 0.
            int mode;
            if (sum == 0.0)
                mode = 0;
            else
                mode = (avg >= 0.5) ? 1 : 0;
            // Count as incoherent any entry that deviates from the mode.
            for (int v : values)
            {
                if (abs(v - mode) > 0)
                    incoherence_counts[col]++;
            }
        }
    }
    return incoherence_counts;
}

//---------------------------------------------------------------------
// Function: computeIncoherenceMatrix
//
// For each row of the full–resolution original binary matrix, groups
// columns by their domain label (from domain_label_matrix). For each
// group, computes the mode (using a threshold of 0.5) and then marks
// as "incoherent" (output 1) every element that deviates from the mode
// (beyond a given tolerance). Returns a matrix (of same shape) of incoherence flags.
//---------------------------------------------------------------------
vector<vector<int>> computeIncoherenceMatrix(const vector<vector<int>> &original_matrix,
                                             const vector<vector<int>> &domain_label_matrix,
                                             double tolerance = 0.0)
{
    int num_rows = original_matrix.size();
    if (num_rows == 0)
        return vector<vector<int>>();
    int num_cols = original_matrix[0].size();
    vector<vector<int>> incoherence_matrix(num_rows, vector<int>(num_cols, 0));

    for (int i = 0; i < num_rows; i++)
    {
        // Group column indices by domain label (including label 0).
        unordered_map<int, vector<int>> groups;
        for (int j = 0; j < num_cols; j++)
        {
            int label = domain_label_matrix[i][j];
            groups[label].push_back(j);
        }
        // For each group, compute the mode and assign incoherence flags.
        for (const auto &kv : groups)
        {
            const vector<int> &indices = kv.second;
            double sum = 0.0;
            for (int j : indices)
                sum += original_matrix[i][j];
            double avg = sum / indices.size();
            // For binary data, this ensures that if the group is all 0s, mode==0.
            int mode = (avg >= 0.5) ? 1 : 0;
            // Now, mark each position in this group as incoherent (1) if its value
            // deviates from the mode (by more than tolerance), otherwise 0.
            for (int j : indices)
            {
                if (abs(original_matrix[i][j] - mode) > tolerance)
                    incoherence_matrix[i][j] = 1;
                else
                    incoherence_matrix[i][j] = 0;
            }
        }
    }
    return incoherence_matrix;
}

//---------------------------------------------------------------------
// Function: writeMetricsToFile
//
// Writes two metrics (displaced changes and incoherence counts) to a file.
// Each line contains the column index, displaced change count, and incoherence count.
//---------------------------------------------------------------------
void writeMetricsToFile(const string& output_file,
                        const vector<int>& displaced_changes,
                        const vector<int>& incoherence_counts) {
    ofstream file(output_file);
    if (!file) {
        cerr << "Error opening output file: " << output_file << endl;
        exit(EXIT_FAILURE);
    }
    int num_columns = static_cast<int>(displaced_changes.size());
    for (int j = 0; j < num_columns; j++) {
        file << j << " " << displaced_changes[j] << " " << incoherence_counts[j] << "\n";
    }
    file.close();
}

//------------------------------------------------------------------------------
// Function: write_domain_label_matrix
//
// Writes an integer matrix to a text file (one row per line, numbers separated by spaces).
//------------------------------------------------------------------------------
void write_incoherece_matrix(const vector<vector<int>> &matrix, const string &filename)
{
    ofstream outfile(filename);
    if (!outfile)
    {
        cerr << "Error opening output file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    for (const auto &row : matrix)
    {
        for (size_t i = 0; i < row.size(); ++i)
        {
            outfile << row[i];
            if (i < row.size() - 1)
                outfile << " ";
        }
        outfile << "\n";
    }
    outfile.close();
}

//Read from temp file for pipe files based approach
vector<vector<int>> readMatrixFromStdin()
{
    vector<vector<int>> matrix;
    string line;
    while (getline(cin, line))
    {
        stringstream ss(line);
        vector<int> row;
        int value;
        while (ss >> value)
        {
            row.push_back(value);
        }
        if (!row.empty())
        {
            matrix.push_back(row);
        }
    }
    return matrix;
}
//---------------------------------------------------------------------
// (Existing) Main function for displaced number changes, now extended.
// 
// Expects three command-line arguments:
//   1. The full-resolution original binary matrix file (e.g., ising_dataR39.dat)
//   2. The domain centroids file (from the reduced matrix via DBSCAN)
//   3. The output file for the metrics (columns: index, displaced changes, incoherence)
//---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <original_file> <centroid_file> <output_file> <incoherence_matrix> <debug_s>" << endl;
        return EXIT_FAILURE;
    }
    
    string original_file = argv[1];
    string centroid_file = argv[2];
    string output_file = argv[3];
    string incoherence_matrix = argv[4];
    string debug_s = argv[5];
    bool debug;
    if (debug_s == "true")
    {
        debug = true;
    }
    else
    {
        debug = false;
    }
    // Read matrices.
    vector<vector<int>> original_matrix = readMatrixFromFile(original_file);
    vector<vector<int>> centroid_matrix_reduced;
    if (centroid_file == "-") {
        centroid_matrix_reduced = readMatrixFromStdin();
    } else {
        centroid_matrix_reduced = readMatrixFromFile(centroid_file);
    }    
    // Replicate the reduced centroid matrix to full resolution.
    // (Assuming each row in the reduced centroid file corresponds to 50 original rows.)
    int replication_factor = 50;
    vector<vector<int>> centroid_matrix_full = replicateRows(centroid_matrix_reduced, replication_factor);
    
    if (centroid_matrix_full.size() != original_matrix.size()) {
        cout << "Warning: The replicated centroid matrix rows (" << centroid_matrix_full.size()
             << ") do not match the original matrix rows (" << original_matrix.size() << ")." << endl;
    }
    
    // Compute displaced column changes.
    vector<int> displaced_changes = countDisplacedColumnChanges(original_matrix, centroid_matrix_full);
    
    // For strength of incoherence, use the domain label matrix.
    // Here we assume that the domain labels are the same as in the centroid matrix.
    // (In your pipeline, you might get the domain labels from your DBSCAN code.)
    vector<vector<int>> domain_label_matrix = centroid_matrix_full; // For demonstration.
    vector<int> incoherence_counts = countDisplacedIncoherencePerColumn(original_matrix, centroid_matrix_full,domain_label_matrix);
    vector<vector<int>> incoh_matrix = computeIncoherenceMatrix(original_matrix, domain_label_matrix, 0.0);

    // Write both metrics to the output file (three columns).
    writeMetricsToFile(output_file, displaced_changes, incoherence_counts);
    if(debug==true){
        write_incoherece_matrix(incoh_matrix, incoherence_matrix);
    }
    
    return EXIT_SUCCESS;
}
