#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <unordered_set>
#include <vector>
#include <utility>
#include <algorithm>

// Use the std namespace for brevity.
using namespace std;

// Define PI if not available.
const double PI = 3.14159265358979323846;

//------------------------------------------------------------------------------
// Structure to hold cluster information (similar to the Python dictionary).
//------------------------------------------------------------------------------
struct ClusterInfo {
    int domain_label;         // A unique label for this accepted cluster.
    double centroid;          // Circular centroid computed from the indices.
    int domain_value;         // The binary value (0 or 1) corresponding to this cluster.
    vector<int> indices;      // The indices (positions) that belong to this cluster.
    int length;               // The size of the cluster (redundant with indices.size()).
};

//------------------------------------------------------------------------------
// Function: cluster_positions_periodic
//
// A DBSCAN–like clustering of 1D positions (indices) on a circle of length n_total.
// The distance between two indices is defined as:
//    d = min(|a - b|, n_total - |a - b|)
// Positions within eps (in index units) are considered neighbors.
// A point must have at least min_samples neighbors (including itself) to be a core point.
//------------------------------------------------------------------------------
vector<int> cluster_positions_periodic(const vector<int>& positions, int n_total, int eps, int min_samples) {
    size_t L = positions.size();
    vector<int> labels(L, -1);
    vector<bool> visited(L, false);
    int cluster_id = 0;

    // Lambda for periodic distance.
    auto pdist = [n_total](int a, int b) -> int {
        int d = abs(a - b);
        return min(d, n_total - d);
    };

    for (size_t i = 0; i < L; ++i) {
        if (visited[i])
            continue;

        // Find neighbors for positions[i]
        vector<int> neighbors;
        for (size_t j = 0; j < L; ++j) {
            if (pdist(positions[i], positions[j]) <= eps)
                neighbors.push_back(j);
        }

        // If there are not enough neighbors, mark as visited (noise) and continue.
        if (neighbors.size() < static_cast<size_t>(min_samples)) {
            visited[i] = true;
            continue;
        }

        // Start a new cluster.
        labels[i] = cluster_id;
        visited[i] = true;

        // Use an unordered_set to hold the "seed" indices.
        unordered_set<int> seed_set(neighbors.begin(), neighbors.end());
        // Remove the current point from the seed set.
        seed_set.erase(i);

        while (!seed_set.empty()) {
            // Get and remove one element from the seed set.
            int j = *seed_set.begin();
            seed_set.erase(seed_set.begin());

            if (!visited[j]) {
                visited[j] = true;
                // Find neighbors of positions[j]
                vector<int> neighbors_j;
                for (size_t k = 0; k < L; ++k) {
                    if (pdist(positions[j], positions[k]) <= eps)
                        neighbors_j.push_back(k);
                }
                if (neighbors_j.size() >= static_cast<size_t>(min_samples)) {
                    // Add all neighbors of j to the seed set.
                    for (int nb : neighbors_j)
                        seed_set.insert(nb);
                }
            }
            // If j was not yet assigned a cluster, assign it.
            if (labels[j] == -1) {
                labels[j] = cluster_id;
            }
        }
        cluster_id++;
    }

    return labels;
}

//------------------------------------------------------------------------------
// Function: process_row_dbscan_periodic
//
// Processes a single 1D row (vector of doubles) using a periodic DBSCAN.
// Steps:
//   1. Binarize the row: value >= threshold becomes 1, else 0.
//   2. For each binary group (0 and 1), extract the indices and run the periodic
//      DBSCAN routine.
//   3. For each cluster found, compute the circular centroid (using angular means).
//   4. Only accept clusters whose size is between n/10 and n/1.5.
//   5. Build a domain label array (same size as the row) where accepted clusters
//      are given unique labels (starting at 1) and all other positions are set to 0.
//   6. Also return a vector of ClusterInfo for each accepted cluster.
//------------------------------------------------------------------------------
pair<vector<int>, vector<ClusterInfo>>
process_row_dbscan_periodic(const vector<double>& row, int eps = 5, int min_samples = 3, double threshold = 0.5) {
    int n = static_cast<int>(row.size());
    vector<int> binary(n, 0);
    // Binarize the row.
    for (int i = 0; i < n; ++i) {
        binary[i] = (row[i] >= threshold) ? 1 : 0;
    }

    vector<int> domain_labels(n, -1);
    vector<ClusterInfo> clusters_info;
    int current_domain_label = 1;  // Domain labels start at 1; 0 means no accepted domain.

    // Process both binary groups: first 0, then 1.
    for (int domain_value = 0; domain_value <= 1; ++domain_value) {
        vector<int> idx;  // Store indices in the row with the current binary value.
        for (int i = 0; i < n; ++i) {
            if (binary[i] == domain_value)
                idx.push_back(i);
        }
        if (idx.empty())
            continue;

        // Cluster the positions in idx.
        vector<int> cluster_labels = cluster_positions_periodic(idx, n, eps, min_samples);

        // Get unique cluster labels (except -1 which denotes noise).
        set<int> unique_labels;
        for (int lab : cluster_labels)
            unique_labels.insert(lab);

        // Process each found cluster.
        for (int cl : unique_labels) {
            if (cl == -1)
                continue;  // Skip noise.

            // Extract the indices corresponding to this cluster.
            vector<int> cluster_idx;
            for (size_t j = 0; j < cluster_labels.size(); ++j) {
                if (cluster_labels[j] == cl)
                    cluster_idx.push_back(idx[j]);
            }
            // Enforce minimum (and maximum) domain size:
            // Accept clusters only if size >= n/10 and <= n/1.5.
            if (cluster_idx.size() < static_cast<size_t>(n / 50.0) ||
                cluster_idx.size() > static_cast<size_t>(n / 1.2))
                continue;

            // Compute the circular centroid.
            double sin_sum = 0.0, cos_sum = 0.0;
            for (int pos : cluster_idx) {
                double angle = 2.0 * PI * pos / n;
                sin_sum += sin(angle);
                cos_sum += cos(angle);
            }
            double sin_mean = sin_sum / cluster_idx.size();
            double cos_mean = cos_sum / cluster_idx.size();
            double mean_angle = atan2(sin_mean, cos_mean);
            if (mean_angle < 0)
                mean_angle += 2.0 * PI;
            double centroid = (mean_angle / (2.0 * PI)) * n;

            // Record the cluster info.
            ClusterInfo info;
            info.domain_label = current_domain_label;
            info.centroid = centroid;
            info.domain_value = domain_value;
            info.indices = cluster_idx;
            info.length = static_cast<int>(cluster_idx.size());
            clusters_info.push_back(info);

            // Mark the corresponding positions in the domain label vector.
            for (int pos : cluster_idx) {
                domain_labels[pos] = current_domain_label;
            }
            current_domain_label++;
        }
    }

    // Replace any remaining -1 (noise or rejected clusters) with 0.
    for (int& lab : domain_labels) {
        if (lab < 0)
            lab = 0;
    }
    return make_pair(domain_labels, clusters_info);
}

//------------------------------------------------------------------------------
// Function: process_matrix_dbscan_periodic
//
// Processes an entire matrix (represented as vector of vector<double>) row–by–row
// using the periodic DBSCAN routine. Returns a pair:
//   - domain_label_matrix: an integer matrix (same shape as the input) with domain labels
//   - centroids_all: a vector (one per row) of vectors of ClusterInfo for each accepted cluster.
//------------------------------------------------------------------------------
pair<vector<vector<int>>, vector<vector<ClusterInfo>>>
process_matrix_dbscan_periodic(const vector<vector<double>>& reduced_matrix,
                               int eps = 5, int min_samples = 3, double threshold = 0.5) {
    int num_rows = static_cast<int>(reduced_matrix.size());
    if (num_rows == 0)
        return make_pair(vector<vector<int>>(), vector<vector<ClusterInfo>>());

    int num_cols = static_cast<int>(reduced_matrix[0].size());

    vector<vector<int>> domain_label_matrix(num_rows, vector<int>(num_cols, 0));
    vector<vector<ClusterInfo>> centroids_all(num_rows);

    for (int i = 0; i < num_rows; ++i) {
        // Process each row.
        auto result = process_row_dbscan_periodic(reduced_matrix[i], eps, min_samples, threshold);
        domain_label_matrix[i] = result.first;
        centroids_all[i] = result.second;
    }
    return make_pair(domain_label_matrix, centroids_all);
}

//------------------------------------------------------------------------------
// Function: build_centroid_matrix
//
// Builds a centroid matrix (one row per time step) that is zero everywhere except
// at the (rounded) centroid positions. At each centroid position (and within a window
// of ±RV), the corresponding domain label is placed.
//------------------------------------------------------------------------------
vector<vector<int>> build_centroid_matrix(const vector<vector<ClusterInfo>>& centroids_all,
                                            int num_cols, int RV) {
    int num_rows = static_cast<int>(centroids_all.size());
    vector<vector<int>> centroid_matrix(num_rows, vector<int>(num_cols, 0));

    for (int i = 0; i < num_rows; ++i) {
        for (const auto& info : centroids_all[i]) {
            // Round the centroid to the nearest column index and wrap using modulo.
            int col = static_cast<int>(round(info.centroid)) % num_cols;
            // Place the domain label in a window around col.
            for (int j = -RV; j < RV; ++j) {
                int pos = (col + j + num_cols) % num_cols;
                centroid_matrix[i][pos] = info.domain_label;
            }
        }
    }
    return centroid_matrix;
}


//------------------------------------------------------------------------------
// Function: write_domain_label_matrix
//
// Writes an integer matrix to a text file (one row per line, numbers separated by spaces).
//------------------------------------------------------------------------------
void write_domain_label_matrix(const vector<vector<int>>& matrix, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error opening output file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    for (const auto& row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i];
            if (i < row.size() - 1)
                outfile << " ";
        }
        outfile << "\n";
    }
    outfile.close();
}

//------------------------------------------------------------------------------
// Function: write_centroid_matrix
//
// Writes an integer matrix (the centroid matrix) to a text file.
//------------------------------------------------------------------------------
void write_centroid_matrix(const vector<vector<int>>& matrix, const string& filename) {
    // This function is essentially the same as write_domain_label_matrix.
    write_domain_label_matrix(matrix, filename);
}

// Function to read a matrix from a file
vector<vector<int>> readMatrixFromFile(const string &filename)
{
    ifstream file(filename);
    vector<vector<int>> matrix;
    if (!file)
    {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(file, line))
    {
        vector<int> row;
        for (char c : line)
        {
            if (c == '0' || c == '1')
            { // Assuming binary matrix
                row.push_back(c - '0');
            }
        }
        if (!row.empty())
        {
            matrix.push_back(row);
        }
    }
    file.close();
    return matrix;
}

vector<vector<double>> reduced_matrix(const vector<vector<int>> &matrix, unsigned rows, unsigned cols, 
    unsigned step, const string output_filename,bool debug)
{
    
    double step_f = static_cast<double>(step);
    unsigned reduced_rows = rows / step;
    vector<vector<double>> matrix_reduced(reduced_rows, vector<double>(cols, 0.0));

    // Compute reduced matrix
    for (unsigned i = 0; i < cols; i++)
    {
        for (unsigned k = 0; k < reduced_rows; k++)
        {
            for (unsigned j = 0; j < step; j++)
            {
                if (k * step + j < rows)
                { // Ensure we don't go out of bounds
                    matrix_reduced[k][i] += static_cast<double>(matrix[k * step + j][i]) / step_f;
                }
            }
        }
    }
    if (debug = true)
    {
        // Open output file
        ofstream outFile(output_filename);
        if (!outFile.is_open())
        {
            cerr << "Error: Could not create output file." << endl;
        }
        // Write transposed matrix directly by swapping loop order

        for (unsigned k = 0; k < reduced_rows; k++)
        {
            for (unsigned i = 0; i < cols; i++)
            {
                outFile << matrix_reduced[k][i] << " ";
            }
            outFile << endl;
        }

        outFile.close();
    }
        return matrix_reduced;
    }

//------------------------------------------------------------------------------
// Main function.
//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <input_file> <reduced_matrix_file> <domain_label_file> <centroid_matrix_file> <debug_s>" << endl;
        return EXIT_FAILURE;
    }
    string input_file = argv[1];
    string reduced_matrix_file = argv[2];
    string domain_label_file = argv[3];
    string centroid_matrix_file = argv[4];
    string debug_s=argv[5];
    bool debug;
    if(debug_s=="true"){
        debug=true;
    }
    else{
        debug=false;
    }
    
    const int step = 50;
    vector<vector<int>> matrix = readMatrixFromFile(input_file);

    // Get correct matrix dimensions
    unsigned rows = static_cast<unsigned>(matrix.size());
    unsigned cols = static_cast<unsigned>(matrix.empty() ? 0 : matrix[0].size());

    if (rows == 0 || cols == 0)
    {
        cerr << "Error: Matrix is empty or not formatted correctly." << endl;
        return EXIT_FAILURE;
    }

    vector<vector<double>> reduced_mat = reduced_matrix(matrix, rows, cols, step, reduced_matrix_file,debug);
    
    // Convert reduced_mat from float to double.
    vector<vector<double>> reduced_mat_double;
    reduced_mat_double.reserve(reduced_mat.size());
    for (const auto& row_f : reduced_mat) {
        vector<double> row_d;
        row_d.reserve(row_f.size());
        for (float value : row_f) {
            row_d.push_back(static_cast<double>(value));
        }
        reduced_mat_double.push_back(row_d);
    }
    
    int RV = 1;  // The window range around the centroid to mark.

    
    int num_rows = static_cast<int>(reduced_mat_double.size());
    int num_cols = static_cast<int>(reduced_mat_double[0].size());

    // Process the matrix row-by-row using periodic DBSCAN.
    // Adjust eps, min_samples, and threshold as needed.
    int eps = 3;
    int min_samples = 5;
    double threshold = 0.5;
    auto processed = process_matrix_dbscan_periodic(reduced_mat_double, eps, min_samples, threshold);
    vector<vector<int>> domain_label_matrix = processed.first;
    vector<vector<ClusterInfo>> centroids_all = processed.second;

    // Build the centroid matrix.
    vector<vector<int>> centroid_matrix = build_centroid_matrix(centroids_all, num_cols, RV);

    // Write the output matrices to files.
    if(debug=true){
        write_domain_label_matrix(domain_label_matrix, domain_label_file);
        write_centroid_matrix(centroid_matrix, centroid_matrix_file);
    }
    

    return 0;
}
