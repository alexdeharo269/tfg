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

const double PI = 3.14159265358979323846;

//---------------------------------------------------------------------
// Structure to hold cluster information.
//---------------------------------------------------------------------
struct ClusterInfo {
    int domain_label;         // Unique label for this accepted cluster.
    double centroid;          // Circular centroid computed from the indices.
    int domain_value;         // The binary value (0 or 1) corresponding to this cluster.
    vector<int> indices;      // The indices (positions) that belong to this cluster.
    int length;               // The size of the cluster.
};

//---------------------------------------------------------------------
// Function: readMatrixFromFile
//
// Reads a whitespace–separated integer matrix from a file.
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
// Function: readReducedMatrix
//
// Reads a whitespace–separated matrix of doubles from a file.
//---------------------------------------------------------------------
vector<vector<double>> readReducedMatrix(const string& filename) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening input file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    vector<vector<double>> matrix;
    string line;
    while (getline(infile, line)) {
        if (line.empty())
            continue;
        istringstream iss(line);
        vector<double> row;
        double value;
        while (iss >> value)
            row.push_back(value);
        if (!row.empty())
            matrix.push_back(row);
    }
    infile.close();
    return matrix;
}

//---------------------------------------------------------------------
// Function: replicateRows
//
// Replicates each row of the given matrix 'replication_factor' times.
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
// Function: cluster_positions_periodic
//
// A DBSCAN–like clustering of 1D positions (indices) on a circle of length n_total.
//---------------------------------------------------------------------
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

        vector<int> neighbors;
        for (size_t j = 0; j < L; ++j) {
            if (pdist(positions[i], positions[j]) <= eps)
                neighbors.push_back(j);
        }
        if (neighbors.size() < static_cast<size_t>(min_samples)) {
            visited[i] = true;
            continue;
        }
        labels[i] = cluster_id;
        visited[i] = true;

        unordered_set<int> seed_set(neighbors.begin(), neighbors.end());
        seed_set.erase(i);

        while (!seed_set.empty()) {
            int j = *seed_set.begin();
            seed_set.erase(seed_set.begin());
            if (!visited[j]) {
                visited[j] = true;
                vector<int> neighbors_j;
                for (size_t k = 0; k < L; ++k) {
                    if (pdist(positions[j], positions[k]) <= eps)
                        neighbors_j.push_back(k);
                }
                if (neighbors_j.size() >= static_cast<size_t>(min_samples)) {
                    for (int nb : neighbors_j)
                        seed_set.insert(nb);
                }
            }
            if (labels[j] == -1) {
                labels[j] = cluster_id;
            }
        }
        cluster_id++;
    }
    return labels;
}

//---------------------------------------------------------------------
// Function: process_row_dbscan_periodic
//
// Processes a single 1D row using periodic DBSCAN. Binarizes the row,
// clusters indices for each binary value, computes circular centroids,
// and returns a domain label vector and a vector of ClusterInfo.
//---------------------------------------------------------------------
pair<vector<int>, vector<ClusterInfo>>
process_row_dbscan_periodic(const vector<double>& row, int eps = 5, int min_samples = 3, double threshold = 0.5) {
    int n = static_cast<int>(row.size());
    vector<int> binary(n, 0);
    for (int i = 0; i < n; ++i) {
        binary[i] = (row[i] >= threshold) ? 1 : 0;
    }

    vector<int> domain_labels(n, -1);
    vector<ClusterInfo> clusters_info;
    int current_domain_label = 1;  // Labels start at 1; 0 means no accepted domain.

    for (int domain_value = 0; domain_value <= 1; ++domain_value) {
        vector<int> idx;
        for (int i = 0; i < n; ++i) {
            if (binary[i] == domain_value)
                idx.push_back(i);
        }
        if (idx.empty())
            continue;
        vector<int> cluster_labels = cluster_positions_periodic(idx, n, eps, min_samples);
        set<int> unique_labels;
        for (int lab : cluster_labels)
            unique_labels.insert(lab);
        for (int cl : unique_labels) {
            if (cl == -1)
                continue;
            vector<int> cluster_idx;
            for (size_t j = 0; j < cluster_labels.size(); ++j) {
                if (cluster_labels[j] == cl)
                    cluster_idx.push_back(idx[j]);
            }
            if (cluster_idx.size() < static_cast<size_t>(n / 10.0) ||
                cluster_idx.size() > static_cast<size_t>(n / 1.5))
                continue;

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

            ClusterInfo info;
            info.domain_label = current_domain_label;
            info.centroid = centroid;
            info.domain_value = domain_value;
            info.indices = cluster_idx;
            info.length = static_cast<int>(cluster_idx.size());
            clusters_info.push_back(info);

            for (int pos : cluster_idx) {
                domain_labels[pos] = current_domain_label;
            }
            current_domain_label++;
        }
    }
    for (int& lab : domain_labels) {
        if (lab < 0)
            lab = 0;
    }
    return make_pair(domain_labels, clusters_info);
}

//---------------------------------------------------------------------
// Function: process_matrix_dbscan_periodic
//
// Processes an entire matrix (row-by-row) using the periodic DBSCAN routine.
// Returns a pair: (domain_label_matrix, centroids_all)
//---------------------------------------------------------------------
pair<vector<vector<int>>, vector<vector<ClusterInfo>>>
process_matrix_dbscan_periodic(const vector<vector<double>>& reduced_matrix,
                               int eps = 5, int min_samples = 3, double threshold = 0.5) {
    int num_rows = static_cast<int>(reduced_matrix.size());
    vector<vector<int>> domain_label_matrix(num_rows, vector<int>());
    vector<vector<ClusterInfo>> centroids_all(num_rows);

    for (int i = 0; i < num_rows; ++i) {
        auto result = process_row_dbscan_periodic(reduced_matrix[i], eps, min_samples, threshold);
        domain_label_matrix[i] = result.first;
        centroids_all[i] = result.second;
    }
    return make_pair(domain_label_matrix, centroids_all);
}

//---------------------------------------------------------------------
// Function: build_centroid_matrix
//
// Builds a centroid matrix (one row per time step) that is zero everywhere
// except at positions around each cluster's centroid (within ±RV).
//---------------------------------------------------------------------
vector<vector<int>> build_centroid_matrix(const vector<vector<ClusterInfo>>& centroids_all,
                                            int num_cols, int RV) {
    int num_rows = static_cast<int>(centroids_all.size());
    vector<vector<int>> centroid_matrix(num_rows, vector<int>(num_cols, 0));

    for (int i = 0; i < num_rows; ++i) {
        for (const auto& info : centroids_all[i]) {
            int col = static_cast<int>(round(info.centroid)) % num_cols;
            for (int j = -RV; j < RV; ++j) {
                int pos = (col + j + num_cols) % num_cols;
                centroid_matrix[i][pos] = info.domain_label;
            }
        }
    }
    return centroid_matrix;
}

//---------------------------------------------------------------------
// Function: writeIntMatrixToFile
//
// Writes an integer matrix to a file (one row per line, space-separated).
//---------------------------------------------------------------------
void writeIntMatrixToFile(const vector<vector<int>>& matrix, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error opening output file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    for (const auto& row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i] << (i < row.size()-1 ? " " : "");
        }
        outfile << "\n";
    }
    outfile.close();
}

//---------------------------------------------------------------------
// Function: computeStrengthOfIncoherence
//
// Computes a strength of incoherence (SI) metric that focuses on intra-domain
// fluctuations. For each row of the original simulation (binary values assumed),
// groups the sites by domain label (from domain_label_matrix), computes the mode
// (using 0.5 threshold), and then determines the fraction of sites in each domain
// whose value deviates from the mode. The SI for the row is the average over all domains,
// and overall SI is the average over rows.
//---------------------------------------------------------------------
double computeStrengthOfIncoherence(const vector<vector<int>>& original_matrix,
                                    const vector<vector<int>>& domain_label_matrix,
                                    double tolerance = 0.0) {
    int num_rows = static_cast<int>(original_matrix.size());
    if (num_rows == 0)
        return 0.0;
    int num_cols = static_cast<int>(original_matrix[0].size());
    double totalSI = 0.0;
    int validRows = 0;

    for (int i = 0; i < num_rows; ++i) {
        unordered_map<int, vector<int>> domainValues;
        for (int j = 0; j < num_cols; ++j) {
            int label = domain_label_matrix[i][j];
            if (label > 0) { // Only consider accepted domains.
                domainValues[label].push_back(original_matrix[i][j]);
            }
        }
        if (domainValues.empty())
            continue;
        double rowSI = 0.0;
        int domainCount = 0;
        for (const auto& kv : domainValues) {
            const vector<int>& vals = kv.second;
            double sum = 0.0;
            for (int v : vals)
                sum += v;
            double avg = sum / vals.size();
            int mode = (avg >= 0.5) ? 1 : 0;
            int countDiff = 0;
            for (int v : vals) {
                if (abs(v - mode) > tolerance)
                    countDiff++;
            }
            double frac = static_cast<double>(countDiff) / vals.size();
            rowSI += frac;
            domainCount++;
        }
        if (domainCount > 0) {
            totalSI += rowSI / domainCount;
            validRows++;
        }
    }
    return (validRows > 0) ? totalSI / validRows : 0.0;
}

//---------------------------------------------------------------------
// Main function
//
// Usage (integrated version):
//   <program> <original_matrix_file> <reduced_matrix_file>
//              <domain_label_output_file> <centroid_matrix_output_file>
//              <si_output_file>
//---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <original_matrix_file> <reduced_matrix_file> "
             << "<domain_label_output_file> <centroid_matrix_output_file> <si_output_file>" << endl;
        return EXIT_FAILURE;
    }

    string original_file         = argv[1];
    string reduced_matrix_file   = argv[2];
    string domain_label_file     = argv[3];
    string centroid_matrix_file  = argv[4];
    string si_output_file        = argv[5];

    cout << "Reading original matrix from " << original_file << " ..." << endl;
    vector<vector<int>> original_matrix = readMatrixFromFile(original_file);
    if (original_matrix.empty()) {
        cerr << "Error: Original matrix is empty." << endl;
        return EXIT_FAILURE;
    }
    int full_rows = original_matrix.size();
    int full_cols = original_matrix[0].size();
    cout << "Original matrix: " << full_rows << " x " << full_cols << endl;

    cout << "Reading reduced matrix from " << reduced_matrix_file << " ..." << endl;
    vector<vector<double>> reduced_matrix = readReducedMatrix(reduced_matrix_file);
    if (reduced_matrix.empty()) {
        cerr << "Error: Reduced matrix is empty." << endl;
        return EXIT_FAILURE;
    }
    int reduced_rows = reduced_matrix.size();
    int reduced_cols = reduced_matrix[0].size();
    cout << "Reduced matrix: " << reduced_rows << " x " << reduced_cols << endl;

    // If the reduced matrix corresponds to a downsampled version (e.g., one row per 50 original rows),
    // then replicate the reduced matrix rows to match the full resolution.
    int replication_factor = full_rows / reduced_rows;
    vector<vector<int>> reduced_centroid_matrix = {}; // We'll obtain domain labels from DBSCAN below.
    // Process the reduced matrix row-by-row.
    auto processed = process_matrix_dbscan_periodic(reduced_matrix, 5, 10, 0.5);
    vector<vector<int>> reduced_domain_labels = processed.first;
    vector<vector<ClusterInfo>> centroids_all = processed.second;

    // Replicate the reduced domain labels to full resolution.
    vector<vector<int>> domain_label_matrix = replicateRows(reduced_domain_labels, replication_factor);
    if (domain_label_matrix.size() != original_matrix.size()) {
        cout << "Warning: Replicated domain label matrix rows (" << domain_label_matrix.size()
             << ") do not match original matrix rows (" << original_matrix.size() << ")." << endl;
    }

    // Build the centroid matrix at full resolution.
    vector<vector<int>> centroid_matrix = build_centroid_matrix(centroids_all, reduced_cols, 5);
    // Replicate centroid matrix rows to full resolution.
    vector<vector<int>> full_centroid_matrix = replicateRows(centroid_matrix, replication_factor);
    if (full_centroid_matrix.size() != original_matrix.size()) {
        cout << "Warning: Replicated centroid matrix rows (" << full_centroid_matrix.size()
             << ") do not match original matrix rows (" << original_matrix.size() << ")." << endl;
    }

    // Write the domain label matrix and centroid matrix to files.
    writeIntMatrixToFile(domain_label_matrix, domain_label_file);
    writeIntMatrixToFile(full_centroid_matrix, centroid_matrix_file);
    cout << "Domain labels written to " << domain_label_file << endl;
    cout << "Centroid matrix written to " << centroid_matrix_file << endl;

    // Compute the Strength of Incoherence (SI) metric.
    double SI = computeStrengthOfIncoherence(original_matrix, domain_label_matrix, 0.0);
    cout << "Strength of Incoherence (SI): " << SI << endl;

    // Write SI value to file.
    ofstream siFile(si_output_file);
    if (!siFile) {
        cerr << "Error opening SI output file: " << si_output_file << endl;
        return EXIT_FAILURE;
    }
    siFile << "Strength of Incoherence (SI): " << SI << "\n";
    siFile.close();
    cout << "SI written to " << si_output_file << endl;

    cout << "Processing complete." << endl;
    return 0;
}
