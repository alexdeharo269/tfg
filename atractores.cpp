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
vector<vector<int>> atractors(vector<vector<int>> &simulation, int min__width)
{
    vector<vector<int>> atractors;
    int time = simulation.size();
    int size = simulation[0].size();
    
    // Number of consecutive time steps to check for stability
    const int STABILITY_CHECK = 50;
    // Helper lambda to check if two attractors overlap in time and space
    auto attractors_overlap = [&](const vector<int>& a1, const vector<int>& a2) {
        // Check if they overlap in time
        if (abs(a1[0] - a2[0]) <= STABILITY_CHECK) {
            // Check if they overlap in space
            int center1 = a1[1];
            int center2 = a2[1];
            return abs(center1 - center2) < min__width/2;
        }
        return false;
    };
    
    for (size_t i = 0; i < time - STABILITY_CHECK; i++)
    {
        int sum = 0;
        for (size_t j = 0; j < size; j++)
        {
            if (simulation[i][j] == 1) {
                sum = 0;
            }
            else {
                sum++;
            }

            // Handle periodic boundary conditions
            if (j == size-1 && simulation[i][j] == 0) {
                int k = 0;
                while (k < size && simulation[i][k] == 0) {
                    sum++;
                    k++;
                }
            }

            if (sum > min__width && j + 1 < size && simulation[i][j+1] == 1) {
                int center = j - sum/2;
                if (center < 0) center += size;

                // Check stability: verify that the region remains zeros for next time steps
                bool is_stable = true;
                int start_pos = j - sum + 1;
                if (start_pos < 0) start_pos += size;

                // Check next STABILITY_CHECK time steps
                for (int t = 1; t <= STABILITY_CHECK && i + t < time; t++) {
                    for (int pos = 0; pos < sum; pos++) {
                        int check_pos = (start_pos + pos) % size;
                        if (simulation[i + t][check_pos] != 0) {
                            is_stable = false;
                            break;
                        }
                    }
                    if (!is_stable) break;
                }

                // Only add if stable and not overlapping with existing attractors
                if (is_stable) {
                    bool found_existing = false;
                    for (auto& existing : atractors) {
                        int existing_center = existing[1];
                        if (abs(existing_center - center) < min__width/2) {
                            if (sum > existing[2]) {
                                existing[2] = sum;
                            }
                            found_existing = true;
                            break;
                        }
                    }

                    if (!found_existing) {
                        vector<int> attractor = {(int)i, center, sum};
                        atractors.push_back(attractor);
                    }
                }
            }
        }
    }

    return atractors;
}

void write_atractors(const vector<vector<int>> &matrix, const string &filename)
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

//---------------------------------------------------------------------
// Function: readMatrixFromFile
//
// Reads a whitespace-separated matrix from a file. Each line is split
// into numbers (of type int). Returns a 2D vector representing the matrix.
//---------------------------------------------------------------------
vector<vector<int>> readMatrixFromFile(const string &input_file)
{
    ifstream file(input_file);
    if (!file)
    {
        cerr << "Error opening file: " << input_file << endl;
        exit(EXIT_FAILURE);
    }
    vector<vector<int>> matrix;
    string line;
    while (getline(file, line))
    {
        if (line.empty())
            continue;
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

// Read from temp file for pipe files based approach
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




int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " <original_file> <min_width> <output_file>" << endl;
        return EXIT_FAILURE;
    }

    string original_file = argv[1];
    string min_width = argv[2];
    string output_file = argv[3];
    
    
    // Read matrices.
    vector<vector<int>> original_matrix = readMatrixFromFile(original_file);
    vector<vector<int>> centroid_matrix_reduced;
    vector<vector<int>> output;
    
    
    output=atractors(original_matrix,stoi(min_width));
    write_atractors(output,output_file);
    

    return EXIT_SUCCESS;
}
