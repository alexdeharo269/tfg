#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

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

vector<vector<float>> reduced_matrix(const vector<vector<int>> &matrix, unsigned rows, unsigned cols, unsigned step, const string output_filename)
{
    // Open output file
    ofstream outFile(output_filename);
    if (!outFile.is_open())
    {
        cerr << "Error: Could not create output file." << endl;
    }
    float step_f = static_cast<float>(step);
    unsigned reduced_rows = rows / step;
    vector<vector<float>> matrix_reduced(reduced_rows, vector<float>(cols, 0.0f));

    // Compute reduced matrix
    for (unsigned i = 0; i < cols; i++)
    {
        for (unsigned k = 0; k < reduced_rows; k++)
        {
            for (unsigned j = 0; j < step; j++)
            {
                if (k * step + j < rows)
                { // Ensure we don't go out of bounds
                    matrix_reduced[k][i] += static_cast<float>(matrix[k * step + j][i]) / step_f;
                }
            }
        }
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
    return matrix_reduced;
}



int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return EXIT_FAILURE;
    }

    string filename = argv[1];
    string output_filename = argv[2];

    const int step = 50;
    vector<vector<int>> matrix = readMatrixFromFile(filename);

    // Get correct matrix dimensions
    unsigned rows = static_cast<unsigned>(matrix.size());
    unsigned cols = static_cast<unsigned>(matrix.empty() ? 0 : matrix[0].size());
    cout << "Time: " << rows << ", Size: " << cols << endl;

    if (rows == 0 || cols == 0)
    {
        cerr << "Error: Matrix is empty or not formatted correctly." << endl;
        return EXIT_FAILURE;
    }

    vector<vector<float>> matrix_reduced = reduced_matrix(matrix, rows, cols, step,output_filename);

    

    return 0;
}
