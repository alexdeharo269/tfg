#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

vector<vector<int>> readMatrixFromFile(const string &filename)
{
    ifstream in(filename);
    if (!in.is_open())
    {
        cerr << "Error: Could not open file '" << filename << "'. Please check if the file exists and you have read permissions." << endl;
        exit(EXIT_FAILURE);
    }

    vector<vector<int>> matrix;
    int value;
    string line;

    while (getline(in, line))
    {
        vector<int> row;
        istringstream iss(line);

        while (iss >> value)
        {
            if (value == 0 || value == 1)
            {
                row.push_back(value);
            }
        }

        if (!row.empty())
        {
            matrix.push_back(row);
        }
    }

    return matrix;
}

vector<vector<float>> reduced_matrix(const vector<vector<int>> &matrix, unsigned rows, unsigned cols, unsigned step)
{
    // Open output file
    ofstream outFile("./ising_data_R_low_temp/reduced_matrix.dat");
    if (!outFile.is_open())
    {
        cerr << "Error: Could not create output file." << endl;
    }

    unsigned reduced_rows = rows / step;
    vector<vector<float>> matrix_reduced(reduced_rows, vector<float>(cols, 0));

    // Compute reduced matrix
    for (unsigned i = 0; i < cols; i++)
    {
        for (unsigned k = 0; k < reduced_rows; k++)
        {
            for (unsigned j = 0; j < step; j++)
            {
                if (k * step + j < rows)
                { // Ensure we don't go out of bounds
                    matrix_reduced[k][i] += static_cast<float>(matrix[k * step + j][i]) / step;
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




void writeMatrixToFile(const vector<vector<int>> &matrix, const string &filename)
{
    ofstream out(filename);
    if (!out.is_open())
    {
        cerr << "Error: Could not create output file '" << filename << "'." << endl;
        exit(EXIT_FAILURE);
    }

    out.close();
}

void write_reduced(const vector<vector<float>> &matrix, const string &filename)
{
    ofstream out(filename);
    if (!out.is_open())
    {
        cerr << "Error: Could not create output file '" << filename << "'." << endl;
        exit(EXIT_FAILURE);
    }

    for (const auto &row : matrix)
    {
        for (int val : row)
        {
            out << val << " ";
        }
        out << endl;
    }

    out.close();
}

int main()
{
    const unsigned step=50;
    string filename = "./ising_data_R_low_temp/ising_dataR40.dat";
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

    vector<vector<float>> matrix_reduced =reduced_matrix(matrix, rows, cols, step);


    return 0;
}

