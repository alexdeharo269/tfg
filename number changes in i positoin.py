import os 
import sys

def read_matrix_from_file(input_file):
    with open(input_file, 'r') as file:
        matrix = [list(map(int, line.strip().split())) for line in file]
    return matrix

def count_column_changes(matrix):
    num_columns = len(matrix[0])
    changes = []
    for col in range(num_columns):
        count = 0
        for row in range(1, len(matrix)):
            if matrix[row][col] != matrix[row - 1][col]:
                count += 1
        changes.append((col, count))
    return changes

def write_changes_to_file(output_file, changes):
    with open(output_file, 'w') as file:
        for col, count in changes:
            file.write(f"{col} {count}\n")

if __name__ == "__main__":

    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    matrix = read_matrix_from_file(input_file)
    changes = count_column_changes(matrix)
    write_changes_to_file(output_file, changes)

