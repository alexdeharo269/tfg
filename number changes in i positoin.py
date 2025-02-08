import os 

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
    input_file = './ising_data_R_low_temp/ising_dataR40.dat'
    output_file = './ising_data_R_low_temp/numberchanges40.txt'
    
    matrix = read_matrix_from_file(input_file)
    changes = count_column_changes(matrix)
    write_changes_to_file(output_file, changes)

'''
def process_all_files_in_folder(folder_path):
        for i in range(1, 41):
            input_file = os.path.join(folder_path, f'ising_dataR{i}.dat')
            output_file = os.path.join(folder_path, f'numberchanges{i}.txt')
            
            matrix = read_matrix_from_file(input_file)
            changes = count_column_changes(matrix)
            write_changes_to_file(output_file, changes)

if __name__ == "__main__":
    folder_path = './ising_data_R_low_temp'
    process_all_files_in_folder(folder_path)'''