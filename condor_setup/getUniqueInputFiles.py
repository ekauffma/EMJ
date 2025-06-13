def filter_unique_lines(file1, file2, output_file):
    with open(file2, 'r', encoding='utf-8') as f2:
        file2_lines = set(line.rstrip() for line in f2)  # Store lines from file2 in a set

    with open(file1, 'r', encoding='utf-8') as f1, open(output_file, 'w', encoding='utf-8') as out:
        for line in f1:
            if line.rstrip() not in file2_lines:
                out.write(line)  # Write only lines not in file2

if __name__ == "__main__":
    file1 = "inputFiles_250402.txt"  # File to filter
    file2 = "inputFiles_250401.txt"  # Reference file
    output_file = "inputFiles_250402_temp.txt"

    filter_unique_lines(file1, file2, output_file)
