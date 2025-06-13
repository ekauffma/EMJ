import sys

def concatenate_and_remove_duplicates(input_files, output_file):
    unique_lines = []
    seen = set()

    for file in input_files:
        with open(file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip()  # Remove trailing newline/whitespace
                if line not in seen:
                    seen.add(line)
                    unique_lines.append(line)

    with open(output_file, 'w', encoding='utf-8') as out:
        out.write("\n".join(unique_lines) + "\n")

if __name__ == "__main__":
    input_files = ["inputFiles5.txt", "inputFiles4.txt", "inputFiles3.txt", "inputFiles2.txt", "inputFiles.txt"]
    output_file = "inputFiles_250401.txt"
    concatenate_and_remove_duplicates(input_files, output_file)
    print(f"Unique lines have been written to {output_file}")

