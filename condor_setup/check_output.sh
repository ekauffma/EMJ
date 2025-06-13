#!/bin/bash

# Loop from 3 to 728
for i in $(seq 3 728); do
    dir="output_inputFiles_250401.txt_$i"
    err_file="$dir/log/scriptcondor_.err"

    # Check if the error file exists
    if [[ -f "$err_file" ]]; then
        # Check if it contains the specific error line
        if grep -q "Run: \[ERROR\] Invalid operation: file exists" "$err_file"; then
            echo "Deleting directory: $dir (contains target error)"
            rm -rf "$dir"
        fi
    fi
done

