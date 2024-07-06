#!/bin/bash

# Get the current working directory
cwd=$(pwd)

# Change to the directory where this script is located
cd "$(dirname "$0")"

# Loop through all Python scripts in the current directory
for script in *.py; do
    # Execute each Python script
    python "$script"
done

# Change back to the original working directory
cd "$cwd"