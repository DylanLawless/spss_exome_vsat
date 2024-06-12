#!/bin/bash

# Generate the directory structure and save it to a file
tree > directory_structure.txt

# Initialize the documentation file with the directory structure
echo "Directory Structure:" > code_documentation.txt
cat directory_structure.txt >> code_documentation.txt
echo "" >> code_documentation.txt

# Use find to list all .R and .sh files recursively
# Loop through each file found by find
find . -type f \( -name "*.R" -o -name "*.sh" \) | while read file; do
    # Append the file name in a readable format to the documentation
    echo -e "\n\nFile: $file\n" >> code_documentation.txt
    # Append the actual file content to the documentation
    cat "$file" >> code_documentation.txt
    echo "" >> code_documentation.txt
done


rm directory_structure.txt
