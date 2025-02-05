#!/bin/bash

# Check if the filename is provided as an argument
if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Get the filename from the first argument
target_file="$1"

# Define the directory to search and the destination directory
search_dir="/mnt/gtklab01/xiaoqing/" # Base directory to start searching
dest_dir="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter" # Destination directory to store renamed files

# echo "Dry Run: Processing files named '$target_file'. No changes will be made."
echo "Processing files named '$target_file'. Files will be copied."


# Create the destination directory if it doesn't exist (for simulation purposes)
mkdir -p "$dest_dir"

# Find all files named "unmapped_d.lst" within directories containing '2025-01-14' or '2025-01-14-mm'
find "$search_dir" -type f -path "*/2025-01-14*" -name "$target_file" | while read -r file; do
    # Extract the subdirectory name from the path
    dir_path=$(dirname "$file")
    ctx=$(echo "$dir_path" | grep -o 'CTX_[0-9]\+')
    
    # Check if the path includes "-mm" and set the appropriate prefix
    if [[ "$dir_path" == *"-mm"* ]]; then
        prefix="${ctx}_mm"
    else
        prefix="${ctx}"
    fi
    
    # Define the new file name with the prefix
    new_file_name="${prefix}-${target_file}"
    
    # # Simulate the move and rename action
    # echo "Would copy or move: $file -> $dest_dir/$new_file_name"

    # Copy the file to the destination directory
    cp "$file" "$dest_dir/$new_file_name"
    echo "Copied: $file -> $dest_dir/$new_file_name"
done

# echo "Dry run complete. No files were modified."
echo "File copy process complete."
