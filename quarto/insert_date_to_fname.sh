#!/bin/bash

# Print a message to the console
echo "Inserting date to filename..."

# Loop over all markdown files in the _posts directory
for file in ../_posts/*.md; do
    # Extract the base name of the file (i.e., the file name without the path)
    filename=$(basename -- "$file")

    # Check if the filename does not start with a date pattern like 2020-12-11
    if [[ ! $filename =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2} ]]; then
        # Search the file for a line starting with "date: " followed by a date
        date_line=$(grep -m 1 "^date: " "$file")

        # If such a line is found and it matches the pattern "date: YYYY-MM-DD"
        if [[ $date_line =~ ^date:\ ([0-9]{4}-[0-9]{2}-[0-9]{2}) ]]; then
            # Extract the date from the matched line
            date=${BASH_REMATCH[1]}

            # Print a message to the console showing the old and new file names
            echo "  $filename to $date-$filename"

            # Rename the file by prepending the date to the filename
            mv "$file" "../_posts/$date-$filename"
        fi
    fi
done
