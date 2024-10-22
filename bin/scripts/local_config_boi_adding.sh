#!/bin/bash

# MN 22/10/24
# Check if a file is provided as an argument
# This code only runned by ML and YMM becuase the proc_dir permission
if [ -z "$1" ]; then
  echo "Usage: $0 <file>, file should includes the LiCS frame id:
173A_05151_131313
005D_04404_131313
034D_04318_131313
 "
  exit 1
fi



# Define a regular expression for the correct format (with 'A' or 'B' only for the fourth character)
regex="^[0-9]{3}[AD]_[0-9]{5}_[0-9]{6}$"

# Loop through the provided file (passed as $1)
while IFS= read -r i; do
  # Check if the line matches the correct format
  if [[ $i =~ $regex ]]; then
    # Remove leading zeros from the first three characters
    prefix=$(echo $i | cut -c 1-3 | sed 's/^0*//')
    path="$prefix/$i"

    # Check if local_config.py exists in the path
    if [ -f "$path/local_config.py" ]; then
      echo "local_config.py found in $path, updating..."
      if ! grep -q "bovl=1" "$path/local_config.py"; then
        echo "bovl=1" >> "$path/local_config.py"
        echo "Added bovl=1 to $path/local_config.py"
      else
        echo "bovl=1 already exists in $path/local_config.py"
      fi
    else
      echo "local_config.py not found in $path, creating..."
      echo "bovl=1" > "$path/local_config.py"
      echo "Created $path/local_config.py with bovl=1"
    fi
  else
    echo "$i # Incorrect format"
  fi
done < "$1"

