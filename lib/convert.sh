#!/bin/bash

# Input file containing bash functions
INPUT_FILE=$1

# Output directory for separate bash scripts
OUTPUT_DIR=$2

# Check if input file and output directory are provided
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
  echo "Usage: $0 <input_file> <output_dir>"
  exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: Input file '$INPUT_FILE' not found."
  exit 1
fi

# Check if output directory exists, create it if not
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# Iterate over each function in the input file
while IFS= read -r line; do
  # Check if the line starts with a function declaration
  #if [[ $line =~ ^([a-zA-Z_][a-zA-Z_0-9]*)\(\) ]]; then
  if [[ $line =~ ^function\ ([a-zA-Z_][a-zA-Z_0-9]*) ]]; then
    func_name=${BASH_REMATCH[1]}
    output_file="$OUTPUT_DIR/$func_name.sh"
    echo "#!/bin/bash" > "$output_file"
    echo "$func_name() {" >> "$output_file"
    
    # Read the function body
    while IFS= read -r line; do
      # Check if the line ends the function body
      if [[ $line =~ \} ]]; then
        echo "}" >> "$output_file"
        break
      fi
      echo "$line" >> "$output_file"
    done
    # Make the output file executable
    chmod +x "$output_file"
  fi
done < "$INPUT_FILE"


cd $OUTPUT_DIR
for x in `ls`; do a=`echo $x | cut -d '.' -f1`; echo $a '$@' >> $x; mv $x $a; done
cd -

echo "Functions have been extracted to separate bash scripts in '$OUTPUT_DIR'."
