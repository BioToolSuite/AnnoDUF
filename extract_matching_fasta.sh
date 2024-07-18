# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fasta_file> <search_string> <output_folder>"
    exit 1
fi

input_file=$1
search_string=$2
output_folder=$3

# Create the output folder if it doesn't exist
mkdir -p "DUF/$output_folder"

# Initialize variables
IFS='.' read -ra parts <<< "$search_string"
newsearch_string="${parts[0]}"
output_file="DUF/$output_folder/$newsearch_string.fasta"
match=false

# Create or clear the output file
> $output_file

# Read the input FASTA file
while IFS= read -r line; do
    if [[ $line == \>* ]]; then
        # If the line starts with '>', it's a header line
        if [[ $line == *"$search_string"* ]]; then
            # If the header contains the search string, mark it as a match
            match=true
            echo $line >> $output_file
        else
            match=false
        fi
    elif $match; then
        # If it's a matching sequence, write the sequence line to the output file
        echo $line >> $output_file
    fi
done < $input_file