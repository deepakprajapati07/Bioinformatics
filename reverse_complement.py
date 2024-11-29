# Program to Reverse Complement a given DNA Sequence
# Run the program by the following command: python reverse_complement.py

import os

# Function to Reverse Complement a Sequence
def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement_dict = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    new_seq = ''.join(complement_dict.get(base.lower(), base) for base in seq)[::-1]
    return new_seq

# Function to Ensure Output Directory Exists
def ensure_output_directory(output_file_path):
    """Ensure the output directory exists; create it if it doesn't."""
    output_dir = os.path.dirname(output_file_path)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except Exception as e:
            raise IOError(f"Error creating output directory: {e}")

# Function to Read Sequence from File
def read_sequence(input_file_path):
    """Reads the DNA sequence from the input file."""
    try:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            description = lines[0].strip()
            sequence = ''.join(line.strip() for line in lines[1:])
        return description, sequence
    except FileNotFoundError:
        raise FileNotFoundError(f"The input file '{input_file_path}' does not exist.")
    except IOError as e:
        raise IOError(f"Error reading the input file: {e}")

# Function to Write Reverse Complement to Output File
def write_reverse_complement(output_file_path, description, reverse_complement_seq):
    """Writes the reverse complement sequence to the output file."""
    try:
        with open(output_file_path, 'w') as outfile:
            outfile.write(f"{description} ~ Reverse Complement \n")
            outfile.write(reverse_complement_seq)
        print(f"Reverse complement saved to: {output_file_path}")
    except IOError as e:
        raise IOError(f"Error writing to output file: {e}")

# Function to Process File: Reverse Complement and Save to Output
def process_file(input_file_path, output_file_path):
    """Main function to process the input file, generate reverse complement, and save it."""
    description, sequence = read_sequence(input_file_path)
    reverse_complement_seq = reverse_complement(sequence)
    ensure_output_directory(output_file_path)
    write_reverse_complement(output_file_path, description, reverse_complement_seq)

# Function to Prompt for File Paths with Validation
def get_valid_file_path(prompt):
    """Prompts the user to enter a valid file path and ensures the file exists."""
    while True:
        file_path = input(prompt)
        if os.path.isfile(file_path):
            return file_path
        else:
            print("Error: File not found. Please enter a valid file path.")

# Function to Handle Output File Path Logic
def handle_output_file_path(input_file_path, output_file_path):
    """Handles the logic for output file path. If only directory is given, generates a default file name."""
    if os.path.isdir(output_file_path):  # If only directory is specified
        # Extract the base file name (without extension) and create a new file name
        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        output_file_path = os.path.join(output_file_path, f"{base_name}_reverse_complement.wgs")
    return output_file_path

# Main Function to Handle User Interaction
def main():
    """Main user interface for interacting with the file processing program."""
    print("\nWelcome to the DNA Sequence Reverse Complement Program!\n")

    # Get input and output file paths from the user
    input_file = get_valid_file_path("Enter the input file path: ")
    output_file = input("Enter the output file path or directory: ")

    # Handle output file path logic
    output_file = handle_output_file_path(input_file, output_file)

    try:
        # Process the file (reverse complement and save it)
        process_file(input_file, output_file)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
