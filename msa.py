# Multiple Sequence Alignment Tool using ClustalW2

# Requirements:
    # 1. Download and Install ClustalW2 based on your OS from the following website: http://www.clustal.org/download/current/
    # 2. Copy the path of the ClustalW2.exe from where it is installed on your system
    # 3. Be ready with two or more nucleotide sequences for performing alignment
    
    # Here are some sample sequences: ACGTTG, ACTTG, ACGTC
    
# Execute this file using the command: python msa.py

import subprocess
import tempfile
import os

def get_user_sequences():
    """
    Prompt the user to input sequences via the terminal.
    """
    print("\n-----Input Instructions-----")
    print("Enter your sequences one by one. Press Enter after each sequence.")
    print("When done, press Enter on an empty line to finish.\n")

    sequences = []
    while True:
        seq = input(f"Enter sequence {len(sequences) + 1} (or press Enter to finish): ").strip()
        if not seq:
            break
        sequences.append(seq)

    if len(sequences) < 2:
        print("At least two sequences are required for alignment.")
        return []

    print(f"\nYou entered {len(sequences)} sequences.")
    return sequences

def msa(sequences, clustalw_path="clustalw2"):
    """
    Align sequences using ClustalW2, print alignment and scores.
    """
    if not sequences:
        print("No sequences provided. Exiting.")
        return

    # Create a temporary file to hold the sequences in FASTA format
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as temp_fasta:
        for i, seq in enumerate(sequences):
            temp_fasta.write(f">seq{i+1}\n{seq}\n")
        temp_fasta_path = temp_fasta.name

    alignment_file_path = temp_fasta_path.replace(".fasta", ".aln")
    guide_tree_path = temp_fasta_path.replace(".fasta", ".dnd")

    print("Aligning sequences with ClustalW2...")

    # Build the ClustalW2 command
    cmd = [
        clustalw_path,
        f"-INFILE={temp_fasta_path}",
        f"-OUTFILE={alignment_file_path}",
        f"-OUTPUT=CLUSTAL"
    ]

    # Run the ClustalW2 command
    process = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if process.returncode != 0:
        print(f"Error: {process.stderr}")
        return

    # Extract alignment output from the alignment file
    if os.path.exists(alignment_file_path):
        with open(alignment_file_path, "r") as aln_file:
            alignment_content = aln_file.read()
        print(alignment_content)
    else:
        print("Error: Alignment file was not created.")

    # Extract scores from the process stdout
    diagnostic_output = process.stdout
    print("\nAlignment Scores and Information:")
    for line in diagnostic_output.splitlines():
        if "Aligned. Score:" in line or "Total Alignment Score" in line:
            print(line)

    # Clean up temporary files
    os.remove(temp_fasta_path)
    if os.path.exists(alignment_file_path):
        os.remove(alignment_file_path)
    if os.path.exists(guide_tree_path):
        os.remove(guide_tree_path)

if __name__ == "__main__":
    print("---------------------------------------------------------")
    print("\n*** Welcome to ClustalW2 Sequence Alignment Tool ***\n")
    print("---------------------------------------------------------")
    clustalw_path = input("Enter the path of ClustalW2 executable: ")
    sequences = get_user_sequences()
    msa(sequences, clustalw_path=clustalw_path)
