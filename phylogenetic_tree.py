# Creating Phylogenetic Tree

# Mandatory Requirements:
    # 1. Email id for NCBI queries
    # 2. Download and Install ClustalW2 based on your OS from the following website: http://www.clustal.org/download/current/
    # 3. Copy the path of the ClustalW2.exe from where it is installed on your system
    
# Optional Requirements (Only if not using the hardcoded information):
    # 1. Gene Name
    # 2. Set of Organisms names
    # 3. Output Sequence filename

# Execute the program with the following command: python phylogenetic_tree.py

from Bio import Entrez, SeqIO, Phylo
import subprocess
import os

# Step 1: Fetch sequences from NCBI Database
def fetch_sequences(gene_name, organism_list, email):
    """
    Fetch nucleotide sequences for a gene from NCBI for a list of organisms.
    """
    Entrez.email = email
    sequences = {}
    
    for organism in organism_list:
        print(f"Fetching sequence for {organism}...")
        query = f"{gene_name}[Gene] AND {organism}[Organism] AND biomol_genomic[PROP]"
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            sequence_id = record["IdList"][0]
            fetch_handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(fetch_handle, "fasta")
            fetch_handle.close()
            sequences[organism] = seq_record
        else:
            print(f"No sequence found for {organism}.")
    
    return sequences

# Step 2: Save fetched sequences to a FASTA file
def save_sequences_to_fasta(sequences, filename):
    """
    Save sequences to a FASTA file.
    """
    with open(filename, "w") as fasta_file:
        SeqIO.write(sequences.values(), fasta_file, "fasta")
    print(f"Sequences saved to {filename}.")

# Step 3: Align sequences using ClustalW2
def align_sequences(fasta_file, clustalw_path):
    """
    Align sequences using ClustalW2.
    """
    alignment_file = fasta_file.replace(".fasta", ".aln")
    print("Aligning sequences with ClustalW2...")

    # Build the ClustalW2 command
    cmd = [
        clustalw_path,
        f"-INFILE={fasta_file}",
        f"-OUTFILE={alignment_file}",
        f"-OUTPUT=CLUSTAL"
    ]

    # Run the ClustalW2 command
    process = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if process.returncode != 0:
        print(f"Error during alignment: {process.stderr}")
        return None

    print("Alignment completed.")
    return alignment_file

# Step 4: Build the phylogenetic tree
def build_phylogenetic_tree(alignment_file):
    """
    Build a phylogenetic tree from the alignment file.
    """
    tree_file = alignment_file.replace(".aln", ".dnd")
    if not os.path.exists(tree_file):
        print(f"Error: Guide tree file {tree_file} not found.")
        return None

    tree = Phylo.read(tree_file, "newick")
    print("Phylogenetic tree generated.")
    return tree

# Step 5: Visualize the tree
def visualize_tree(tree):
    """
    Visualize a phylogenetic tree.
    """
    print("Visualizing phylogenetic tree...")
    Phylo.draw(tree)

# Main program
if __name__ == "__main__":
    print("---------------------------------------------------------")
    print("\n*** Welcome to Phylogenetic Tree Visualization Tool ***\n")
    print("---------------------------------------------------------")
    # Mandatory Inputs
    email = input("Enter your email address for NCBI queries: ")
    clustalw_path = input("Enter the full path to the ClustalW2 executable: ")

    while True:
        print("\n*** Menu ***")
        print("1. Use hardcoded values")
        print("2. Input custom values")
        print("3. Exit")
        
        choice = input("Enter your choice: ")
        
        if choice == "1":
            # Hardcoded values
            gene_name = "COX1"
            organism_list = ["Homo sapiens", "Pan troglodytes", "Mus musculus", "Gallus gallus"]
            fasta_file = "sequences.fasta"
            
        elif choice == "2":
            # Custom input
            gene_name = input("Enter the gene name: ")
            organism_list = input("Enter the organism names (comma-separated): ").split(",")
            fasta_file = input("Enter the output FASTA file name (like sequences.fasta): ")
            
        elif choice == "3":
            print("Exiting program.")
            break
        
        else:
            print("Invalid choice. Please try again.")
            continue

        # Fetch sequences
        sequences = fetch_sequences(gene_name, organism_list, email)
        if not sequences:
            print("No sequences fetched. Returning to menu.")
            continue
        
        save_sequences_to_fasta(sequences, fasta_file)
        
        # Align sequences
        alignment_file = align_sequences(fasta_file, clustalw_path)
        if not alignment_file:
            print("Alignment failed. Returning to menu.")
            continue
        
        # Build and visualize the phylogenetic tree
        tree = build_phylogenetic_tree(alignment_file)
        if tree:
            visualize_tree(tree)
