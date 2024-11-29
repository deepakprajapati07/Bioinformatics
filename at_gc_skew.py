# Program to plot AT & GC Skew and Cumulative Sum
# Run the file using the following command: python at_gc_skew.py

import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Helper Functions

# Process nucleotide sequences from a file (can handle both single and multiple sequence formats)
def process_nucleotide_sequences(file_path, substring_length=120):
    nested_dict = {}
    
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Validate file structure
    if len(lines) < 2:
        raise ValueError("File format is incorrect. Expected a description and sequence.")
    
    # Parse the description (first line)
    description = lines[0].strip()[1:]

    # Process sequence (handle multiple lines or single sequence)
    sequence = "".join(line.strip().lower() for line in lines[1:]).replace(" ", "")
    substrings = [sequence[i:i + substring_length] for i in range(0, len(sequence), substring_length)]

    for line_number, substring in enumerate(substrings, start=1):
        d = Counter(substring)
        nested_dict[line_number] = {key: d.get(key, 0) for key in 'atgc'}

    df = pd.DataFrame.from_dict(nested_dict, orient='index').reset_index(drop=True)
    return description, df

# Add cumulative sum columns to the DataFrame
def add_cumulative_sum_columns(df):
    nucleotides = ['a', 't', 'g', 'c']
    for nucleotide in nucleotides:
        df[f'{nucleotide}_cum_sum'] = df[nucleotide].cumsum()
    return df

# Calculate AT and GC skew and add to DataFrame
def add_skew_columns(df):
    df['AT_skew'] = (df['a'] - df['t']) / (df['a'] + df['t'])
    df['GC_skew'] = (df['g'] - df['c']) / (df['g'] + df['c'])
    return df

# Add cumulative sum of AT and GC skew
def add_skew_cumulative_sums(df):
    df['AT_skew_cum_sum'] = df['AT_skew'].cumsum()
    df['GC_skew_cum_sum'] = df['GC_skew'].cumsum()
    return df

# Plot cumulative sum of nucleotides
def cumulative_sum_plot(df, subtitle='Cumulative Sum of Nucleotides', title='Cumulative Sum of Nucleotides', figsize=(10, 6)):
    cumulative_sum_columns = ['a_cum_sum', 't_cum_sum', 'g_cum_sum', 'c_cum_sum']
    plt.figure(figsize=figsize)
    for column in cumulative_sum_columns:
        plt.plot(df.index, df[column], label=column, marker='o')
    plt.title(title, fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel('Index')
    plt.ylabel('Cumulative Sum')
    plt.legend(title='Nucleotides')
    plt.grid(True)
    plt.show()

# Plot AT and GC skew
def AT_GC_Skew_plot(df, subtitle='Cumulative Skew of DNA Sequence', title='AT and GC Skew Analysis', figsize=(10, 6)):
    plt.figure(figsize=figsize)
    plt.plot(df.index, df['AT_skew_cum_sum'], label='AT Skew Cumulative Sum', marker='o')
    plt.plot(df.index, df['GC_skew_cum_sum'], label='GC Skew Cumulative Sum', marker='o')
    plt.title(title, fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel('Index')
    plt.ylabel('AT & GC Cumulative Sum')
    plt.legend(title='AT & GC Skew')
    plt.grid(True)
    plt.show()

# Visualize Cumulative Sum of Nucleotides
def visualize_cumulative_sum(file_path):
    description, df = process_nucleotide_sequences(file_path)
    df = add_cumulative_sum_columns(df)
    cumulative_sum_plot(df, subtitle=description)

# Visualize AT and GC Skew
def visualize_at_gc_skew(file_path):
    description, df = process_nucleotide_sequences(file_path)
    df = add_skew_columns(df)
    df = add_skew_cumulative_sums(df)
    AT_GC_Skew_plot(df, subtitle=description)

# Main Program Interface
def get_valid_file_path():
    while True:
        file_path = input("Enter file path of DNA sequence: ")
        if os.path.isfile(file_path):
            return file_path
        else:
            print("File not found. Please enter a valid file path.")

def menu():
    print("\n--- Menu ---")
    print("1. Visualize AT & GC Skew of a DNA Sequence")
    print("2. Visualize Cumulative Sum of Nucleotides in a DNA Sequence")
    print("3. Change file path")
    print("4. Exit")

# Main Program Loop
def main():
    print("Welcome to the DNA Sequence Analysis Program!")
    file_path = get_valid_file_path()

    while True:
        menu()
        choice = input("Choose an option (1-4): ")

        if choice == '4':
            print("\nExiting program. Goodbye!")
            break

        if choice == '3':
            print("\nYou can now enter a new file path.")
            file_path = get_valid_file_path()
            print(f"New file path selected: {file_path}")
        elif choice in ('1', '2'):
            try:
                if choice == '1':
                    visualize_at_gc_skew(file_path)
                elif choice == '2':
                    visualize_cumulative_sum(file_path)
            except ValueError as e:
                print(f"Error: {e}")
        else:
            print("Invalid choice. Please select a valid option.")

if __name__ == "__main__":
    main()
