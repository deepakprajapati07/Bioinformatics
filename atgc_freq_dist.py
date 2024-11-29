# Visualization of Cumulative Sum of Nucleotides & Various Frequency Distribution of Nucleotides in a DNA Sequence

# Run the program using the command : python atgc_freq_dist.py

import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ThreadPoolExecutor

# Global cache to store processed data
processed_data = None

# Helper Functions

# Function to count nucleotides in a substring
def count_nucleotides(substring):
    return Counter(substring)

# Function to process nucleotide sequences
def process_file_once(file_path, substring_length=120):
    """
    Reads and processes the DNA sequence file.
    Uses cached processed data if the file was already processed.
    """
    global processed_data
    if processed_data is None or processed_data[0] != file_path:
        nested_dict = {}
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            description = first_line[1:] if first_line.startswith(">") else ""
            sequence = "".join(line.strip().lower() for line in file)

        if not description or not sequence:
            raise Exception("Invalid file format or empty sequence")

        substrings = [sequence[i:i + substring_length] for i in range(0, len(sequence), substring_length)]
        with ThreadPoolExecutor() as executor:
            counts = list(executor.map(count_nucleotides, substrings))

        nested_dict = {i + 1: {key: counts[i].get(key, 0) for key in 'atgc'} for i in range(len(counts))}
        df = pd.DataFrame.from_dict(nested_dict, orient='index')
        df = add_cumulative_sum_columns(df)
        processed_data = (file_path, description, df)

    return processed_data[1], processed_data[2]

# Add cumulative sum columns

def add_cumulative_sum_columns(df):
    """
    Adds cumulative sum columns for each nucleotide in the DataFrame.
    """
    for nucleotide in 'atgc':
        if nucleotide in df.columns:
            df[f'{nucleotide}_cum_sum'] = df[nucleotide].cumsum()
    return df

# Plotting Functions

def plot_nucleotide_data(df, columns, kind='line', title='', subtitle='', xlabel='', ylabel='', figsize=(10, 6)):
    """
    Generic function to plot nucleotide data. , marker='o', label=column
    """
    df[columns].plot(kind=kind, figsize=figsize, marker='o')
    plt.title(title, fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(title="Nucleotide")
    plt.grid(True)
    plt.show()

def cumulative_sum_plot(df, subtitle='Cumulative Sum of Nucleotides in DNA Sequence'):
    """
    Plots the cumulative sum of nucleotides.
    """
    columns = ['a_cum_sum', 't_cum_sum', 'g_cum_sum', 'c_cum_sum']
    plot_nucleotide_data(
        df, columns,
        title='Cumulative Sum of Nucleotides',
        subtitle=subtitle,
        xlabel='Index',
        ylabel='Cumulative Sum'
    )

def nucleotides_frequencies(df, subtitle='Frequency of Nucleotides in DNA Sequence'):
    """
    Plots the frequency of nucleotides.
    """
    columns = ['a', 't', 'g', 'c']
    plot_nucleotide_data(
        df, columns,
        title='Nucleotides Frequency Distribution',
        subtitle=subtitle,
        xlabel='Index',
        ylabel='Frequency'
    )

def individual_nucleotide_frequencies(df, subtitle='Frequency of Nucleotides in DNA Sequence'):
    """
    Plots individual nucleotide frequencies in subplots.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    nucleotides = ['a', 't', 'g', 'c']
    axes = axes.flatten()

    for i, nucleotide in enumerate(nucleotides):
        axes[i].plot(df.index, df[nucleotide], marker='o', linestyle='-', label=f'{nucleotide.upper()} Frequency')
        axes[i].set_title(f'{nucleotide.upper()} Frequency Distribution')
        axes[i].set_xlabel('Index')
        axes[i].set_ylabel('Frequency')
        axes[i].grid(True)
        axes[i].legend()

    plt.tight_layout()
    plt.show()

def nucleotides_boxplot(df, subtitle):
    """
    Plots a boxplot of nucleotide frequencies.
    """
    df[['a', 't', 'g', 'c']].boxplot()
    plt.title('Nucleotides Frequency Distribution (Box Plot)', fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel("Nucleotide")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.show()

def nucleotides_histograms(df, subtitle):
    """
    Plots histograms of nucleotide frequencies.
    """
    df[['a', 't', 'g', 'c']].plot(kind='hist', bins=15, alpha=0.7, figsize=(10, 6), stacked=True)
    plt.title('Nucleotides Frequency Distribution (Histogram)', fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel("Frequency")
    plt.ylabel("Count")
    plt.legend(title="Nucleotide")
    plt.grid(True)
    plt.show()

def nucleotides_violin(df, subtitle='Frequency of Nucleotides in DNA Sequence', figsize=(10, 6)):
    """
    Plots a violin plot of nucleotide frequencies, excluding cumulative sum columns.
    """
    nucleotide_columns = ['a', 't', 'g', 'c']
    df_melted = df[nucleotide_columns].melt(var_name="Nucleotide", value_name="Frequency")
    plt.figure(figsize=figsize)
    sns.violinplot(x="Nucleotide", y="Frequency", data=df_melted, inner="quartile")
    plt.title('Nucleotides Frequency Distribution (Violin Plot)', fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel("Nucleotide")
    plt.ylabel("Frequency")
    plt.show()

def rolling_average(df, subtitle, window_size=10):
    """
    Plots a rolling average of nucleotide frequencies.
    """
    for nucleotide in 'atgc':
        df[nucleotide].rolling(window=window_size).mean().plot(label=nucleotide)
    plt.title('Rolling Average of Nucleotides Frequency', fontsize=14, y=1.08)
    plt.suptitle(subtitle, fontsize=10, y=0.92)
    plt.xlabel("Index")
    plt.ylabel("Rolling Average Frequency")
    plt.legend(title="Nucleotide")
    plt.grid(True)
    plt.show()

# Menu Functions

def cumulative_sum_visualization(file_path):
    description, df = process_file_once(file_path)
    cumulative_sum_plot(df, subtitle=description)

def all_nucleotides_frequency_visualization(file_path):
    description, df = process_file_once(file_path)
    nucleotides_frequencies(df, subtitle=description)

def individual_nucleotides_frequency_visualization(file_path):
    description, df = process_file_once(file_path)
    individual_nucleotide_frequencies(df, subtitle=description)

def nucleotides_boxplot_visualization(file_path):
    description, df = process_file_once(file_path)
    nucleotides_boxplot(df, subtitle=description)

def nucleotides_histogram_visualization(file_path):
    description, df = process_file_once(file_path)
    nucleotides_histograms(df, subtitle=description)

def nucleotides_violinplot_visualization(file_path):
    description, df = process_file_once(file_path)
    nucleotides_violin(df, subtitle=description)

def nucleotides_rolling_avg_visualization(file_path):
    description, df = process_file_once(file_path)
    rolling_average(df, subtitle=description)

# Menu and Main Loop

menu_actions = {
    '1': cumulative_sum_visualization,
    '2': all_nucleotides_frequency_visualization,
    '3': individual_nucleotides_frequency_visualization,
    '4': nucleotides_boxplot_visualization,
    '5': nucleotides_histogram_visualization,
    '6': nucleotides_violinplot_visualization,
    '7': nucleotides_rolling_avg_visualization,
}

def get_valid_file_path():
    """
    Prompts the user for a valid file path.
    """
    while True:
        file_path = input("Enter file path of DNA Sequence: ")
        if os.path.isfile(file_path):
            return file_path
        else:
            print("Error: File not found. Please enter a valid file path.\n")

def menu():
    """
    Displays the menu options.
    """
    print("\n--- Menu ---")
    for key, action in menu_actions.items():
        print(f"{key}. {action.__name__.replace('_', ' ').capitalize()}")
    print("8. Change file path")
    print("9. Exit")

def main():
    """
    Main function to handle user input and program flow.
    """
    print("Welcome to the DNA Sequence Analysis Program!")
    file_path = get_valid_file_path()

    while True:
        menu()
        choice = input("Choose an option (1-9): ")
        
        if choice == '9':
            print("\nExiting program. Goodbye!")
            break

        if choice == '8':
            print("\nYou can now enter a new file path.")
            file_path = get_valid_file_path()
            print(f"New file path selected: {file_path}")
        elif choice in menu_actions:
            try:
                menu_actions[choice](file_path)
            except Exception as e:
                print(f"Error: {e}")
        else:
            print("Invalid choice. Please select a valid option.")

if __name__ == "__main__":
    main()
