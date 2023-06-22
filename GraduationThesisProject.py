# region Working Directory
import os
# Get current working directory
cwd = os.getcwd()
print("Current working directory:", cwd)
# Change working directory
new_dir = "C:\DEPO\Bitirme Tezi\Yenitez\Allcode"
os.chdir(new_dir)
# Get new working directory
new_cwd = os.getcwd()
print("New working directory:", new_cwd)
# endregion

# region diff_finder
import pandas as pd

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Filter the data based on the GeneName
filtered_df = df[df['GeneName'] == 'rpoB']  # Change 'rpoB' to any desired GeneName or remove this line to compare all sequences

# Extract the DNA sequences
sequences = filtered_df['Pre-Sequence'].tolist()

# Compare the DNA sequences
for i in range(len(sequences)):
    seq1 = sequences[i]
    for j in range(i + 1, len(sequences)):
        seq2 = sequences[j]
        if seq1 == seq2:
            print(f"Sequences {i+1} and {j+1} are identical: {seq1}")
        else:
            print(f"Sequences {i+1} and {j+1} are different.")

# endregion

# region phylogenetic tree
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment

# Read the Excel file
data = pd.read_excel("sequences.xlsx")

# Filter the data based on gene names
gene_names = ["rpoB"]  # Specify the gene names you want to compare
filtered_data = data[data["GeneName"].isin(gene_names)]

# Create a list to store the DNA sequences
sequences = []
names = []

# Iterate over the filtered data and extract the sequences and names
for _, row in filtered_data.iterrows():
    sequence = row["Pre-Sequence"]
    seq_name = row["SpecName"]
    sequences.append(sequence)
    names.append(seq_name)

# Create a temporary FASTA file for sequence input
temp_fasta = "temp.fasta"

# Write the sequences to a temporary FASTA file
with open(temp_fasta, "w") as fasta_file:
    for i, sequence in enumerate(sequences):
        fasta_file.write(f">Sequence{i}\n")
        fasta_file.write(f"{sequence}\n")

# Read the sequences from the temporary FASTA file
records = SeqIO.parse(temp_fasta, "fasta")

# Convert records to MultipleSeqAlignment object
alignment = MultipleSeqAlignment(records)

# Calculate the distance matrix
calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(distance_matrix)

# Draw the phylogenetic tree
Phylo.draw(tree)

# Clean up the temporary FASTA file
import os
os.remove(temp_fasta)

# endregion

# region Nucleotide Frequency Line
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Filter rows based on GeneName (optional)
filtered_df = df[df['GeneName'] == 'rpoB']  # Change 'rpoB' to filter for specific genes

# Extract the DNA sequences
sequences = filtered_df['Pre-Sequence']
#sequences = df['Pre-Sequence']

# Get the length of the sequences
sequence_length = len(sequences.iloc[0])

# Calculate the nucleotide frequencies
nucleotide_counts = {
    'A': [0] * sequence_length,
    'C': [0] * sequence_length,
    'G': [0] * sequence_length,
    'T': [0] * sequence_length
}

for sequence in sequences:
    for i, nucleotide in enumerate(sequence):
        nucleotide_counts[nucleotide][i] += 1

total_sequences = len(sequences)

# Plot the nucleotide frequencies
positions = list(range(1, sequence_length + 1))
nucleotides = ['A', 'C', 'G', 'T']

plt.figure(figsize=(10, 6))
for nucleotide in nucleotides:
    plt.plot(positions, nucleotide_counts[nucleotide], label=nucleotide)

plt.xticks(positions)  # Set x-axis ticks to positions
plt.xlabel('Position')
plt.ylabel('Number of Sequences (Total: {})'.format(total_sequences))
plt.title('Nucleotide Frequency Analysis')
plt.legend()
plt.tight_layout()
plt.show()
# endregion

# region Nucleotide Frequency Bar
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Filter rows based on GeneName (optional)
filtered_df = df[df['GeneName'] == 'infB']  # Change 'rpoB' to filter for specific genes

# Extract the DNA sequences
sequences = filtered_df['Pre-Sequence']
#sequences = df['Pre-Sequence']

# Get the length of the sequences
sequence_length = len(sequences.iloc[0])

# Calculate the nucleotide frequencies
nucleotide_counts = {
    'A': [0] * sequence_length,
    'C': [0] * sequence_length,
    'G': [0] * sequence_length,
    'T': [0] * sequence_length
}

for sequence in sequences:
    for i, nucleotide in enumerate(sequence):
        nucleotide_counts[nucleotide][i] += 1

# Get the total number of sequences
total_sequences = len(sequences)

# Plot the nucleotide frequencies as bar plots
positions = np.arange(sequence_length) + 1
nucleotides = ['A', 'C', 'G', 'T']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

plt.figure(figsize=(10, 6))
bottom = np.zeros(sequence_length)
for i, nucleotide in enumerate(nucleotides):
    counts = nucleotide_counts[nucleotide]
    plt.bar(positions, counts, label=nucleotide, color=colors[i], bottom=bottom)
    bottom += counts

plt.xlabel('Position')
plt.ylabel('Number of Sequences (Total: {})'.format(total_sequences))
plt.title('Nucleotide Frequency Analysis')
plt.legend(loc="upper right")
plt.tight_layout()
plt.show()

# endregion

# region Fasta Downloader
import os
import pandas as pd
import requests

# Load the Excel file
df = pd.read_excel('sequences.xlsx')

# Create the main folder
main_folder = 'Fasta_files'
os.makedirs(main_folder, exist_ok=True)

# Iterate over the rows in the Excel file
for _, row in df.iterrows():
    gene_id = str(row['ID'])
    gene_name = str(row['GeneName'])
    url = 'https://www.example.com/api/get_fasta?id=' + gene_id  # Replace with the appropriate URL to fetch the FASTA file
    
    # Create the gene-specific subfolder
    gene_folder = os.path.join(main_folder, gene_name)
    os.makedirs(gene_folder, exist_ok=True)
    
    # Download the FASTA file
    response = requests.get(url)
    fasta_data = response.text
    
    # Save the FASTA file
    file_path = os.path.join(gene_folder, gene_id + '.fasta')
    with open(file_path, 'w') as f:
        f.write(fasta_data)

# endregion

# region Fasta sequences
import os
import pandas as pd

# Read the Excel file
df = pd.read_excel("sequences.xlsx")

# Create a main folder for fasta sequences
main_folder = "Fasta_sequences"
os.makedirs(main_folder, exist_ok=True)

# Iterate over the rows in the Excel file
for _, row in df.iterrows():
    gene_id = str(row['ID'])
    gene_name = str(row['GeneName'])
    spec_name = str(row['SpecName'])
    
    # Replace invalid characters in the subfolder name
    invalid_chars = [ '/', ':', '*', '?', '"', '<', '>', '|']
    for char in invalid_chars:
        gene_name = gene_name.replace(char, '_')
    
    # Create the gene-specific subfolder
    gene_folder = os.path.join(main_folder, gene_name)
    os.makedirs(gene_folder, exist_ok=True)
    
    # Create the file path for the empty text file
    file_name = f"{gene_id}_{spec_name}.txt"
    file_path = os.path.join(gene_folder, file_name)
    
    # Create an empty text file
    open(file_path, 'w').close()


# endregion

# region Gene Phylogenetic Tree

import os
# Get current working directory
cwd = os.getcwd()
print("Current working directory:", cwd)
# Change working directory
new_dir = "C:\DEPO\Bitirme Tezi\Yenitez\Allcode\Fasta_files"
os.chdir(new_dir)
# Get new working directory
new_cwd = os.getcwd()
print("New working directory:", new_cwd)

from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
import os

# Define the folder containing the fasta files
folder_path = "rpoB"

# Create a list to store the DNA sequences
sequences = []
names = []

# Iterate over the fasta files and extract the sequences and names
for file_name in os.listdir(folder_path):
    if file_name.endswith(".fasta"):
        file_path = os.path.join(folder_path, file_name)
        with open(file_path) as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append(record.seq)
                names.append(record.name)

# Convert sequences to MultipleSeqAlignment object
alignment = MultipleSeqAlignment(sequences)

# Calculate the distance matrix
calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.upgma(distance_matrix, names=names)  # Provide sequence names for labeling the tree

# Draw the phylogenetic tree
Phylo.draw(tree)




# endregion

# region Partial Sequence 
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Extract the pre-sequences column
pre_sequences = df['Pre-Sequence']

# Function to find the most seen sequence partial in a given length
def find_most_seen_partial(length):
    partials_count = {}
    for sequence in pre_sequences:
        for i in range(len(sequence) - length + 1):
            partial = sequence[i:i+length]
            partials_count[partial] = partials_count.get(partial, 0) + 1

    most_seen_partial = max(partials_count, key=partials_count.get)
    times_seen = partials_count[most_seen_partial]
    return most_seen_partial, times_seen

# Find the most seen sequence partials for different lengths
results = []
max_length = max(df['Pre-Length'])
min_length = 1
current_length = max_length

while current_length >= min_length:
    most_seen_partial, times_seen = find_most_seen_partial(current_length)
    results.append((current_length, most_seen_partial, times_seen))
    current_length -= 1

# Create a DataFrame from the results
result_df = pd.DataFrame(results, columns=['Length', 'RepeatingSequence', 'TimesSeen'])

# Display the result table
print(result_df)

plt.figure(figsize=(10, 6))
plt.plot(result_df['Length'], result_df['TimesSeen'], marker='o', linestyle='-', color='b')
plt.xlabel('Length')
plt.ylabel('Times Seen')
plt.title('Most Seen Sequence Partial')
plt.grid(True)
# Add labels on top of each data point
for x, y in zip(result_df['Length'], result_df['TimesSeen']):
    plt.text(x, y, str(y), ha='center', va='bottom')
plt.show()


# Plotting the graph
plt.figure(figsize=(10, 6))
plt.bar(result_df['Length'], result_df['TimesSeen'], color='b')
plt.xlabel('Length')
plt.ylabel('Times Seen')
plt.title('Most Seen Sequence Partial')
plt.grid(True)
plt.show()
# endregion

# region Partial Sequence Gene Specific
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Specify the GeneName for which you want to draw the graph
target_gene = 'rpoB'

# Filter the dataframe based on the GeneName
filtered_df = df[df['GeneName'] == target_gene]

# Extract the pre-sequences column
pre_sequences = filtered_df['Pre-Sequence']

# Function to find the most seen sequence partial in a given length
def find_most_seen_partial(length):
    partials_count = {}
    for sequence in pre_sequences:
        for i in range(len(sequence) - length + 1):
            partial = sequence[i:i+length]
            partials_count[partial] = partials_count.get(partial, 0) + 1

    most_seen_partial = max(partials_count, key=partials_count.get)
    times_seen = partials_count[most_seen_partial]
    return most_seen_partial, times_seen

# Find the most seen sequence partials for different lengths
results = []
max_length = max(filtered_df['Pre-Length'])
min_length = 1
current_length = max_length

while current_length >= min_length:
    most_seen_partial, times_seen = find_most_seen_partial(current_length)
    results.append((current_length, most_seen_partial, times_seen))
    current_length -= 1

# Create a DataFrame from the results
result_df = pd.DataFrame(results, columns=['Length', 'RepeatingSequence', 'TimesSeen'])

# Display the result table
print(result_df)

plt.figure(figsize=(10, 6))
plt.plot(result_df['Length'], result_df['TimesSeen'], marker='o', linestyle='-', color='b')
plt.xlabel('Length')
plt.ylabel('Times Seen')
plt.title(f"Most Seen Sequence Partial for Gene: {target_gene}")
plt.grid(True)
# Add labels on top of each data point
for x, y in zip(result_df['Length'], result_df['TimesSeen']):
    plt.text(x, y, str(y), ha='center', va='bottom')
plt.show()

# endregion

# region RBS Table
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel("sequences.xlsx")

# Define the RBS motifs to search for
rbs_motifs = ["AGGAGG", "GGAGG", "GGGAGG", "AGGA", "AGG"]
# Create an empty list to store the results
results = []

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Retrieve the necessary columns from the DataFrame
    id_val = row["ID"]
    gene_name = row["GeneName"]
    spec_name = row["SpecName"]
    sequence = row["Pre-Sequence"]
    pre_length = row["Pre-Length"]
    
    # Initialize variables to track the best RBS location
    best_rbs_location = -1
    best_rbs_motif = ""
    
    # Find the most suitable RBS location for each motif
    for motif in rbs_motifs:
        rbs_location = sequence.find(motif)
        
        if rbs_location != -1 and (best_rbs_location == -1 or rbs_location < best_rbs_location):
            best_rbs_location = rbs_location
            best_rbs_motif = motif
    
    # Append the results to the list
    results.append({
        "ID": id_val,
        "GeneName": gene_name,
        "SpecName": spec_name,
        "Pre-Sequence": sequence,
        "Pre-Length": pre_length,
        "RBS_Sequence": best_rbs_motif,
        "RBS_Location": best_rbs_location
    })

# Create a new DataFrame from the results list
results_df = pd.DataFrame(results, columns=["ID", "GeneName", "SpecName", "Pre-Sequence", "Pre-Length", "RBS_Sequence", "RBS_Location"])


# Print the table
print(results_df)

# Save the DataFrame to an Excel file
results_df.to_excel("RBS_results.xlsx", index=False)

# Count the frequency of each RBS motif
rbs_counts = results_df["RBS_Sequence"].value_counts()
# Plotting the graph
plt.bar(rbs_motifs, rbs_counts[rbs_motifs])
plt.xlabel("RBS Motif")
plt.ylabel("Frequency")
plt.title("Frequency of RBS Motifs")
plt.xticks(rotation=45)
# Add labels to the top of each bar
for i, count in enumerate(rbs_counts[rbs_motifs]):
    plt.text(i, count, str(count), ha='center', va='bottom')
plt.tight_layout()
plt.show()


# Count the frequency of each RBS motif
rbs_counts = results_df["RBS_Location"].value_counts()
# Plotting the graph
plt.bar(results_df["RBS_Location"], rbs_counts[results_df["RBS_Location"]])
plt.xlabel("RBS Motif Locations")
plt.ylabel("Frequency")
plt.title("Frequency of RBS Motifs Locations")
plt.xticks(rotation=45)
# Add labels on top of each bar
for x, y in zip(rbs_counts.index, rbs_counts):
    plt.text(x, y, str(y), ha='center', va='bottom')
plt.tight_layout()
plt.show()

# endregion

# region Selected RBS Table
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel("sequences.xlsx")

# Define the RBS motifs to search for
#rbs_motifs = ["AGGAGG", "GGAGG", "GGGAGG", "AGGA", "AGG"]
rbs_motifs = ["AAGGA","GGAA", "AATG", "AGGA", "GGAG", "TGAGG"]
# Create an empty list to store the results
results = []

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Retrieve the necessary columns from the DataFrame
    id_val = row["ID"]
    gene_name = row["GeneName"]
    spec_name = row["SpecName"]
    sequence = row["Pre-Sequence"]
    pre_length = row["Pre-Length"]
    
    # Initialize variables to track the best RBS location
    best_rbs_location = -1
    best_rbs_motif = ""
    
    # Find the most suitable RBS location for each motif
    for motif in rbs_motifs:
        rbs_location = sequence.find(motif)
        
        if rbs_location != -1 and (best_rbs_location == -1 or rbs_location < best_rbs_location):
            best_rbs_location = rbs_location
            best_rbs_motif = motif
    
    # Append the results to the list
    results.append({
        "ID": id_val,
        "GeneName": gene_name,
        "SpecName": spec_name,
        "Pre-Sequence": sequence,
        "Pre-Length": pre_length,
        "RBS_Sequence": best_rbs_motif,
        "RBS_Location": best_rbs_location
    })

# Create a new DataFrame from the results list
results_df = pd.DataFrame(results, columns=["ID", "GeneName", "SpecName", "Pre-Sequence", "Pre-Length", "RBS_Sequence", "RBS_Location"])


# Print the table
print(results_df)

# Save the DataFrame to an Excel file
results_df.to_excel("RBS_results.xlsx", index=False)

# Count the frequency of each RBS motif
rbs_counts = results_df["RBS_Sequence"].value_counts()
# Plotting the graph
plt.bar(rbs_motifs, rbs_counts[rbs_motifs])
plt.xlabel("RBS Motif")
plt.ylabel("Frequency")
plt.title("Frequency of Selected RBS Motifs")
plt.xticks(rotation=45)
# Add labels to the top of each bar
for i, count in enumerate(rbs_counts[rbs_motifs]):
    plt.text(i, count, str(count), ha='center', va='bottom')
plt.tight_layout()
plt.show()


# Count the frequency of each RBS motif
rbs_counts = results_df["RBS_Location"].value_counts()
# Plotting the graph
plt.bar(results_df["RBS_Location"], rbs_counts[results_df["RBS_Location"]])
plt.xlabel("RBS Motif Locations")
plt.ylabel("Frequency")
plt.title("Frequency of Selected RBS Motifs Locations")
plt.xticks(rotation=45)
# Add labels on top of each bar
for x, y in zip(rbs_counts.index, rbs_counts):
    plt.text(x, y, str(y), ha='center', va='bottom')
plt.tight_layout()
plt.show()

# endregion

# region most seen partial sequence graph and table
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Specify the GeneName for which you want to draw the graph
target_gene = 'infB'

# Filter the dataframe based on the GeneName
filtered_df = df[df['GeneName'] == target_gene]

# Extract the pre-sequences column
pre_sequences = filtered_df['Pre-Sequence']

# Function to find the most seen sequence partial in a given length
def find_most_seen_partial(length):
    partials_count = {}
    for sequence in pre_sequences:
        unique_partials = set()
        for i in range(len(sequence) - length + 1):
            partial = sequence[i:i+length]
            unique_partials.add(partial)
        
        for partial in unique_partials:
            partials_count[partial] = partials_count.get(partial, 0) + 1

    most_seen_partial = max(partials_count, key=partials_count.get)
    sequences_seen = partials_count[most_seen_partial]
    return most_seen_partial, sequences_seen

# Find the most seen sequence partials for different lengths
results = []
max_length = max(filtered_df['Pre-Length'])
min_length = 1
current_length = max_length

while current_length >= min_length:
    most_seen_partial, sequences_seen = find_most_seen_partial(current_length)
    results.append((current_length, most_seen_partial, sequences_seen))
    current_length -= 1

# Create a DataFrame from the results
result_df = pd.DataFrame(results, columns=['Length', 'RepeatingSequence', 'SequencesSeen'])

# Display the result table
print(result_df)

plt.figure(figsize=(10, 6))
plt.plot(result_df['Length'], result_df['SequencesSeen'], marker='o', linestyle='-', color='b')
plt.xlabel('Length')
plt.ylabel('Sequences Seen')
plt.title(f"Most Seen Sequence Partial for Gene: {target_gene}")
plt.grid(True)
# Add labels on top of each data point
for x, y in zip(result_df['Length'], result_df['SequencesSeen']):
    plt.text(x, y, str(y), ha='center', va='bottom')

# Display the table graph
fig, ax = plt.subplots(figsize=(10, 6))
ax.axis('off')  # Hide the axis
ax.table(cellText=result_df.values, colLabels=result_df.columns, loc='center')

plt.show()
# endregion

# region automated most seen partial sequence graph and table
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('sequences.xlsx')

# Get unique gene names
unique_genes = df['GeneName'].unique()

# Iterate over gene names
for target_gene in unique_genes:
    # Filter the dataframe based on the current gene name
    filtered_df = df[df['GeneName'] == target_gene]

    # Extract the pre-sequences column
    pre_sequences = filtered_df['Pre-Sequence']

    # Function to find the most seen sequence partial in a given length
    def find_most_seen_partial(length):
        partials_count = {}
        for sequence in pre_sequences:
            unique_partials = set()
            for i in range(len(sequence) - length + 1):
                partial = sequence[i:i+length]
                unique_partials.add(partial)

            for partial in unique_partials:
                partials_count[partial] = partials_count.get(partial, 0) + 1

        most_seen_partial = max(partials_count, key=partials_count.get)
        sequences_seen = partials_count[most_seen_partial]
        return most_seen_partial, sequences_seen

    # Find the most seen sequence partials for different lengths
    results = []
    max_length = max(filtered_df['Pre-Length'])
    min_length = 1
    current_length = max_length

    while current_length >= min_length:
        most_seen_partial, sequences_seen = find_most_seen_partial(current_length)
        results.append((current_length, most_seen_partial, sequences_seen))
        current_length -= 1

    # Create a DataFrame from the results
    result_df = pd.DataFrame(results, columns=['Length', 'RepeatingSequence', 'SequencesSeen'])

    # Add gene name column
    result_df['GeneName'] = target_gene

    # Display the result table
    print(f"Gene: {target_gene}")
    print(result_df)
    print()

    plt.figure(figsize=(10, 6))
    plt.plot(result_df['Length'], result_df['SequencesSeen'], marker='o', linestyle='-', color='b')
    plt.xlabel('Length')
    plt.ylabel('Sequences Seen')
    plt.title(f"Most Seen Sequence Partial for Gene: {target_gene}")
    plt.grid(True)
    # Add labels on top of each data point
    for x, y in zip(result_df['Length'], result_df['SequencesSeen']):
        plt.text(x, y, str(y), ha='center', va='bottom')

    # Display the table graph
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axis('off')  # Hide the axis
    ax.table(cellText=result_df.values, colLabels=result_df.columns, loc='center')

    plt.show()
# endregion

# region Escherichia coli most seen partial sequence graph and table
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('Escherichia_coli_sequences.xlsx')

# Extract the pre-sequences column
pre_sequences = df['Pre-Sequence']

# Function to find the most seen sequence partial in a given length
def find_most_seen_partial(length):
    partials_count = {}
    for sequence in pre_sequences:
        unique_partials = set()
        for i in range(len(sequence) - length + 1):
            partial = sequence[i:i+length]
            unique_partials.add(partial)
        
        for partial in unique_partials:
            partials_count[partial] = partials_count.get(partial, 0) + 1

    most_seen_partial = max(partials_count, key=partials_count.get)
    sequences_seen = partials_count[most_seen_partial]
    return most_seen_partial, sequences_seen

# Find the most seen sequence partials for different lengths
results = []
max_length = max(df['Pre-Length'])
min_length = 1
current_length = max_length

while current_length >= min_length:
    most_seen_partial, sequences_seen = find_most_seen_partial(current_length)
    results.append((current_length, most_seen_partial, sequences_seen))
    current_length -= 1

# Create a DataFrame from the results
result_df = pd.DataFrame(results, columns=['Length', 'RepeatingSequence', 'SequencesSeen'])

# Display the result table
print(result_df)

plt.figure(figsize=(10, 6))
plt.plot(result_df['Length'], result_df['SequencesSeen'], marker='o', linestyle='-', color='b')
plt.xlabel('Length')
plt.ylabel('Sequences Seen')
plt.title("Most Seen Sequence Partial for Escherichia coli")
plt.grid(True)
# Add labels on top of each data point
for x, y in zip(result_df['Length'], result_df['SequencesSeen']):
    plt.text(x, y, str(y), ha='center', va='bottom')

# Display the table graph
fig, ax = plt.subplots(figsize=(10, 6))
ax.axis('off')  # Hide the axis
ax.table(cellText=result_df.values, colLabels=result_df.columns, loc='center')

plt.show()
# endregion

# region Escherichia coli RBS Motifs
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('Escherichia_coli_sequences.xlsx')

# Define the RBS motifs to search for
#rbs_motifs = ["AGGAGG", "GGAGG", "GGGAGG", "AGGA", "AGG"]
rbs_motifs = ["AAGGA","GGAA", "AATG", "AGGA", "GGAG", "TGAGG"] # selected
# Create an empty list to store the results
results = []

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Retrieve the necessary columns from the DataFrame
    id_val = row["ID"]
    gene_name = row["GeneName"]
    spec_name = row["SpecName"]
    sequence = row["Pre-Sequence"]
    pre_length = row["Pre-Length"]
    
    # Initialize variables to track the best RBS location
    best_rbs_location = -1
    best_rbs_motif = ""
    
    # Find the most suitable RBS location for each motif
    for motif in rbs_motifs:
        rbs_location = sequence.find(motif)
        
        if rbs_location != -1 and (best_rbs_location == -1 or rbs_location < best_rbs_location):
            best_rbs_location = rbs_location
            best_rbs_motif = motif
    
    # Append the results to the list
    results.append({
        "ID": id_val,
        "GeneName": gene_name,
        "SpecName": spec_name,
        "Pre-Sequence": sequence,
        "Pre-Length": pre_length,
        "RBS_Sequence": best_rbs_motif,
        "RBS_Location": best_rbs_location
    })

# Create a new DataFrame from the results list
results_df = pd.DataFrame(results, columns=["ID", "GeneName", "SpecName", "Pre-Sequence", "Pre-Length", "RBS_Sequence", "RBS_Location"])


# Print the table
print(results_df)

# Save the DataFrame to an Excel file
results_df.to_excel("RBS_results.xlsx", index=False)

# Count the frequency of each RBS motif
rbs_counts = results_df["RBS_Sequence"].value_counts()
# Plotting the graph
plt.bar(rbs_motifs, rbs_counts[rbs_motifs])
plt.xlabel("RBS Motif")
plt.ylabel("Frequency")
plt.title("Frequency of Selected RBS Motifs - Escherichia coli")
plt.xticks(rotation=45)
# Add labels to the top of each bar
for i, count in enumerate(rbs_counts[rbs_motifs]):
    plt.text(i, count, str(count), ha='center', va='bottom')
plt.tight_layout()
plt.show()


# Count the frequency of each RBS motif
rbs_counts = results_df["RBS_Location"].value_counts()
# Plotting the graph
plt.bar(results_df["RBS_Location"], rbs_counts[results_df["RBS_Location"]])
plt.xlabel("RBS Motif Locations")
plt.ylabel("Frequency")
plt.title("Frequency of Selected RBS Motifs Locations - Escherichia coli")
plt.xticks(rotation=45)
# Add labels on top of each bar
for x, y in zip(rbs_counts.index, rbs_counts):
    plt.text(x, y, str(y), ha='center', va='bottom')
plt.tight_layout()
plt.show()

# endregion