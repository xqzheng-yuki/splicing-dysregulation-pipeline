import matplotlib.pyplot as plt
import pandas as pd

# Load data from Mashmap output file
cd $outfolder
data_file = 'mashmap.out'

# Define column names based on the provided format
columns = [
    'query_name', 'length', '0_based_start', '0_based_end', 
    'strand', 'chr_id', 'length', 'start', 'end', 
    'mapping_quality', 'identity', 'kc', 'fc'
]

# Read the file into a pandas DataFrame
df = pd.read_csv(data_file, sep='\s+', names=columns)

# Plotting
plt.figure(figsize=(10, 6))

# Create a scatter plot
plt.scatter(df['ref_start'], df['query_start'], 
            c=df['identity'], cmap='viridis', s=5)

# Adding labels and titles
plt.xlabel('Reference Start Position')
plt.ylabel('Query Start Position')
plt.title('Mashmap Genome Mapping Dot Plot')
plt.colorbar(label='Identity')

# Show the plot
plt.show()
