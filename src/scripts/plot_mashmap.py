import matplotlib.pyplot as plt
import pandas as pd

# Load data from Mashmap output file
data_file = '/mnt/gtklab01/xiaoqing/salmon/index/decoy2/mashmap.out'


# Define column names based on the provided format
columns = [
    'query_id', 'query_len', 'query_start', 'query_end', 
    'strand', 'ref_id', 'ref_len', 'ref_start', 'ref_end',
    'mapping_quality', 'identity', 'kc', 'fc'
]

# Read the file into a pandas DataFrame
df = pd.read_csv(data_file, sep='\\s+', names=columns)

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
# plt.show()

# Save the plot as a .jpg file
plt.savefig("mashmap_dotplot.jpg", format='jpg', dpi=300)