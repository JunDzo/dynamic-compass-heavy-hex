import pandas as pd

# Read the two CSV files
path1 = 'data/stats/stats_r_0.csv'
path2 = 'data/stats/stats_r_4.csv'
df1 = pd.read_csv(path1)
df2 = pd.read_csv(path2)
# df3 = pd.read_csv('figures/manuscript/no-reset-finer/stats.csv')

# Combine them
combined_df = pd.concat([df1, df2], ignore_index=True)
# combined_df = pd.concat([combined_df, df3], ignore_index=True)
# Save to a new CSV
combined_df.to_csv('data/stats/stats_r_04.csv', index=False)
print(f"Combined CSV saved as {path1}")
