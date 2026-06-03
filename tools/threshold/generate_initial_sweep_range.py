import pandas as pd
import numpy as np

# Parameters to generate sweep_range.csv
# b1p (x): from 0 to 0.004, 20 points
x_values = np.linspace(0, 0.004, 10)

# b2p (y): from 0 to 0.1, 20 points
y_values = np.linspace(0, 0.1, 10)

# z: from 0 to 0.003, 20 points (same list for all rows)
z_values = np.linspace(0, 0.003, 10)
z_str = ' '.join(f'{z:.6f}' for z in z_values)

# Basis: X and Z
basis_values = ['X', 'Z']

# Fixed values
diameters = '17 23 29 35'
base_prob = 0.0002

rows = []
for x in x_values:
    for y in y_values:
        for basis in basis_values:
            rows.append({
                'x': f'{x:.6f}',
                'y': f'{y:.6f}',
                'z': z_str,
                'basis': basis,
                'diameters': diameters,
                'base_prob': f'{base_prob:.4f}'
            })

df = pd.DataFrame(rows)
df.to_csv('data/analysis/swp_sug/sweep_range.csv', index=False)
print("Generated sweep_range.csv")