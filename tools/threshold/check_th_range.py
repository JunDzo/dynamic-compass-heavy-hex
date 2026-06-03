from pathlib import Path
import numpy as np
import pandas as pd

in_path = Path("data/approx_thresholds.csv")
df = pd.read_csv(in_path)
th_hi = df["p_high"].to_numpy(float)
th_lo = df["p_low"].to_numpy(float)
th_hi = np.nan_to_num(th_hi, nan=0.0)
th_lo = np.nan_to_num(th_lo, nan=0.0)

range_th = abs(th_hi - th_lo)
range_th_max = np.max(range_th)
print(f"Maximum threshold range across approximate thresholds: {range_th_max}")