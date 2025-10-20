#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Count threshold results by basis.")
    parser.add_argument("--csv", required=True, help="Path to common_thresholds_*.csv")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    if "basis" not in df.columns or "status" not in df.columns:
        raise SystemExit("CSV must contain 'basis' and 'status' columns.")

    df["basis"] = df["basis"].astype(str).str.upper()
    df["status"] = df["status"].astype(str)

    out = []
    for b in ["X", "Z"]:
        sub = df[df["basis"] == b]
        total = len(sub)
        ok = (sub["status"] == "ok").sum()
        no_thr = (sub["status"] == "no_common_threshold").sum()
        other = total - ok - no_thr
        out.append({
            "basis": b,
            "attempted": total,
            "with_threshold": ok,
            "no_threshold": no_thr,
            "other": other,
        })

    res = pd.DataFrame(out, columns=["basis", "attempted", "with_threshold", "no_threshold", "other"])
    print(res.to_string(index=False))

if __name__ == "__main__":
    main()