#!/usr/bin/env python3
import os
import re
import argparse
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from multiprocessing import Pool

# =========================================================
# Helpers
# =========================================================
def infer_freq_mag_columns(df):
    cols = [c.strip() for c in df.columns]
    low = [c.lower() for c in cols]
    freq_candidates = ['frequency','freq','f']
    mag_candidates = ['magnitude','mag','amplitude']
    freq_col = None
    mag_col = None

    for i,name in enumerate(low):
        if any(name == c for c in freq_candidates):
            freq_col = cols[i]
        if any(name == c for c in mag_candidates):
            mag_col = cols[i]

    if not freq_col:
        for i,name in enumerate(low):
            if 'freq' in name:
                freq_col = cols[i]
                break

    if not mag_col:
        for i,name in enumerate(low):
            if 'mag' in name:
                mag_col = cols[i]
                break

    return freq_col, mag_col


def load_fft_two_col(path):
    df = pd.read_csv(path, comment="#", engine="python")
    if df.shape[1] >= 2:
        freq_col, mag_col = infer_freq_mag_columns(df)
        if freq_col and mag_col:
            freqs = pd.to_numeric(df[freq_col], errors="coerce").to_numpy()
            mags  = pd.to_numeric(df[mag_col], errors="coerce").to_numpy()
            mask  = ~np.isnan(freqs)
            return freqs[mask], mags[mask]

        # fallback
        numeric_cols = []
        for c in df.columns:
            try:
                pd.to_numeric(df[c], errors="raise")
                numeric_cols.append(c)
            except:
                pass
        if len(numeric_cols) >= 2:
            freqs = pd.to_numeric(df[numeric_cols[0]], errors="coerce").to_numpy()
            mags  = pd.to_numeric(df[numeric_cols[1]], errors="coerce").to_numpy()
            mask  = ~np.isnan(freqs)
            return freqs[mask], mags[mask]

    raise RuntimeError("Could not parse Freq/Magnitude columns")


def parse_start_end_from_filename(fname):
    base = os.path.basename(fname)
    nums = re.findall(r"(\d{4,12})", base)
    if len(nums) < 2:
        raise ValueError(f"No start/end coords found in {base}")

    start = int(nums[-2])
    end   = int(nums[-1])
    if end < start:
        start, end = end, start
    return start, end


# =========================================================
# Peak extraction — now returns Top_Peak_Count
# =========================================================
def compute_top_peaks_from_fft(freqs, mags, N, top_n=3):
    half = N // 2
    use_len = min(half, len(mags))
    if use_len <= 1:
        return [], [], [], 0, 0

    mean_mag = np.mean(mags[:use_len])
    std_mag  = np.std(mags[:use_len])
    threshold = mean_mag + 2 * std_mag

    peaks, props = find_peaks(mags[:use_len], height=threshold)
    heights = props.get("peak_heights", np.array([]))

    raw_peak_count = len(peaks)

    if raw_peak_count == 0:
        return [], [], [], raw_peak_count, 0

    top_idx = np.argsort(heights)[::-1][:top_n]
    selected = [peaks[i] for i in top_idx]

    sel_f = []
    sel_m = []
    sel_p = []

    for idx in selected:
        if idx == 0 or idx >= len(freqs):
            continue
        f = float(freqs[idx])
        m = float(mags[idx])
        period = int(round(N / idx))
        sel_f.append(f)
        sel_m.append(m)
        sel_p.append(period)

    top_peak_count = len(sel_f)

    return sel_f, sel_m, sel_p, raw_peak_count, top_peak_count


# =========================================================
# Parallel worker
# =========================================================
def process_one_csv(args):
    csvf, top_n = args
    try:
        freqs, mags = load_fft_two_col(csvf)
    except Exception:
        return None

    try:
        start, end = parse_start_end_from_filename(csvf)
        N = end - start + 1
    except:
        N = len(freqs) * 2

    sel_freqs, sel_mags, sel_periods, raw_peaks, top_peaks = \
        compute_top_peaks_from_fft(freqs, mags, N, top_n)

    scfr_name = os.path.basename(csvf).replace(".csv", "")
    freqs_str   = ";".join(f"{v:.6f}" for v in sel_freqs)
    mags_str    = ";".join(f"{v:.3f}" for v in sel_mags)
    periods_str = ";".join(str(v) for v in sel_periods)

    return [
        scfr_name,
        raw_peaks,
        top_peaks,
        freqs_str,
        mags_str,
        periods_str
    ]


# =========================================================
# Main Walker
# =========================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species_folder")
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--cores", type=int, default=4)
    args = ap.parse_args()

    species = args.species_folder
    top_n   = args.top
    cores   = args.cores

    window_dirs = []
    for root, dirs, files in os.walk(species):
        if root == species:
            continue
        if root.count(os.sep) == species.count(os.sep) + 1:
            window_dirs.append(root)

    for wdir in window_dirs:
        print(f"\nProcessing window size folder: {wdir}")

        spectrum_dirs = []
        for root, dirs, files in os.walk(wdir):
            if root.endswith("spectrum_data"):
                spectrum_dirs.append(root)

        for spectrum_path in spectrum_dirs:
            chrom_dir = os.path.dirname(spectrum_path)
            summary_dir = os.path.join(chrom_dir, "chromosome_wise_summary")
            os.makedirs(summary_dir, exist_ok=True)

            summary_file = os.path.join(summary_dir, "summary.tsv")

            csv_files = sorted(
                [os.path.join(spectrum_path, f)
                 for f in os.listdir(spectrum_path)
                 if f.endswith(".csv")]
            )
            if not csv_files:
                continue

            work_items = [(csvf, top_n) for csvf in csv_files]
            rows = []

            with Pool(processes=cores) as pool:
                for result in tqdm(pool.imap_unordered(process_one_csv, work_items),
                                   total=len(work_items),
                                   desc=f"{os.path.basename(wdir)}",
                                   unit="csv"):
                    if result is not None:
                        rows.append(result)

            df = pd.DataFrame(rows, columns=[
                "SCFR_Name",
                "Num_Raw_Peaks",
                "Top_Peak_Count",
                "Frequencies",
                "Magnitudes",
                "Periods"
            ])
            df.to_csv(summary_file, sep="\t", index=False)
            print(f"  → Wrote {summary_file}")

    print("\nDone.")


if __name__ == "__main__":
    main()
