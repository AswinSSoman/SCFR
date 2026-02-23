import os
import re
import argparse
import numpy as np
from collections import defaultdict, Counter
from scipy.fft import fft
from Bio import SeqIO
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# ----------------------------
# Original functions unchanged
# ----------------------------

def safe_filename(name, max_length=100):
    return re.sub(r'[^A-Za-z0-9_]+', '_', name)[:max_length]

def escape_latex(s):
    specials = {
        '\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#',
        '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\^{}', '|': r'\textbar{}'
    }
    for char, replacement in specials.items():
        s = s.replace(char, replacement)
    return s

def parse_header(header):
    """
    Parse two possible FASTA header types:
    
    1. Original CDS headers:
       >lcl|NC_060925.1_cds_NP_001005277.1_1 [gene=OR4F16] ...

    2. SCFR headers:
       >-2::NC_060925.1:2194196-2200379
    """

    # ===============================
    # Case 1 → SCFR HEADER
    # ===============================
    scfr = re.match(r'([+-]?\d+)::([A-Za-z0-9_.]+):(\d+)-(\d+)', header)
    if scfr:
        frame = scfr.group(1)
        acc = scfr.group(2)
        start = scfr.group(3)
        end = scfr.group(4)

        return {
            "gene": f"SCFR_{acc}_{start}_{end}",
            "gene_id": "SCFR",
            "prot_id": f"frame{frame}",
            "prot_name": f"Stop-codon-free region ({start}-{end})",
            "acc": acc
        }

    # ===============================
    # Case 2 → Original CDS HEADER
    # ===============================
    gene = re.search(r'\[gene=(.*?)\]', header)
    gene_id = re.search(r'\[db_xref=GeneID:(\d+)\]', header)
    protein_id = re.search(r'\[protein_id=(.*?)\]', header)
    protein_name = re.search(r'\[protein=(.*?)\]', header)
    accession = re.search(r'lcl\|([A-Z0-9_.]+)', header)

    return {
        "gene": gene.group(1) if gene else "Unknown",
        "gene_id": gene_id.group(1) if gene_id else "Unknown",
        "prot_id": protein_id.group(1) if protein_id else "Unknown",
        "prot_name": protein_name.group(1) if protein_name else "Unknown",
        "acc": accession.group(1) if accession else "Unknown"
    }

def encode_sequence(seq):
    seq = seq.upper()
    return (
        np.array([1 if b == 'A' else 0 for b in seq]),
        np.array([1 if b == 'C' else 0 for b in seq]),
        np.array([1 if b == 'G' else 0 for b in seq]),
        np.array([1 if b == 'T' else 0 for b in seq])
    )

def compute_fft_magnitude(signal):
    return np.abs(fft(signal))

def extract_motifs(seq, period, top_n=5):
    motifs = [seq[i:i+period] for i in range(0, len(seq)-period+1)]
    return Counter(motifs).most_common(top_n)

def classify_peaks(freqs):
    if len(freqs) == 1 and abs(freqs[0] - 0.33) < 0.02:
        return 1
    elif len(freqs) == 1:
        return 2
    else:
        return 3

# -------------------------------------------
# Modified: analyze_sequence is multiprocessing-friendly
# -------------------------------------------
def analyze_sequence_wrapper(args):
    record_dict, output_dir, spectrum_dir = args
    # Rebuild SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    record = SeqRecord(Seq(record_dict["seq"]), id=record_dict["id"], description=record_dict["description"])
    return analyze_sequence(record, output_dir, spectrum_dir)

def analyze_sequence(record, output_dir, spectrum_dir):
    try:
        header = record.description
        meta = parse_header(header)
        seq = str(record.seq).upper()
        N = len(seq)
        meta["length"] = N
        if N < 20:
            raise ValueError("Sequence too short for FFT")

        A, C, G, T = encode_sequence(seq)
        freqs = np.fft.fftfreq(N)
        mag = sum(compute_fft_magnitude(x) for x in [A, C, G, T])

        peaks, props = find_peaks(mag[:N//2], height=np.mean(mag[:N//2]) + 2*np.std(mag[:N//2]))
        heights = props["peak_heights"]
        top = np.argsort(heights)[::-1][:3]
        selected = [(peaks[i], freqs[peaks[i]]) for i in top if peaks[i] > 0]

        motifs = []
        for peak, freq in selected:
            period = int(N / peak)
            if 2 <= period <= 30:
                motifs.append((freq, extract_motifs(seq, period)))

        safe = safe_filename(meta["gene"] + "_" + meta["prot_id"])
        csv_path = os.path.join(spectrum_dir, f"{safe}_fft.csv")
        np.savetxt(csv_path, np.column_stack((freqs, mag)), delimiter=",",
                   header="Freq,Magnitude", comments="")

        fig_path = os.path.join(output_dir, f"{safe}_fft.png")
        plt.figure(figsize=(10, 4))
        plt.plot(freqs[:N//2], mag[:N//2])

        for peak, freq in selected:
            plt.axvline(x=freq, color='red', linestyle='--')
            label = ""
            for f, mlist in motifs:
                if abs(f - freq) < 1e-4 and mlist:
                    label = mlist[0][0]
                    break
            if label:
                plt.text(freq, max(mag)*0.7, label, rotation=90, ha='center', fontsize=6, color='blue')

        plt.title(f"{meta['gene']} ({meta['prot_id']})")
        plt.xlabel("Frequency")
        plt.ylabel("Magnitude")
        plt.tight_layout()
        plt.savefig(fig_path)
        plt.close()

        return {**meta, "fft": fig_path, "motifs": motifs,
                "peak_freqs": [round(f, 3) for _, f in selected]}

    except Exception as e:
        return {"error": (record.id, str(e))}

# ------------------------------------------
# write_latex_report unchanged
# ------------------------------------------

def write_latex_report(grouped, errors, output_dir):
    tex = os.path.join(output_dir, "fft_report.tex")
    with open(tex, "w") as f:
        f.write(r"""\documentclass{article}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{geometry}
\usepackage{float}
\usepackage{hyperref}
\geometry{margin=1in}
\title{Fourier Spectral Analysis of Coding Sequences}
\author{Automated Report}
\date{\today}
\begin{document}
\maketitle
\tableofcontents
\newpage
""")
        for group in [1, 2, 3]:
            group_title = {
                1: "Group 1: Single Peak at 0.33",
                2: "Group 2: Single Peak Elsewhere",
                3: "Group 3: Multiple Peaks"
            }[group]
            f.write(f"\\section{{{escape_latex(group_title)}}}\n")

            genes = defaultdict(list)
            for r in grouped[group]:
                genes[r["gene"]].append(r)

            for gene, records in genes.items():
                label = safe_filename(f"{group}_{gene}")
                f.write(f"\\subsection{{{escape_latex(gene)}}}\\label{{sec:{label}}}\n")
                for r in records:
                    f.write(f"\\subsubsection{{Protein: {escape_latex(r['prot_id'])}}}\n")
                    f.write(f"\\textbf{{Gene ID}}: {escape_latex(r['gene_id'])}\\\\\n")
                    f.write(f"\\textbf{{Nucleotide Accession}}: {escape_latex(r['acc'])}\\\\\n")
                    f.write(f"\\textbf{{Sequence Length}}: {r['length']} bp\\\\\n")
                    f.write(f"\\textbf{{Protein Description}}: \\\\ \\begin{{minipage}}[t]{{\\linewidth}} {escape_latex(r['prot_name'])} \\end{{minipage}} \\\\ \n")
                    f.write(f"\\begin{{figure}}[H]\\centering\n")
                    f.write(f"\\includegraphics[width=0.9\\linewidth]{{{os.path.basename(r['fft'])}}}\n")
                    f.write(f"\\caption{{FFT Spectrum for {escape_latex(r['prot_id'])}}}\n\\end{{figure}}\n")
                    if r["motifs"]:
                        f.write("\\begin{longtable}{|c|c|}\n\\hline\nFrequency & Motifs\\\\\\hline\n")
                        for freq, motif_list in r["motifs"]:
                            mtext = ", ".join([f"{escape_latex(m[0])} ({m[1]})" for m in motif_list])
                            f.write(f"{round(freq,3)} & {mtext}\\\\\\hline\n")
                        f.write("\\end{longtable}\n")
                f.write("\\newpage\n")

        if errors:
            f.write("\\section{Errors Encountered}\n")
            for seq_id, msg in errors:
                f.write(f"\\textbf{{{escape_latex(seq_id)}}}: {escape_latex(msg)}\\\\\n")
        f.write("\\end{document}\n")
    return tex

# ------------------------------------------
# Main: now uses parallel execution + tqdm
# ------------------------------------------

def main(fasta, output_dir, threads=8):
    os.makedirs(output_dir, exist_ok=True)
    spectrum_dir = os.path.join(output_dir, "spectrum_data")
    os.makedirs(spectrum_dir, exist_ok=True)

    # Prepare small dicts to pass through ProcessPool
    tasks = []
    for record in SeqIO.parse(fasta, "fasta"):
        tasks.append({
            "seq": str(record.seq),
            "id": record.id,
            "description": record.description
        })

    results = []
    errors = []

    with ProcessPoolExecutor(max_workers=threads) as exe:
        futures = {exe.submit(analyze_sequence_wrapper, (task, output_dir, spectrum_dir)): task for task in tasks}
        for fut in tqdm(as_completed(futures), total=len(tasks), desc="Processing sequences"):
            res = fut.result()
            if res is None:
                continue
            if "error" in res:
                errors.append(res["error"])
            else:
                results.append(res)

    grouped = defaultdict(list)
    for r in results:
        group = classify_peaks(r["peak_freqs"])
        grouped[group].append(r)

    tex = write_latex_report(grouped, errors, output_dir)
    os.chdir(output_dir)
    os.system("pdflatex -interaction=nonstopmode fft_report.tex > /dev/null")
    os.system("pdflatex -interaction=nonstopmode fft_report.tex > /dev/null")

    print(f"\n PDF generated: {os.path.abspath('fft_report.pdf')}")
    print(f" Spectrum data: {spectrum_dir}")
    if errors:
        print("⚠️ Some sequences failed. See end of PDF.")

# CLI
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Input multi-FASTA file")
    parser.add_argument("-o", "--output", default="fft_output", help="Output folder")
    parser.add_argument("-t", "--threads", default=8, type=int, help="Number of CPU cores")
    args = parser.parse_args()
    main(args.fasta, args.output, args.threads)
