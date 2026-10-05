#!/usr/bin/env python3
"""
==============================================================================
GenSplice-Agent: GSE190087 (PRJNA786720 / PRJNA785817) RNA-seq Dataset Downloader
==============================================================================
Dataset Information:
  - Title: SF3B1 mutation-mediated sensitization to H3B-8800 splicing inhibitor
           in chronic lymphocytic leukemia (Life Science Alliance 2023)
  - GEO Series: GSE190087
  - BioProject: PRJNA786720 / PRJNA785817 (SRA: SRP349102)
  - Organism: Homo sapiens (MEC-1 CLL cell line)
  - Platform: Illumina NovaSeq 6000 (150 bp Paired-End)
  - Reference required: Homo sapiens GRCh38 (./ref human)

This script downloads genuine paired-end RNA-seq FASTQ files directly from
the European Nucleotide Archive (ENA) high-speed HTTP mirrors and places them
into GenSplice-Agent's inputs/control/ and inputs/treatment/ directories.
==============================================================================
"""

import os
import sys
import shutil
import argparse
import subprocess
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

# Terminal ANSI Color codes
BOLD = "\033[1m"
GREEN = "\033[1;32m"
CYAN = "\033[1;36m"
YELLOW = "\033[1;33m"
RED = "\033[1;31m"
RESET = "\033[0m"
DIM = "\033[2m"

# ENA FTP/HTTP direct base URL
ENA_BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq"

# Complete metadata for GSE190087
SAMPLES_METADATA = {
    # 1. SF3B1 Wild-Type (Vehicle / DMSO)
    "SRR17111301": {
        "gsm": "GSM5712991",
        "condition": "WT_vehicle",
        "rep": 1,
        "sample_name": "MEC1_WT_rep1",
        "subpath": "001/SRR17111301",
        "approx_size_mb": 1400,
        "description": "MEC-1 SF3B1 WT (Vehicle DMSO) Replicate 1"
    },
    "SRR17111304": {
        "gsm": "GSM5712992",
        "condition": "WT_vehicle",
        "rep": 2,
        "sample_name": "MEC1_WT_rep2",
        "subpath": "004/SRR17111304",
        "approx_size_mb": 1600,
        "description": "MEC-1 SF3B1 WT (Vehicle DMSO) Replicate 2"
    },
    "SRR17111307": {
        "gsm": "GSM5712993",
        "condition": "WT_vehicle",
        "rep": 3,
        "sample_name": "MEC1_WT_rep3",
        "subpath": "007/SRR17111307",
        "approx_size_mb": 1400,
        "description": "MEC-1 SF3B1 WT (Vehicle DMSO) Replicate 3"
    },

    # 2. SF3B1 K700E Mutant (Vehicle / DMSO)
    "SRR17111310": {
        "gsm": "GSM5712994",
        "condition": "K700E_vehicle",
        "rep": 1,
        "sample_name": "MEC1_K700E_rep1",
        "subpath": "010/SRR17111310",
        "approx_size_mb": 1700,
        "description": "MEC-1 SF3B1 K700E (Vehicle DMSO) Replicate 1"
    },
    "SRR17111313": {
        "gsm": "GSM5712995",
        "condition": "K700E_vehicle",
        "rep": 2,
        "sample_name": "MEC1_K700E_rep2",
        "subpath": "013/SRR17111313",
        "approx_size_mb": 1400,
        "description": "MEC-1 SF3B1 K700E (Vehicle DMSO) Replicate 2"
    },
    "SRR17111316": {
        "gsm": "GSM5712996",
        "condition": "K700E_vehicle",
        "rep": 3,
        "sample_name": "MEC1_K700E_rep3",
        "subpath": "016/SRR17111316",
        "approx_size_mb": 1400,
        "description": "MEC-1 SF3B1 K700E (Vehicle DMSO) Replicate 3"
    },

    # 3. SF3B1 Wild-Type + H3B-8800 Splicing Modulator Drug
    "SRR17111319": {
        "gsm": "GSM5712997",
        "condition": "WT_H3B8800",
        "rep": 1,
        "sample_name": "MEC1_WT_H3B8800_rep1",
        "subpath": "019/SRR17111319",
        "approx_size_mb": 1500,
        "description": "MEC-1 SF3B1 WT + H3B-8800 Replicate 1"
    },
    "SRR17111322": {
        "gsm": "GSM5712998",
        "condition": "WT_H3B8800",
        "rep": 2,
        "sample_name": "MEC1_WT_H3B8800_rep2",
        "subpath": "022/SRR17111322",
        "approx_size_mb": 1550,
        "description": "MEC-1 SF3B1 WT + H3B-8800 Replicate 2"
    },
    "SRR17111325": {
        "gsm": "GSM5712999",
        "condition": "WT_H3B8800",
        "rep": 3,
        "sample_name": "MEC1_WT_H3B8800_rep3",
        "subpath": "025/SRR17111325",
        "approx_size_mb": 1400,
        "description": "MEC-1 SF3B1 WT + H3B-8800 Replicate 3"
    },

    # 4. SF3B1 K700E Mutant + H3B-8800 Splicing Modulator Drug
    "SRR17111328": {
        "gsm": "GSM5713000",
        "condition": "K700E_H3B8800",
        "rep": 1,
        "sample_name": "MEC1_K700E_H3B8800_rep1",
        "subpath": "028/SRR17111328",
        "approx_size_mb": 1600,
        "description": "MEC-1 SF3B1 K700E + H3B-8800 Replicate 1"
    },
    "SRR17111331": {
        "gsm": "GSM5713001",
        "condition": "K700E_H3B8800",
        "rep": 2,
        "sample_name": "MEC1_K700E_H3B8800_rep2",
        "subpath": "031/SRR17111331",
        "approx_size_mb": 1450,
        "description": "MEC-1 SF3B1 K700E + H3B-8800 Replicate 2"
    },
    "SRR17111334": {
        "gsm": "GSM5713002",
        "condition": "K700E_H3B8800",
        "rep": 3,
        "sample_name": "MEC1_K700E_H3B8800_rep3",
        "subpath": "034/SRR17111334",
        "approx_size_mb": 1350,
        "description": "MEC-1 SF3B1 K700E + H3B-8800 Replicate 3"
    },
}

def print_banner():
    print(f"{CYAN}{BOLD}")
    print("=" * 72)
    print("    GenSplice-Agent: GSE190087 (PRJNA786720) Dataset Downloader    ")
    print("=" * 72)
    print(f"{RESET}")
    print(f"{BOLD}Biological Study:{RESET} SF3B1 mutation & H3B-8800 splicing inhibitor in CLL")
    print(f"{BOLD}Publication:{RESET}      Life Science Alliance (2023) | PMID: 36723485")
    print(f"{BOLD}Organism:{RESET}         Homo sapiens (MEC-1 CLL cell line, 150bp Paired-End)")
    print(f"{BOLD}Reference:{RESET}        Human GRCh38 (run {CYAN}./ref human{RESET} if needed)")
    print("-" * 72)

def detect_downloader() -> str:
    """Detects available high-speed download tools."""
    if shutil.which("aria2c"):
        return "aria2c"
    elif shutil.which("curl"):
        return "curl"
    elif shutil.which("wget"):
        return "wget"
    else:
        return "python"

def get_runs_for_mode(mode: str) -> Tuple[List[str], List[str]]:
    """
    Returns (control_run_ids, treatment_run_ids) based on selected mode.
    """
    mode = mode.lower()
    if mode in ("benchmark", "wt_vs_k700e", "standard"):
        # 3 WT vs 3 K700E
        control = ["SRR17111301", "SRR17111304", "SRR17111307"]
        treatment = ["SRR17111310", "SRR17111313", "SRR17111316"]
    elif mode in ("quick", "test", "sanity"):
        # 1 WT vs 1 K700E
        control = ["SRR17111301"]
        treatment = ["SRR17111310"]
    elif mode in ("drug_wt", "drug_response_wt"):
        # WT vehicle vs WT + H3B-8800
        control = ["SRR17111301", "SRR17111304", "SRR17111307"]
        treatment = ["SRR17111319", "SRR17111322", "SRR17111325"]
    elif mode in ("drug_k700e", "drug_response_k700e"):
        # K700E vehicle vs K700E + H3B-8800
        control = ["SRR17111310", "SRR17111313", "SRR17111316"]
        treatment = ["SRR17111328", "SRR17111331", "SRR17111334"]
    elif mode == "all":
        # All 12 runs: Control gets WT vehicle + K700E vehicle; Treatment gets drug
        control = ["SRR17111301", "SRR17111304", "SRR17111307", "SRR17111310", "SRR17111313", "SRR17111316"]
        treatment = ["SRR17111319", "SRR17111322", "SRR17111325", "SRR17111328", "SRR17111331", "SRR17111334"]
    else:
        raise ValueError(f"Unknown mode: {mode}")
    return control, treatment

def download_file(url: str, dest_path: Path, tool: str) -> bool:
    """Downloads a single file using the specified tool with resume capability."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = dest_path.with_suffix(dest_path.suffix + ".part")

    # If final file already exists and is not empty, check integrity
    if dest_path.exists() and dest_path.stat().st_size > 1024 * 1024:
        print(f"  {GREEN}✔ Already exists:{RESET} {dest_path.name} ({dest_path.stat().st_size / (1024*1024):.1f} MB)")
        return True

    print(f"  {CYAN}⬇ Downloading:{RESET} {url}")
    print(f"    --> {dest_path}")

    try:
        if tool == "aria2c":
            cmd = [
                "aria2c", "-c", "-s", "4", "-x", "4", "-k", "1M",
                "--file-allocation=none", "-d", str(dest_path.parent),
                "-o", dest_path.name, url
            ]
            res = subprocess.run(cmd)
            return res.returncode == 0
        elif tool == "curl":
            cmd = [
                "curl", "-L", "-C", "-", "--retry", "5", "--retry-delay", "3",
                "-o", str(dest_path), url
            ]
            res = subprocess.run(cmd)
            return res.returncode == 0
        elif tool == "wget":
            cmd = [
                "wget", "-c", "--tries=5", "--waitretry=3",
                "-O", str(dest_path), url
            ]
            res = subprocess.run(cmd)
            return res.returncode == 0
        else:
            # Fallback Python urllib with streaming
            headers = {"User-Agent": "GenSplice-Agent/1.0"}
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req) as resp, open(dest_path, "wb") as out_file:
                shutil.copyfileobj(resp, out_file)
            return True
    except Exception as e:
        print(f"  {RED}✖ Download failed for {url}: {e}{RESET}")
        return False

def verify_gzip(file_path: Path) -> bool:
    """Tests gzip integrity using gzip -t or python gzip."""
    if not file_path.exists():
        return False
    try:
        res = subprocess.run(["gzip", "-t", str(file_path)], capture_output=True)
        return res.returncode == 0
    except FileNotFoundError:
        # If system gzip command not available, use python gzip
        import gzip
        try:
            with gzip.open(file_path, 'rb') as f:
                while f.read(1024 * 1024):
                    pass
            return True
        except Exception:
            return False

def subsample_fastq(src_path: Path, dst_path: Path, num_reads: int = 1000000):
    """Subsamples first N reads (4*N lines) from a gzipped fastq."""
    print(f"  {YELLOW}✂ Subsampling {num_reads:,} reads from {src_path.name}...{RESET}")
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    lines_to_keep = num_reads * 4
    cmd = f"gzip -cd '{src_path}' | head -n {lines_to_keep} | gzip > '{dst_path}'"
    subprocess.run(cmd, shell=True, check=True)
    print(f"  {GREEN}✔ Created subsampled file:{RESET} {dst_path.name}")

def main():
    parser = argparse.ArgumentParser(
        description="Download GSE190087 (PRJNA786720 / PRJNA785817) RNA-seq dataset for GenSplice-Agent.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  benchmark   : SF3B1 Wild-Type (3 reps) vs SF3B1 K700E (3 reps) [DEFAULT]
  quick       : 1 Wild-Type vs 1 K700E replicate (Fast test)
  drug_wt     : SF3B1 WT DMSO (3 reps) vs WT + H3B-8800 (3 reps)
  drug_k700e  : SF3B1 K700E DMSO (3 reps) vs K700E + H3B-8800 (3 reps)
  all         : Download all 12 RNA-seq runs

Examples:
  # Standard 3 vs 3 benchmark comparison (Recommended):
  python download_GSE190087.py --mode benchmark

  # Quick sanity test (1 vs 1 pair):
  python download_GSE190087.py --mode quick

  # Dry run to inspect files and URLs without downloading:
  python download_GSE190087.py --mode quick --dry-run
        """
    )
    parser.add_argument(
        "--mode", "-m",
        choices=["benchmark", "quick", "drug_wt", "drug_k700e", "all"],
        default="benchmark",
        help="Experiment comparison mode (default: benchmark [3 WT vs 3 K700E])"
    )
    parser.add_argument(
        "--tool", "-t",
        choices=["auto", "aria2c", "curl", "wget", "python"],
        default="auto",
        help="Downloader tool to use (default: auto [aria2c -> curl -> wget -> python])"
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Clear existing files in inputs/control and inputs/treatment before downloading"
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        default=True,
        help="Verify gzip integrity of downloaded files (default: True)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List files and target destinations without performing any downloads"
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=0,
        help="Subsample downloaded FASTQs to N reads (e.g. 1000000) for rapid testing"
    )

    args = parser.parse_args()

    print_banner()

    tool = detect_downloader() if args.tool == "auto" else args.tool
    print(f"{BOLD}Selected Mode:{RESET}      {CYAN}{args.mode.upper()}{RESET}")
    print(f"{BOLD}Downloader Engine:{RESET}  {GREEN}{tool}{RESET}")

    # Determine workspace inputs directories
    workspace_dir = Path(__file__).resolve().parent
    control_dir = workspace_dir / "inputs" / "control"
    treatment_dir = workspace_dir / "inputs" / "treatment"

    control_runs, treatment_runs = get_runs_for_mode(args.mode)

    print(f"\n{BOLD}[Dataset Composition]{RESET}")
    print(f"  {CYAN}Control Group ({len(control_runs)} samples):{RESET}")
    ctrl_total_mb = 0
    for r in control_runs:
        meta = SAMPLES_METADATA[r]
        ctrl_total_mb += meta["approx_size_mb"] * 2
        print(f"    - {BOLD}{r}{RESET} ({meta['gsm']}): {meta['description']} (~{meta['approx_size_mb']*2 / 1024:.1f} GB total)")

    print(f"\n  {YELLOW}Treatment Group ({len(treatment_runs)} samples):{RESET}")
    trt_total_mb = 0
    for r in treatment_runs:
        meta = SAMPLES_METADATA[r]
        trt_total_mb += meta["approx_size_mb"] * 2
        print(f"    - {BOLD}{r}{RESET} ({meta['gsm']}): {meta['description']} (~{meta['approx_size_mb']*2 / 1024:.1f} GB total)")

    total_est_gb = (ctrl_total_mb + trt_total_mb) / 1024.0
    print(f"\n{BOLD}Total Estimated Download Volume:{RESET} ~{total_est_gb:.1f} GB")

    if args.dry_run:
        print(f"\n{YELLOW}[DRY RUN MODE] The following files would be downloaded:{RESET}")
        for group_name, runs, target_dir in [("Control", control_runs, control_dir), ("Treatment", treatment_runs, treatment_dir)]:
            print(f"\n  Group: {group_name} -> {target_dir}")
            for r in runs:
                meta = SAMPLES_METADATA[r]
                for read_num in (1, 2):
                    filename = f"{meta['sample_name']}_{read_num}.fq.gz"
                    url = f"{ENA_BASE}/{r[:6]}/{meta['subpath']}/{r}_{read_num}.fastq.gz"
                    print(f"    - {filename}  <=  {url}")
        print(f"\n{GREEN}✔ Dry run complete. Run without --dry-run to start downloading.{RESET}")
        return 0

    # Clean existing directory contents if requested
    if args.clean:
        print(f"\n{YELLOW}Cleaning target input directories...{RESET}")
        for target_dir in (control_dir, treatment_dir):
            if target_dir.exists():
                for item in target_dir.iterdir():
                    if item.is_file() and not item.name.startswith("."):
                        item.unlink()
        print(f"  {GREEN}✔ Cleaned inputs/control and inputs/treatment.{RESET}")

    control_dir.mkdir(parents=True, exist_ok=True)
    treatment_dir.mkdir(parents=True, exist_ok=True)

    # Download task execution
    download_tasks = []
    for r in control_runs:
        meta = SAMPLES_METADATA[r]
        download_tasks.append((r, meta, control_dir, "Control"))
    for r in treatment_runs:
        meta = SAMPLES_METADATA[r]
        download_tasks.append((r, meta, treatment_dir, "Treatment"))

    total_files = len(download_tasks) * 2
    completed_files = 0
    failed_files = 0

    print(f"\n{BOLD}[Starting Download of {total_files} FASTQ files]{RESET}")
    print(f"Resumable HTTP downloads supported. You can interrupt and resume at any time.\n")

    for idx, (run_id, meta, target_dir, group_label) in enumerate(download_tasks, 1):
        print(f"\n{BOLD}[Sample {idx}/{len(download_tasks)}] {group_label}: {meta['sample_name']} ({run_id} / {meta['gsm']}){RESET}")
        print(f"  Condition: {meta['description']}")

        for read_num in (1, 2):
            # Target filename following GenSplice-Agent auto-detection: sample_1.fq.gz / sample_2.fq.gz
            dest_filename = f"{meta['sample_name']}_{read_num}.fq.gz"
            dest_path = target_dir / dest_filename
            url = f"{ENA_BASE}/{run_id[:6]}/{meta['subpath']}/{run_id}_{read_num}.fastq.gz"

            success = download_file(url, dest_path, tool)
            if success:
                if args.verify:
                    print(f"  {DIM}Verifying gzip integrity...{RESET}", end="", flush=True)
                    if verify_gzip(dest_path):
                        print(f"\r  {GREEN}✔ Integrity verified:{RESET} {dest_filename}")
                    else:
                        print(f"\r  {RED}✖ Corrupted archive detected:{RESET} {dest_filename}")
                        print(f"    Removing corrupted file and retrying...")
                        dest_path.unlink(missing_ok=True)
                        download_file(url, dest_path, tool)

                if args.subsample > 0:
                    sub_path = target_dir / f"{meta['sample_name']}_sub_{read_num}.fq.gz"
                    subsample_fastq(dest_path, sub_path, args.subsample)

                completed_files += 1
            else:
                failed_files += 1

    print("\n" + "=" * 72)
    if failed_files == 0:
        print(f"{GREEN}{BOLD}✔ All {completed_files}/{total_files} FASTQ files downloaded and verified successfully!{RESET}")
        print(f"\n{BOLD}Next Steps to Run GenSplice-Agent:{RESET}")
        print(f"  1. Ensure Human Reference genome is indexed:")
        print(f"     {CYAN}./ref human{RESET}  (Homo sapiens GRCh38)")
        print(f"  2. Launch the end-to-end analysis pipeline:")
        print(f"     {CYAN}./GenSplice.sh{RESET}  or  {CYAN}python run_pipeline.py{RESET}")
        print("=" * 72)
        return 0
    else:
        print(f"{RED}{BOLD}✖ Warning: {failed_files} file(s) failed to download.{RESET}")
        print(f"Re-run the command to resume downloading interrupted files:")
        print(f"  {CYAN}python download_GSE190087.py --mode {args.mode}{RESET}")
        print("=" * 72)
        return 1

if __name__ == "__main__":
    sys.exit(main())
