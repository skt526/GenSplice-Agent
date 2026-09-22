import os
import re
from pathlib import Path
from typing import List, Dict, Any

FASTQ_EXTENSIONS = ('.fastq', '.fastq.gz', '.fq', '.fq.gz')

def get_sample_key_and_read(filename: str):
    """
    Extracts sample base name and read indicator (R1 or R2) from filename.
    Returns: (sample_base, read_number) or (filename, None)
    
    Supported patterns:
    - sample_1.fq.gz / sample_2.fq.gz
    - sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
    - sample_R1.fastq.gz / sample_R2.fastq.gz
    """
    # Pattern 1: Illumina style _R1_001.fastq.gz / _R2_001.fastq.gz
    m1 = re.search(r'^(?P<base>.+?)_R(?P<read>[12])_(?P<lane>\d+)\.(?:fastq|fq)(?:\.gz)?$', filename, re.IGNORECASE)
    if m1:
        base = f"{m1.group('base')}_{m1.group('lane')}"
        return base, f"R{m1.group('read')}"

    # Pattern 2: Standard _R1.fastq.gz / _R2.fastq.gz or _R1_*.fastq.gz / _R2_*.fastq.gz
    m2 = re.search(r'^(?P<base>.+?)_R(?P<read>[12])(?:_.*)?\.(?:fastq|fq)(?:\.gz)?$', filename, re.IGNORECASE)
    if m2:
        return m2.group('base'), f"R{m2.group('read')}"

    # Pattern 3: Simple _1.fq.gz / _2.fq.gz or _1.fastq.gz / _2.fastq.gz
    m3 = re.search(r'^(?P<base>.+?)_(?P<read>[12])\.(?:fastq|fq)(?:\.gz)?$', filename, re.IGNORECASE)
    if m3:
        return m3.group('base'), f"R{m3.group('read')}"

    return None, None

def detect_fastq_samples_in_dir(directory_path: str, group_name: str) -> List[Dict[str, Any]]:
    """
    Scans a directory for FASTQ files and automatically pairs R1 and R2 reads.
    Returns a list of sample dictionaries.
    """
    dir_path = Path(directory_path)
    if not dir_path.exists() or not dir_path.is_dir():
        return []

    files = [f for f in dir_path.iterdir() if f.is_file() and f.name.endswith(FASTQ_EXTENSIONS)]
    
    r1_files: Dict[str, Path] = {}
    r2_files: Dict[str, Path] = {}
    single_files: List[Path] = []

    for file in sorted(files):
        sample_base, read_num = get_sample_key_and_read(file.name)
        if sample_base and read_num == "R1":
            r1_files[sample_base] = file
        elif sample_base and read_num == "R2":
            r2_files[sample_base] = file
        else:
            single_files.append(file)

    samples = []
    all_sample_bases = sorted(list(set(list(r1_files.keys()) + list(r2_files.keys()))))

    for sample_base in all_sample_bases:
        r1 = r1_files.get(sample_base)
        r2 = r2_files.get(sample_base)

        if r1 and r2:
            read_type = "paired"
            status = "paired"
        elif r1:
            read_type = "single"
            status = "unpaired_r1_only"
        else:
            read_type = "single"
            status = "unpaired_r2_only"

        samples.append({
            "group": group_name,
            "sample_name": sample_base,
            "read1": str(r1) if r1 else "",
            "read1_filename": r1.name if r1 else "",
            "read2": str(r2) if r2 else "",
            "read2_filename": r2.name if r2 else "",
            "read_type": read_type,
            "status": status
        })

    for file in single_files:
        sample_name = re.sub(r'\.(fastq|fq)(\.gz)?$', '', file.name, flags=re.IGNORECASE)
        samples.append({
            "group": group_name,
            "sample_name": sample_name,
            "read1": str(file),
            "read1_filename": file.name,
            "read2": "",
            "read2_filename": "",
            "read_type": "single",
            "status": "single_end"
        })

    return samples

def scan_all_inputs(inputs_base_dir: str = "./inputs") -> Dict[str, Any]:
    """
    Scans inputs/control and inputs/treatment directories.
    Returns structured results and summary counts.
    """
    control_dir = os.path.join(inputs_base_dir, "control")
    treatment_dir = os.path.join(inputs_base_dir, "treatment")

    control_samples = detect_fastq_samples_in_dir(control_dir, "control")
    treatment_samples = detect_fastq_samples_in_dir(treatment_dir, "treatment")

    all_samples = control_samples + treatment_samples

    total_samples = len(all_samples)
    paired_samples = sum(1 for s in all_samples if s["read_type"] == "paired")

    return {
        "control_samples": control_samples,
        "treatment_samples": treatment_samples,
        "all_samples": all_samples,
        "summary": {
            "total_samples": total_samples,
            "control_count": len(control_samples),
            "treatment_count": len(treatment_samples),
            "paired_count": paired_samples,
            "single_count": total_samples - paired_samples
        }
    }
