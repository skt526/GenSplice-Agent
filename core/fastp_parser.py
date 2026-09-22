import os
import json
from pathlib import Path
from typing import List

def get_mean_read_length_from_fastp(fastp_json_paths: List[str]) -> int:
    """
    Parses fastp.json report files and calculates the overall average read length.
    Returns integer read length (e.g., 150, 100, 75, 50).
    Defaults to 150 if JSON reports are missing.
    """
    lengths = []
    
    for json_path in fastp_json_paths:
        if os.path.exists(json_path):
            try:
                with open(json_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                # Check summary -> before_filtering -> read1_mean_length
                summary = data.get("summary", {}).get("before_filtering", {})
                r1_len = summary.get("read1_mean_length")
                r2_len = summary.get("read2_mean_length")

                if r1_len:
                    lengths.append(float(r1_len))
                if r2_len:
                    lengths.append(float(r2_len))
            except Exception:
                pass

    if lengths:
        avg_len = sum(lengths) / len(lengths)
        return int(round(avg_len))
    
    return 150
