import os
import re
from pathlib import Path
from typing import Dict, Any

try:
    import yaml
except ImportError:
    yaml = None

DEFAULT_CONFIG_PATH = "config.yaml"

def read_config_file(config_path: str = DEFAULT_CONFIG_PATH) -> Dict[str, Any]:
    if not os.path.exists(config_path):
        return {}

    if yaml is not None:
        with open(config_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
            return data if data else {}

    # Simple line parser for config.yaml if pyyaml is missing
    config = {
        "threads": 8,
        "reference": {},
        "inputs": {},
        "outputs": {}
    }
    cur_sec = None
    with open(config_path, "r", encoding="utf-8") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            if line_str.endswith(":") and not ":" in line_str[:-1]:
                sec_name = line_str[:-1].strip()
                config[sec_name] = {}
                cur_sec = sec_name
            elif ":" in line_str:
                k, v = line_str.split(":", 1)
                k = k.strip()
                v = v.strip().strip('"').strip("'")
                if "#" in v:
                    v = v.split("#")[0].strip()
                if v.isdigit():
                    v = int(v)
                if cur_sec and cur_sec in config and isinstance(config[cur_sec], dict):
                    config[cur_sec][k] = v
                else:
                    config[k] = v
    return config

def load_config(config_path: str = DEFAULT_CONFIG_PATH) -> Dict[str, Any]:
    return read_config_file(config_path)
