"""
GenSplice-Agent Configuration Loader Engine
Loads YAML configuration files cleanly using PyYAML.
"""

import os
import yaml
from typing import Dict, Any

DEFAULT_CONFIG_PATH = "config.yaml"

def load_config(config_path: str = DEFAULT_CONFIG_PATH) -> Dict[str, Any]:
    """
    Reads and parses a YAML configuration file. Returns an empty dict if file does not exist.
    """
    if not os.path.exists(config_path):
        return {}

    with open(config_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
        return data if isinstance(data, dict) else {}

read_config_file = load_config
