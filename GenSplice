#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Main Execution Wrapper
# ==============================================================================

set -e

# Find conda python
PYTHON_BIN="python3"

for candidate in \
    "$HOME/miniforge3/envs/gensplice-agent/bin/python" \
    "$HOME/miniconda3/envs/gensplice-agent/bin/python" \
    "/root/miniconda3/envs/gensplice-agent/bin/python" \
    "/opt/conda/envs/gensplice-agent/bin/python" \
    "$CONDA_PREFIX/bin/python" \
    ".venv/bin/python"; do
    if [ -x "$candidate" ]; then
        PYTHON_BIN="$candidate"
        break
    fi
done

"$PYTHON_BIN" run_pipeline.py "$@"
