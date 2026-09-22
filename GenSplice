#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Main Execution Wrapper
# ==============================================================================

set -e

# Find conda python
PYTHON_BIN="python3"
CONDA_ENV_PY="$HOME/miniconda3/envs/gensplice-agent/bin/python"
MINIFORGE_ENV_PY="$HOME/miniforge3/envs/gensplice-agent/bin/python"

if [ -n "$CONDA_PREFIX" ] && [ -x "$CONDA_PREFIX/bin/python" ]; then
    PYTHON_BIN="$CONDA_PREFIX/bin/python"
elif [ -x "$CONDA_ENV_PY" ]; then
    PYTHON_BIN="$CONDA_ENV_PY"
elif [ -x "$MINIFORGE_ENV_PY" ]; then
    PYTHON_BIN="$MINIFORGE_ENV_PY"
fi

"$PYTHON_BIN" run_pipeline.py "$@"
