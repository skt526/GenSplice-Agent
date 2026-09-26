#!/usr/bin/env bash
# ==============================================================================
# GenSplice-Agent Main Execution Wrapper
# ==============================================================================

set -e

# Find conda python
PYTHON_BIN="python3"

for candidate in \
    "$HOME/miniforge3/envs/gensplice-agent/bin/python3" \
    "$HOME/miniforge3/envs/gensplice-agent/bin/python" \
    "$HOME/miniconda3/envs/gensplice-agent/bin/python3" \
    "$HOME/miniconda3/envs/gensplice-agent/bin/python" \
    "/root/miniconda3/envs/gensplice-agent/bin/python3" \
    "/root/miniconda3/envs/gensplice-agent/bin/python" \
    "/opt/conda/envs/gensplice-agent/bin/python3" \
    "/opt/conda/envs/gensplice-agent/bin/python" \
    "$CONDA_PREFIX/bin/python3" \
    "$CONDA_PREFIX/bin/python" \
    ".venv/bin/python3" \
    ".venv/bin/python"; do
    if [ -x "$candidate" ]; then
        PYTHON_BIN="$candidate"
        ENV_BIN="$(dirname "$PYTHON_BIN")"
        export PATH="$ENV_BIN:$PATH"
        export CONDA_PREFIX="$(dirname "$ENV_BIN")"
        export CONDA_DEFAULT_ENV="gensplice-agent"
        break
    fi
done

"$PYTHON_BIN" run_pipeline.py "$@"
