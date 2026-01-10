#!/bin/bash
source ~/.bashrc
conda activate torchgpu

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

set -euo pipefail

cd "$HOME/LLM-TS/Weibo-final"

mkdir -p nf_results

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "GPU not detected (nvidia-smi missing). Exiting." >&2
  exit 1
fi

PYTHON_BIN="$CONDA_PREFIX/bin/python"
CSV_PATH="./lda_final.csv"

for h in {8..15}; do
  "$PYTHON_BIN" run_nf_gpu.py --csv "$CSV_PATH" --h "$h" --out "nf_results/no_unem_h${h}.json"
done
