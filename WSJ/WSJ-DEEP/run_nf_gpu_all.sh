#!/bin/bash
source ~/.bashrc
conda activate torchgpu

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p nf_results

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "GPU not detected (nvidia-smi missing). Exiting." >&2
  exit 1
fi

PYTHON_BIN="$CONDA_PREFIX/bin/python"
CSV_PATH="$SCRIPT_DIR/2025-10-MD.csv"
INDICES=(CPIAUCSL CPIAPPSL CPITRNSL CPIMEDSL CPIULFSL)

for idx in "${INDICES[@]}"; do
  for h in {5..12}; do
    "$PYTHON_BIN" run_nf_gpu.py --csv "$CSV_PATH" --index "$idx" --h "$h" \
      --out "nf_results/${idx}_h${h}.json"
  done
done
