#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT=./dataset/ETT-small
SUR_LAG=48
PRED_LEN=96

run_dataset() {
  local dataset="$1"
  local data_path="$2"
  shift 2
  local base_args=("$@")

  local data_file="${DATA_ROOT%/}/${data_path}"
  if [[ ! -f "${data_file}" ]]; then
    echo "Missing data file: ${data_file}" >&2
    exit 1
  fi

  mapfile -t TARGET_COLUMNS < <(python3 - "${data_file}" <<'PY'
import csv, sys
path = sys.argv[1]
with open(path, newline='') as fh:
    reader = csv.reader(fh)
    header = next(reader)
cols = [c.strip() for c in header if c.strip().lower() != 'date']
for col in cols:
    print(col)
PY
)

  if [[ ${#TARGET_COLUMNS[@]} -eq 0 ]]; then
    echo "No target columns found in ${data_file}" >&2
    exit 1
  fi

  for target in "${TARGET_COLUMNS[@]}"; do
    local model_id="${dataset}_${target}_pl${PRED_LEN}_PatchTSTDual"
    echo ">>> Running PatchTSTDual | dataset=${dataset} | target=${target}"
    python -u run.py "${base_args[@]}" \
      --model_id "${model_id}" \
      --model PatchTSTDual \
      --data "${dataset}" \
      --data_path "${data_path}" \
      --target "${target}" \
      --pred_len "${PRED_LEN}" \
      --dlinear_surrogate_lag "${SUR_LAG}"
  done
}

# ETTh1 configuration
run_dataset "ETTh1" "ETTh1.csv" \
  --task_name long_term_forecast \
  --is_training 1 \
  --root_path "${DATA_ROOT}/" \
  --features M \
  --seq_len 512 --label_len 128 \
  --enc_in 7 --dec_in 7 --c_out 7 \
  --d_model 320 --d_ff 1024 --n_heads 8 \
  --batch_size 128 --train_epochs 50 --patience 10 \
  --learning_rate 0.0001 \
  --embedding_root ./embedding \
  --dlinear_embedding_topk 4 \
  --dlinear_embedding_lag 96 \
  --split_strategy 16_2_2 \
  --itr 1 \
  --loss_type mse

# ETTh2 configuration
run_dataset "ETTh2" "ETTh2.csv" \
  --task_name long_term_forecast \
  --is_training 1 \
  --root_path "${DATA_ROOT}/" \
  --features M \
  --seq_len 512 --label_len 192 \
  --enc_in 7 --dec_in 7 --c_out 7 \
  --d_model 64 --d_ff 128 --n_heads 2 \
  --batch_size 64 --train_epochs 10 --patience 3 \
  --learning_rate 0.0001 \
  --split_strategy 12_4_4 \
  --embedding_root ./embedding \
  --dlinear_embedding_topk 4 \
  --dlinear_embedding_lag 24 \
  --itr 1 \
  --loss_type mse

# ETTm1 configuration
run_dataset "ETTm1" "ETTm1.csv" \
  --task_name long_term_forecast \
  --is_training 1 \
  --root_path "${DATA_ROOT}/" \
  --features M \
  --seq_len 336 --label_len 192 \
  --enc_in 7 --dec_in 7 --c_out 7 \
  --d_model 32 --d_ff 70 --n_heads 2 \
  --batch_size 32 --train_epochs 10 --patience 3 \
  --learning_rate 0.001 \
  --split_strategy 12_4_4 \
  --embedding_root ./embedding \
  --dlinear_embedding_topk 4 \
  --dlinear_embedding_lag 24 \
  --itr 1 \
  --loss_type mse

# ETTm2 configuration
run_dataset "ETTm2" "ETTm2.csv" \
  --task_name long_term_forecast \
  --is_training 1 \
  --root_path "${DATA_ROOT}/" \
  --features M \
  --seq_len 332 --label_len 192 \
  --enc_in 7 --dec_in 7 --c_out 7 \
  --d_model 64 --d_ff 256 --n_heads 2 \
  --batch_size 64 --train_epochs 5 --patience 1 \
  --learning_rate 0.00005 \
  --split_strategy 12_4_4 \
  --embedding_root ./embedding \
  --dlinear_embedding_topk 4 \
  --dlinear_embedding_lag 24 \
  --itr 1 \
  --loss_type mse
