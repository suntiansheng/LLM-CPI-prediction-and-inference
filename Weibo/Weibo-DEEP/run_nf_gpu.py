import argparse
import gc
import inspect
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import torch
from neuralforecast import NeuralForecast
from neuralforecast.models import TimeLLM, PatchTST


                                                              
               
                                                              
FREQ = "MS"
UNIQUE_ID = "cpi"


                                                              
         
                                                              
def _scale_series(values):
    values = np.asarray(values, dtype=float)
    mean = np.mean(values)
    std = np.std(values, ddof=1)
    if std == 0:
        return np.zeros_like(values), mean, std
    return (values - mean) / std, mean, std


def _prepare_monthly(csv_path):
    df = pd.read_csv(csv_path)
    df["time"] = pd.to_datetime(df["time"])
    df["month"] = df["time"].dt.to_period("M").astype(str)

    cpi_month = df.groupby("month", sort=True)["cpi_lastmonth_100"].mean()

    months = cpi_month.index.tolist()
    ds = pd.to_datetime([f"{m}-01" for m in months])

    y_scaled, _, _ = _scale_series(cpi_month.values - 100.0)
    return ds, y_scaled


def _input_size(train_len, h, cap=12):
    max_context = min(cap, train_len - 1)
    target = max(4, h * 2)
    return min(max_context, target)


def _build_model(
    model_cls,
    h,
    input_size,
    max_steps,
    seed,
    use_gpu,
    hist_exog=None,
    futr_exog=None,
    model_size=None,
    llm=None,
):
    """Generic model builder with optional exogenous + LLM support."""
    sig = inspect.signature(model_cls.__init__)
    has_var_kw = any(
        p.kind == inspect.Parameter.VAR_KEYWORD
        for p in sig.parameters.values()
    )
    kwargs = {}

    def set_if(name, value):
        if name in sig.parameters or has_var_kw:
            kwargs[name] = value

    def set_if_strict(name, value):
        if name in sig.parameters:
            kwargs[name] = value

                 
    set_if("h", h)
    set_if("input_size", input_size)
    set_if("max_steps", max_steps)
    set_if("random_seed", seed)

                     
    if hist_exog is not None:
        set_if("hist_exog_list", hist_exog)
    if futr_exog is not None:
        set_if("futr_exog_list", futr_exog)

                          
    if model_size:
        for key, val in model_size.items():
            set_if_strict(key, val)

                                                                 
    if llm is not None:
        set_if_strict("llm", llm)

                     
    if use_gpu:
        set_if("accelerator", "gpu")
        set_if("devices", 1)
        set_if("precision", 16)               
    else:
        set_if("accelerator", "cpu")

    return model_cls(**kwargs)


def _metrics(pred, actual):
    pred = np.asarray(pred, dtype=float)
    actual = np.asarray(actual, dtype=float)
    rmse = np.sqrt(np.mean((pred - actual) ** 2))
    sign = np.mean(np.sign(pred) != np.sign(actual))
    return rmse, sign


def _cleanup(use_gpu, *objs):
    for obj in objs:
        try:
            del obj
        except NameError:
            pass
    if use_gpu:
        torch.cuda.empty_cache()
    gc.collect()


                                                              
      
                                                              
def main():
    parser = argparse.ArgumentParser(
        description="Tiny GPT2 TimeLLM + PatchTST CPI Forecasting"
    )
    parser.add_argument("--csv", default="./lda_final.csv")
    parser.add_argument("--h", type=int, default=3)
    parser.add_argument("--max-steps", type=int, default=5)
    parser.add_argument("--seed", type=int, default=2025)
    parser.add_argument("--no-gpu", action="store_true")
    parser.add_argument("--input-cap", type=int, default=6)
    parser.add_argument("--llm-model", default="openai-community/gpt2")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

               
    ds, y_scaled = _prepare_monthly(args.csv)
    n = len(y_scaled)
    h = args.h

    train_end = n - h
    if train_end < h + 4:
        raise ValueError("Not enough training data for given horizon.")

    df_train = pd.DataFrame({
        "unique_id": UNIQUE_ID,
        "ds": ds[:train_end],
        "y": y_scaled[:train_end],
    })

    futr_df = pd.DataFrame({
        "unique_id": UNIQUE_ID,
        "ds": ds[train_end:],
    })

    hist_exog = None
    futr_exog = None

                                         
    use_gpu = (not args.no_gpu) and torch.cuda.is_available()
    time_llm_gpu = False

    input_size = _input_size(train_end, h, cap=args.input_cap)

                                                                           
    if "llm" not in inspect.signature(TimeLLM.__init__).parameters:
        raise RuntimeError("TimeLLM does not accept an llm model name in this version.")
    llm = args.llm_model

                              
    top_k = 1
    patch_len = 2
    min_input = patch_len + top_k
    if input_size < min_input:
        input_size = min_input

    time_llm_size = {
        "d_llm": 768,
        "d_model": 16,
        "n_heads": 2,
        "top_k": top_k,
        "patch_len": patch_len,
        "stride": 1,
        "batch_size": 1,
        "windows_batch_size": 1,
        "llm_output_hidden_states": True,
        "llm_output_attention": False,
        "prompt_prefix": "The dataset contains monthly CPI data.",
    }

    actual = y_scaled[train_end:]
    results = {}

                                                
    model_time_llm = _build_model(
        TimeLLM,
        h,
        input_size,
        args.max_steps,
        args.seed,
        use_gpu=time_llm_gpu,
        hist_exog=None,
        futr_exog=None,
        model_size=time_llm_size,
        llm=llm,
    )

    nf = NeuralForecast(models=[model_time_llm], freq=FREQ)
    nf.fit(df_train)
    preds = nf.predict(futr_df=futr_df)

    if "TimeLLM" in preds.columns:
        pred = preds["TimeLLM"].values
    else:
        pred = preds.filter(like="TimeLLM").values.flatten()
        if pred.size == 0:
            raise KeyError("TimeLLM predictions not found in output columns.")

    rmse, sign = _metrics(pred, actual)
    results["TimeLLM"] = {
        "rmse": float(rmse),
        "sign_mismatch": float(sign),
        "pred": pred.tolist(),
    }

    _cleanup(time_llm_gpu, nf, model_time_llm, preds, llm)

                                                
    model_patchtst = _build_model(
        PatchTST,
        h,
        input_size,
        args.max_steps,
        args.seed,
        use_gpu=use_gpu,
        hist_exog=hist_exog,
        futr_exog=futr_exog,
        model_size={"d_model": 8, "n_layers": 1, "n_heads": 1},
    )

    nf = NeuralForecast(models=[model_patchtst], freq=FREQ)
    nf.fit(df_train)
    preds = nf.predict(futr_df=futr_df)

    if "PatchTST" in preds.columns:
        pred = preds["PatchTST"].values
    else:
        pred = preds.filter(like="PatchTST").values.flatten()
        if pred.size == 0:
            raise KeyError("PatchTST predictions not found in output columns.")

    rmse, sign = _metrics(pred, actual)
    results["PatchTST"] = {
        "rmse": float(rmse),
        "sign_mismatch": float(sign),
        "pred": pred.tolist(),
    }

    _cleanup(use_gpu, nf, model_patchtst, preds)

    payload = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "h": h,
        "use_gpu": bool(use_gpu),
        "results": results,
    }

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
    else:
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
