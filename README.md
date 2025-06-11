# Project Code for [LLM-Powered CPI Prediction Inference with Online Text Time Series]

This repository contains the R scripts and RMarkdown files used for the figures, real data analysis, and simulations in the paper **[LLM-Powered CPI Prediction Inference with Online Text Time Series]**. Below is an overview of the files and their purposes in `LLM-CPI code.zip'.

---

## **Contents**

### **1. Plot**
This folder contains the scripts for generating Figure 1 and Figure 3 from the paper.
- `introduction_plot.R`: R script for generating **Figure 1**.
- `error_density_plot.R`: R script for generating **Figure 3**.

---

### **2. Real Data**
This folder includes the scripts for real data analysis.
- `read_data.Rmd`: RMarkdown file for all real data analyses, excluding Section 5.4.
- `COV9.R`: R script for all real data analyses, excluding Section 5.4.

---

### **3. Simulation**
This folder contains all R scripts used for the simulations described in **Section 4** of the paper.
- `interval_miss_model` : Prediction interval simulations under Omitted relevant predictor scenario
- `interval_over_model` : Prediction interval simulations under model overfitting scenario
- `interval_t_model` :  Prediction interval simulations under t-distribution scenario
- `interval_true_model` :  Prediction interval simulations under true scenario
- `miss_model_simulation.R` : Point prediction under Omitted relevant predictor scenario
- `interval_over_model.R` : Point prediction under under model overfitting scenario
- `interval_t_model.R` : Point prediction under t-distribution scenario
- `interval_true_model.R` :  Point prediction under true scenario
---

### **4. Data**
- `meanTopicDistri_mergeCPI100_Unemployment_df19-23_20250301.csv`: LDA embedding
- `meanVector_mergeCPI_df19-23_20250201.csv`: BERT embedding


---
