# 🧬 Pathogen-Host Protein Interaction Analysis Pipeline

This repository provides a comprehensive pipeline for analyzing pathogen-host protein interactions using a combination of bioinformatics tools and AlphaFold3. It includes scripts for preprocessing protein sequences, preparing input JSONs for AlphaFold3, running predictions, and summarizing the results.

---

## 📁 Repository Contents

- `prepare_json_from_fa.py`: Prepares JSON files from FASTA sequences for AlphaFold3 multimer predictions.
- `post_analysis.py`: Summarizes AlphaFold3 output by extracting and ranking confidence scores.
- `all_in_one_final.sh`: A full pipeline script that integrates SignalP, OrthoFinder, CD-HIT, AlphaFold3, and post-analysis into a Slurm job.

---

## ⚙️ Prerequisites

- Python 3.x
- Slurm Workload Manager
- Singularity
- Required tools:
  - [SignalP 6.0](https://services.healthtech- AlphaFold3

Python packages:
```bash
pip install pandas
