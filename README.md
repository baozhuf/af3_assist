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
  - SignalP 6.0
  - OrthoFinder
  - CD-HIT
  - AlphaFold3

Python packages:

pip install pandas


## 🚀 Installation
Clone the repository:

git clone https://github.com/baozhuf/af3_assist.git

cd af3_assist

## 🧪 Usage
### Scenario 1. Run the Full Pipeline with Slurm
bash all_in_one_final.sh \
  -a your_slurm_account \
  -e your_email@example.com \
  -p ./pathogen_fasta_dir \
  -i ./host_fasta_dir \
  -l /path/to/af3_model_parameters \
  -f 0.5 \
  -o ./AF3_out

### Scenario 2. You only want to Prepare JSONs for AlphaFold3 on alphafoldserver.com
python prepare_json_from_fa.py \
  --fa1_path path/to/pathogen.fa \
  --fa2_path path/to/host.fa \
  --protein1_cnt 1 \
  --protein2_cnt 1 \
  --num 30 \
  --today 20250501 \
  --out_dir ./output_jsons

### Scenario 3. You only want to Run Post-Analysis to extract AlphaFold3 metrics (pTM, ipTM)
python post_analysis.py \
  --af3_out_dir ./AF3_out \
  --summary_path ./AF3_out/af3_results_summary.csv


## 🧾 Script Argument Descriptions (all_in_one_final.sh)
| Flag | Description |
|------|-------------|
| `-a` | Slurm account name (required) |
| `-e` | Email for job notifications (required) |
| `-c` | Number of CPUs (default: 4) |
| `-m` | Memory in GB (default: 62) |
| `-d` | Number of days for job runtime (default: 1) |
| `-p` | Pathogen FASTA directory (required) |
| `-i` | Host FASTA directory (required) |
| `-l` | AlphaFold3 model parameter directory (required) |
| `-f` | CD-HIT cutoff (range: 0.4–1.0, default: 0.5) |
| `-o` | Output directory (default: `./AF3_out`) |



## 🔄Pipeline Stages
SignalP: Predicts secreted proteins from pathogen sequences.
OrthoFinder: Identifies orthologous groups.
CD-HIT: Clusters proteins to reduce redundancy.
JSON Preparation: Generates input files for AlphaFold3.
AlphaFold3: Predicts protein-protein interactions.
Post-Analysis: Extracts and ranks interaction confidence scores.

## 📬Contact
For questions or contributions, please contact:
Zhenghong Bao
📧 z.bao@ufl.edu


## 📄 License
This project is licensed under the MIT License.

