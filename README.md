# Splicing-Analysis

A reproducible RNA-seq splicing analysis pipeline built for viral infection models in rhesus macaque.  
This repository contains:

- **HPC/SLURM bash scripts** to run the full pipeline from raw FASTQ files to rMATS outputs.
- **Python utilities** for post-splicing analysis (filtering, gene mapping, ranking).
- **Transcript-level tools** to infer and summarize transcript consequences of alternative splicing.

The code was originally developed and run on the UMass Lowell Unity SLURM cluster for multiple viral infection datasets.

---

## Repository structure

Splicing-Analysis/
├─ scripts/
│  ├─ download_sra.sh
│  ├─ trim_qc.sh
│  ├─ readLength.sh
│  ├─ readLengthMode.sh
│  ├─ build_star_index.sh
│  ├─ run_star_PE.sh
│  ├─ align_star_SE.sh
│  ├─ rMATS_run.sh
│  └─ (other SLURM/bash helpers)
│
├─ post-splicing-analysis/
│  ├─ filter_rmats_events.py
│  ├─ gene_name_map.py
│  ├─ rank_genes_by_event_count.py
│  └─ (additional notebooks / scripts for exploratory analysis)
│
├─ transcripts-mapping/
   ├─ (scripts to aggregate rMATS events at the transcript level)
   └─ (utilities to summarize per-transcript consequences)


---

# 1. Overview

This pipeline supports:

- **Adapter trimming & QC** (Trimmomatic, FastQC, MultiQC)
- **Genome alignment** (STAR)
- **Alternative splicing quantification** (rMATS-turbo)
- **Post-analysis filtering, annotation, and ranking**
- **Transcript-level interpretation**

The workflow was used in analyzing:
- SARS-CoV-2 (GSE156701)
- Mpox (GSE234118)
- H5N1 (GSE287709)

Samples were grouped into **baseline**, **peak**, **resolution**, and **early** (Mpox only).

---

#  2. Methods Summary

### **Reference Genome**
- **Species:** *Macaca mulatta*  
- **Genome:** Mmul_10  
- **Annotation:** Ensembl 108 GTF  
- **STAR Index:** `sjdbOverhang = readLength - 1`  

### **QC + Trimming**
- **Trimmomatic v0.39** for adapter removal  
- **FastQC v0.11.9** for per-sample QC  
- **MultiQC** for aggregated summaries  

### **Alignment**
- **STAR v2.7.x**, outputs sorted BAM files  
- Key params:
  - `--outFilterMultimapNmax 1`
  - `--outSAMtype BAM SortedByCoordinate`
  - `--alignEndsType EndToEnd`

### **rMATS Analysis**
- **rMATS-turbo v4.1.2**
- Uses **JCEC** counts
- Event types:
  - SE, RI, A3SS, A5SS, MXE
- Significance:
  - `FDR < 0.05`
  - `|ΔPSI| > 0.1`

---

#  3. HPC Pipeline Scripts (`scripts/`)

Below is a  summary of what each script does.

### **`download_sra.sh`**
Download raw FASTQ files (via `prefetch` / `fasterq-dump`).  
for personal use, edit:
- SRA accessions  
- Output paths  
---

### **`trim_qc.sh`**
Runs:
- **Trimmomatic** for adapter trimming  
- **FastQC** for raw & trimmed reads  
- Generates MultiQC report  

Supports **paired-end and single-end** modes.
---

### **`readLength.sh` & `readLengthMode.sh`**
Utility scripts to:
- Sample reads  
- Report observed read lengths  
- Extract **mode read length**  
Used for STAR and rMATS parameters.
---

### **`build_star_index.sh`**
Builds STAR index for **Mmul_10** using Ensembl GTF.  
Sets correct `sjdbOverhang`.

---

### **`run_star_PE.sh` / `align_star_SE.sh`**
Aligns FASTQ files to STAR index.

- **PE version:** paired-end libraries  
- **SE version:** single-end libraries  

Outputs: sorted BAM files in designated folder.

---

### **`rMATS_run.sh`**
Runs **rMATS-turbo** using:

- `--b1`: BAM list (group 1)
- `--b2`: BAM list (group 2)
- `--readLength`
- `--libType fr-unstranded`

Outputs event tables for **SE, RI, A3SS, A5SS, MXE**.

---

#  4. Post-Splicing Analysis (`post-splicing-analysis/`)

### **`filter_rmats_events.py`**
Filters rMATS results based on:
- FDR threshold  
- ΔPSI threshold  
- Optional read support filters  

Example:

```bash
python filter_rmats_events.py \
  --input SE.MATS.JCEC.txt \
  --fdr 0.05 --delta_psi 0.1 \
  --output SE.filtered.tsv
