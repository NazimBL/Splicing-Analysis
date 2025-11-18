#  Splicing-Analysis  
*A Modular RNA-seq Alternative Splicing Pipeline *

This repository provides a **modular**, and **reproducible** pipeline for performing alternative splicing analysis on RNA-seq data using **rMATS-turbo**, followed by **post-processing**, **gene annotation**, and **transcript-level interpretation**.

This workflow was originally developed for the splicing analysis of multiple viral infection models in *Macaca mulatta* (rhesus macaque), but the scripts can adapted to *any species or dataset.

---

#  Features

### ** End-to-end RNA-seq splicing pipeline**
1. Download raw FASTQ data  
2. Perform trimming + QC  
3. Build STAR index  
4. Align reads (Pair-Ended (PE) or Single-Ended (SE))  
5. Run rMATS for differential splicing  
6. Post-process events (filter, annotate, rank)  
7. Map events to transcript-level consequences  

### ** HPC-friendly (SLURM)**
- All steps come with ready-to-run **SLURM scripts**  
- Modular design → you can run steps independently  
- Works with any genome FASTA + GTF  , just need to plug the reference genome file and the appropriate path

### ** downstream analysis tools**
- Filter significant SE/RI/A5SS/A3SS/MXE events  
- Map Ensembl IDs → gene symbols  
- Rank genes by splicing burden  
- Summarize transcript-level impact  
 
---
#  Repository Structure

Splicing-Analysis/
│
├── scripts/ # Main HPC pipeline (FASTQ → QC → STAR → rMATS)
│ ├── download_sra.sh
│ ├── trim_qc.sh
│ ├── readLength.sh
│ ├── readLengthMode.sh
│ ├── build_star_index.sh
│ ├── run_star_PE.sh
│ ├── align_star_SE.sh
│ └── rMATS_run.sh
│
├── post-splicing-analysis/ # Event filtering + annotation + ranking
│ ├── filter_rmats_events.py
│ ├── gene_name_map.py
│ └── rank_genes_by_event_count.py
│
├── transcripts-mapping/ # Transcript-level impact assessment
  └── (scripts summarizing isoform-level consequences)

