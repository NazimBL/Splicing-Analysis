#  Splicing-Analysis  
*A Reproducible and Modular RNA-seq Alternative Splicing Pipeline for HPC Environments*

This repository provides a **modular**, and **reproducible** pipeline for performing alternative splicing analysis on RNA-seq data using **STAR** and **rMATS-turbo**, followed by **post-processing**, **gene annotation**, and **transcript-level interpretation**.

This workflow was originally developed for the splicing analysis of multiple viral infection models in *Macaca mulatta* (rhesus macaque), but it is written in a **dataset-agnostic**, **genome-agnostic**, and **reusable** way.  Anyone can adapt it to *any species*, *any dataset*.

---

#  Features

### ** End-to-end RNA-seq splicing pipeline**
1. Download raw FASTQ data  
2. Perform trimming + QC  
3. Build STAR index  
4. Align reads (PE or SE)  
5. Run rMATS for differential splicing  
6. Post-process events (filter, annotate, rank)  
7. Map events to transcript-level consequences  

### ** HPC-friendly (SLURM)**
- All steps come with ready-to-run **SLURM scripts**  
- Modular design → you can run steps independently  
- Works with any genome FASTA + GTF  , just need to edit the reference genome file and the appropriate path

### ** downstream analysis tools**
- Filter significant SE/RI/A5SS/A3SS/MXE events  
- Map Ensembl IDs → gene symbols  
- Rank genes by splicing burden  
- Summarize transcript-level impact  
 

