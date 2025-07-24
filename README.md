## cefaclor-anaphylaxis
WES-based association study of cefaclor-induced anaphylaxis


This repository contains the resources used for performing association analysis and pre-ranked Gene Set Enrichment Analysis (GSEA) in our study. 
The repository includes analysis scripts, metadata submitted to EGA, and reproducible workflows.


## Exome-wide association analysis #변경필
- Input: Genotype and phenotype data (not shared due to controlled access)
- Usage: See `scripts/run_analysis.sh`
- `scripts/run_analysis.sh`: Shell pipeline for filtering
- `scripts/analysis.R`: Prepares inputs and visualizations (e.g., Manhattan plots)
  - Assumes Linux (Ubuntu) environment with R and Java installed


## Pre-ranked GSEA  #변경필
- Input file: `preranked_gsea_input.txt`
  - Format: tab-delimited, with gene symbol and ranking score (e.g., t-statistic or beta)
- Tool: GSEA (Broad Institute)


## EGA Submission
- Metadata associated with the dataset is submitted to the **European Genome-Phenome Archive (EGA)**
- EGA accession: **[EGAD50000001660](https://ega-archive.org/datasets/EGAD50000001660)**
- Access requires approval via the corresponding Data Access Committee (DAC)


## License
This repository is released under the MIT License.  
Please cite our associated publication if you use any part of this repository.


## 🧪 Contact
For questions or access requests, please contact:
**SUNGRYEOL KIM**
Division of Pulmonology, Allergy and Critical Care Medicine, Department of Internal Medicine, Yongin Severance Hospital, Yonsei University College of Medicine
sungryeol@yuhs.ac

**DAEUN LEE**
College of Pharmacy, Seoul National University, Seoul 08826, Republic of Korea
ldaeun2@snu.ac.kr
