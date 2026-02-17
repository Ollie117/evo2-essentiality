# Evo2 Gene Essentiality Prediction
This repository includes the scripts for replicating and extending the methodology employed by the Arc Institute in their paper [Evo2 paper](https://arcinstitute.org/manuscripts/Evo2), implementing the Evo2 genome foundation model to predict gene essentiality in *Mycobacterium tuberculosis*.

## Overview

The workflow consists of three main stages:

**1. Data preparation (Python)**: The *M. tuberculosis* H37Rv whole genome was loaded from UniProt, and gene essentiality classifications were obtained from DeJesus et al. Genes were matched to classifications, then sequences were perturbed according to each experimental strategy (mutation size variation, genomic context variation, and alternative scoring approaches). The [Evo2 python package](https://pypi.org/project/evo2/#setup) was used implementing the forward module to obtain likelihood scores. Full list of packages listed in requirements below.

**2. Likelihood scoring (Evo2)**: The [Evo2 python package](https://pypi.org/project/evo2/#setup) was used implementing the forward module to obtain likelihood scores for wild-type perturbated sequences. Evo2 analyses were performed on a HPC cluster with a 48GB Nvidia GPU. Full list of hardware and software requrements listed in requirements below.

**3. Statistical analysis and visualisation (R)**: Full AUROC analysis and performance metrics (sensitivity, specificitywere computed using the likelihood scores. Visualizations were generated (ROC curves, boxplots, distribution plots).

## What's included

- **Code**: Python scripts for sequence preparation and perturbation generation; SLURM job submission scripts for Evo2 scoring; R scripts for statistical analysis and visualisation
- **Results**: Performance metrics, figures, and processed outputs from all experiments

**Note**: Raw and processed data files are not included in this repository due to file size limitations. Source data is available from UniProt (H37Rv genome) and DeJesus et al. (2017) (essentiality classifications).

## Requirements
- Python 3.11+ (see `src/requirements.txt` for package dependencies)
- R 4.3+ (see `analysis/requirements.txt` for package dependencies)
- GPU access for Evo2 scoring (NVIDIA L40 / L40S or similar)
- SLURM job scheduler (for HPC submission scripts)
- evo2 0.4.0

## References
- DeJesus MA, Gerrick ER, Xu W, Park SW, Long JE, Boutte CC, Rubin EJ, Schnappinger D, Ehrt S, Fortune SM, Sassetti CM, Ioerger TR. Comprehensive Essentiality Analysis of the Mycobacterium tuberculosis Genome via Saturating Transposon Mutagenesis. mBio. 2017;8(1):e02133-16.
- Brixi, G., Durrant, M.G., Ku, J., Poli, M., Brockman, G., Chang, D., Gonzalez, G.A., King, S.H., Li, D.B., Merchant, A.T., Naghipourfar, M., Nguyen, E., Ricci-Tam, C., Romero, D.W., Sun, G., Taghibakshi, A., Vorontsov, A., Yang, B., Deng, M. and Gorton, L. (2025). Genome modeling and design across all domains of life with Evo 2. Biorxiv.
- M. tuberculosis H37Rv genome: NCBI RefSeq GCF_000195955.2
- [Evo2 python package](https://pypi.org/project/evo2/#setup)
