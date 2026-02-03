# Evo2 Gene Essentiality Prediction

This project uses Evo2, a large language model trained on DNA sequences, to predict gene essentiality in *Mycobacterium tuberculosis*

## Key Results

82.5% AUC
81% sensitivity
96.2% NPV 

## Installation

### Requirements
- Python 3.8+
- R 4.0+
- A 48GB minimum NVIDIA GPU

### Setup

```bash
# Python dependencies
pip install -r scripts/model_setup/requirements.txt

# R dependencies
install.packages(c("readxl", "dplyr", "tidyr", "pROC", "ggplot2", "readr"))
```

## References
DeJesus MA, Gerrick ER, Xu W, Park SW, Long JE, Boutte CC, Rubin EJ, Schnappinger D, Ehrt S, Fortune SM, Sassetti CM, Ioerger TR. Comprehensive Essentiality Analysis of the Mycobacterium tuberculosis Genome via Saturating Transposon Mutagenesis. mBio. 2017;8(1):e02133-16.
G. Brixi, M. G. Durrant, J. Ku, M. Poli, G. Brockman, D. Chang, G. A. Gonzalez, S. H. King, D. B. Li, A. T. Merchant, M. Naghipourfar, E. Nguyen, C. Ricci-Tam, D. W. Romero, G. Sun, A. Taghibakshi, A. Vorontsov, B. Yang, M. Deng, L. Gorton, N. Nguyen, N. K. Wang, E. Adams, S. A. Baccus, S. Dillmann, S. Ermon, D. Guo, R. Ilango, K. Janik, A. X. Lu, R. Mehta, M. R. K. Mofrad, M. Y. Ng, J. Pannu, C. Ré, J. C. Schmok, J. St. John, J. Sullivan, K. Zhu, G. Zynda, D. Balsam, P. Collison, A. B. Costa, T. Hernandez-Boussard, E. Ho, M.-Y. Liu, T. McGrath, K. Powell, D. P. Burke, H. Goodarzi, P. D. Hsu and B. L. Hie, Genomics, 2025, preprint
