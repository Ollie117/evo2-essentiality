# Evo2 Gene Essentiality Prediction

This project validates the use of the Evo2 genome foundation model's ability to predict gene essentiality in *Mycobacterium tuberculosis*
The methodology used by the authors in first replicated to validate the results obtained in the paper. https://arcinstitute.org/manuscripts/Evo2
Then, in the interest of understanding the model's biological understanding, further analyses are performed including:
  1. Changing the approach taken by the authors to obtain the likelihood scores from comparing the wildtype sequences to perturbed sequences to calculating likelihood from wiltype alone.
  2. Varying the size of the mutations introduced into the mutated sequences.
  3. Varying the size of the genomic context given to the model.

## References
DeJesus MA, Gerrick ER, Xu W, Park SW, Long JE, Boutte CC, Rubin EJ, Schnappinger D, Ehrt S, Fortune SM, Sassetti CM, Ioerger TR. Comprehensive Essentiality Analysis of the Mycobacterium tuberculosis Genome via Saturating Transposon Mutagenesis. mBio. 2017;8(1):e02133-16.

G. Brixi, M. G. Durrant, J. Ku, M. Poli, G. Brockman, D. Chang, G. A. Gonzalez, S. H. King, D. B. Li, A. T. Merchant, M. Naghipourfar, E. Nguyen, C. Ricci-Tam, D. W. Romero, G. Sun, A. Taghibakshi, A. Vorontsov, B. Yang, M. Deng, L. Gorton, N. Nguyen, N. K. Wang, E. Adams, S. A. Baccus, S. Dillmann, S. Ermon, D. Guo, R. Ilango, K. Janik, A. X. Lu, R. Mehta, M. R. K. Mofrad, M. Y. Ng, J. Pannu, C. Ré, J. C. Schmok, J. St. John, J. Sullivan, K. Zhu, G. Zynda, D. Balsam, P. Collison, A. B. Costa, T. Hernandez-Boussard, E. Ho, M.-Y. Liu, T. McGrath, K. Powell, D. P. Burke, H. Goodarzi, P. D. Hsu and B. L. Hie, Genomics, 2025, preprint
