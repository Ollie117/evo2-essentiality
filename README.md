# Evo2 Gene Essentiality Prediction
This project validates the use of the Evo2 genome foundation model's ability to predict gene essentiality in *Mycobacterium tuberculosis*.

The methodology used by the authors is first replicated to validate the results obtained in the paper [https://arcinstitute.org/manuscripts/Evo2](https://arcinstitute.org/manuscripts/Evo2). Then, an alternative strategy to obtaining likelihood through absolute likelihood scores was tested. Then we employed various perturbation strategies in order to glean further information about the model's understanding of the biology. The strategies being: varying the size of the mutations introduced into the mutated sequences; and varying the size of the genomic context given to the model.

## Findings
  When replicating the authors methods for obtianing likelihood scores, results were very close to what the authors reported in the paper: ~ 0.825 AUROC. 
  Absolute likelihood scores came out lower than the authors: ~ 0.664 AUROC.
  Varying the mutation size resulted in the AUROC increasing monotonically (AUROC: 0.760 at 3 bp → 0.825 at 15 bp)
  Varying flanking sequence from 512 bp to 8,192 bp produced stable performance (AUROC: 0.819–0.826), indicating the model does not require distal genomic information for gene scoring. Contrary to what what the paper states about the model using distant genomic elements for predicting gene essentiality.
  
## References
DeJesus MA, Gerrick ER, Xu W, Park SW, Long JE, Boutte CC, Rubin EJ, Schnappinger D, Ehrt S, Fortune SM, Sassetti CM, Ioerger TR. Comprehensive Essentiality Analysis of the Mycobacterium tuberculosis Genome via Saturating Transposon Mutagenesis. mBio. 2017;8(1):e02133-16.

Brixi, G., Durrant, M.G., Ku, J., Poli, M., Brockman, G., Chang, D., Gonzalez, G.A., King, S.H., Li, D.B., Merchant, A.T., Naghipourfar, M., Nguyen, E., Ricci-Tam, C., Romero, D.W., Sun, G., Taghibakshi, A., Vorontsov, A., Yang, B., Deng, M. and Gorton, L. (2025). Genome modeling and design across all domains of life with Evo 2. Biorxiv.
