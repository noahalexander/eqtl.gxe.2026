# eqtl.gxe.2026
Analyses and figures for Alexander et al. (2026)

In order of use through the analysis:

fastq.procesing.sh: initial data processing with fastq files as input (not necessarily a shell script to be run at once but a record of commands I ran. I could make multiple smaller scripts to reflect processing of each dataset)

monocle.fxns.R: we initially used monocle3 for analyses and conducted ESR activity quantification using mononcle3 objects. These functions are related to this phase of analysis.

byxrm.sp.monocle.aucell.R: processing of cellranger outputs using the above functions as well as monocle3 functions. We only used the ESR activity values from this phase of analysis in the later mapping. We switched to seurat for analyses but kept these numbers. This script focuses on the BYxRM stationary phase perturbation experiment/data and following scripts correspond to BYxRM adn CBSxYJM salt perturbation experiments/data.

cbsxyjm.salt.monocle.aucell.R: same as the above but with data for the CBSxYJM cross.

byxrm.salt.monocle.aucell.R: same as the above but with data for the BYxRM cros during the salt perturbation.

pca.segs.R: to conduct pca and generate plots of PC1 and PC2 coordinates with ESR activity values. 

seurat.processing.segregants.R: initial read in and clustering/state assignment of each cross:timepoint combination.

seurat.object.postprocessing.R: generation of output files reflecting the cluster/state assignments.

enrichment.testing.fxn.R: a function to conduct enrichment testing using gprofiler.

combined.level.enrichment.testing.R: enrichment testing for hotspots identified at the combined level.

byxrm.state.enrichment.testing.R: enrichment testing for hotspots identified at the state level in byxrm experiments.

3004.enrichment.testing.R: same as above but for the CBSxYJM cross.

cbsxyjm.state.qtl.granges.R: modify ESR activity QTL confidence intervals to be at most 15kb. note: add same for BYxRM salt and sp experiments.

hotspot.bins.vis.R: generate .bed files of 50kb bins classified as eQTL hotspots to upload to a genome browser for visualization.

esr.qtl.lod.traces.R: generate LOD traces for ESR activity QTL mapping results conducted in all six cross:enviornment combinations.
