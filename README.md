
                        ANALYSIS FOR ARS2 PAPER


GENERAL INFO:

* Analysis involves
  - ChIPseq data presented first in this paper
  - RNAseq data for siARS2, siZC3H18 and siEGFP presented first in this paper
  - RNAseq data from Meola et al., Mol. Cell 2016
  - TIFseq data from Chen et al., Nature Genetics 2016

* ChIPseq mapping as desribed in Materials and Methods
* RNAseq mapping as desribed in Materials and Methods


ABOUT THE ANALYSIS

* attempts to collect all processing for Figures in paper.

* Essentially all analysis presented is based counting from bigwigs files.
I did this with the deeptools python package. However, this is intended for ChIPseq mostly, to be useful for strand-specific bigwigs, ie RNAseq data requires some home-made scripts. This is not part of this package but can be found under **https://github.com/manschmi/MS_Metagene_Tools**

* output from deeptools is fed into R and then saved as '.RData' files.

* all RMarkdown files should be reproducible from starting with those .RData files.



ABOUT THE R package arsRtools

* Everything to run the analysis scripts should be within the package.

* Structure:
  /R  ....  common functions
  /data  ....  saved data
  /analysis  ....  html output from RMarkdown analysis files

* R code:
  everything is heavily dependent on the "tidyverse" ie perhaps more commonly known as dplyr+ggplot2 but uses also magrittr and perhaps some other parts of tidyverse too.

* Intermediate data is not saved @Github. Please contact ms@mbg.au.dk if you need access to intermediate files.




###############################################################################
############################## Analysis #######################################

All .html output files are in Folder: /analysis...



############# ChIPseq ############

Pasha.R
	example R script applied to all .bam files (ie bam to bigwig)

FigureS0a_ChIPseq_quintiles_and_scaling.html:
	defines quintiles and protein-coding gene body (2nd half) scaling
	(no Figure panels for paper)

FigureS1B_ChIPseq_correlations.html:
	based on ChIPseq signal  from TSS to TES (mean values used for analysis)
	correlations for Figure S1

FigureS1E_BGSub_Scaling_Procedure_scaleregions.html:
	ChIPseq scaling procedure visual explanation.

FigureS1FG_ChIP_all_quintiles_relTSS.html:
	metaprofiles and heatmaps for all quintile genes.

FigureS1H_ChIPseq_vs_RNAseq.html:
	based on mean ChIPseq signal from TSS to TES
	using RNAseq_DESeq2_introns for RNAseq part

FigureS1H_ChIPseq_vs_RNAseq.html:
	based on mean ChIPseq signal from TSS to TES
	using RNAseq_DESeq2_introns for RNAseq part

Figures showing metaprofiles and heatmaps:
  all very similar, essentially same procedure though some details for snRNA and RDH (ie gene names included in heatmap)

  . Figure1A_ChIP_pc_RefseqAllq1_TSSTES_final.html
  . Figure2AB_3AB_S2A_S3A_ChIPseq_snRNA_RDH.html
  . Figure4AB_S4AB_ChIPseq_PROMPT_eRNA.html



############# RNAseq ############

FigureS1D_RNAseq_PCA:
  . based on all Refseq exons

FigureS0b_RNAseq_DESeq2_introns:
  . RNAseq differential expression data used for Figure S1E

Figures metagene and heat maps:
  . Figure2CD_3CD_RNAseq_snRNA_RDH.html
  . Figure4CD_S4CDE_RNAseq_PROMPT_eRNA.html
  . Figure5BC_S5BC_RNAseq_intron1and2_rel_5SS.html
  . Figure5A_S5A_RNAseq_rel_pc_multiexon_TSS.html



############# TIFseq ############

  . Figure5C_S5FG_TIFseq.html
  . Figure5DE_RNAseq_ChIPseq_by_TIFseq_intron1.html
