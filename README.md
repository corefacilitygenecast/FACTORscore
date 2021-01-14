# FACTORscore

## 1. Introduction
-----------------
This script can process any continuous data(eg. RNA expressing data and CNV log2ratio values) to establish the FACTORscore, refering to this literature: Tumor microenvironment characterization in gastric cancer identifies prognostic and immunotherapeutically relevant gene signatures.

FACTORscore is an R script. 

Main steps: 
1. Normalized the data for unsupervised clustering. Cluster analysis was conducted using the ConsensusClusterPlus package in R. 
2. Sample groups are identified by Kaplan-Meier (KM) survival analysis among cluster groups and unsupervised clustering methods in item-consensus and cluster-consensus plots. 
3. Subsequently, the wilcoxon rank-sum test or Kruskal-Wallis test was utilized to perform differential analysis for each gene among the sample groups, and the genes with p-value < 0.05 were identified as differential genes. The consensus clustering algorithm was also applied to define the cluster of differential genes. 
4. the Principal Component Analysis (PCA) algorithm was used to perform dimension reduction in order to obtain principal component 1 of each gene group. The PC1 value of gene groups was performed cox regression analysis to validate positive or negative of coefficient. Then we applied a method to define the FACTORscore of each sample:
    CNVscore=∑PC1i – ∑PC1j
where i is the PC1 of groups whose HR is larger than 1, and j is the PC1 of groups whose HR less than 1.

## 2. R packages version:
   R                    : 3.5.1<br>
   getopt               : 1.20.2<br>
   ConsensusClusterPlus : 1.46.0<br>
   ggsci                : 2.9<br>
   ggplot2              : 3.3.2<br>
   ComplexHeatmap       : 1.18.1<br>
   circlize             : 0.4.9<br>
   ggpubr               : 0.2<br>
   tidyverse            : 1.3.0<br>
   FactoMineR           : 1.41<br>
   factoextra           : 1.0.5<br>
   survminer            : 0.4.4.999<br>
   survival             : 3.1.8<br>
   broom                : 0.5.6<br>
   stringr              : 1.4.0<br>
   
## 3. Usage
  Rscript FactorScore.R [-[-help|h]] [-[-input|i] <character>] [-[-time|t] <character>] [-[-fpkm|f] <character>] [-[-outdir|o] <character>] [-[-prefix|p] <character>]<br>
    -h|--help      this help<br>
    -i|--input     input matrix<br>
    -t|--time      survival time and survival state<br>
    -f|--fpkm      T/TR/C（T: fpkm/tpm by gene; TR: fpkm/tpm by transcript; C: other data)<br>
    -o|--outdir    pathway of output files<br>
    -p|--prefix    prefix of the results<br>

format of input matrix file：<br>
        top-left cell： must be "sample" character.<br>
        CNV data     ： row names are samples, columns names are genes;<br>
        RNA data     ： （1）fpkm/tpm by gene：row names are samples, column names are genes;<br>
                       （2）fpkm/tpm by transcript：column names are samples, row names are genes.<br>
format of time file：the title of this file must be "sample    OS  re" separated by tab.<br>

## Citation

If you use FACTORscore in published research, please cite: Tumor microenvironment characterization in gastric cancer identifies prognostic and imunotherapeutically relevant gene signatures. Cancer Immunology Research, 2019, 7(5), 737-750. DOI: 10.1158/2326-6066.CIR-18-0436, PMID: 30842092<br>

Contact<br>
E-mail any questions to shuijingnvhaifl@163.com<br>
        
