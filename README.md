# cSupB
VARI-cSupB is an application of the succinct colored de Bruijn graph constructed by VARI-merge. It overcomes the complexity of graphs and organize a set of species- or population-specific haplotype sequences of interest. Based on this model, a tri-tuple coordinate system that combines an offset value, topological structure and sample information. Additionally, VARI-cSupB provides a novel method that utilizes complete topological information and efficiently detects variants for highly similar samples, which can be validated by simulated datasets. Moreover, VARI-cSupB can adapt to a complex cycle structure.

# Installation
VARI-cSupB is based on the colored de Bruijn graph constructed by [VARI-merge](https://github.com/cosmo-team/cosmo/tree/VARI-merge),so the installation processing is the same as the VARI-merge. We only need to add the files (Makefile and cosmo-cSupB.cpp) to VARI-merge.

# Usage
The construction of graph and color matrix are seen  [VARI-merge](https://github.com/cosmo-team/cosmo/tree/VARI-merge).To obtain the coordinate system or variants, simply use the command line as follows:<tab>
* cosmo-cSupB primates_mtDNAs_kmc_list.dbg primates_mtDNAs_kmc_list.colors.sd_vector >log.txt
