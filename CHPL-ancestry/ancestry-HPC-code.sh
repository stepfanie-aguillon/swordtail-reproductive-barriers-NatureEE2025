# This includes scripts used on the HPC for ancestry analyses
#
# Authors: Aguillon SM, et al.
# Year: 2025
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: Nature Ecology and Evolution
# DOI: https://doi.org/10.1038/s41559-025-02669-9
#
# Edited date: 14 Nov 2023
#
# Please cite the paper if you use these scripts
#

## calculate read counts for sequencing files
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/count_reads_fastq_list.pl sample-list.txt


## calculate hybrid indices and heterozygosity for ancestry HMM results
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/parsetsv_ancestry_v2.pl ancestry-probs-par1_allchrs.tsv ancestry-probs-par2_allchrs.tsv > my-hybrid-indices.txt
