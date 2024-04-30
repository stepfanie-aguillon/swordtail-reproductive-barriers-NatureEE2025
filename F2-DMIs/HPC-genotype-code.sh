# This includes scripts used on the HPC for genotype/DMI analyses in F2s
#
# Authors: Aguillon SM, et al.
# Year: 2024
# Title: Pervasive gene flow despite strong and varied reproductive barriers in swordtails
#        
# Journal Info: TBD
# bioRxiv DOI: https://doi.org/10.1101/2024.04.16.589374
#
# Edited date: 14 Nov 2023
#
# Please cite the paper if you use these scripts
#

## approach to call genotypes from ancestry HMM results
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par1_allchrs.tsv ancestry-probs-par2_allchrs.tsv mygenotypes.txt

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/transpose_nameout.pl mygenotypes.txt

perl -pi -e 's/id/chr\tpos/g' mygenotypes.txt_transposed

perl -pi -e 's/:/\t/g' mygenotypes.txt_transposed


## calculate ancestry by site from a genotypes file
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/parse_genotypes_ancestry_bysite.pl mygenotypes_file.txt /home/groups/schumer/shared_bin/Lab_shared_scripts/


## pulling particular genotypes to assess the DMI
head -n 1 mygenotypes.txt_transposed > selected_DMI_mygenotypes.txt

grep chr-13 mygenotypes.txt_transposed > chr13-subset
grep 2112372 chr13-subset >> selected_DMI_mygenotypes.txt

grep chr-06 mygenotypes.txt_transposed > chr06-subset
grep 12520420 chr06-subset >> selected_DMI_mygenotypes.txt

grep mitochondria mygenotypes.txt_transposed > mtDNA-subset
grep 5159 mtDNA-subset >> selected_DMI_mygenotypes.txt

perl /home/groups/schumer/shared_bin/Lab_shared_scripts/transpose_nameout.pl selected_DMI_mygenotypes.txt


# pull just a subset of individuals
perl /home/groups/schumer/shared_bin/Lab_shared_scripts/grep_list_keep_header.pl sample_list.txt mygenotypes.txt mygenotypes_subset.txt 1
