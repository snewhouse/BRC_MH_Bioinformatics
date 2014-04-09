#!/bin/bash/

#----------------------------------------------------------------------------------------------#
# --- convertGen2RDose_v1.0.sh
# --- Stephen Newhouse
# --- stephen.newhouse@kcl.ac.uk
#----------------------------------------------------------------------------------------------#

# http://sourceforge.net/projects/fcgene/files/

#----------------------------------------------------------------------------------------------#
# Set PLINK exe                                                                                #
#----------------------------------------------------------------------------------------------#
FCGENE_EXE="/scratch/project/pipelines/Hamel/adn_hg18_hg19_impute/impute/imputation/anm_batch1_and_batch2_merged/anm_imputed1000G_31032014/fcgene-1.0.7/fcgene"

#----------------------------------------------------------------------------------------------#
# get options                                                                                  #
#----------------------------------------------------------------------------------------------#
gen_file=${1}
out_file=${2}

#----------------------------------------------------------------------------------------------#
# STANDARD GEN FORMAT:-                                                                        #
#----------------------------------------------------------------------------------------------#

echo ""
echo "Converting to dosage using FCGENE (http://sourceforge.net/projects/fcgene/files/)"
echo ""
echo "calling:- fcgene --gen ${1} --oformat r-dose --out ${2}"
echo ""

${FCGENE_EXE} --gens ${1} --oformat r-dose --out ${2};
echo ""


echo "rename as .dose"
mv ${2}_dose.txt ${2}_dose
echo ""

echo "gzipping" ${2}_dose
gzip ${2}_dose


# test_affstat.txt
# test_dose.txt
# test_fcgene.log
# test_pedinfo.txt
# test_snp.frq
# test_snpinfo.txt
