Illumina Gene Expression Analysis Workflow
===========================================

***Authors:***  hamel.patel@kcl.ac.uk, stephen.newhouse@kcl.ac.uk

Illumina Chip based Gene Expression Pipeline (QC, Normalization, Differential Expression)  

This is an example of a typical workflow we apply at the NIHR BRC-MH Bioinformatics unit.  

***Please note that this is under constant development and may need tweaking for your needs***  

***NB:*** This is one way of doing it - that works well in our hands, and ensures you arent testing for association and building classification models with noise and unreliably detected probes!

*********

### Some Notes, references and papers : Please read
This isn’t an exhaustive list, but a place to start and get an overview of some of the important stuff when it comes to gene expression analysis, qc, pre-processing...ignore the "learn R" stuff if you already useR!..read the papers and try the tutorials...  

Quick R :- http://www.statmethods.net/  
Online learning RStudio: http://www.rstudio.com/resources/training/online-learning/  
Using Bioconductor with Microarray Analysis: http://www.bioconductor.org/help/workflows/arrays/  
Du, P., Kibbe, W.A., Lin and S.M. (2008). “lumi: a pipeline for processing Illumina microarray.” Bioinformatics.  
Lin, S.M., Du, P., Kibbe and W.A. (2008). “Model-based Variance-stabilizing Transformation for Illumina Microarray Data.” Nucleic Acids Res.  
Du, P., Kibbe, W.A., Lin and S.M. (2007). “nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays.” Biology Direct.  
Smyth GK (2005). “Limma: linear models for microarray data.” In Gentleman R, Carey V, Dudoit S, Irizarry R and Huber W (eds.), Bioinformatics and Computational Biology Solutions Using R and Bioconductor, pp. 397–420. Springer, New York.
Lumi:- http://bioconductor.org/packages/2.0/bioc/vignettes/lumi/inst/doc/lumi.pdf  
Xie Y (2010). MBCB: MBCB (Model-based Background Correction for Beadarray). R package version 1.18.0, http://www.utsouthwestern.edu.  
WGCNA: an R package for weighted correlation network analysis : http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/index.html  
Tutorials for the WGCNA package : http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/  
WGCNA Background and glossary: http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-00-Background.pdf  
Integrating Genetic and Network Analysis to Characterize Genes Related to Mouse Weight : http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.0020130  
Getting Started in Gene Expression Microarray Analysis:http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000543  
Network methods for describing sample relationships in genomic datasets: application to Huntington’s disease: http://ccforum.com/1752-0509/6/63  
Discovering statistically significant pathways in expression profiling studies: http://www.pnas.org/content/102/38/13544.full?maxtoshow=&HITS=10&hits=10&RESULTFORMAT=1&author1=tian&andorexacttitle=and&andorexacttitleabs=and&andorexactfulltext=and&searchid=1&FIRSTINDEX=0&sortspec=relevance&resourcetype=HWCIT
Removing Batch Effects in Analysis of Expression Microarray Data: An Evaluation of Six Batch Adjustment Methods: http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0017238  
Adjusting batch effects in microarray expression data using empirical Bayes methods. : http://biostatistics.oxfordjournals.org/content/8/1/118.abstract  
Tackling the widespread and critical impact of batch effects in high-throughput data: http://www.nature.com/nrg/journal/v11/n10/abs/nrg2825.html  
Comparison of normalization methods for Illumina BeadChip HumanHT-12 v3: http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3091625&tool=pmcentrez&rendertype=abstract  
RIN: http://en.wikipedia.org/wiki/RNA_integrity_number  
PCA: http://en.wikipedia.org/wiki/Principle_components_analysis  
XIST: http://en.wikipedia.org/wiki/XIST_%28gene%29  

**NB1**:our workflow..This is one way of doing it, that works well in my hands...it ensures good clean data and picks out potential problematic samples/issues....those with a little R skill are free to edit and adjust as needed. This isnt a "one shoe fits all" workflow and folks are free to edit/adjust all steps/thresholds as needed...just document the hows and whys and edits when you come to sharing the data and any and all changes and steps made when processing your data.

**NB2**: Small data sets where you have <10 samples per group of interest tend to break the workflow., same with missing data....so....Preaching to the converted  ( I hope ) :- explore your data before blindly stepping through any workflow/pipeline/analysis...box plots, histograms, missing rates, skewness, outliers, correlations, pca, means, mins, max, frequency tables, counts per group, variance, zero variance and small (<10) sample groups/predictors...and the list goes on and on...then adjust your analyses accordingly 


*********

illumina_gene_expression_workflow_preProcessing.Rmd
-----------------------------------------------------

- Use this [illumina_gene_expression_workflow_preProcessing.Rmd](https://github.com/snewhouse/BRC_MH_Bioinformatics/tree/master/Illumina_expression_workflow/illumina_gene_expression_workflow_preProcessing.Rmd)
- check for updates/changes on a regular basis
- Do not blindly apply this workflow. 
- Step through it manually and edit to your needs  
- This is set up for best use through [RStudio](http://www.rstudio.com/)  
- The script file format is written in [R Markdown (*.Rmd)](https://www.rstudio.com/ide/docs/r_markdown)  
- Be smart with the input [data](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#dos-and-donts-for-data)  
- This works best for studies with at least 10 samples per group/phenotype
- For small sample sizes ie < 10, then manually step through this script using [RStudio](http://www.rstudio.com/), and edit where needed  

Required R Pacakges
-----------------------
```
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)
library(gplots)
library(limma)
library(vsn)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)
library(plyr)
library(ggplot2)
library(gdata)
library(sva)
library(pamr)
library(glmnet)
library(snow)
library(parallel)
library(foreach)
library(doParallel)
library(caret)
```

Required Input Files
----------------------

1. [Group_Probe_Profile_FinalReport.txt](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#4-creating-final-reports)  
2. [Control_Probe_Profile_FinalReport.txt](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#4-creating-final-reports)   
3. [Sample_Table_FinalReport.txt](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#4-creating-final-reports)    
4. [Phenotype file](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#phenotype-file-format)    
5. [Batch information file](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#batch-file-format)      

Files [1-3](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#making-lumi-input-files-from-genomestudio) are generated in Genometudio.   
Files [4-5](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#additional-files) are provided by the user


******

### Basic routines included in workflow 

- Genomestudio Final Reports read into ***LumiBatch***  using `lumi`
- ***QC Plots***: Boxplots, PCA Plots, hierarchical clustering, Coloured Dendograms...
- ***Background Correction***: [MBCB (Model-based Background Correction for Beadarray)](http://www.bioconductor.org/packages/release/bioc/html/MBCB.html)
- ***Probe Detection*** : based on expression level greater than the mean on the [NEGATIVE beads](https://github.com/snewhouse/Misc/blob/master/Dont_Trust_Illumina_Detection_P_Values.md#expression-arrays--dont-trust-illumina-detection-p-values-at-th001)
- ***Gender Checks*** : based on XIST Probe expression
- ***SampleNetwork Outliers*** : (based on : BMC Syst Biol. 2012 Jun 12;6:63. doi: 10.1186/1752-0509-6-63.
Network methods for describing sample relationships in genomic datasets: application to Huntington's disease.
Oldham MC1, Langfelder P, Horvath S. URL http://ccforum.com/1752-0509/6/63)
- Analysis of ***Batch Effects***: PCA `prcomp` using non-detected probes (optional:plus Housekeeping probes)
- Removal of ***Batch Effects*** using `sva`, `ComBat` or `lm` PCA Batch regressions or a combination of all

**NB:*** This is one way of doing it - that works well in our hands, and ensures you arent testing for association and building classification models with noise and unreliably detected probes!

******

### Summary of basic worklfow

1. raw data
2. background correct
3. probe detection rates 
4. gender checks
5. transformation and normalisation
6. sample group outliers
7. Identify bacth effects
8. Remove batch effects
9. final clean data

******

### Additional Files

Along with the output from Genomestudio (see below: [Making Lumi input files from Genomestudio](https://github.com/snewhouse/BRC_MH_Bioinformatics/blob/master/Illumina_expression_workflow/README.md#making-lumi-input-files-from-genomestudio)), the worklfow requires the user to provide two additional files.  

1) Phenotype file  
2) Batch information file  

#### Phenotype File Format

See [pheno_info.txt](https://raw.githubusercontent.com/snewhouse/BRC_MH_Bioinformatics/master/Illumina_expression_workflow/pheno_info.txt) as an example.

- Required fields   
-- Sample.ID : Must match the chip id from the illumina array eg 	9020374058_A  
-- GROUPS  : Biological groups of interest. May be identical to phenotype codes. Examples below.  
-- SEX  : Gender (MALE or FEMALE)  
-- TISSUE : Tissue type (BLOOD, BRAIN etc)  
-- PHENOTYPE : Primary Phenotype of interest. CASE, CONTROL etc  
-- Study_ID : Study ID  

******

##### GROUPS

- How you set up the ***GROUPS*** field is very important.
- Make sure you hav at least 5 samples per group or the stats is a bit weird and may break
- GROUPS should represent some sensible biology or natural groupings of data
- GROUPS is used in the sampleNetwork outliers detection - this is used to qc your sample data/phenotype and tissue groups, with the aim of making a set of *homogenous* samples per GROUP/PENOTYPE
- This helps prevent problems that sample heterogenity can cause
- The end result is that for e.g. your CASES are more like each other and *funny* outlier samples are removed

*** Some examples :- ***  
*Mulitple Tissuese and Multiple phenotypes*
- 2 tissues, 2 disease groups.
- GROUPS:- tissue_1_disease_group_1, tissue_2_disease_group_1, tissue_1_disease_group_2, tissue_2_disease_group_2  

*Case/Control* 
- GROUPS:- cases, controls

*Multiple Phenotypes*
- GROUPS:- phenotype_1, phenotype_2, phenotype_3

******

#### Batch File Format

See [batch_info.txt](https://raw.githubusercontent.com/snewhouse/BRC_MH_Bioinformatics/master/Illumina_expression_workflow/batch_info.txt) as an example.  
This is a file with any and all known *batch* or *technical* data related to sample processing eg  
- RIN
- Date of RNA extraction
- Date Chip Run
- Processing Batch
- RNA concentrations
- .......

******

### Do's and Dont's for data
- ***Stick to these formats or the workflow will break!***  
- ***USE UPPER CASE FOR GROUPS, SEX, PHENOTYPE****  
- ***DONT MIX NUMBERS WITH TEXT IN THE SAME FIELD EG 1 and 1 year***  
- ***MISSING DATA : USE "NA" OR "" IE LEAVE BLANK FOR MISSING DATA***  


******


Making Lumi input files from Genomestudio
==============================================

## 1. Creating a new gene expression project in GenomeStudio

- Select “File” > “New Project” > “Gene Expression” and follow the GenomeStudio Project Wizard.
- Select “Direct Hyb”.
- Locate the “GenomeStudio_project” folder created in step one under “Project Repository”.
- Create a “project name” as [PROJECT_NAME]_expression_[DATE]_01.
- In the next window, under “Repository” specify the location of the “Data” folder created in step one. Once selected this should show all available chips within the “Sentrix Array Products” sub-window.
- Select and highlight the chip barcodes which require processing. Use the shift key to select multiple barcodes.
- Use the “Add Selected Samples to Project” icon to add the chip barcodes into the “Project Data” window. Click “Next”.
- Under “Groupset” repeat the project name followed by “raw” e.g project_expression_131101_01_raw.
- In the “Sentrix Array Products” sub-window select and highlight the chip barcodes which require processing.
- Use the “Add Selected Samples to Project” icon to add the chip barcodes into the “Project Data” window. Click “Next”.
- Under the “Name” tab select “Default” and click “Finish”.
- A message window will appear stating: “GenomeStudio detected that some samples have missing bead types. Would you like to impute missing data?” select “No”. Poor performing samples will be removed at a later step, after which the samples can be imputed.

***NOTE:*** The data will take a long time to load, check using Windows Task Manager (Performance panel) that at least one core is maxed out (otherwise it could mean GS has hung unexpectedly, restart the whole process if that's the case

*******

## 2.  Initial quality control

- The screen will present you next with the various tables of data. Select "Sample Table"
- Click on the Scatter Graph, then look at the following plots
- x=index, y=P05, labels=Sample Id
- x=index, y=Average Signal, labels=Sample Id
- x=index, y=Detected Genes, labels=Sample Id
- In each case note any samples which deviate substantially from the main chord of data points, some discretion is required here. Refer back to the gene expression log spread sheet to investigate deviated samples. This spread sheet is created by the geneticist conducting the gene expression assay, and can indicate samples of low amplification, gene expression failure or any technical problems associated with the assay.

*******

## 3. Excluding samples and imputing

- Select “Analysis” > “Manage Group Sets...”
- Under “Groupset” create a project name in the format of: [PROJECT_NAME]_expression_[DATE]_02 e.g “project_expression_130111_02”
- From the “Sentrix Array Products” sub window select and highlight the chip barcodes which require processing.
- Use the “Create a group for each selected sample” icon to add the chip barcodes into the “Project Groups” window. This will show each individual sample from the chips selected. Using the “Remove selected groups and samples from the project” icon, remove all samples which have failed the initial quality control step. Select “Next”
Under the “Name” tab select “Default”. Click “Finish”.
- A message window will appear stating: “GenomeStudio detected that some samples have missing bead types. Would you like to impute missing data?” select “Yes”.

*******

## 4. Creating Final Reports

***Group_Probe_Profile_Final_Report.txt***

- Select “Analysis” > “Reports” > “Final Report”
- Browse” and locate the “GenomeStudio_Final_Report” folder created in step 1 and name the file to be created as: “Group_Probe_Profile_Final_Report.txt”. Click “OK”
- From the “Tables” sub-window check the “Group Probe Profile” box
- From the “Columns” sub window check all boxes. Use the shift key to check multiple boxes.
- From the “Subcolumns” sub window check all boxes. Click “OK”.
- You could use "Sample Probe Profile" depending on how your data us set up

***Control_Probe_Profile_Final_Report.txt***

- Select “Analysis” > “Reports” > “Final Report”
- Browse” and locate the “GenomeStudio_Final_Report” folder created in step 1 and name the file to be created as: “Control_Probe_Profile_Final_Report.txt”. Click “OK”
- From the “Tables” sub-window check the “Control Probe Profile” box
- From the “Columns” sub window check all boxes. Use the shift key to check multiple boxes.
- From the “Subcolumns” sub window check all boxes. Click “OK”.

***Samples_Table_Final_Report.txt***

- Select “Analysis” > “Reports” > “Final Report”
- Browse” and locate the “GenomeStudio_Final_Report” folder created in step 1 and name the file to be created as: - “Samples_Table_Final_Report.txt”. Click “OK”
- From the “Tables” sub-window check the “Samples Table” box
- From the “Columns” sub window check all boxes except those beginning with “ERCC”. Click “OK”.

*******

## 5 Make Lumi Input

To make the input file for `lumiR()` the three final report file generated in step 4 above need to be merged.

### Input files

1 Group_Probe_Profile_FinalReport.txt  
2 Control_Probe_Profile_FinalReport.txt  
3 Sample_Table_FinalReport.txt  

***example script for generating lumiR input***

```
# Get Control Data
awk 'NR >7' Control_Probe_Profile_FinalReport.txt > control.probe;

# Get Sample Data
awk 'NR >7' Sample_Table_FinalReport.txt > sample.table;

# Concatenate All Data 
cat \
Group_Probe_Profile_FinalReport.txt \
control.probe \
sample.table > lumiR_input.txt;
```
