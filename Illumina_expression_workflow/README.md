Illumina Gene Expression Analysis Workflow
===========================================

***Authors:***  amos.folarin@kcl.ac.uk, hamel.patel@kcl.ac.uk, stephen.newhouse@kcl.ac.uk

Illumina Chip based Gene Expression Pipeline (QC, Normalization, Differential Expression)  

This is an example of a typical workflow we apply at the NIHR BRC-MH Bioinformatics unit.  

***Please note that this is under constant development and may need tweaking for your needs***

### Basic routines included in workflow 

- Genomestudio Final Reports read into ***LumiBatch***  using `lumi`
- ***QC Plots***: Boxplots, PCA Plots, hierarchical clustering, Coloured Dendograms...
- ***Background Correction***: [MBCB (Model-based Background Correction for Beadarray)](http://www.bioconductor.org/packages/release/bioc/html/MBCB.html)
- ***Probe Detection*** : based on expression level greater than the mean on the NEGATIVE beads
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

Along with the output from Genomestudio (see below: Making Lumi input files from Genomestudio), the worklfow requires the user to provide two additional files.  

1) Phenotype file
2) Batch information file

#### Phenotype File Format

See [pheno_info.txt] as an example.

- Required fields   
-- Sample.ID Must match the chip id from the illumina array eg 	9020374058_A  
-- GROUPS  Biological groups of interest. May be identical to phenotype codes. Examples below.  
-- SEX  Gender (MALE or FEMALE)
-- TISSUE  	Tissue type (BLOOD, BRAIN etc)
-- PHENOTYPE  	Primary Phenotype of interest. CASE, CONTROL etc
-- Study_ID  Study ID 

***GROUPS***
How you set up the ***GROUPS*** field is very important.

#### Batch File Format

See [batch_info.txt] as an example.  
This is a file with any and all known *batch* or *technical* data related to sample processing eg  
- RIN
- Date of RNA extraction
- Date Chip Run
- Processing Batch
- RNA concentrations
- .......


### Do's and Dont's
***Stick to these formats or the workflow will break!***
***USE UPPER CASE FOR GROUPS, SEX, PHENOTYPE****
***DONT MIX NUMBERS WITH TEXT IN THE SAME FIELD***
***USE "NA" OR "" IE LEAVE BLANK FOR MISSING DATA***


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

***Sample_Probe_Profile_Final_Report.txt***

- Select “Analysis” > “Reports” > “Final Report”
- Browse” and locate the “GenomeStudio_Final_Report” folder created in step 1 and name the file to be created as: “Sample_Probe_Profile_Final_Report.txt”. Click “OK”
- From the “Tables” sub-window check the “Sample Probe Profile” box
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
