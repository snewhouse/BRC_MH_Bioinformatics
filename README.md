Illumina Gene Expression Analysis Workflow
===========================================
Illumina Chip based Gene Expression Pipeline (QC, Normalization, Differential Expression)  

This is an example of a typical workflow we apply at the NIHR BRC-MH Bioinformatics unit.  

***Please note that this is under constant development and may need tweaking for your needs***

### Basic routines included in workflow 

- Genomestudio Final Reports read into ***LumiBatch***  using `lumi`
- ***QC Plots***: Boxplots, PCA Plots, hierarchical clustering, Coloured Dendograms...
- ***Probe Detection*** : based on expression level greater than the mean on the NEGATIVE beads
- ***Gender Checks*** : XIST Probe expression
- ***SampleNetwork Outliers*** : (based on : BMC Syst Biol. 2012 Jun 12;6:63. doi: 10.1186/1752-0509-6-63.
Network methods for describing sample relationships in genomic datasets: application to Huntington's disease.
Oldham MC1, Langfelder P, Horvath S. URL http://ccforum.com/1752-0509/6/63)
- Analysis of ***Batch Effects***
- Removal of ***Batch Effects*** using `sva`, `ComBat` or `rlm` PCA Batch regressions

More soon...
