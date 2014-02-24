#!/bin/bash

RMDFILE=GAP_gx_illumina_pipeline_preProcessing_step_001;

Rscript -e "require(knitr); require(markdown); knit('$RMDFILE.Rmd','$RMDFILE.md','$RMDFILE.tex'); markdownToHTML('$RMDFILE.md', '$RMDFILE.html', options=c('use_xhml'));"

pandoc -s GAP_gx_illumina_pipeline_preProcessing_step_001.html -o GAP_gx_illumina_pipeline_preProcessing_step_001.pdf



knit2html('$RMDFILE.Rmd')