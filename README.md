## msgbsR: an R package for analysing methyla-tion sensitive genotyping by sequencing data  ###

Current data analysis tools do not fulfil all experimental designs. For example, GBS experiments using methylation sensitive restriction enzymes (REs), which is also known as methylation sensitive genotyping by sequencing (msGBS). msgbsR is an R package for the data analysis of msGBS experiments. Read counts and cut sites from a msGBS experiment can be read directly into the R environment from a sorted and indexed BAM file(s) to perform differential methylation analyses.

## Installation ##
To install in R:
```
install.packages("devtools")
library("devtools")
devtools::install_git("https://github.com/BenjaminAdelaide/msgbsR")
```

## Uses ##
msgbsR can be used to read data generated from a msGBS experiment directly into the R environment. The data can be used to identify methyalted sites within the genome and to conduct differential methylation analyses.

## Example ##
Please refer to the vignette within the msgbsR package containing a work flow of an example msGBS data set.
