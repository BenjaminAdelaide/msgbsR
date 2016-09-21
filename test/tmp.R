dat <- function (x)
{
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

tmp <- list.files('inst', pattern = '.bai')
length(dat(x = tmp))

dat(bamFiles)

# Get a list of all the bam files
bamFiles <- list.files(path=bamFilepath,
                       full.names=TRUE, pattern=".bam$")


library(R.utils)

tmpFun <- function(fasta){
print('zipping fasta file');
suppressWarnings(gzip(fasta))
print('Finished zipping');

}
tmpFun(fasta = 'inst/extdata/chr1.fa')


library(Rsamtools)
library(plyr)

datCounts <- rawCounts(bamFilepath = 'inst/extdata/', threads = 1)






