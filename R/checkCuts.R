library(bedr)

checkCuts <- function(cutSites, fastaPath, seq){
  
  # Firstly turn the cutSites into the appropiate format to be sorted
  cutSites <- paste(cutSites[ ,1], ':', cutSites[ ,2], '-', cutSites[ ,3], sep='')
  
  # Sort the cutSites  
  cutSites_sorted <- bedr.sort.region(cutSites)
  
  # Use the cutSites_sorted to get the sequence of each cut site along with the supplied fasta file
  sequences <- get.fasta(cutSites_sorted, fasta=fastaPath)
  
  # Filter out cutSites that do not match the seq given 
  sequences <- sequences[which(sequences[ ,2] == seq | sequences[ ,2] == toupper(seq)), ]
  
  # Return the correct sequences
  return(sequences)
  
}









