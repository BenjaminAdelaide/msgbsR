countMatrix = x; phenotype1 = y; phenotype2 = z

x <- data.frame(matrix(sample(0:100, size = 10000*10, replace = TRUE), nrow = 10000, ncol = 10))
y <- c(rep('A', 5), rep('B', 5))
z <- c(rep('C', 3), rep('D', 2), rep('E', 2), rep('F', 3))



plotCounts(countMatrix, phenotype1, phenotype2 = NULL){
    
  # Determine the total number of cuts sites per sample
  # cut sites with > 1 read produced for each sample
  cuts <- data.frame(t(data.frame(lapply(1:ncol(countMatrix),function(x) 
                      length(countMatrix[,x][which(countMatrix[,x] > 0)])))))
  
  # Determine the library size (total number of reads)
  libSize <- data.frame(t(data.frame(lapply(1:ncol(countMatrix),
                                       function(x) sum(countMatrix[,x])))))
  
  # Make a data frame out of the cuts and libSize
  datPlot <- as.data.frame(cbind(cuts[,1], libSize[,1]))
  datPlot[,3] <- phenotype1
  colnames(datPlot) <- c('cuts', 'libSize', 'phenotype1')
  
  # If else statment depending on whether or not there is a second phenotype characteristic 
  if(phenotype2 == NULL){
    
    qplot(x = libSize, y = cuts, colour = phenotype1, data = datPlot,
          xlab = 'Total number of Reads per sample',
          ylab = 'Total number of cut sites per sample',
          main = '')
    
  } else {
    
  datPlot[,4] <- phenotype2
  colnames(datPlot) <- c('cuts', 'libSize', 'phenotype1', 'phenotype2')
    
  qplot(x = libSize, y = cuts, colour = phenotype1, shape = phenotype2, data = datPlot,
        xlab = 'Total number of Reads per sample',
        ylab = 'Total number of cut sites per sample',
        main = '')
  colnames(datPlot) <- c('cuts', 'libSize')
  
  }

}

