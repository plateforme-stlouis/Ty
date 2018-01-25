#Script Josh original
ProcessDataset = function(inputSet){
  AB15_sub = cbind(inputSet[seq(1,nrow(inputSet), 2),c(-1,-6)], inputSet[seq(2,nrow(inputSet), 2),c(-1,-6)]) #Merge paired reads in to a single line, ABN: Takes 1 line out of two in "inputset" and removes the first and last columns
  AB15_sub = AB15_sub[which(t(sapply(AB15_sub[,1], intToBits))[,3] != 1),] #Remove unmapped reads
  AB15_dim = t(apply(AB15_sub, 1, function(x){return(c(chr1=x[2],chr2=x[6], minPos=min(x[c(3,4,7,8)]), maxPos=max(x[c(3,4,7,8)])))})) #Get dataframe with full range that the DNA fragments mapped to
  AB15_sub = AB15_sub[!duplicated(AB15_dim),] #Remove duplicate fragments mapping to exactly the same chromosomal sequence
  AB15_sub = AB15_sub[,1:4] #Keep only info for the R1 read
  AB15_sub = cbind(AB15_sub, data.frame(pos=AB15_sub[,3], orient='+', Weight=1, stringsAsFactors =FALSE)) # Set placeholder values
  AB15_sub_flags = t(sapply(AB15_sub[,1], intToBits)) #Create a binary matrix of all FLAG info
  AB15_sub[which(AB15_sub_flags[,5] == 1), c('pos','orient')] = data.frame(pos = AB15_sub[which(AB15_sub_flags[,5] == 1), 4], orient='-', stringsAsFactors =FALSE) #For all alignments in reverse orientation, mark the orientation as "-"
  names(AB15_sub) = c("flags", "chr", "minPos", "maxPos", "pos", "orient", "weight") # Assign names to columns of data.frame
  print(unique(AB15_sub$flags))
  print(dim(AB15_sub)) # print dimensions for verification and integration count
  AB15_sub = summaryBy(weight~pos+orient+chr, AB15_sub, FUN=sum, keep.names=TRUE) #Summarize Weight based on unique Position,Orientation,Chromosome combinations.
  levels(AB15_sub[,3]) = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) # Convert Chromosome IDs to Chromosome Numbers (This works for these files, this won't necessarially work on other files)
  return(AB15_sub)
}


#Les valeurs calculées pour minPos et maxPos ne sont pas correctes, pas besoin pour le calcul du paramètre "weight" après apparement.


#Script Josh modif Amna

new_JDR4 <- subset(mini_JDR4, X2 %in% c(83,163,99,147))
JDR4_sub = cbind(new_JDR4[seq(1,nrow(new_JDR4), 2),c(-1,-6)], new_JDR4[seq(2,nrow(new_JDR4), 2),c(-1,-6)]) #Merge paired reads in to a single line,  Takes 1 line out of two in "inputset" and removes the first and last columns
JDR4_dim = t(apply(JDR4_sub, 1, function(x){return(c(chr1=x[2],chr2=x[6], minPos=min(x[c(3,4,7,8)]), maxPos=max(x[c(3,4,7,8)])))})) #Get dataframe with full range that the DNA fragments mapped to
JDR4_sub = JDR4_sub[!duplicated(JDR4_dim),] #Remove duplicate fragments mapping to exactly the same chromosomal sequence
JDR4_sub = JDR4_sub[,1:6] #Keep only info for the R1 read
JDR4_sub = JDR4_sub$X5 <- NULL
JDR4_sub = JDR4_sub$X7 <- NULL
JDR4_sub = cbind(JDR4_sub, data.frame(pos=JDR4_sub[,3], orient='+', Weight=1, stringsAsFactors =FALSE)) # Set placeholder values
JDR4_sub_flags = t(sapply(JDR4_sub[,1], intToBits)) #Create a binary matrix of all FLAG info
JDR4_sub[which(JDR4_sub_flags[,5] == 1), c('pos','orient')] = data.frame(pos = JDR4_sub[which(JDR4_sub_flags[,5] == 1), 4], orient='-', stringsAsFactors =FALSE) #For all alignments in reverse orientation, mark the orientation as "-"
names(JDR4_sub) = c("flags", "chr", "minPos", "maxPos", "pos", "orient", "weight") # Assign names to columns of data.frame
print(unique(JDR4_sub$flags))
print(dim(JDR4_sub)) # print dimensions for verification and integration count
library(doBy)
JDR4_sub = summaryBy(weight~pos+orient+chr, JDR4_sub, FUN=sum, keep.names=TRUE) #Summarize Weight based on unique Position,Orientation,Chromosome combinations.




#Script Amna

goodFlag <- 99
reverseFlag <- 83
new_JDR4 <- subset(mini_JDR4, X2 %in% c(83,99))
new_JDR4 <- new_JDR4[,2:4]
names(new_JDR4) <- c("flag", "chr", "pos")
new_JDR4 <- cbind(new_JDR4,data.frame(orient="+",weight=1,stringsAsFactors = FALSE))
new_JDR4[bitwAnd(new_JDR4$flag,reverseFlag)==reverseFlag,"orient"]<- "-"
new_JDR4 <- aggregate(formula = weight ~ chr+pos+orient , data = new_JDR4, FUN = sum)


#Je retrouve (presque) les mêmes positions que Josh mais pas exactement les même valeurs pour weight....
