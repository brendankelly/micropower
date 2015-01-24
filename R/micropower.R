#' @title Simulate a sample OTU vector
#' @description Generates a labelled uniform vector from specified number of OTUs.
#' @details For weighted analysis, specify a sequence depth greater than 1; for unweighted analyses, sequence depth is 1 (default).
#' @param otu_number number of OTUs
#' @param sequence_depth number of sequence counts per OTU bin
#' @return numeric vector labelled with OTU names
#' @seealso \code{\link{simSampU}}
#' @export
#' @examples
#' simSamp(100,10)
simSamp <- function(otu_number=1000,sequence_depth=1) {
  s <- structure(.Data=rep(sequence_depth,otu_number),.Names=paste0("OTU",as.character(seq(otu_number))))
  return(s)
}


#' @title Rarefy a simulated OTU vector
#' @description Performs random subsampling without replacement from simulated OTU vector.
#' @param otu_vector numeric vector labelled with OTU names
#' @param rare_depth proportion of sequence counts to retain after subsampling
#' @return numeric vector labelled with OTU names
#' @seealso \code{\link{simSampU}}, \code{\link{simSampW}}
#' @export
#' @examples
#' rareSamp(simSamp(100,10),0.6)
rareSamp <- function(otu_vector,rare_depth=0.5) {
  r <- as.vector(rrarefy(otu_vector,rare_depth*sum(otu_vector)))
  r <- structure(.Data=r,.Names=names(otu_vector))
  return(r)
}


#' @title Calculate unweighted Jaccard distance between two OTU vectors
#' @description Wrapper for \code{\link{vegan}} package distance function, to calculate unweighted Jaccard distance.
#' @import vegan
#' @param sampleA numeric vector labelled with OTU names
#' @param sampleB numeric vector labelled with OTU names
#' @return numeric unweighted Jaccard distance
#' @seealso \code{\link{simSamp}}, \code{\link{rareSamp}}
#' @export
#' @examples
#' calcUJsample(rareSamp(simSamp(100),0.5),rareSamp(simSamp(100),0.5))
calcUJsample <- function(sampleA=rareSamp(simSamp(),0.5),sampleB=rareSamp(simSamp(),0.5)) {
  d <- vegdist(as.data.frame(rbind(sampleA,sampleB)),method="jaccard",binary=T)
  return(as.numeric(d))
}


#' @title Calculate weighted Jaccard distance between two OTU vectors
#' @description Wrapper for \code{\link{vegan}} package distance function, to calculate weighted Jaccard distance.
#' @import vegan
#' @param sampleA numeric vector labelled with OTU names
#' @param sampleB numeric vector labelled with OTU names
#' @return numeric weighted Jaccard distance
#' @seealso \code{\link{simSamp}}, \code{\link{rareSamp}}
#' @export
#' @examples
#' calcWJsample(rareSamp(simSamp(100,10),0.5),rareSamp(simSamp(100,10),0.5))
calcWJsample <- function(sampleA=rareSamp(simSamp(,10),0.5),sampleB=rareSamp(simSamp(,10),0.5)) {
  d <- vegdist(as.data.frame(rbind(sampleA,sampleB)),method="jaccard")
  return(as.numeric(d))
}


#' @title Simulate OTU table for presence-absence or abundance-weighted analysis
#' @description Incorporates within-group distance simulation from OTU number and depth of subsampling (rarefying) with group-level differences modeled by segregating OTU membership between groups.
#' @details For weighted analysis, specify a sequence depth greater than 1; for unweighted analyses, sequence depth is 1 (default).
#' @param group_size_vector numeric vector representing subjects per exposure/intervention group
#' @param otu_number number of simulated OTUs
#' @param sequence_depth number of sequence counts per OTU bin
#' @param rare_depth proportion of sequence counts to retain after subsampling
#' @param effect proportion of unique community membership in affected group of subjects
#' @return two-dimensional matrix OTU table, with row and column names to suit downstream analysis
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' simStudy(c(16,16,16),100,10,0.8,0.1)
simStudy <- function(group_size_vector=c(100,100,100),otu_number=1000,sequence_depth=1,rare_depth=0.4,effect=0.2) {
  otus <- simSamp(otu_number/(1+effect),sequence_depth)
  rare <- structure(lapply(as.list(seq(sum(group_size_vector))),FUN=function(x) {rareSamp(otus,rare_depth)}),.Names=paste0("s",as.character(seq(sum(group_size_vector)))))
  groups <- unlist(mapply(FUN=function(a,b) {rep(paste0("g",a),b)},a=as.character(seq_along(group_size_vector)),b=group_size_vector))
  names(rare) <- paste0(groups,names(rare))
  effect_otus <- sample(seq(length(otus)),round(effect*length(otus)),replace=FALSE)
  effect_group <- sample(seq_along(group_size_vector),1)
  effect_samples <- grep(paste0("g",as.character(effect_group)),sapply(strsplit(names(rare),"s"),FUN=function(x) {x[1]}))
  effected <- lapply(rare[effect_samples],FUN=function(x) {names(x)[effect_otus] <- paste0("X",names(x)[effect_otus]);x})
  rare[effect_samples] <- effected
  all_otus <- unique(do.call(c,lapply(rare,names)))
  empty <- structure(.Data=lapply(names(rare),FUN=function(x) {structure(.Data=rep(0,length(all_otus)),.Names=all_otus)}),.Names=names(rare))
  study <- mapply(FUN=function(e,r) {e[names(r)] <- r;e},e=empty,r=rare)
  return(study)
}


#' @title Simulate OTU table with null group-level effect for presence-absence or abundance-weighted analysis
#' @description Within-group distances simulated from OTU number and depth of subsampling (rarefying); no group-level differences.
#' @details For weighted analysis, specify a sequence depth greater than 1; for unweighted analyses, sequence depth is 1 (default).
#' @param group_size_vector numeric vector representing subjects per exposure/intervention group
#' @param otu_number number of simulated OTUs
#' @param sequence_depth number of sequence counts per OTU bin
#' @param rare_depth proportion of sequence counts to retain after subsampling
#' @return two-dimensional matrix OTU table, with row and column names to suit downstream analysis
#' @seealso \code{\link{simStudy}}, \code{\link{simPower}}
#' @export
#' @examples
#' simNull(c(16,16,16),100,10,0.8)
simNull <- function(group_size_vector=c(100,100,100),otu_number=1000,sequence_depth=1,rare_depth=0.5) {
  n <- simStudy(group_size_vector,otu_number,sequence_depth,rare_depth,0)
  n <- n[,grep("g1",colnames(n))]
  n <- lapply(seq_along(group_size_vector),FUN=function(x) {x <- n})
  n <- do.call(cbind,n)
  colnames(n) <- paste0("g",as.character(rep(seq_along(group_size_vector),each=group_size_vector[1])))
  colnames(n) <- paste0(colnames(n),"s",as.character(rep(seq(group_size_vector[1]),length(group_size_vector))))
  return(n)
}


#' @title Simulate a list of OTU tables encoding a range of effect sizes for presence-absence or abundance-weighted analysis
#' @description Extends the \code{\link{simStudy}} function to generate a list of OTU tables according to a specified range of group-level effects.
#' @details For weighted analysis, specify a sequence depth greater than 1; for unweighted analyses, sequence depth is 1 (default).
#' @param group_size_vector numeric vector representing subjects per exposure/intervention group
#' @param otu_number number of simulated OTUs
#' @param sequence_depth number of sequence counts per OTU bin
#' @param rare_depth proportion of sequence counts to retain after subsampling
#' @param effect_range range of proportions of unique community membership in affected group of subjects
#' @return list of two-dimensional-matrix OTU tables, with row and column names to suit downstream analysis
#' @seealso \code{\link{simStudy}}, \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' simPower(c(16,16,16),100,10,0.8,seq(0,0.3,length.out=100))
#' sapply(simPower(c(16,16,16),100,10,0.8,seq(0,0.3,length.out=100)),FUN=function(x) {calcOmega2(calcWJstudy(x))})
simPower <- function(group_size_vector=c(100,100,100), otu_number=1000, sequence_depth=1, rare_depth=0.5, effect_range=seq(0,0.3,length.out=100)) {
  p <- structure(.Data=as.list(effect_range),.Names=as.character(effect_range))
  p <- lapply(p,FUN=function(x) {simStudy(group_size_vector,otu_number,sequence_depth,rare_depth,x)})
  n <- structure(.Data=list(simNull(group_size_vector,otu_number,sequence_depth,rare_depth)),.Names="null")
  p <- c(n,p)
  return(p)
}


#' @title Simulate a list of OTU tables to create a hash table relating depth of subsampling (rarefying) to mean pairwise distance
#' @description Extends the \code{\link{simSamp}} function to generate a list of OTU tables representing different levels of subsampling.
#' @details For weighted analysis, specify a sequence depth greater than 1; for unweighted analyses, sequence depth is 1 (default).
#' @param rare_levels the levels of random subsampling (rarefying), as the proportion of sequence counts to retain after subsampling, to be tested
#' @param rep_per_level number of different levels (proportion of sequence counts to retain after subsampling) tested
#' @param otu_number number of simulated OTUs
#' @param sequence_depth number of sequence counts per OTU bin
#' @return list of two-dimensional-matrix OTU tables
#' @seealso \code{\link{simSamp}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' hashMean(,10,1000,10)
#' sapply(hashMean(runif(100,0,1),1,1000,10),FUN=function(x) {mean(lowerTriDM(calcWJstudy(x)))})
hashMean <- function(rare_levels=runif(1000,0,1),rep_per_level=1,otu_number=1000,sequence_depth=1) {
  r <- rep(rare_levels,each=rep_per_level)
  r <- structure(.Data=lapply(r,FUN=function(x) {cbind(sampleA=rareSamp(simSamp(otu_number,sequence_depth),x),sampleB=rareSamp(simSamp(otu_number,sequence_depth),x))}),.Names=as.character(r))
  names(r) <- paste0(names(r),"_index",as.character(seq_along(r)))
  return(r)
}


#' @title Simulate a list of OTU tables to create a hash table relating OTU number to standard deviation of pairwise distances
#' @description Extends the \code{\link{simSamp}} function to generate a list of OTU tables representing different OTU numbers.
#' @details For weighted analysis, specify a sequence depth greater than 1; for unweighted analyses, sequence depth is 1 (default).
#' @param rare_depth proportion of sequence counts to retain after subsampling
#' @param otu_number_range numeric vector, including the range of numbers of simulated OTUs to test
#' @param sim_number number of simulated subjects per OTU number
#' @param sequence_depth number of sequence counts per OTU bin
#' @return list of two-dimensional-matrix OTU tables
#' @seealso \code{\link{simSampW}}, \code{\link{calcDistW}}
#' @export
#' @examples
#' hashSD(0.5,c(10,100,1000),1000,10)
#' sapply(hashSD(),FUN=function(x) {sd(lowerTriDM(calcWJstudy(x)))})
hashSD <- function(rare_depth=0.5,otu_number_range=c(10,100,1000,10000),sim_number=100,sequence_depth=10) {
  o <- structure(.Data=as.list(otu_number_range),.Names=paste0("otu",as.character(otu_number_range)))
  o <- lapply(o,FUN=function(x) {replicate(sim_number,rareSamp(simSamp(x,sequence_depth),rare_depth))})
  o <- lapply(o,FUN=function(x) {colnames(x) <- paste0("s",as.character(seq(ncol(x))));x})
  return(o)
}


#' @title Calculate pairwise unweighted Jaccard distances from an OTU table
#' @description Generates a square matrix of pairwise unweighted Jaccard distances from an OTU table.
#' @import reshape2
#' @param otu_table a matrix with one column per subject and one row per OTU
#' @return square matrix of pairwise distances
#' @seealso \code{\link{simStudyU}}, \code{\link{calcDistW}}
#' @export
#' @examples
#' calcUJstudy(simStudy(c(16,16,16),100,,0.8,0.1))
calcUJstudy <- function(otu_table) {
  dm <- expand.grid(r=colnames(otu_table),c=colnames(otu_table),stringsAsFactors=FALSE)
  dm$d <- unlist(Map(function(r,c) {calcUJsample(otu_table[,r],otu_table[,c])},dm$r,dm$c),use.names=FALSE)
  dm <- acast(dm,r~c,value.var="d")
  return(dm)
}


#' @title Calculate pairwise weighted Jaccard distances from an OTU table
#' @description Generates a square matrix of pairwise weighted Jaccard distances from an OTU table.
#' @import reshape2
#' @param  otu_table a matrix with one column per subject and one row per OTU
#' @return square matrix of pairwise distances
#' @seealso \code{\link{simStudyU}}, \code{\link{calcDistW}}
#' @export
#' @examples
#' calcWJstudy(simStudy(c(16,16,16),100,10,0.8,0.1))
calcWJstudy <- function(otu_table) {
  dm <- expand.grid(r=colnames(otu_table),c=colnames(otu_table),stringsAsFactors=FALSE)
  dm$d <- unlist(Map(function(r,c) {calcWJsample(otu_table[,r],otu_table[,c])},dm$r,dm$c),use.names=FALSE)
  dm <- acast(dm,r~c,value.var="d")
  return(dm)
}


#' @title Simulate a phylogenetic tree with specified OTU names
#' @description Generates a random phylogenetic tree based upon the OTUs included in a single OTU table.
#' @details A wrapper for the \code{\link{rtree}} function included with the \code{\link{ape}} package.
#' @import ape
#' @param  otu_table a matrix with one column per subject and one row per OTU
#' @return An object of the class \code{\link{phylo}}.
#' @seealso \code{\link{simStudyU}}, \code{\link{simStudyW}}, \code{\link{rtree}}, \code{\link{write.tree}}
#' @export
#' @examples
#' simTreeTable(simStudy(c(16,16,16),100,,0.8,0.1))
simTreeTable <- function(otu_table) {
  if(is.vector(otu_table)) {otus <- names(otu_table)} else
  {otus <- unique(rownames(otu_table))}
  tree <- rtree(n=length(otus),rooted=F,tip.label=otus,br=rlnorm)
  return(tree)
}


#' @title Simulate a phylogenetic tree with specified OTU names
#' @description Generates a random phylogenetic tree based upon the OTUs included in a list of OTU tables.
#' @details A wrapper for the \code{\link{rtree}} function included with the \code{\link{ape}} package.
#' @import ape
#' @param  otu_table_list a list of OTU table matrices with one column per subject and one row per OTU
#' @return An object of the class \code{\link{phylo}}.
#' @seealso \code{\link{simPowerU}}, \code{\link{simPowerW}}, \code{\link{rtree}}, \code{\link{write.tree}}
#' @export
#' @examples
#' simTreeList(simPower(c(100,100,100),100,10,0.5,seq(0,0.3,length.out=100)))
simTreeList <- function(otu_table_list) {
  otus <- unique(do.call(c,lapply(otu_table_list,rownames)))
  tree <- rtree(n=length(otus),rooted=F,tip.label=otus,br=rlnorm)
  return(tree)
}


#' @title Export a simulated OTU table to permit analysis of non-Jaccard distances
#' @description Exports a simulated OTU table as a tab-separated file to permit processing for pairwise distances by bioinformatic pipelines.
#' @param  otu_table a matrix with one column per subject and one row per OTU
#' @return Writes a tab-seperated file to the working directory.
#' @seealso \code{\link{simStudyU}}, \code{\link{simStudyW}}
#' @export
#' @examples
#' writeOTUtable(simStudy(c(16,16,16),100,1,0.8,0.1),"otu_table_export")
writeOTUtable <- function(otu_table, otu_table_name) {
  otu_table <- cbind("#OTU ID" = rownames(otu_table), otu_table)
  filename <- paste0(otu_table_name, ".txt")
  write.table(otu_table, file=filename, append=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}


#' @title Export a list of simulated OTU tables to permit analysis of non-Jaccard distances
#' @description Exports a list of simulated OTU tables as tab-separated files to permit processing for pairwise distances by bioinformatic pipelines.
#' @param  otu_table_list a list of OTU table matrices with one column per subject and one row per OTU
#' @return Writes tab-seperated files to the working directory.
#' @seealso \code{\link{simPowerU}}, \code{\link{simPowerW}}
#' @export
#' @examples
#' writeOTUlist(simPower(c(100,100,100),100,1,0.5,seq(0,0.3,length.out=100)))
writeOTUlist <- function(otu_table_list) {
  Map(writeOTUtable, otu_table_list, as.list(names(otu_table_list)))
}


#' @title Read QIIME-formatted distance matrices
#' @description Imports a tab-separated distance matrix file, as output by the QIIME pipeline.
#' @details A wrapper for the \code{\link{read.delim}} function, tailored to the output from QIIME's beta_diversity.py script.
#' @param  filepath the path to the tab-separated distance matrix file
#' @return A square distance matrix.
#' @seealso \code{\link{writeOTUtable}}, \code{\link{readDMdir}}
#' @export
#' @examples
#' readDM()
readDM <- function(filepath=file.choose()) {
  dm <- as.matrix(read.delim(filepath,header=T,row.names=1))
  return(dm)
}


#' @title Read all QIIME-formatted distance matrices in working directory
#' @description Imports all tab-separated distance matrix files in the working directory, as output by the QIIME pipeline.
#' @details A wrapper for the \code{\link{read.delim}} function, tailored to the output from QIIME's beta_diversity.py script.
#' @return A list of square distance matrices.
#' @seealso \code{\link{writeOTUlist}}, \code{\link{readDM}}
#' @export
#' @examples
#' readDMdir()
readDMdir <- function() {
  files <- structure(.Data=as.list(list.files()),.Names=list.files())
  files <- lapply(files,readDM)
  return(files)
}


#' @title Read all QIIME-formatted distance matrices in working directory that are labelled with a specified distance metric
#' @description Imports all tab-separated distance matrix files in the working directory with names that include the specified distance metric, as output by the QIIME pipeline.
#' @details A wrapper for the \code{\link{read.delim}} function, tailored to the output from QIIME's beta_diversity.py script.
#' @return A list of square distance matrices.
#' @seealso \code{\link{writeOTUlist}}, \code{\link{readDMdir}}
#' @export
#' @examples
#' readDMmetric("weighted_normalized_unifrac")
readDMmetric <- function(metric="weighted_normalized_unifrac") {
  files <- structure(.Data=as.list(grep(metric,list.files(),value=T)),.Names=grep(metric,list.files(),value=T))
  files <- lapply(files,readDM)
  return(files)
}


#' @title Convert distance matrix subject names to group-level names.
#' @description Convert distance matrix subject names to group-level names to permit PERMANOVA testing and analysis of group-level effects.
#' @details Input distance matrix colnames and rownames must be formatted "g2s12" for group 2, subject 12.
#' @param dm a square distance matrix
#' @return A square distance matrix.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' groupNames(calcUJstudy(simStudy()))
groupNames <- function(dm) {
  dimnames(dm) <- lapply(dimnames(dm),FUN=function(x) {gsub("(s.*)","",x,perl=T)})
  return(dm)
}


#' @title Extract the pairwise distances from the lower triangle of a square distance matrix.
#' @description Produces a vector of pairwise distances from the lower triangle of a symmetric, square distance matrix.
#' @details A wrapper for the \code{\link{lower.tri}} function.
#' @param dm a square distance matrix
#' @return A vector of pairwise distances.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}, \code{\link{readDM}}
#' @export
#' @examples
#' lowerTriDM(calcUJstudy(simStudy()))
lowerTriDM <- function(dm) {
  l <- as.vector(dm[lower.tri(dm,diag=F)])
  return(l)
}


#' @title Calculate omega-squared.
#' @description The proportion of distance accounted for by the grouping factor, corrected for the mean-squared error.
#' @param dm a square distance matrix
#' @return A numeric value.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' calcOmega2(calcUJstudy(simStudy()))
calcOmega2 <- function(dm) {
  dm <- groupNames(dm)
  groups <- unique(rownames(dm))
  df <- length(groups) - 1
  dm_within <- structure(.Data=lapply(as.list(groups),FUN=function(x) {dm[rownames(dm)==x,colnames(dm)==x]}),.Names=groups)
  sst <- sum(lowerTriDM(dm)^2)/nrow(dm)
  ssw <- sum(sapply(dm_within,FUN=function(x) {sum(lowerTriDM(x)^2)/nrow(x)}))
  ssa <- sst-ssw
  group_ms <- ssa/df
  error_ms <- ssw/(nrow(dm) - length(groups))
  omega2 <- (ssa-(df*error_ms))/(sst+error_ms)
  return(omega2)
}


#' @title Perform PERMANOVA testing.
#' @description A wrapper for the \code{\link{adonis}} function included in the \code{\link{vegan}} package, restricted to a single grouping factor.
#' @details Does not permit specifying strata within which to constrain permutations; \code{\link{adonis}} should be used for this purpose.
#' @param dm a square distance matrix
#' @return Typical output for analysis of variance.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' PERMANOVA(calcUJstudy(simStudy()))
PERMANOVA <- function(dm) {
  dm <- groupNames(dm)
  dm <- adonis(as.dist(dm)~colnames(dm),permutations=1000)
  return(dm)
}


#' @title Calculate the coefficient of determination (R-squared).
#' @description The proportion of distance accounted for by the grouping factor.
#' @details A wrapper for the \code{\link{PERMANOVA}} function, which itself utilizes the \code{\link{adonis}} function included in the \code{\link{vegan}} package.
#' @param perm output from the \code{\link{PERMANOVA}} function
#' @return A numeric value.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' calcR2(calcUJstudy(simStudy()))
calcR2 <- function(perm=PERMANOVA(dm)) {
  cod <- perm$aov.tab$R2[1]
  return(cod)
}


#' @title Calculate the PERMANOVA p-value.
#' @description The probability of the observed or greater effect under the null hypothesis, calculated by label permutation.
#' @details A wrapper for the \code{\link{PERMANOVA}} function, which itself utilizes the \code{\link{adonis}} function included in the \code{\link{vegan}} package.
#' @param perm output from the \code{\link{PERMANOVA}} function
#' @return A numeric value.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' calcPERMANOVAp(calcUJstudy(simStudy()))
calcPERMANOVAp <- function(perm=PERMANOVA(dm)) {
  p <- perm$aov.tab$'Pr(>F)'[1]
  return(p)
}


#' @title Take a bootstrap sample from a square distance matrix.
#' @description Random sample with replacement from a distance matrix, with resulting matrix specified by number of subjects per group.
#' @param dm a square distance matrix
#' @param subjects_group_vector vector with number of subjects in each group to sample.
#' @return A square distance matrix.
#' @seealso \code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
#' @export
#' @examples
#' bootDM(calcUJstudy(simStudy()),c(3,4,5))
bootDM <- function(dm,subject_group_vector) {
  groups <- sort(unique(rownames(groupNames(dm))))
  s <- Map(function(g,v) sample(grep(g,rownames(dm),value=TRUE),v,replace=TRUE),groups,subject_group_vector)
  s <- do.call(c,s)
  s <- dm[s,s]
  return(s)
}


#' @title Perform bootstrap power analysis on a list of square distance matrices.
#' @description Estimates the statistical power of PERMANOVA testing to detect the group-level effect in the given distance matrices, based upon bootstrap sampling.
#' @param dm_list a list of square distance matrices, with names
#' @param boot_number the number of bootstrap samples to perform on each distance matrix in the list
#' @param subjects_per_group number of subjects per group to sample.
#' @param alpha the threshold for PERMANOVA type I error
#' @return A data frame relating PERMANOVA power to effect size quantified by the coefficient of determination (R^2) and omega-squared.
#' @seealso \code{\link{simPower}}, \code{\link{bootDM}}
#' @export
#' @examples
#' bootPower(lapply(simPower(),calcUJstudy))
bootPower <- function(dm_list,boot_number=100,subjects_per_group=10,alpha=0.05) {
  e <- rep(names(dm_list),each=boot_number)
  simulated_omega2 <- rep(sapply(dm_list,calcOmega2),each=boot_number)
  dm <- lapply(dm_list,FUN=function(x) {lapply(seq(boot_number),FUN=function(y) {bootDM(x,subjects_per_group)})})
  o <- lapply(dm,FUN=function(x) {sapply(x,calcOmega2)})
  p <- lapply(dm,FUN=function(x) {lapply(x,PERMANOVA)})
  r <- lapply(p,FUN=function(x) {sapply(x,calcR2)})
  p <- lapply(p,FUN=function(x) {sapply(x,calcPERMANOVAp)})
  dm <- data.frame(effect=e,simulated_omega2=simulated_omega2,observed_omega2=do.call(c,o),observed_R2=do.call(c,r),p=do.call(c,p))
  dm <- ddply(dm,.(effect),here(transform),power=length(p[p<alpha])/length(p))
  dm$subjects_per_group <- rep(subjects_per_group,nrow(dm))
  dm$observed_omega2[dm$observed_omega2<0] <- 0
  dm$simulated_omega2[dm$simulated_omega2<0] <- 0
  return(dm)
}
