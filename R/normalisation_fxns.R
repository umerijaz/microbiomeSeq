#===== normalisation functions =========================#
#Each of the normalisation functions takes up a phloseq object and returns
#returns a physeq object whose otu-table is transformed.
#1. edgernorm
edgeRnorm = function(physeq, ...){
  
  abund_table <- otu_table(physeq)
  
  # Enforce orientation.
  if(!taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }
  x = as(abund_table, "matrix")
  # See if adding a single observation, 1, 
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify 
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems. 
  # Can the 1 be reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
  y = edgeR::DGEList(counts=x, remove.zeros=TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  otu_table(physeq) <- otu_table(z$counts, taxa_are_rows=TRUE)
  return(physeq)
}
#2. variance stabilisation
deseq_varstab = function(physeq, sampleConditions=rep("A", nsamples(physeq)), ...){
  
  abund_table <- otu_table(physeq)
  
  # Enforce orientation.
  if(!taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }
  
  x = as(abund_table, "matrix")
  # The same tweak as for edgeR to avoid NaN problems
  # that cause the workflow to stall/crash.
  x = x + 1
  
  cds = DESeqDataSetFromMatrix(x, DataFrame(sampleConditions), ~ 1)
  # First estimate library size factors
  cds = estimateSizeFactors(cds)
  # Variance estimation, passing along additional options
  cds = estimateDispersions(cds, ...)
  #vsmat = varianceStabilizingTransformation(cds)
  vsmat <- varianceStabilizingTransformation(cds)
  otu_table(physeq) <- otu_table(assay(vsmat), taxa_are_rows=TRUE)
  return(physeq)
}
#3.proportion
# Normalize total sequences represented 

# Scale by dividing each variable by its standard deviation.
#physeq = transform_sample_counts(physeq, function(x) x/sd(x))
# Center by subtracting the median
#physeq = transform_sample_counts(physeq, function(x) (x-median(x)))
proportion = function(physeq){
  normf = function(x, tot=max(sample_sums(physeq))){ tot*x/sum(x) }
  physeq = transform_sample_counts(physeq, normf)
  return(physeq)
}

#4.random sampling 
randomsubsample = function(physeq, smalltrim=0.15, replace=TRUE,meta=F){
  # Set the minimum value as the smallest library quantile, n`smalltrim` 
  samplemin = sort(sample_sums(physeq))[-(1:floor(smalltrim*nsamples(physeq)))][1]
  physeqr = rarefy_even_depth(physeq, samplemin, rngseed=TRUE,replace=replace, trimOTUs=TRUE)
  return(physeqr)
}

#5.relative transformation
relative <- function(physeq,norm.meta=F,select.variables=NULL){
  if(norm.meta){
    get.vars <- get.num.variables(physeq)
    
    norm.variables <- get.vars$num.variables/rowSums(get.vars$num.variables)
    
    meta_table <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
    
    sample_data(physeq) <- meta_table
  }
  else if(!norm.meta){
    abund_table <- otu_table(physeq)
    otu_table(physeq) <- abund_table/rowSums(abund_table)
  }
  return(physeq)
}
#6. log relative transformation
log_relative <- function(physeq, norm.meta=F, select.variables=NULL){
  if(norm.meta){
    get.vars <- get.num.variables(physeq)
    norm.variables <- log(get.vars$num.variables/rowSums(get.vars$num.variables))
    meta_table <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
    sample_data(physeq) <- meta_table
  }
  else if(!norm.meta){
    abund_table <- otu_table(physeq)
    otu_table(physeq) <- log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
  }
    return(physeq)
}

#7.
scale.meta <- function(physeq,type="scale",select.variables=NULL){
  
  get.vars <- get.num.variables(physeq)
  num.variables <- get.vars$num.variables
  
  type <- match.arg(type, c("scale","log","sqrt"))
  
  if(type=="scale"){
    norm.variables <- data.frame(scale(num.variables))
  }
  else if(type=="log"){
    norm.variables <- log2(num.variables)
  }
  else if(type=="sqrt"){
   norm.variables <- sqrt(num.variables)
  }
  
  sample_data(physeq) <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
  
  return(physeq)
}
