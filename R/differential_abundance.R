#'Differential abundance analysis
#'
#'This function finds the features that are significantly differentially abundant in the provided data,
#'using DESeq implementation which models taxa abundance as a negative binomial distribution. See \link[DESeq2]{DESeq}  for more details. The significance
#'of differentially abundant taxa is defined by log2 fold change and pvalue thresholds. These features are
#'then assigned importance using random forest classifer. The measure of importance used in this implementation
#'is mean decrease in accuracy.
#'
#' @param physeq (Required). A link[phyloseq]{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the.
#'        information, this should be one of the components of grouping vector.
#' @param output_norm (optional). A character string specifying method to be used for transforming abundance data to be used for plotting.
#'        note that, this normalisation occurs after DESeq analysis. Therefore, it is strictly for purposes of the output data
#'        and not for differential expression analysis.
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#' @param  pvalue.threshold. Cut off p-value for significance of differentially abundance taxa, default is 0.05.
#' @param  lfc.threshold. Threshold for log2 fold change over which significance of differentially expressed taxa is considered.
#' @return Returns a list of three items: \itemize{
#'         \item  SignfeaturesTable: A \code{data.frame} of taxa with raw and adjusted pvalues, basemean, log2 fold change and significance labels.
#'         computed by using raw p-values.
#'         \item  importance: A \code{data.frame} of taxa mean decrease in accuracy as obtained by random forest classifier.
#'         \item plotdata: A \code{data.frame} of taxa and corresponding corrected p-values, importance rank organised
#'         in form accepted by ggplot.
#'        }
#' @examples
#' data(pitlatrine)
#' physeq <- taxa_level(pitlatrine,"Phylum")
#' deseq_sig  <- differential_abundance(physeq, grouping_column = "Country")
#' #plot the significant features
#' plot_signif(deseq_sig$plotdata) #see function \link[microbiomeSeq]{plot_signif}
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export differential_abundance
#'
differential_abundance <- function(physeq, grouping_column,pvalue.threshold=0.05,lfc.threshold=0,
                          filename = "NB_significant",output_norm=NULL){

  abund_table <- otu_table(physeq)
  meta_table <-data.frame(sample_data(physeq))

  #==get count data to be used in deseq ==========#
  countData = round(as(abund_table, "matrix"), digits = 0)+1
  if(!taxa_are_rows(physeq)){
    countData = t(countData)
  }
  #==add a dummy column corresponding to grouping variable ==#
  meta_table$Groups <- meta_table[,grouping_column]
  #== create deseq compatible matrix and then test differential abundance in the groups=#
  dds <- DESeq2::DESeqDataSetFromMatrix(countData, meta_table, as.formula(~Groups))
  data_deseq_test = DESeq2::DESeq(dds)
  #==Extract results of the test =======#
  res = DESeq2::results(data_deseq_test, cooksCutoff = FALSE)
  res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))
  res_tax_sig = subset(res_tax, padj < pvalue.threshold & lfc.threshold < abs(log2FoldChange))
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
  res_tax$Significant[is.na(res_tax$Significant)] <- "No"

  #==normalising abundance data for output ==#
  if(!is.null(output_norm)){
    physeq <- normalise_data(physeq,output_norm)
    data <- as.data.frame(otu_table(physeq))
  }
  else{
    data <- abund_table
  }

  #==get importance measure of significant features using random forest
  subset.data<-data.frame(data[,as.character(res_tax[rownames(res_tax_sig),"OTU"])])
  rownames(res_tax_sig) <- colnames(subset.data) #enforce rownames of res_tax_sig to be the same as colnames of subset data for easy indexing during plotting.

  rf_res <- randomforest_res(subset.data, meta_table$Groups)
  df_accuracy <- rf_res$importance

  #organise data for output
  df <- NULL
  for(i in df_accuracy$Sample){
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," padj = ", sprintf("%.5g",res_tax_sig[i,"padj"]), sep=""), dim(data)[1]))
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    df <- na.omit(df)
    }

  data_to_write=NULL
  if(!is.null(filename)){
    data_to_write<-res_tax_sig[,c("baseMean","log2FoldChange","pvalue","padj")]
    data_to_write$Upregulated<-levels(meta_table[,grouping_column])[as.numeric(data_to_write$log2FoldChange>0)+1]
    write.csv(data_to_write,paste(filename,paste(levels(meta_table$Groups),collapse="_vs_"),".csv",sep=""))
  }

  out_put <- list("SignFeaturesTable"=res_tax,"plotdata"=df,"importance"=df_accuracy)

  return(out_put)

}


