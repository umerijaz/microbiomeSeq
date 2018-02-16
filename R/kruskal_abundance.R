#'Kruskal-Wallis differential abundance  analysis
#'
#'This function finds the features that are significantly
#'differentially abundant  in the provided taxa abundance data under different conditions
#'using Kruskal-Wallis test. The p-values values generated are corrected for multiple testing
#'using  family wise error rate. Significance is based on the corrected pvalue threshold. Significant features
#'are assigned importance using random forest classifier. The measure of importance used in this case is mean decrese accuracy.
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information, this should be one of the components of grouping vector.
#' @param  pvalue.threshold. Cut off p-value for significance of differentially abundant taxa, default is 0.05.
#' @return Returns a list of three items: \itemize{
#'         \item  SignfeaturesTable: A \code{data.frame} of taxa with coresponding raw p-values, corrected p-values, family wise error rate and expected abundance
#'         computed by using raw p-values.
#'         \item  importance: A \code{data.frame} of taxa mean decrease in accuracy as obtained by random forest classifier.
#'         \item plotdata: A \code{data.frame} of taxa and corresponding corrected p-values, importance rank organised
#'         in form accepted by ggplot.
#'        }
#' @examples
#' data(pitlatrine)
#' physeq <- pitlatrine
#' physeq <- taxa_level(physeq,"Phylum")
#' kw_sig <-  kruskal_abundance(physeq, "Country")
#' #plot the significant features
#' plot_signif(kw_sig$plotdata) #see function \link[microbiomeSeq]{plot_signif}
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references \url{http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics/practicals/microarrays_berry_2010/berry_feature_selection.html}
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export kruskal_abundance
#'
kruskal_abundance <- function(physeq, grouping_column,pvalue.threshold=0.05)
{
  abund_table <- otu_table(physeq)
  meta_table <-data.frame(sample_data(physeq))
  meta_table$Groups <- meta_table[,grouping_column]

  kruskal.wallis.table <- data.frame()
  data <- as.data.frame(abund_table)
  for (i in 1:dim(data)[2]){
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
  }
  kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
  kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
  kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
  rownames(kruskal.wallis.table) <- kruskal.wallis.table$id

  #==================significant feature selection =================================#
  last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
  selected <- 1:last.significant.element
  sig_res <-kruskal.wallis.table$id[selected]

  #==random forest classifier ==#
  subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
  kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
  subset.data <- subset.data[,sig_res]
  rf_res <- randomforest_res(subset.data, meta_table$Groups)
  df_accuracy <- rf_res$importance

  df <- NULL
  for(i in df_accuracy$Sample){
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," p.adj = ",sprintf("%.10g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    df <- na.omit(df)
    }

  out <- list("SignfeaturesTable"=kruskal.wallis.table, "plotdata"=df, "importance"=df_accuracy)

  return(out)

}
