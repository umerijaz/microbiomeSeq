#'Relationship between sub-communities and environment
#'
#' This function calculates correlation between modules (sub-communities) and environmental
#' variables. The most connected node(taxon) in a given module is considered a good representation
#' of the module. This is implemented in accordance to (Deng et al. 2012) where
#' correlations between module-based eigengenes and environmental factors are
#' used to detect the modules' response to environmental change.
#'
#' @param  co_occur_res  an object returned by `co_occurence_network` function.
#' @param  select.variables  environmental variables of interest.
#'
#' @examples
#' mod.env.cor <- module_env_correlation(co_occur)
#' p <- plot_taxa_env(mod.env.cor)
#' print(p)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references Deng, Ye, Yi-Huei Jiang, Yunfeng Yang, Zhili He, Feng Luo, and Jizhong Zhou. 2012.
#' “Molecular Ecological Network Analyses.” BMC Bioinformatics 13 (1). BioMed Central: 113.
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export module_env_correlation
#'

module_env_correlation <- function(co_occur_res, select.variables=NULL, ...){

  comdat <- co_occur_res$comm_data

  comm_cor <- co_occur_res$net$comm_results

  comm_graph <- co_occur_res$net$graph

  meta_table <- co_occur_res$meta_table[sapply(co_occur_res$meta_table, is.numeric)]

  net.comm <- data.frame(comdat[, colnames(comdat)%in%comm_cor$taxon])

  net.comm <- net.comm[rownames(net.comm)%in%rownames(meta_table),]

  if(!is.null(select.variables)){ meta_table<- meta_table[, select.variables]}

  comm.taxa.cor <- tables.correlate(net.comm, meta_table, method= "pearson")

  colnames(comm.taxa.cor) <- c("Taxa", "Env", "Correlation", "Pvalue")

  comm.taxa.cor$AdjPvalue<-rep(0,dim(comm.taxa.cor)[1])

  # correct pvalues for multiple testing
  comm.taxa.cor$Pvalue <- as.numeric(as.character(comm.taxa.cor$Pvalue))

  comm.taxa.cor <- p.adjust.cor(comm.taxa.cor)

  comm_cor$Taxa <- comm_cor$taxon

  comm.cor.merge <- merge.data.frame(comm_cor, comm.taxa.cor, by="Taxa")

  graph.btn <- betweenness(comm_graph)

  graph.btn <- data.frame(row.names = names(graph.btn), Taxa = names(graph.btn) ,betweenness=graph.btn)

  graph.btn <- graph.btn[rownames(graph.btn)%in%colnames(net.comm),]

  comm.cor.merge  <- merge.data.frame(comm.cor.merge, graph.btn, by="Taxa")

  df <- NULL

  for(i in levels(comm.cor.merge$module)){

    modi <- subset(comm.cor.merge, module==i)

    if(dim(modi)[1]>2){

      modi <- modi[which(modi$betweenness == max(modi$betweenness)),]
      tax <- modi$taxon
      pvalue <- as.numeric(as.character(modi$AdjPvalue))
      corr <- as.numeric(as.character(modi$Correlation))
      env <- modi$Env
      trt <- modi$trt

      tmp <- data.frame("Taxa"=paste("mod",i,"-", tax) , "Env"=env, "Type"=trt, "tax"=tax, "AdjPvalue"=pvalue, "Correlation"=corr)

      if(is.null(df)){df<-tmp} else {df <- rbind(df, tmp)}

    }

  }

  df$Significance<-cut(df$AdjPvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

  df <- na.omit(df)

  return(df)

}
