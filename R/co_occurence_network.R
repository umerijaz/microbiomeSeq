#' Co-occurence network
#'
#'This function uses co-occurence pattern analysis to identify co-occuring features/taxa in community
#'data under specified environmental conditions. Co-occurence is measured as positive correlation whose threshold(s)
#'can be specified as indicated in arguments section. Amongst these features, pairwise co-occurences which are outstanding within sub communities
#'are detected. p-values generated during pairwise correlations are adjusted for multiple comparisons
#'by false discovery rate. The network statistics used to assign importance of taxa/features include betweenness, closeness,
#'and eigenvector centrality.
#'
#' @param physeq(Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method. A character string for the correlation method to be used; options include: "cor" and "bicor"
#' @param rhos (required). A vector specifying thresholds for correlation between co-occuring pairs of features.
#' @param select.condition (optional). A character string speifying name of a desired condition/group. If not supplied, co-occuence analysis
#'       is performed amongst all conditions present in the grouping variable.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information.
#' @param qval_threshold Cut off for "fdr" adjusted p-values.
#'
#' @return Files of visual representation of the network showing subcommunities (identified by colors),
#' network statistics and file containing pairwise corrrelations of taxa in under different conditions.
#' @examples
#' data(pitlatrine)
#' physeq <- taxa_level(pitlatrine, "Phylum")
#' co_occr <- co_occurence_network(physeq, grouping_column = "Country", rhos = 0.35, select.condition = "V", scale.vertex.size=3, scale.edge.width=15)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references [Ryan J. Williams, Adina Howe and Kirsten S. Hofmockel.
#'            Demonstrating microbial co-occurrence pattern analyses within and between ecosystems,
#'            Frontiers in Microbial Ecology, 5:358, 2014].
#'
#' @author  Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export co_occurence_network
#'

co_occurence_network <- function(physeq ,qval_threshold=0.05, grouping_column,rhos=c(-0.75,-0.5,0.5,0.75),select.condition=NULL,method="cor", filename=NULL, ...){

  corr_results <- network_correlation(physeq, grouping_column, select.condition, method, filename)

  comm_results<-comm_stats(corr_results$results, rhos)

  edge_lists<-co_occur_pairs(corr_results$results, rhos)

  edge_lists<-edge_lists[edge_lists$qval<=qval_threshold,]

  net <- construct_network(edge_lists, comm_results, ...)

  out <- list(net=net, comm_data=corr_results$comm.data, meta_table=corr_results$meta_table)

  return(out)
}

#computing correlation of taxa/features
network_correlation <- function(physeq, grouping_column, select.condition, method, filename){
  abund_table <- as.data.frame(otu_table(physeq))
  meta_table <- data.frame(sample_data(physeq))
  meta_table$Groups <- meta_table[,grouping_column]

  if(!is.null(select.condition)){meta_table <- subset(meta_table, Groups==select.condition)}

  # === rarefy
  comm.data<-vegan::rrarefy(abund_table,min(rowSums(abund_table)))
  trts<-unique(meta_table$Groups)

  results<-matrix(nrow=0,ncol=7)
  for(a in 1:length(trts)){
    trt.temp<-trts[a]
    #subset the dataset for those treatments
    temp<-subset(comm.data, meta_table$Groups==trt.temp)
    #remove empty species
    temp<-temp[,colSums(temp)>0]

    new.row<-NULL
    if(method=="bicor"){
      #Biweight midcorrelation
      res<- WGCNA::bicorAndPvalue(temp,use="p")

      res$bicor[upper.tri(res$bicor)] <- NA
      diag(res$bicor)<-NA
      res$p[upper.tri(res$p)]<-NA
      diag(res$p)<-NA

      new.row<-data.frame(rep(trts[a],dim(temp)[2]),reshape2::melt(as.matrix(res$bicor)),reshape2::melt(as.matrix(res$p))[,3])

    } else if (method=="cor"){
      #Pearson correlation
      res<- WGCNA::corAndPvalue(temp,use="p")

      res$bicor[upper.tri(res$cor)] <- NA
      diag(res$cor)<-NA
      res$p[upper.tri(res$p)]<-NA
      diag(res$p)<-NA

      new.row<-data.frame(rep(trts[a],dim(temp)[2]), reshape2::melt(as.matrix(res$cor)),reshape2::melt(as.matrix(res$p))[,3])
    }
    colnames(new.row)<-c("trt","taxa1","taxa2","rho","p.value")
    new.row<-data.frame(new.row,ab1=as.vector(sapply(as.character(new.row$taxa1),function(x){colSums(temp)[x]})),ab2=as.vector(sapply(as.character(new.row$taxa2),function(x){colSums(temp)[x]})))
    results<-rbind(results,new.row)
  }

  results$taxa1<-as.character(results$taxa1)
  results$taxa2<-as.character(results$taxa2)
  #We remove the upper triangle and diagonal from the correlation just calculated using NA assignments done before
  results<-results[complete.cases(results),]
  results[(results$ab1 <= 1) | (results$ab2 <= 1),"rho"]<-0
  results[(results$ab1 <= 1) | (results$ab2 <= 1),"p.value"]<-1

  if(!is.null(filename)){write.csv(results,paste(filename,"_co-occurence","_",paste(unique(meta_table$Groups),collapse="_vs_"),"_correlations.csv",sep=""))}

  out <- list(comm.data=comm.data, results=results, meta_table=meta_table)
  return(out)
}
