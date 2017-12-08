#'Taxa level collating
#'
#'This function takes a physeq object and returns a physeq object with taxa at a specified taxonomic level.
#'It extends the \link[phyloseq]{tax_glom}  to include names of corresponding taxa
#'at a specified taxonomy level and  updating tree tip labels accordingly.
#'
#' @param physeq (Required). A \link[phyloseq]{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param which_level (Required). Character string specifying taxonomic level.
#' @return physeq object at specified taxonomic level.
#'
#' @examples
#'
#' data(pitlatrine)
#' physeq<-data(physeq)
#' physeq <- taxa_level(physeq = physeq,which_level = "Family")
#'
#'@export taxa_level
#'
taxa_level <- function(physeq,which_level){
    OTU <- otu_table(physeq)
    SAM <- sample_data(physeq)
    OTU_taxonomy <- tax_table(physeq)
    new_abund_table<-NULL
    if(which_level=="Otus"){
      OTU_tree <- phy_tree(physeq)
      new_abund_table<-OTU
    } else {
      list<-unique(OTU_taxonomy[,which_level])
      new_abund_table<-NULL
      for(i in list){
        tmp<-data.frame(rowSums(OTU[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i]]))
        if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
        if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
      }
    }
    OTU<-as.data.frame(as(new_abund_table,"matrix"))
    #Convert the data to phyloseq format
    OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
    TAX = tax_table(as.matrix(OTU_taxonomy))
    SAM = sample_data(SAM)

    physeq<-NULL
    if(which_level=="Otus"){
      physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
    } else {
      physeq<-merge_phyloseq(phyloseq(OTU),SAM)
    }
  return(physeq)
}
