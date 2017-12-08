#'Local Contribution to Beta Diversity (LCBD)
#'
#'This function finds the most abundant taxa  for each sample. It also calculates local
#'contribution to beta  diversity using a selected dissimilarity coefficient.
#'It returns a ggplot object which a visual representation of the most abundant
#'taxa for each of the samples. The number of top taxa can be suggested as an argument. The visual
#'representation is limited to 21 top taxa, more information can be availed by supplying a file name to which
#'details of LCBD can be written.
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method. A \code{method} to used to calculate dissimilarity coefficents.
#'        Available methods are: "hellinger", "chord", "chisquare", "profiles",
#'        "percentdiff", "ruzicka", "divergence", "canberra",
#'        "whittaker", "wishart", "kulczynski", "jaccard", "sorensen", "ochiai",
#'        "ab.jaccard", "ab.sorensen", "ab.ochiai", "ab.simpson", "euclidean". the function uses
#'        "hellinger" as a default.
#' @param grouping_column(Required). Name of a categorical variable that is preffered for grouping the.
#'        information.
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#'
#' @examples plot_taxa(physeq, "Country")
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#'
#' @export plot_taxa
#'

plot_taxa <- function(physeq,grouping_column,method="hellinger",number.taxa=21,filename=NULL){

  #==extract components of the phyloseq object
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))
  OTU_taxonomy <- data.frame(tax_table(physeq))

  #Enforce orientation of the phyloseq object
  if(taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }

  #===Calculate beta diversity and extract measure for local contribution to beta diversity
  beta_div<-beta.div(abund_table,method=method,sqrt.D=F,samp=T,nperm=999)
  df_LCBD<-data.frame(Sample=names(beta_div$LCBD),LCBD=beta_div$LCBD,p.LCBD=beta_div$p.LCBD)

  #=== add grouping information to the LCBD results
  df_LCBD<-data.frame(df_LCBD,Groups=meta_table[rownames(df_LCBD),grouping_column])

  if(!is.null(filename)){
    write.csv(df_LCBD,paste(filename,"_LCBD",".csv",sep=""))
  }

  select.top.taxa <- top.taxa(abund_table, number.taxa)
  new_x <- select.top.taxa$abund_table
  number.taxa <- select.top.taxa$number.taxa

  #arrange data for plotting in a format compatible to ggplot
  df<-NULL
  for (i in 1:dim(new_x)[2]){
    tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Groups=meta_table[,grouping_column])
    if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
  }
  df<-data.frame(df,df_LCBD[as.character(df$Sample),c("LCBD","p.LCBD","Groups")])
  #==plot the data
  colours <- microbiomeseq_cols()
  p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Groups, drop=TRUE,scale="free",space="free_x")
  p<- p+ guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=colours[1:(number.taxa+1)])+theme_bw()+xlab("Samples")
  p<-p+ scale_y_continuous(expand = c(0.02,0))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  p<-p+geom_point(aes(Sample,-0.02,size=LCBD))+theme(strip.background = element_rect(fill = "white"))
  return(p)
}

top.taxa <- function(abund_table, number.taxa){
  #==== sort the abundance table by total abundance of each taxa  in decreasing order
  abund_table<-abund_table[,order(colSums(abund_table),decreasing=TRUE)]
  #Extract list of top number.taxa Taxa
  taxa_list<-colnames(abund_table)[1:number.taxa]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
  number.taxa<-length(taxa_list)
  #Generate a new table with everything added to Others
  new_x<-data.frame(abund_table[,colnames(abund_table) %in% taxa_list],Others=rowSums(abund_table[,!colnames(abund_table) %in% taxa_list]))

  out <- list("abund_table"=new_x, "number.taxa"=number.taxa)
  return(out)
}

microbiomeseq_cols <- function(){
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00",
               "#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000",
               "#FFFF00",grey.colors(1000));
  return(colours)
}
