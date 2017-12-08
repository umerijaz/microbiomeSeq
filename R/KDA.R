#'kernel-based differential analysis
#'
#'This function performs differential analysis using a distance-based kernel score test.
#'It obtains set(s) of differentially abundant features between a specified pair of conditions/groups.
#'The sets are obtained by grouping based on correlation of abundance among features. The strength of the relationship that must exist for
#'features/variables to belong to the same set can be specified. A similar approach applies to environmental
#'variables which can also be grouped into sets depending on a desired level of correlation amongst member variables.
#'
#' @param physeq (Required). A \link[phyloseq]{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the.
#'        information, this should be one of the components of grouping vector.
#' @param analyse (optional). A character string specifying whether to analyse taxa abundance ("abundance") or sample data ("meta"). Default is set to analyse taxa abundance.
#' @param method (optional). A character string specifying the kernel based method to be used for differential analysis,
#'        two options are available for this; "dscore" and "sscore". See \link[KMDA]{dscore} and  \link[KMDA]{sscore} for details.
#' @param  adjusted.p.value.threshold Cut off p-value for significance of differentially abundant/correlated taxa/environmental variables, default is 0.05.
#' @param p.adjust.method (optional). method to adjust raw pvalues obtained by the distance based score tests.
#' @param  corr (optional). Threshold correlation on which grouping of features/ variables is based.
#'
#' @return Returns a list of two items: \itemize{
#'         \item  kscore_table: A \code{data.frame} of taxa/ environmental variable with kernel based score stats, raw and adjusted pvalues.
#'         \item plotdata: A \code{data.frame} similar to kcore_table but organised in a way compatible to ggplot.
#'        }
#'
#' @examples
#' data(pitlatrine)
#' physeq<-data(pitlatrine)
#' physeq<- taxa_level(physeq, "Phylum")
#' kda_sig  <- KDA(physeq,grouping_column = "Country")
#' plot_kda(kda_sig$plotdata)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export KDA
#'
KDA <- function(physeq, grouping_column, analyse="abundance", method="dscore",
                p.adjust.method="BH",corr=0.99,adjusted.p.value.threshold= 0.05, select.variables=NULL, ...){

  meta_table <- data.frame(sample_data(physeq))

  # choose whether to analyse taxa differential abundance or differential variation
  # among measured environmental variables
  analyse <- match.arg(analyse,c("abundance","meta"))

  if(analyse=="abundance"){
    #enforce orientation of Otu_table , this is beacause dscore/sscore requiresamples to be rows.
    if(!taxa_are_rows(physeq)){

      abund_table <- t(abund_table)
    }
    analyse.data <- data.frame(abund_table)
  }
  else if(analyse=="meta"){

    analyse.data <- (get.num.variables(physeq))$num.variables

    if(!is.null(select.variables)){
      analyse.data <- analyse.data[,colnames(analyse.data)%in%select.variables]
    }

    analyse.data <- t(analyse.data) #transpose such that samples are rows
  }
  # pick out grouping variable
  meta_table$Groups <- meta_table[,grouping_column]
  # get levels of grouping variable
  group_levels<-levels(meta_table$Groups)
  # change the variables to 1 or 0 as required by the kmda procedure
  meta_table$Groups <- sapply(meta_table$Groups, function(x) ifelse(x == group_levels[1], 1, 0))
  # group taxa/environmental variables into sets using pearson correlation coefficients
  Sets<- KMDA::pearson.group(as.matrix(analyse.data),corr)
  Sets_mapping<-data.frame(row.names=rownames(analyse.data),Set=Sets)
  #obtain distance based kernel scores of each otu set
  #specify the lower and upper bounds of kernel parameter and number of gridpoints
  #selected in the interval
  KMtable<-NULL
  for(i in unique(Sets)){
    x=analyse.data[Sets==i,,drop=F]
    if(method=="dscore"){
      ks=KMDA:: dscore(x,meta_table$Groups,lower=1,upper=10,m=3)
    }
    if(method=="sscore"){
      ks=KMDA::sscore(x,meta_table$Groups,lower=10^-3,upper=10^3,m=10)
    }
    tmp<-data.frame(row.names=i,kscore=ks)
    ifelse(is.null(KMtable),KMtable<-tmp,KMtable<-rbind(KMtable,tmp))
  }
  #Adjust p-values
  if(is.null(p.adjust.method)){
    p.adjust.method <- "BH"
  }
  KMtable$kscore.adjusted<-p.adjust(KMtable$kscore,method=p.adjust.method,n=dim(KMtable)[1])

  #Select the sets of otus/env_variables that are differentially expressed
  kscore_selected<-rownames(KMtable)[KMtable$kscore.adjusted<=adjusted.p.value.threshold]
  #organise the data in a format accepted by ggplot2
  df<-NULL
  for(i in kscore_selected){
    tmp<- reshape2::melt(as.matrix(analyse.data[Sets==i,,drop=F]))
    colnames(tmp)<-c("Set_variable","Sample","Value")
    tmp$Groups<-meta_table[as.character(tmp$Sample),"Groups"]
    tmp$Set<-Sets_mapping[as.character(tmp$Set_variable),]
    tmp$Set<-paste("Set",tmp$Set,sprintf("Adj.p = %0.5g",KMtable[as.character(tmp$Set),"kscore.adjusted"]),
                   cut(KMtable[as.character(tmp$Set),"kscore.adjusted"],breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))
    ifelse(is.null(df),df<-tmp,df<-rbind(df,tmp))
  }
  #change the grouping values back to original representation
  df$Groups<-sapply(df$Groups, function(x) ifelse(x == 1, group_levels[1], group_levels[2]))

  out <- list("plotdata"=df, "kscore_table"=KMtable)
  return(out)
}

plot_kda <- function(df){
  p<-ggplot2::ggplot(aes(x=Set_variable,y=Value,color=Groups),data=df)+theme_bw()+geom_boxplot(outlier.size = NA)+facet_grid(. ~ Set, drop=TRUE,scales="free",space="free_x")
  p<-p+ ggplot2::theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 12, colour = "black", angle = 90,vjust = 0))
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))+ylab("")+xlab("")
  return(p)
}
