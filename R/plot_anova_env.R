#'ANOVA of environmental variables
#'
#'This functions applies analysis of variance on the measured environmental variables
#'in specified groups and constructs a plot showing how the environmental variables vary in the groups.
#'The plot is annotated with significance level as obtained from \code{anova}.
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Name of a categorical variable that is preffered for grouping the.
#'        information, this should be one of the components of grouping vector.
#' @param pValueCutoff. p-value threshold for significance of anova results, default set to 0.05.
#' @param select.variables. A vector of character strings(s) specifying environmental variable(s) to be analysed. If
#' not supplied, all numeric variables are analysed.
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#' @examples
#' data(pitlatrine)
#' physeq<-pitlatrine
#' p1<-plot_anova_env(physeq,grouping_column =  "Country",select.variables=c("Temp","pH"))
#' print(p1)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export plot_anova_env
#'
plot_anova_env <- function(physeq, grouping_column, pValueCutoff=0.05,select.variables=NULL){

  #get meta data from phyloseq object
  meta_table <- as.data.frame(sample_data(physeq))
  #pick numerical variables of environmental data
  env_table <- meta_table[,sapply(meta_table,is.numeric)]
  df<- reshape2::melt(as.matrix(env_table))
  names(df)<-c("sample","measure","value")
  #Incorporate categorical data in df
  df<-data.frame(df,(meta_table[, grouping_column])[as.character(df$sample),])

  #do anova of environmental variables between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  df <- anova_res$df #get updated environmental measure information

  #pick selected variables
  if(!is.null(select.variables)){
    df <- df[which(df$measure%in%select.variables),]
    df_pw<-df_pw[which(df_pw$measure%in%select.variables),]
  }
  #Draw the boxplots
  p<-ggplot2::ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+ggplot2::geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+ggplot2::theme_bw()+geom_point(size=1,alpha=0.2)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+ggplot2::facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Groups")
  p<-p+ggplot2::theme(strip.background = element_rect(fill = "white"))
  #This loop will generate the lines and signficances
  for(i in 1:dim(df_pw)[1]){
    p<-p+ ggplot2::geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
    p<-p+ ggplot2::geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
  }
  return(p)
}
