#'Alpha diversity with ANOVA
#'
#'This function calculates alpha diversity of provided community data using
#'selected indices/method(s). It performs pair-wise ANOVA of diversity measures between groups
#'and outputs a plot for each of the selected methods(indices) annotated with significance levels.
#' @importFrom phyloseq t
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "richness", "fisher", "simpson", "shannon" and "evenness".
#' @param grouping_column (Required). A character string specifying the name of a categorical variable containing  grouping information.
#' @return Returns a ggplot object which can further be manipulated further.
#'
#' @examples
#' data(pitlatrine)
#' physeq <- pitlatrine
#' p<-plot_anova_diversity(physeq, method = c("richness","simpson"),grouping_column =  "Country",pValueCutoff=0.05)
#' plot_anova_diversity(physeq, method = c("richness","shannon"), grouping_column = "Depth")
#' print(p)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export plot_anova_diversity
#'

plot_anova_diversity <- function(physeq, method, grouping_column,pValueCutoff=0.05, outfile="anova_diversity.csv")
{
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)

  #get diversity measure using selected methods
  div.df <- alpha_div(physeq,method)

  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  df[,grouping_column]<-as.factor(df[,grouping_column])
  #perform anova of diversity measure between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  write.csv(df_pw, file = outfile)
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Samples")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Groups")

  #This loop will generate the lines and signficances
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
    }
  }
  return(p)
}
