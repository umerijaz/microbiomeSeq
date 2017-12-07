#' Canonical Correspondence Analysis plot
#'
#' This function finds a set of best environmental variables that describe community structure
#' superimposes them on the ordination plot.
#'
#' @param `physeq` is a required phyloseq object containing taxa abundance and meta data.
#' @param `grouping_column` is the variable in the meta data with respect to which the data should be grouped,
#' @param `pvalueCutoff` the threshold p-value in `anova` of distance matrices, default set to `0.05`.
#' @param `env.variables` is a list of variables prefered to be on the cca plot.
#' @param `exclude.variables` a list of variables to be excluded from the cca plot.
#' @param `num.env.variables` is an integer specifying the number of variables to show on the cca plot.
#'
#' @import vegan
#' @import gridBase
#' @import ggplot2
#'
#' @return A ggplot object
#'
#' @examples
#' data(pitlatrine)
#' physeq<-data(physeq)
#' ccaplot  <- plot_cca(physeq, grouping_column = "Country")
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}

#'
#' @export plot_cca

plot_cca <- function(physeq, grouping_column, pvalueCutoff=0.01,norm_method=NULL, env.variables=NULL,
                     num.env.variables=NULL, exclude.variables=NULL, draw_species=F){
  #extract abundance and meta data from supplied phyloseq object
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))

  #Use adonis to find significant environmental variables
  abund_table.adonis <- vegan::adonis(abund_table ~ ., data=meta_table)
  #pick significant features
  bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=pvalueCutoff]
  #throw out na in case any exist
  bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
  #pick out the selected variables if at all they are part of the significant ones
  if(!is.null(env.variables)&&(env.variables%in%bestEnvVariables)){
    bestEnvVariables <- env.variables
  }
  #provide number of variables to display
  if(!is.null(num.env.variables)){
    if(num.env.variables>length(bestEnvVariables)){
      stop(cat(paste("Choose a number less than",length(bestEnvVariables))))
    }else{
      bestEnvVariables <- bestEnvVariables[1:num.env.variables]
    }
  }
  #exclude selected variables from appearing on the plot
  if(!is.null(exclude.variables)&&(exclude.variables%in%bestEnvVariables)){
    bestEnvVariables <- bestEnvVariables[!(bestEnvVariables%in%exclude.variables)]
  }
  #We are now going to use only those environmental variables in cca that were found significant
  eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=meta_table)",sep="")))

  scrs<- vegan::scores(sol,display=c("sp","wa","lc","bp","cn"))
  #Extract site data first
  df_sites<-data.frame(scrs$sites,meta_table[,grouping_column])
  colnames(df_sites)<-c("x","y","Groups")

  #Draw sites
  p<-ggplot2::ggplot()
  p<-p+ggplot2::geom_point(data=df_sites,aes(x,y,colour=Groups))
  #Draw biplots
  multiplier <- vegan:::ordiArrowMul(scrs$biplot)
  df_arrows<- scrs$biplot*multiplier
  colnames(df_arrows)<-c("x","y")
  df_arrows=as.data.frame(df_arrows)
  p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)
  p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)
  # Draw species
  df_species<- as.data.frame(scrs$species)
  colnames(df_species)<-c("x","y")
  # Either choose text or points
  #p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species)))
  if(draw_species){
    p<-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2)
  }
  p<-p+theme_bw()+xlab("CCA1")+ylab("CCA2")
  return(p)
}
