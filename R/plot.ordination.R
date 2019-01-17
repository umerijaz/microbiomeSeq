#'Ordination plots
#'
#'This function produces visualisation of ordination and beta dispersion results.
#'
#' @param physeq(Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param ordination.res A solution of ordination results returned from \link[microbiomeSeq]{ordination}
#' @param  meta_table: A \code{data.frame} of environmental variables bearing atleast the grouping variable.
#' @param method: A character string specifying ordination method used to generate the ordination to plot.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#' @param adn_res: \code{adonis} results to be annotated on the plot.
#' @param betadisper_res: \code{data.frame} of betadispersion results to be annotated on the ordination plot.
#' @param pvalue.cutoff: pvalue threshold for significant dispersion results.
#' @param show.pvalues: A logical value specifying whether to annotate p-values or only significant labels for dispersion results.
#' @param N: Integer specifying number of group pairs to show on plot.
#' @param extra_marginspace: units to add extra space in margins to accomodate off plot area annotations.
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#' @examples
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references  \url{http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplot}
#'
#' @author  Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export plot.ordination
#' @export plot_ordisurf

plot.ordination <- function(ordination.res, method, pvalue.cutoff=0.05, show.pvalues=T, N=5,
                            extra_marginspace=0.35){
  sol <- ordination.res$solution
  adn_res <- ordination.res$adonis_res
  betadisper_res <- ordination.res$betadispersion
  groups <- ordination.res$groups
  #This function is applied to each level of NMDS (group) and it uses also function
  #cov.wt to calculate covariance matrix.
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  ord_res<-data.frame(x=sol$points[,1],y=sol$points[,2],Groups=groups)
  ord_res.mean=aggregate(ord_res[,1:2],list(group=ord_res$Groups),mean)

  plot.new()
  ord<-ordiellipse(sol, groups,display = "sites", kind ="se", conf = 0.95, label = T)
  dev.off()
  #Generate ellipse points
  df_ell <- data.frame()
  #include a condition to check positive definitenessas required by
  #cholskey decomposition to avoid errors in vegancovellipse
  for(g in levels(ord_res$Groups)){
    if(g!="" && (g %in% names(ord)) && all(eigen(ord[[g]]$cov)$values>0)){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(ord_res[ord_res$Groups==g,],
                                                       veganCovEllipse(ord[[g]]$cov ,ord[[g]]$center,ord[[g]]$scale))),Groups=g))
    }
  }
  colnames(df_ell)<-c("x","y","Groups")

  #coloring function
  gg_color_hue<-function(n){
    hues=seq(15,375,length=n+1)
    hcl(h=hues,l=65,c=100)[1:n]
  }

  cols=gg_color_hue(length(unique(ord_res$Groups)))

  p<-ggplot2::ggplot(data=ord_res,aes(x,y,colour=Groups))
  p<-p + ggplot2::geom_point(alpha=0.5,size = 2)
  p<-p+ggplot2::theme_bw()
  p<-p+ ggplot2::annotate("text",x=ord_res.mean$x,y=ord_res.mean$y,label=ord_res.mean$group,size=6,colour=cols,family="Courier",fontface="bold",alpha=0.8,vjust=0.3)
  p<-p+ ggplot2::geom_path(data=df_ell, aes(x=x, y=y), size=1, linetype=1,alpha=0.3)

  #axis labels
  if(method=="NMDS"){
    #annotate plot with stress value from NMDS results
    stress.value <- sol$stress
    stress.label <- paste("STRESS=",round(stress.value,4))
    p <- p + ggplot2::annotation_custom(grob = textGrob(label = stress.label, hjust = 0, gp = gpar(cex = 1.5,fontsize=8)),
                               ymin = max(ord_res$y), ymax = max(ord_res$y),
                               xmin = extra_marginspace+max(ord_res$x),xmax = extra_marginspace+max(ord_res$x))
    p<-p+xlab("NMDS1")+ylab("NMDS2")
  }
  else if(method=="PCoA"){
    p<-p+xlab(paste("Dim1 (",sprintf("%.4g",sol$eig[1]),"%)",sep=""))+ylab(paste("Dim2 (",sprintf("%.4g",sol$eig[2]),"%)",sep=""))
  }

  #add the adonis results on the plot using custom annoatation
  #this only happens if adonis results turn out significant
  gt<-NULL #ggtable table to accomodate ordination plot and adonis results
  if(!is.null(adn_res)){
    adn_pvalue<-adn_res[[1]][["Pr(>F)"]][1]
    adn_rsquared<-round(adn_res[[1]][["R2"]][1],3)
    #use the bquote function to format adonis results to be annotated on the ordination plot.
    signi_label <- paste(cut(adn_pvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ".")))
    adn_res_format <- bquote(atop(atop("PERMANOVA",R^2==~.(adn_rsquared)), atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))
    if(adn_pvalue<=pvalue.cutoff){
      gt <- plot_adonis_res(p,ord_res, adn_res_format,extra_marginspace)
    }
  }

  # add a table of beta dispersion results
  anova_table<-NULL
  if(!is.null(betadisper_res)){
    anova_table <- plot_betadisper(betadisper_res, show.pvalues,pvalue.cutoff, N)
    th <- sum(anova_table$heights)
  }

  out<-p
  if(!is.null(gt) && !is.null(anova_table)){
    out <-gridExtra::grid.arrange(gt,anova_table,heights = unit.c(unit(1, "null"), th))
  }
  else if(is.null(gt) && !is.null(anova_table)){
    out <- gridExtra::grid.arrange(p,anova_table,heights = unit.c(unit(1, "null"), th))
  }
  return(out)
}

# ==========function to annotate adonis (PERMANOVA) results onto an ordination plot
plot_adonis_res <- function(p, ord_res, adn_res,extra_marginspace){
  p<-p+theme(legend.position = "none") #get rid of the legend
  p<-p+theme(plot.margin = unit(c(1,8,1,1), "lines")) #allow extra space in margins
  #annotate plot with adonis results
  p <- p + annotation_custom(grob = textGrob(label = adn_res, hjust = 0, gp = gpar(cex = 1.5,fontsize=12)),
                             ymin = median(ord_res$y), ymax = median(ord_res$y), xmin =extra_marginspace+ max(ord_res$x), xmax = extra_marginspace+max(ord_res$x))
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  return(gt)
}

#=====================A function to plot betadispersion results onto an ordination plot
plot_betadisper <- function(betadisper_res, show.pvalues=T, pvalue.cutoff, N){
  #order the significantly dispersed groups in incresing size of p-value
  if(dim(betadisper_res)[1]>=2){
    betadisper_res<-betadisper_res[order(betadisper_res$p.value),]
  }

  colnames(betadisper_res)<-c("groups","p_value","label")
  #display significantly dispersed groups,select the number of groups to display
  #and whether or not to show p-values select groups to display on plot
  anova_table_display <- subset(betadisper_res, p_value<=pvalue.cutoff)
  print(anova_table_display)
  if(dim(anova_table_display)[1]<N){
    N <- dim(anova_table_display)[1]
  }
  if(show.pvalues){
    anova_table_display$p_label <- paste(rep("p-value="),sprintf("%.e",anova_table_display$p_value)," ",anova_table_display$label, sep="")
    anova_table_display<-anova_table_display[1:N, c("groups","p_label")]
  }
  else{
    anova_table_display<-anova_table_display[1:N,c("groups","label")]
  }
  #theme_minimal in tableGrob provides a white background
  anova_table <- gridExtra::tableGrob(anova_table_display,rows = NULL,cols = NULL,theme =ttheme_minimal()) #set rows/cols to null to get rid of the row numbering
  title <- grid::textGrob("BETA-DISPERSION",gp=gpar(fontsize=12))
  anova_table <- gtable::gtable_add_rows( anova_table, heights = grobHeight(title) + unit(5,"mm"), pos = 0)
  anova_table <- gtable::gtable_add_grob( anova_table,  title, 1, 1, 1, ncol(anova_table))
  ##calculate height of table and pass it to glob to avoid white space between plot and table
  return(anova_table)
}

plot_ordisurf <- function(sol, meta_table, env.variable, grouping_column){

  groups <- meta_table[,grouping_column] #get grouping information from meta data

  df=data.frame(x=sol$point[,1],y=sol$point[,2],Groups=groups)
  #Add a dummy variable corrresponding to the selected variable
  meta_table$var <- meta_table[,env.variable]

  #fit a surface for a selected variable onto ordination stats
  ordi<- vegan::ordisurf(sol,meta_table$var ,plot = FALSE, bs="ds")
  ordi.grid <- ordi$grid #extracts the ordisurf object
  #str(ordi.grid) #it's a list though - cannot be plotted as is
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas

  #make the plot
  p<-ggplot2::ggplot()+stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = ..level..),positon="identity") #can change the binwidth depending on how many contours you want
  p<-p+ ggplot2::geom_point(data=df,aes(x,y,fill=Groups),pch=21,size=3)
  p<-p+ ggplot2::scale_colour_continuous(high = "darkgreen", low = "darkolivegreen1") #here we set the high and low of the colour scale.  Can delete to go back to the standard blue, or specify others
  p<-p+ ggplot2::labs(colour = paste(env.variable)) #another way to set the labels, in this case, for the colour legend
  p<-p+ ggplot2::theme_bw()
  return(p)
}
