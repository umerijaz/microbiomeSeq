#'Taxa-Environment correlation
#'
#'This function finds the relationship between most abundant taxa
#'and numerical environmental variables based on correlation. The abundance of each feature/taxa is
#'correlated with each of the environmental variables. A correlation test is performed and associated
#'p-values are adjusted for multiple testing. The scheme of adjustment is elaborated in the arguments section. The function returns a data.frame with raw p-values, corrected p-values,
#'and correlation results.
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information, this should be one of the components of grouping vector.
#' @param method A character string indicating which correlation coefficient is to be computed, available options are "pearson" which is also the default, "kendall" and "spearman".
#' @param pvalue.threshold. Cut off p-value for significance of correlation between taxa abundance and environmental variables, default is 0.05.
#' @param padjust.method. Method for adjusting p-values.
#' @param adjustment (optional). An integer to specify how the p-values should be adjusted;
#'        \itemize{
#'        \item 1 - donot adjust
#'        \item 2 - adjust environmental variables + Groups (column on the correlation plot)
#'        \item 3 - adjust Taxa + Groups (row on the correlation plot for each Groups)
#'        \item 4 - adjust Taxa (row on the correlation plot)
#'        \item 5 - adjust environmental variables (panel on the correlation plot)
#'        }
#' @param select.variables (optional). Character string for the names(s) of environmental variables to be considered in correlation compution.
#'        if not given, all numerical variables are considered.
#' @return A \code{data.frame} of correlation results.
#'
#' @examples
#' data(pitlatrine)
#' physeq <- pitlatrine
#' physeq <- taxa_level(physeq, "Phylum")
#' tax_env_cor <- taxa.env.correlation(physeq,"Country")
#' plot_taxa_env(tax_env_cor)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export taxa.env.correlation
#' @export plot_taxa_env

taxa.env.correlation <- function(physeq, grouping_column, method="pearson", pvalue.threshold=0.05,
                                 padjust.method="BH", adjustment=1, num.taxa=50, select.variables=NULL){

  method<- match.arg(method,c("pearson", "kendall", "spearman"),several.ok = F)
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- as.data.frame(t(otu_table(physeq)))
  #abund_table <- t(otu_table(physeq))
  meta_table <- data.frame(sample_data(physeq))
  #get grouping information
  groups<-meta_table[,grouping_column]
  #select variables to show (exclude) in (from) correlation plot.
  if(!is.null(select.variables)){
    meta_table <- subset(meta_table,select=select.variables)
  }
  #pick the numerical environmental variables since correlation function only accepts numerical variable
  mt_env <- meta_table[,sapply(meta_table,is.numeric)]

  #Now get a filtered abundance table based on selected variables
  abund_table_filt<-abund_table[rownames(mt_env),]

  # #pick top most  num.taxa taxa
  # select.top.taxa <- top.taxa(abund_table, num.taxa)
  # abund_table_filt <- select.top.taxa$abund_table
  abund_table_filt<-abund_table_filt[,order(colSums(abund_table_filt),decreasing=TRUE)]
  #Extract list of top N Taxa
  taxa_list<-colnames(abund_table_filt)[1:num.taxa]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
  abund_table_filt<-data.frame(abund_table_filt[,colnames(abund_table_filt) %in% taxa_list])

  #Now calculate the correlation between individual Taxa and the environmental data
  df <- tables.correlate(abund_table_filt, mt_env, groups, method)
  colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$Correlation<-as.numeric(as.character(df$Correlation))
  # add column for adjusted p-values
  df$AdjPvalue<-rep(0,dim(df)[1])

  # correct pvalues for multiple testing
  df <- p.adjust.cor(df, adjustment, padjust.method)
  #Now we generate the labels for signifant values
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  #We ignore NAs
  df<-df[complete.cases(df),]
  return(df)
}

tables.correlate<-function(table1, table2, groups=NULL, method){
  df<-NULL
  for(i in colnames(table1)){
    for(j in colnames(table2)){

      if(!is.null(groups)){
        for(k in unique(groups)){
          a<-table1[groups==k,i,drop=F]
          b<-table2[groups==k,j,drop=F]
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)

          if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}
         }
      }
      else{

        a<-table1[,i,drop=F]
        b<-table2[,j,drop=F]
        tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value)

        if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}

        }

      }
    }

  df<-data.frame(row.names=NULL,df)
  return(df)
}


plot_taxa_env <- function(df){
  p <-ggplot2::ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
  p <- p + ggplot2::geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")
  p<-p+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  p<-p+ ggplot2::geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)
  p<-p+ ggplot2::facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")
  p<-p+ ggplot2::xlab("Groups")
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))
  return(p)
}

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
#adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")

# df is a data frame
p.adjust.cor <- function(df,adjustment=1,padjust.method="BH"){
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$Env)){
      for(j in unique(df$Type)){
        sel<-df$Env==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$Taxa)){
      for(j in unique(df$Type)){
        sel<-df$Taxa==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$Taxa)){
      sel<-df$Taxa==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
    }
  } else if (adjustment==5){
    for(i in unique(df$Env)){
      sel<-df$Env==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
    }
  }
  return(df)
}
