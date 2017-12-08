#'pairwise ANOVA
#'
#'This function performs ANOVA of a given measure in specified groups, in addition it
#'also performs pairwise ANOVA of the measure between possible pairs of levels in the grouping variable. It returns
#'p-values obtained.
#'
#' @param df (Required). A \code{data.frame} containg measures a measure to be analysed.
#' @param meta_table (Required). A data frame containing atleast one variable (grouping variable).
#' @param grouping_column (Required). A character string specifying name of the grouping variable in the supplied meta table.
#' @return It returns a list of two data frames: One being p-values obtained from the pairwise ANOVA of
#' the measure and levels of grouping variable , the other is containing updated measure.
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export perform_anova
#'

perform_anova <- function(df,meta_table,grouping_column,pValueCutoff){

  dt<-data.table::data.table(data.frame(df,.group.=meta_table[,grouping_column]))
  #specifying a p-value cutoff for the ggplot2 strips
  pval<-dt[, list(pvalue = sprintf("%.2g",
                                   tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
           by=list(measure)]
  #Filter out pvals that we are not significant
  pval<-pval[!pval$pvalue=="",]
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]

  #using sapply to generate significances for pval$pvalue using the cut function.
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, pValueCutoff, Inf),label=c("***", "**", "*", "")))})

  #Update df$measure to change the measure names if the grouping_column has more than three classes
  if(length(unique(as.character(meta_table[,grouping_column])))>2){
    df$measure<-as.character(df$measure)
    if(dim(pval)[1]>0){
      for(i in seq(1:dim(pval)[1])){
        df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
      }
    }
    df$measure<-as.factor(df$measure)
  }
  #Get all possible pairwise combination of values in the grouping_column
  s<-combn(unique(as.character(df[,grouping_column])),2)

  #df_pw will store the pair-wise p-values
  df_pw<-NULL
  for(k in unique(as.character(df$measure))){
    #We need to calculate the coordinate to draw pair-wise significance lines
    #for this we calculate bas as the maximum value
    bas<-max(df[(df$measure==k),"value"])
    #Calculate increments as 10% of the maximum values
    inc<-0.1*bas
    #Give an initial increment
    bas<-bas+inc
    for(l in 1:dim(s)[2]){
      #Do a pair-wise anova
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      #Ignore if anova fails
      if(!is.na(as.numeric(tmp[length(tmp)]))){
        #Only retain those pairs where the p-values are significant
        if(as.numeric(tmp[length(tmp)])<pValueCutoff){
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
          #Generate the next position
          bas<-bas+inc
        }
      }
    }
  }
  df_pw<-data.frame(row.names=NULL,df_pw)
  names(df_pw)<-c("measure","from","to","y","p")
  out <- list("df_pw"=df_pw, "df"=df)
  return(out)
}
