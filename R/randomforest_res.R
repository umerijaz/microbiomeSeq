#'Importance of significant features.
#'
#'This function uses Breiman's random forest algorithm that can be used in unsupervised mode for assessing proximities among data points.
#'Random forest classifier is used to determine the importance of differentially expressed features/taxa to the microbial community.
#'The measure used in this case is Mean Descrease in Accuracy aso know as permutation importance which is
#'reported for each of the features. This is obtained by removing the relationship of a feature  and measuring increase in error.
#'Consequently, the feature with high mean decrease in accuracy is considered most important.
#'
#' @param data (Required). A \code{data.frame} of features to be assigned importance.
#' @param groups (Required). Vector corresponding to groups/conditions or source of variation in data.
#' @return Returns a \code{data.frame} of importance measure (mean decrease in accuracy) and ranks.
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},  Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export randomforest_res
#'

#it returns the measure of importance as mean decrease accuaracy
randomforest_res <- function(data, groups){
  IDs_map<-data.frame(row.names=colnames(data),"taxa"=colnames(data))
  val<-randomForest::randomForest(groups~ ., data=data, importance=T, proximity=T,ntree=1500,keep.forest=F)
  imp<- randomForest::importance(val)
  df_accuracy<-data.frame(row.names=NULL,Sample=rownames(imp),Value=abs(as.numeric(imp[,"MeanDecreaseAccuracy"])),Index=rep("Mean Decrease Accuracy",dim(imp)[1]))

  #Rearrange the features in terms of importance for ggplot2 by changing factor levels
  df_accuracy$Sample<-IDs_map[as.character(df_accuracy$Sample),"taxa"]
  df_accuracy_order<-as.character(IDs_map[rownames(imp),"taxa"][order(abs(as.numeric(imp[,"MeanDecreaseAccuracy"])),decreasing=T)])
  df_accuracy$Sample<-factor(as.character(df_accuracy$Sample),levels=df_accuracy_order)
  df_accuracy$rank <- base::rank(df_accuracy$Value, ties.method = "min")
  df_accuracy$rank <- max(df_accuracy$rank)-df_accuracy$rank+1

  out<-list("importance"=df_accuracy)
  return(out)
}


