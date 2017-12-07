identify.num.variables <- function(physeq){
  meta_table <- sample_data(physeq)
  num.variables <- meta_table[,sapply(meta_table,is.numeric)]
  names.num.variables <- names(num.variables)
  return(names.num.variables)
}

get.num.variables <- function(physeq){
  meta_table <- sample_data(physeq)
  num.variables <- meta_table[,identify.num.variables(physeq)]
  not.num.variables <- meta_table[,!colnames(meta_table)%in%colnames(num.variables)]
  out <- list(num.variables=num.variables, notnum.variables=not.num.variables)
  return(out)
}

select.vars <- function(df1, df2, select.variables=NULL){
  if(!is.null(select.variables)){
    df1 <- df1[,colnames(df1)%in%select.variables]
  }
  out<-cbind(df1,df2)
  return(out)
}
