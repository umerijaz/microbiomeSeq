#' Fuzzy Set Ordination
#'
#' This function uses fuzzy set ordination to test effects of pertubation in
#' environmental variables to community structure.
#' It returns  an ordination plot which is annotated with a correlation between the original values of
#' the variable and the fuzzy values together with corresponding significance label.
#'
#' @param physeq(Required) A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method An integer specifying method for computing similarity indices, 1 = Baroni-Urbani & Buser, 2 = Horn and 3 = Yule .
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information.
#' @param step_across (optional) logical variable setting it to TRUE is for step-across correction which might improve the ordination
#' @param type (optional) An integer (1 or 2), for fso or mfso respectively.
#'
#' @return  A ggplot object of fuzzy set against original values which is
#' annotated with a correlation between them and a significance label.
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references \url{http://www.nku.edu/~boycer/fso/}
#'
#' @author  Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @examples
#' data(pitlatrine)
#' physeq <- taxa_level(pitlatrine, "Phylum")
#' p <- generateFSO(physeq, grouping_column="Country", method=1, type=1)
#' print(p)
#'
#'
#' @export generateFSO
#'
generateFSO<-function(physeq,grouping_column,method=1,indices=NULL,filename=NULL, type=1,step_across=F){
  data <- as.data.frame(otu_table(physeq))
  meta_table <- data.frame(sample_data(physeq))
  #pick only the numeric variables of metadata
  Env <- meta_table[,sapply(meta_table,is.numeric)]
  #create dummy column for grouping information
  meta_table$Groups <- meta_table[,grouping_column]
  #get indices of variables to investigate
  if(is.null(indices)){
    indices<-seq(1:dim(Env)[2])
  }
  sim <- sim.abund(data,method=method)
  dis.ho <- 1 - sim
  df<-NULL
  if(type==1){
    for(i in names(Env)[indices]){
      param.fso<-fso(Env[,i],dis.ho,permute=1000)
      tmp<-data.frame(mu=param.fso$mu,param=param.fso$data, Groups=as.factor(meta_table$Groups),
                      label=rep(paste(i,"(",round(param.fso$r,2)," ",formatPvalues(param.fso$p),")",sep=""),dim(data)[1]))
      if(is.null(df)){df<-tmp} else {df<-rbind(df,tmp)}
    }
  }
  else if(type==2){
    if(step_across){
      sim <- st.acr(sim)
    }
    whole.mfso<-mfso(as.formula(paste("~",paste(lapply(indices,function(x) paste("Env[,",x,"]",sep="")),collapse="+"))),dis=dis.ho,data=Env,scaling=2,permute=1000)
    whole.mfso$var<-names(Env)[indices]
    names(whole.mfso$data)<-names(Env)[indices]
    #print(whole.mfso)
    for (i in 1:(ncol(whole.mfso$mu) - 1)) {
      for (j in (i + 1):ncol(whole.mfso$mu)) {
        cat(paste("Processing",names(whole.mfso$data)[i]," and ",names(whole.mfso$data)[j],"\n"))
        tmp<-data.frame(x=whole.mfso$mu[, i],y=whole.mfso$mu[, j], Groups=as.factor(meta_table$Groups),label=rep(paste("x=mu(",names(whole.mfso$data)[i],")",", y=mu(",names(whole.mfso$data)[j],")",sep=""),dim(data)[1]))
        if(is.null(df)){df<-tmp} else {df<-rbind(df,tmp)}
      }
    }
  }
  p <- ggplot(df, aes(param, mu))+geom_point(aes(colour = Groups)) +geom_smooth(method="lm", size=1, se=T) +theme_bw()+facet_wrap( ~ label , scales="free", ncol=3)
  p<-p+theme(strip.background = element_rect(fill = "white"))
  return(p)
}
