#'Alpha diversity measure
#'
#'This function calculates alpha diversity of provided community data using
#'selected indices/method.
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "richness", "fisher", "simpson", "shannon" and "evenness".
#' @return It returns a data frame of diversity measure and corresponding indices/methods
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export alpha_div
#'
alpha_div <- function(physeq,method){
  #==check for validity of selected methods
  method<- match.arg(method,c("richness", "fisher", "simpson", "shannon", "evenness","pd"), several.ok = TRUE)

  abund_table <- otu_table(physeq)
  df <- NULL
  if("richness"%in%method){
    R<- vegan::rarefy(abund_table,min(rowSums(abund_table)))
    df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))
    if(is.null(df)){
      df<-df_R}
    else {
      df<-rbind(df,df_R)}
  }
  if("fisher"%in%method){
    alpha <- vegan::fisher.alpha(abund_table)
    df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))
    if(is.null(df)){
      df<-df_alpha}
    else {
      df<-rbind(df,df_alpha)}
  }
  if("simpson"%in%method){
    simp <- vegan::diversity(abund_table, "simpson")
    df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))
    if(is.null(df)){
      df<-df_simp}
    else {
      df<-rbind(df,df_simp)}
  }
  if("shannon"%in%method){
    H<- vegan::diversity(abund_table)
    df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
    if(is.null(df)){
      df<-df_H}
    else {
      df<-rbind(df,df_H)}
  }
  if("evenness"%in%method){
    H<-vegan::diversity(abund_table)
    S <- specnumber(abund_table)
    J <- H/log(S)
    df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))
    if(is.null(df)){
      df<-df_J}
    else {
      df<-rbind(df,df_J)}
  }
  if("pd"%in%method){
    otu_tree <- phyloseq::phy_tree(physeq)
    PD <- pd(abund_table, otu_tree ,include.root = TRUE)
    df_PD<-data.frame(sample=names(PD),value=PD,measure=rep("PD",length(PD)))
    if(is.null(df)){
      df<-df_PD}
    else {
      df<-rbind(df,df_PD)}
  }
  return(df)
}
