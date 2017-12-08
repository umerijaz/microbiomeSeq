#'Data normalisation
#'
#'This function normalises taxa abundance
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Required). A \code{method} to used to normalise the data. Available methods include:
#'        "edgeRnorm", "varstab", "randomsubsample", "proportion" .
#' @return Returns a phyloseq object whose abundance data is normalised.
#' @examples
#' data(pitlatrine)
#' physeq<-data(pitlatrine)
#' physeq <- normalise_data(physeq,method = "randomsubsample")
#' physeq <- normalise_data(physeq,method = "edgeRnorm")
#' physeq <- normalise_data(physeq,method = "proportion")
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export normalise_data
#'
normalise_data <- function(physeq, norm.method, ...){
  norm.method = match.arg(norm.method,c("edgernorm","varstab","randomsubsample",
                                        "proportion","relative","log-relative","scale"))
  switch(norm.method,
        "randomsubsample"=randomsubsample(physeq),
        "proportion"=proportion(physeq),
        "varstab"=deseq_varstab(physeq, ...),
        "edgernorm"=edgeRnorm(physeq, ...),
        "log-relative"=log_relative(physeq, ...),
        "relative"=relative(physeq, ...),
        "scale"=scale.meta(physeq, ...)
        )
}
