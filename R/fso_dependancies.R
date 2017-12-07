#===fso dependencies===#

#We need a Pvalue formatter
formatPvalues <- function(pvalue) {
  ra<-""
  if(pvalue <= 0.1) ra<-"."
  if(pvalue <= 0.05) ra<-"*"
  if(pvalue <= 0.01) ra<-"**"
  if(pvalue <= 0.001) ra<-"***"
  return(ra)
}
#these are dependance functions for fso function as presented b
#Reference: http://www.nku.edu/~boycer/fso/
sim.binary <- function (df, method = NULL, diag=FALSE, upper=FALSE){
  df <- as.matrix(df)
  a <- df %*% t(df)
  b <- df %*% (1 - t(df))
  c <- (1 - df) %*% t(df)
  d <- ncol(df) - a - b - c
  #1=Jaccard, 2=Baroni-Urbani & Buser, 3=Kulczynski, 4=Ochiai, 5=Sorensen
  if (method == 1) {
    sim <- a/(a + b + c)
  }
  else if (method == 2) {
    sim <- (a + sqrt(a*d))/(a + b + c + sqrt(a*d))
  }
  else if (method == 3) {
    sim <- 0.5* (a/(a + b) + a/(a + c))
  }
  else if (method == 4) {
    sim <- a/sqrt((a + b) * (a + c))
  }
  else if (method == 5) {
    sim <- 2 * a/(2 * a + b + c)
  }
  sim2 <- sim[row(sim) > col(sim)]
  class(sim2) <- "dist"
  attr(sim2, "Labels") <- dimnames(df)[[1]]
  attr(sim2, "Diag") <- diag
  attr(sim2, "Upper") <- upper
  attr(sim2, "Size") <- nrow(df)
  attr(sim2, "call") <- match.call()
  
  return(sim2)
}

sim.abund <- function (df, method = NULL, diag = FALSE, upper = FALSE) 
{
  METHODS <- c("Baroni-Urbani & Buser", "Horn", "Yule (Modified)")
  if (!inherits(df, "data.frame")) 
    stop("df is not a data.frame")
  if (any(df < 0)) 
    stop("non negative value expected in df")
  if (is.null(method)) {
    cat("1 = Baroni-Urbani & Buser\n")
    cat("2 = Horn\n")
    cat("3 = Yule\n")
    cat("Select an integer (1-3): ")
    method <- as.integer(readLines(n = 1))
  }
  df <- as.matrix(df)
  sites <- nrow(df)
  species <- ncol(df)
  sim <- array(0, c(as.integer(sites),as.integer(sites)))
  spmax <- apply(df,2,max)
  
  if (method == 1) {
    #compute similarities (Baroni-Urbani & Buser)
    for (x in 1:sites) {
      for (y in 1:sites) {
        h1 <- 0
        h2 <- 0
        h3 <- 0
        for (i in 1:species) {
          h1 <- h1 + min(df[x,i],df[y,i])
          h2 <- h2 + max(df[x,i],df[y,i])
          h3 <- h3 + spmax[i] - max(df[x,i],df[y,i])
        }
        numer <- h1 + sqrt(h1*h3)
        denom <- h2 + sqrt(h1*h3)
        sim[x,y] <- ifelse(identical(denom,0), 0, numer/denom)
      }
    }
  }
  
  else if (method == 2) {
    #compute similarities (Horn)
    for (x in 1:sites) {
      for (y in 1:sites) {
        h1 <- 0
        h2 <- 0
        h3 <- 0
        for (i in 1:species) {
          if((df[x,i] + df[y,i]) > 0) h1 <- h1 + (df[x,i] + df[y,i]) * log10(df[x,i] + df[y,i])
          if(df[x,i] > 0) h2 <- h2 + df[x,i] * log10(df[x,i])
          if(df[y,i] > 0) h3 <- h3 + df[y,i] * log10(df[y,i])
        }
        x.sum <- sum(df[x,])
        y.sum <- sum(df[y,])
        xy.sum <- x.sum + y.sum
        if (identical(xy.sum, 0)) (sim[x,y] <- 0) else (sim[x,y] <- (h1 - h2 - h3)/(xy.sum * log10(xy.sum)-x.sum * log10(x.sum) - y.sum * log10(y.sum)))
        if (sim[x,y] < 1.0e-10) sim[x,y] <- 0
      }
    }
    
  }
  
  else if (method == 3) {
    #compute similarities (Yule)
    for (x in 1:sites) {
      for (y in 1:sites) {
        h1 <- 0
        h2 <- 0
        h3 <- 0
        h4 <- 0
        for (i in 1:species) {
          h1 <- h1 + min(df[x,i], df[y,i])
          h2 <- h2 + max(df[x,i], df[y,i]) - df[y,i]
          h3 <- h3 + max(df[x,i], df[y,i]) - df[x,i]
          h4 <- h4 + spmax[i] - max(df[x,i], df[y,i])
        }
        numer <- sqrt(h1*h4)
        denom <- sqrt(h1*h4) + sqrt(h2*h3)
        sim[x,y] <- ifelse(identical(denom,0), 0, numer/denom)
      }
    }
  }  
  
  else stop("Non convenient method")
  sim >- t(sim)
  sim2 <- sim[row(sim) > col(sim)]
  attr(sim2, "Size") <- sites
  attr(sim2, "Labels") <- dimnames(df)[[1]]
  attr(sim2, "Diag") <- diag
  attr(sim2, "Upper") <- upper
  attr(sim2, "method") <- METHODS[method]
  attr(sim2, "call") <- match.call()
  class(sim2) <- "dist"
  return(sim2)
}

st.acr <- function(sim, diag=FALSE, upper = FALSE)
{
  dis <- 1 - sim
  #use "shortest" or "extended"
  edis <- as.matrix(stepacross(dis, path = "shortest", toolong = 1))
  sim <- 1 - edis
  amax <- max(sim)
  amin <- min(sim)
  sim <- (sim-amin)/(amax-amin)
  sim2 <- sim[row(sim) > col(sim)]
  attr(sim2, "Size") <- nrow(sim)
  attr(sim2, "Labels") <- dimnames(sim)[[1]]
  attr(sim2, "Diag") <- diag
  attr(sim2, "Upper") <- upper
  attr(sim2, "call") <- match.call()
  class(sim2) <- "dist"
  return(sim2)
}
#/===dependencies===#
