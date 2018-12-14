
MyBUGSOutput <- function(xx,vars){
  x <- xx$sims.matrix
  OUT <- matrix(nrow = length(vars), ncol=6)
  j<-1
  for(i in vars){
    xi <- x[,i]
    OUT[j,3:4] <- quantile(xi, probs = c(0.025, 0.975))
    OUT[j,5:6] <- quantile(xi, probs = c(0.25, 0.75))
    OUT[j,1] <- mean(xi)
    OUT[j,2] <- sd(xi)
    j <- j + 1
  }
  colnames(OUT) <- c("mean", "se", "2.5%", "97.5%", "25%", "75%")
  rownames(OUT) <- vars
  OUT
}

################################################################

uNames <- function(k,Q){
  #Function to make a string of variables names of the form:
  #c("u[1]","u[2]", etc, "u[50]")
  #Q=50 knots were used
  String<-NULL
  for (j in 1:Q){String <- c(String, paste(k,"[",j,"]",sep = ""))}
  String
}