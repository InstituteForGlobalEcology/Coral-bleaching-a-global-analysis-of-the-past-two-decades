GetCIs <- function(x) {
	#Support file to calculate the posterior
	#mean, se and 2.5 and 97.5 quantiles for 
	#fitted values. The object x needs to be 
	#an observation - by - MCMC iteration object.
    OUT <- matrix(nrow = nrow(x), ncol = 4) 
	for(i in 1:nrow(x)){
	  xi <- x[i,]	
      OUT[i,3:4] <- quantile(xi, probs = c(0.025, 0.975))
      OUT[i,1] <- mean(xi)
      OUT[i,2] <- sd(xi)
	}
	colnames(OUT) <- c("mean", "se", "lo", "up")
	OUT
}


MyBUGSChains <- function(xx, vars, PanelNames = NULL){
#Small function to make an xyplot of the iterations per chain,
#for each variable 
  x <- xx$sims.array
  idchain.All <- NULL
  x1.All <- NULL
  ChainLength.All <- NULL
  id.All <- NULL

  NumBerChains <- ncol(x[,,vars[1]])
 
  for (i in vars){
	x1          <- as.vector(x[,,i])
	id          <- rep(rep(i, length = nrow(x[,,i])),NumBerChains)
	idchain     <- rep(1:NumBerChains, each = nrow(x[,,i]))
    ChainLength <- rep(1: nrow(x[,,i]), NumBerChains)

    x1.All <- c(x1.All, x1)
    ChainLength.All <- c(ChainLength.All, ChainLength)
    id.All <- c(id.All, id)
    idchain.All <- c(idchain.All, idchain)
   }


  if (!is.null(PanelNames)) { 
     if (length(unique(id.All)) != length(PanelNames)) {stop("Wrong number of panel names")}
  	  AllNames <- unique(id.All)
  	  for (i in 1:length(AllNames)) {
  	  	id.All[id.All == AllNames[i]] <- PanelNames[i] 
  	  }
  	  id.All <- factor(id.All, levels = PanelNames)
  	  }


  Z <- xyplot(x1.All ~ ChainLength.All | factor(id.All) ,
       type = "l",
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       scales = list(x = list(relation = "same", draw = TRUE),
                     y = list(relation = "free", draw = TRUE)),
       groups = idchain.All,  col = 1:NumBerChains,
       xlab = list(label = "MCMC iterations", cex = 1.5),
       ylab = list(label = "Sampled values", cex = 1.5))
  print(Z)
}





##################################################################
MyBUGSACF <- function(Output  = Output, SelectedVar = SelectedVar, PanelNames = NULL){
#Small function to make an xyplot of the ACF per chain,
#for each variable 
#Output  = is the out object from JAGS
#SelectedVar is a character string of variables in xx
#PanelNames are matching names for the panels

  if (is.null(PanelNames)) { PanelNames <- SelectedVar }
  if (length(SelectedVar) != length(PanelNames)) {stop("Wrong number of panel names")}
  xx <- Output

  x <- xx$sims.array
  idchain.All <- NULL
  acf.Var.All <- NULL
  Lag.Var.All <- NULL
  id.All <- NULL

  NumBerChains <- ncol(x[,, SelectedVar[1]])
 
  for (i in SelectedVar){

    #Extract data from variable i and calculate ACF for each chain
	x1       <- x[,,i]
    acf.Vari <- Lag.Vari <- NULL
    for (j in 1:NumBerChains){
      acf.Varij <- acf(x1[,j], plot = FALSE)
      acf.Vari <- c(acf.Vari, acf.Varij$acf)
      Lag.Vari <- c(Lag.Vari, acf.Varij$lag)
    }
    NumLags <- length(acf.Varij$lag)

    ii <- PanelNames[SelectedVar==i]
    #print(i)
    #print(ii)
	id          <- rep(rep(ii, length = NumLags),NumBerChains)
	idchain     <- rep(1:NumBerChains, each = NumLags)
    

    acf.Var.All <- c(acf.Var.All, acf.Vari)
    Lag.Var.All <- c(Lag.Var.All, Lag.Vari)
    id.All <- c(id.All, id)
    idchain.All <- c(idchain.All, idchain)
   }

  Z <- xyplot(acf.Var.All ~ Lag.Var.All | factor(id.All) ,
       type = "l",
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       scales = list(x = list(relation = "same", draw = TRUE),
                     y = list(relation = "same", draw = TRUE)),
       groups = idchain.All,  col = 1:NumBerChains,
       xlab = list(label = "Lags",cex = 1.5),
       ylab = list(label = "ACF", cex = 1.5))
  print(Z)
}





##############################################################
MyBUGSHist <- function(Output  = Output, SelectedVar = SelectedVar, PanelNames = NULL){
#Small function to make an histogram of the ACF per chain.
#xx$BUGSoutput is the out object from JAGS
#vars is a character string of variables in xx
#PanelNames are matching names for the panels

#for each variable 
  x <- Output$sims.matrix
  AllParams <- NULL
  
  if (is.null(PanelNames)) { PanelNames <- SelectedVar }
  if (length(SelectedVar) != length(PanelNames)) {stop("Wrong number of panel names")}

 
  for (i in SelectedVar){
    #Extract data from variable i
	Paramsi <- x[,i]
    AllParams <- c(AllParams, Paramsi)	
    }

   #AllID <- rep(vars, each = nrow(x))
   AllID2 <- rep(PanelNames, each = nrow(x))
   AllID2 <- factor(AllID2, levels = PanelNames)


MyPanelCex <- 1.2
Z <- histogram( ~ AllParams | factor(AllID2),
           strip = strip.custom(bg = 'white',
           par.strip.text = list(cex = MyPanelCex)),
           type = "count" ,
           nint = 100,
           xlab = list(label = "Posterior distribution", cex = 1.5),
           col = gray(0.5), 
           ylab = list(label = "Frequencies", cex = 1.5),
           scales = list(alternating = FALSE, 
                         x = list(relation = "free"),
                         y = list(relation = "free")),
           breaks=NULL,              
           panel = function(x, ...) {
             panel.histogram(x, ...)
             panel.abline(v = 0, lwd = 3, col =2)
             CI <- quantile(x, probs = c(0.025,0.975))
             panel.arrows (CI[1],-2, CI[2],-2, col = 2, lwd= 7, length=0)
             })
  print(Z)
}

##################################################
MyBUGSOutput <- function(Output  = Output, SelectedVar = SelectedVar, VarNames = NULL){
	xx   <- Output
	vars <- SelectedVar
	
	if (is.null(VarNames)) { VarNames <- SelectedVar }
	if (length(SelectedVar) != length(VarNames)) {stop("Wrong number of variable names")}

	x <- xx$sims.matrix
    OUT <- matrix(nrow = length(vars), ncol=4) 
    j<-1
	for(i in vars){
	  xi <- x[,i]	
   	  OUT[j,3:4] <- quantile(xi, probs = c(0.025, 0.975))
   	  OUT[j,1]   <- mean(xi)
   	  OUT[j,2]   <- sd(xi)
   	  j          <- j + 1
	}
	colnames(OUT) <- c("mean", "se", "2.5%", "97.5%")
	rownames(OUT) <- VarNames
	OUT
}






#Plot 3 curve
#####################################################
My3DPlot <- function(BirdsS, K, KnotPositionsS, u, b, MyTitle ){
	
	#I am defining a mesh of MxM for the smoother f(X,Y). Use M = 40?
    M <- 40
    xM <- seq(from = min(BirdsS$Xkms),
              to  = max(BirdsS$Xkms),
              length = M)
    yM <- seq(from = min(BirdsS$Ykms),
              to  = max(BirdsS$Ykms),
              length = M)
  
    XX.MM <- expand.grid(xM, yM)           
    #On this grid I will predict the smoother for the output           

    #Also calculate the X and Z equivalents for my grid
    XM <- cbind(XX.MM[,1], XX.MM[,2])
    XZM <- GetXandZ(XM[,1], XM[,2], K, KnotPositionsS)
    #XZM$X is the X matrix
    #XZM$Z is the Z matrix

    NumMCMCIterations <- nrow(u)

    #Now calculate the smoother F for each MCMC iteration
    SmoothFMmcmc <- matrix(nrow = M*M, ncol = NumMCMCIterations)
    for (k in 1:NumMCMCIterations){
      SmoothFMmcmc[,k] <- #XZM$X %*% b[k,] + 
                          XZM$Z %*% u[k,]
	}

    #And calculate the posterior mean, and the 95% credible interval for F
    SF.Info <- matrix(nrow = M*M, ncol = 3)
    MtM <- M * M
    for (i in 1:MtM) {
	  SF.Info[i,2:3] <- quantile(SmoothFMmcmc[i,], probs=c(0.025, 0.975))
	  SF.Info[i,1] <- mean(SmoothFMmcmc[i,])	
	} 
	
   colnames(SF.Info) <- c("Mean", "Lower", "Upper")
   SF.Info

   #Skecth the smoother
   #MyCol <- rep("lightblue", M * M)
   #persp(x=xM, y=yM,
   #   z=matrix(SF.Info[,1], nrow = M, byrow = FALSE),
   #   scale = TRUE,
   #   theta = 20, phi = 30, expand = 1, 
   #   col = MyCol,
   #   ltheta = 120, shade = 0.75, ticktype = "detailed",
   #   xlab = "X", ylab = "Y", zlab = "Smoother") -> res
   #round(res, 3)
   
   #IPData <- BirdsS[BirdsS$inpark01==1, c("Xkms", "Ykms")]
   #Close <- rep(0,1600)
   #NI <- nrow(IPData)
   #for (i in 1:1600){
	#a <- as.numeric(XX.MM[i,])
	#dist.aANDPark <- sqrt((IPData[,1]-rep(a[1], NI))^2 + 
	#                      (IPData[,2]-rep(a[2], NI))^2)	
  	#Close[i] <- sum(dist.aANDPark < 0.02)   
   #}
   #MyCol2 <- rep("lightblue", M * M)
   #MyCol2[Close > 0] <- "yellow"
   #Close01 <- Close
   #Close01[Close01>1]<-1

   #nrz <- 40
   #ncz <- 40
   #Z <-matrix(Close01, ncol = M, byrow = FALSE)

   #jet.colors <- colorRampPalette( c("lightblue", "red") ) 
   # Generate the desired number of colors from this palette
   #nbcol <- 100
   #color <- jet.colors(nbcol)
   # Compute the z-value at the facet centres
   #zfacet <- Z[-1, -1] + Z[-1, -ncz] + Z[-nrz, -1] + Z[-nrz, -ncz]
   # Recode facet z-values into color indices
   #facetcol <- cut(zfacet, nbcol)

#   Graph <- 
   persp(x=xM, y=yM, z=matrix(SF.Info[,1], nrow = M, byrow = FALSE),
                 scale = TRUE,
                 theta = 30, phi = 40, expand = 1, 
                 #col = color[facetcol],
                 ltheta = 120, shade = 0.75, ticktype = "detailed",
                 xlab = "X", ylab = "Y", zlab = "Smoother",
                 main = MyTitle) 
                 #-> res
   #round(res, 3)
   #print(Graph)
   
#   pdf(file = paste(MyTitle,".pdf", sep = ""))
#   persp(x=xM, y=yM, z=matrix(SF.Info[,1], nrow = M, byrow = FALSE),
#                 scale = TRUE,
#                 theta = 30, phi = 40, expand = 1, 
#                 col = color[facetcol],
#                 ltheta = 120, shade = 0.75, ticktype = "detailed",
#                 xlab = "X", ylab = "Y", zlab = "Smoother",
#                 main = MyTitle) 
#   dev.off()
   
}
#################################################
Plot3DFit <- function(BirdsS,u,Whichb,XZX,XZZ,aBMCMC,
                      beta1,beta2,beta3,XPoisM){
	
	b1 <- out1$sims.matrix[,Whichb[1]]
    b2 <- out1$sims.matrix[,Whichb[2]]
    b3 <- out1$sims.matrix[,Whichb[3]]


    NumMCMCIterations <- nrow(u)
	SmoothFMmcmc <- matrix(nrow = nrow(BirdsS), ncol = NumMCMCIterations)
    for (k in 1:NumMCMCIterations){
    SmoothFMmcmc[,k] <- b1[k] * XZX[,1] + 
                        b2[k] * XZX[,2] + 
                        b3[k] * XZX[,3] + 
                        XZZ %*% u[k,]
		}

    ZIPFit <- matrix(nrow = nrow(BirdsS), ncol = NumMCMCIterations)
    for (k in 1:NumMCMCIterations){
     aB <- aBMCMC[k]
     ProbFalseZero <- exp(aB)/(1+exp(aB))	
     Park <-  beta1[k] * XPoisM[,1] + 
              beta2[k] * XPoisM[,2] +
              beta3[k] * XPoisM[,3]
    EtaPoisson <- Park + SmoothFMmcmc[,k]
    ZIPFit[,k] <- (1-ProbFalseZero) * exp(EtaPoisson )    
   }

   #And calculate the posterior mean, 
   #and the 95% credible interval for ZIPFit
   ZIP.Info <- matrix(nrow = nrow(BirdsS), ncol = 3)
   for (i in 1:nrow(BirdsS)) {
	ZIP.Info[i,2:3] <- quantile(ZIPFit[i,], probs=c(0.025, 0.975))
	ZIP.Info[i,1]   <- mean(ZIPFit[i,])
	} 	
   colnames(ZIP.Info) <- c("Mean", "Lower", "Upper")
   ZIP.Info

   library(MBA)
   #Sketch ZIPFit
   tmp.x<-seq(min(BirdsS$Xkms),max(BirdsS$Xkms),length.out=50)
   tmp.z<-seq(min(BirdsS$Ykms),max(BirdsS$Ykms),length.out=50)
   tmp.grd<-as.data.frame(expand.grid(tmp.x,tmp.z))
   names(tmp.grd)<-c("x","z")

   WhatToPlot <- ZIP.Info[,1]   #Fitted ZIP values

   pd <- data.frame(Fit1 = WhatToPlot, 
                 x = BirdsS$Xkms, 
                 z = BirdsS$Ykms)
   tmp.mba <-mba.points(pd[,c(2,3,1)], 
                        tmp.grd, n = 1, m = 1, h = 8, extend=FALSE)

   out.pd<-data.frame(Fit1 = tmp.mba$xyz.est[,3],
                     x=tmp.mba$xyz.est[,1],
                     z=tmp.mba$xyz.est[,2])

   xx <- tmp.x
   yy <- tmp.z
   zz <- matrix(out.pd$Fit1, nrow =50 )
   Graph <- drape.plot(xx, yy, zz, col=topo.colors(128),
                       theta = 50, phi = 10)-> pm 

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "OWEZ" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "OWEZ" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "OWEZ",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="magenta", pch=16, cex=1)

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "PAWP" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "PAWP" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "PAWP",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="red", pch=16, cex=1)

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "IJmuiden" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "IJmuiden" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "IJmuiden",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="yellow", pch=16, cex=1)
 
   print(Graph)
   
   quartz()
   MyX <- 240
   SS <- wireframe(Fit1 ~ x + z, data = out.pd,    shade = TRUE,
          aspect = c(61/87, 0.4),
          screen = list(z = 40, x = i, y = 3),
          light.source = c(0,0,10), distance = .5,
          shade.colors.palette = function(irr, ref, height, w = .5)
          grey(w * irr + (1 - w) * (1 - (1-ref)^.4)),
          zlab = list(label = "Fitted values", rot = 90))
   print(SS)
}
#################################################


#################################################
Plot3DFit2 <- function(BirdsS,u,Whichb,XZX,XZZ,aBMCMC,
                      beta1,beta2,beta3,XPoisM, PE, MyTitle){
	
	b1 <- PE[,Whichb[1]]
    b2 <- PE[,Whichb[2]]
    b3 <- PE[,Whichb[3]]


    NumMCMCIterations <- nrow(u)
	SmoothFMmcmc <- matrix(nrow = nrow(BirdsS), ncol = NumMCMCIterations)
    for (k in 1:NumMCMCIterations){
    SmoothFMmcmc[,k] <- b1[k] * XZX[,1] + 
                        b2[k] * XZX[,2] + 
                        b3[k] * XZX[,3] + 
                        XZZ %*% u[k,]
		}

    ZIPFit <- matrix(nrow = nrow(BirdsS), ncol = NumMCMCIterations)
    for (k in 1:NumMCMCIterations){
     aB <- aBMCMC[k]
     ProbFalseZero <- exp(aB)/(1+exp(aB))	
     Park <-  beta1[k] * XPoisM[,1] + 
              beta2[k] * XPoisM[,2] +
              beta3[k] * XPoisM[,3]
    EtaPoisson <- Park + SmoothFMmcmc[,k]
    ZIPFit[,k] <- (1-ProbFalseZero) * exp(EtaPoisson )    
   }

   #And calculate the posterior mean, 
   #and the 95% credible interval for ZIPFit
   ZIP.Info <- matrix(nrow = nrow(BirdsS), ncol = 3)
   for (i in 1:nrow(BirdsS)) {
	ZIP.Info[i,2:3] <- quantile(ZIPFit[i,], probs=c(0.025, 0.975))
	ZIP.Info[i,1]   <- mean(ZIPFit[i,])
	} 	
   colnames(ZIP.Info) <- c("Mean", "Lower", "Upper")
   ZIP.Info

   library(MBA)
   #Sketch ZIPFit
   tmp.x<-seq(min(BirdsS$Xkms),max(BirdsS$Xkms),length.out=50)
   tmp.z<-seq(min(BirdsS$Ykms),max(BirdsS$Ykms),length.out=50)
   tmp.grd<-as.data.frame(expand.grid(tmp.x,tmp.z))
   names(tmp.grd)<-c("x","z")

   WhatToPlot <- ZIP.Info[,1]   #Fitted ZIP values

   pd <- data.frame(Fit1 = WhatToPlot, 
                 x = BirdsS$Xkms, 
                 z = BirdsS$Ykms)
   tmp.mba <-mba.points(pd[,c(2,3,1)], 
                        tmp.grd, n = 1, m = 1, h = 8, extend=FALSE)

   out.pd<-data.frame(Fit1 = tmp.mba$xyz.est[,3],
                     x=tmp.mba$xyz.est[,1],
                     z=tmp.mba$xyz.est[,2])

   xx <- tmp.x
   yy <- tmp.z
   zz <- matrix(out.pd$Fit1, nrow =50 )
   Graph <- drape.plot(xx, yy, zz, col=topo.colors(128),
                       theta = 50, phi = 10, main = MyTitle)-> pm 

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "OWEZ" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "OWEZ" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "OWEZ",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="magenta", pch=16, cex=1)

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "PAWP" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "PAWP" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "PAWP",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="red", pch=16, cex=1)

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "IJmuiden" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "IJmuiden" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "IJmuiden",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="yellow", pch=16, cex=1)
 
   print(Graph)
   
   

   #Same graph...but now saved as a pdf   
    pdf(file = paste("Fit3D ",MyTitle,".pdf", sep = ""))
    Graph <- drape.plot(xx, yy, zz, col=topo.colors(128),
                       theta = 50, phi = 10)-> pm 

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "OWEZ" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "OWEZ" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "OWEZ",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="magenta", pch=16, cex=1)

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "PAWP" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "PAWP" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "PAWP",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="red", pch=16, cex=1)

   xxi <- BirdsS$Xkms[BirdsS$fInPark == "IJmuiden" ]
   yyi <- BirdsS$Ykms[BirdsS$fInPark == "IJmuiden" ]
   zzi <- ZIP.Info[BirdsS$fInPark == "IJmuiden",1]
   trans3d( xxi, yyi,zzi,pm)-> uv
   points( uv, col="yellow", pch=16, cex=1)
 
   print(Graph) 
   dev.off()

   
   quartz()   
   MyX <- 240
   SS <- wireframe(Fit1 ~ x + z, data = out.pd,    shade = TRUE,
          aspect = c(61/87, 0.4),
          screen = list(z = 40, x = i, y = 3),
          light.source = c(0,0,10), distance = .5,
          shade.colors.palette = function(irr, ref, height, w = .5)
          grey(w * irr + (1 - w) * (1 - (1-ref)^.4)),
          zlab = list(label = "Fitted values", rot = 90),
          main = MyTitle)
   print(SS)
}
#################################################






uNames <- function(k,Q){
  #Function to make a string of variables names of the form:
  #c("u[1]","u[2]", etc, "u[50]")	
  #Q=50 knots were used	
  String<-NULL
  for (j in 1:Q){String <- c(String, paste(k,"[",j,"]",sep = ""))}
  String
}


########################################################################
#Code below is based on:
#On semiparametric regression with O'Sullivan Penalized Splines
#Wand and Ormerod (2008)

OSullivan <- function(x, numIntKnots = 20, AddIntercept = FALSE, intKnots) {
  library(splines)  	
  #x is the covariate for the smoother
  #y is response variable
  a <- min(x) 
  b <- max(x)
  xg <- seq(a, b, length = 101)
   
  names(intKnots) <- NULL
  B <- bs(x, 
          knots = intKnots, 
          degree = 3,
          Boundary.knots = c(a, b),
          intercept = TRUE)


  #Create the Omega matrix:
  formOmega <- function(a,b,intKnots) {
    allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
    K <- length(intKnots) ; L <- 3*(K+8)
    xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
    wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
    Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
    Omega     <- t(Bdd*wts)%*%Bdd     
    return(Omega)
  }

  Omega <- formOmega(a,b,intKnots)
  eigOmega <- eigen(Omega)

  #Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$:
  indsZ <- 1:(numIntKnots+2)
  UZ    <- eigOmega$vectors[,indsZ]
  LZ    <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
   
  #Perform stability check:  
  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX    <- eigOmega$vectors[,indsX]   
  L     <- cbind( UX, LZ )
  stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
  if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
    print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")

  # Form the X and Z matrices:
  X <- cbind(rep(1,length(x)),x)
  Z <- B %*%LZ

  if (!AddIntercept) { 
  	  X <- X[,-1]
  	  X <- as.matrix(X, ncol = 1) }
 
  list(X = X, Z = Z)
}
###########################################################################

default.knots <- function(x,num.knots) {
  if (missing(num.knots)) num.knots <- max(5,min(floor(length(unique(x))/4),35))
  return(quantile(unique(x), seq(0,1, length= (num.knots+2))[-c(1,(num.knots+2))]))
}


#Function to calculate the low rank thin plate regression spline
#Run as:
#Z <- GetZ_LRTP(YourStandardizedCovariates, Knots)

GetZ_LRTP <- function(x, Knots) {
  Z_K            <- (abs(outer(x, Knots,"-")))^3
  OMEGA_all      <- (abs(outer(Knots, Knots,"-")))^3
  svd.OMEGA_all  <- svd(OMEGA_all)
  sqrt.OMEGA_all <-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
  Z <- t(solve(sqrt.OMEGA_all,t(Z_K)))
  return(Z)
  }

MySmoother <- function(f) {
 #Get the posterior mean, and 95% credible interval
 Sm1 <- matrix(nrow = nrow(f), ncol = 5)
 for (i in 1:nrow(f)) {
	Sm1[i,1:3] <- quantile(f[i,], 
                        probs = c(0.025, 0.5, 0.975))
    Sm1[i,4]  <- mean(f[i,])
	Sm1[i,5]  <- sd(f[i,])
	} 	
 colnames(Sm1) <- c("2.5%", "Median", "97.5%", "Mean", "SE")
 Sm1
 }


#####################################################################
#Standardize covariates
Mystd <- function(x) {(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)}


#####################################################################

MySmootherV2 <- function(f) {
 #Get the posterior mean, and 95% credible interval
 Sm1 <- matrix(nrow = nrow(f), ncol = 4)
 for (i in 1:nrow(f)) {
	Sm1[i,3:4] <- quantile(f[i,], 
                        probs = c(0.025, 0.975))
    Sm1[i,1]  <- mean(f[i,])
	Sm1[i,2]  <- sd(f[i,])
	} 	
 colnames(Sm1) <- c("Mean", "SE", "SeLo", "SeUp")
 Sm1
 }


