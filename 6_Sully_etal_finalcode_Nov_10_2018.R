# ------------------------------------------------------------------------
# Global bleaching patterns (Sully et al. Nature Climate change - submitted)
# ------------------------------------------------------------------------
# October 2018

# initialization ----------------------------------------------------------
library(rgdal)
library(R2jags)
library(ggplot2)
library(RColorBrewer)
library(plotrix)
library(GISTools)
library(R2OpenBUGS)

# import and format data --------------------------------------------------
# bleaching data
graphs_directory<-"C:/Users/Shannon/Desktop/Ecoregions/output"
#graphs_directory="C:/RobsR/Shannon/Bleaching data"
home <- getwd()
home<-"C:/Users/Shannon/Desktop/Ecoregions"
#home="C:/RobsR/Shannon/Bleaching data"
data <- read.csv(file=file.path(home,"Reef_Check_with_cortad_variables_with_annual_rate_of_SST_change.csv"), header=TRUE, sep=",")
data <- subset(data, Organism.Code=="Bleaching (% of population)")
data$Lat <- data$Latitude.Degrees; data$Long <- data$Longitude.Degrees

# calculate bleaching metrics
data$bleach.prez <- 0
data$bleach.prez[rowSums(data[c('S1','S2','S3','S4')],na.rm=T) >= 1] <- 1
data$avg.bleach <- rowSums(data[c('S1','S2','S3','S4')],na.rm=T)/rowSums(!is.na(data[c('S1','S2','S3','S4')]))
data$avg.bleach[data$avg.bleach < 1] <-0

# temperature
data$Temp <- data$Temperature_Kelvin - 273.15
data <- data[!is.na(data$Temperature_Kelvin),]

# shapefiles
###ecor.sg <- readOGR(file.path(home,'data','shapefiles','ecoregion_FixGeom_Simp1km_PC150.shp')) # ecoregions
###wlrd.p <- readOGR(file.path(home,'data','shapefiles','TM_WORLD_BORDERS_SIMPL_PC150.shp')) # from NASA (https://github.com/nasa/World-Wind-Java) and reprojected with Pacific in center
ecor.sg <- readOGR(file.path(home,'shapefiles','ecoregion_exportPolygon.shp')) # ecoregions
wlrd.p <- readOGR(file.path(home,'shapefiles','TM_WORLD_BORDERS_SIMPL_PC150.shp'))


# predictors by ecoregion
ecor_data <- read.csv(file.path(home,"EcoRegions_mean_variables.csv"))

setwd(home)
source(file= "MyBUGSOutput.R")

# Figure 1 ----------------------------------------------------------------
# maptools version
pal <- c( "#3B9AB2","#6BB0C0", "#78B7C5", "gold", "#EBCC2A", "#E1AF00","darkorange", "#F21A00", "darkred")

tiff(file=file.path(home,'output','Figure1.tif'),height=800,width=3000,res=300)
par(mgp=c(0.5,0.6,0), mar=c(1,1,1,1))
plot(wlrd.p,ylim=c(-4400000,4400000),xlim=c(-2000000,2000000), col='grey90',border='grey70')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60° '),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
box()

xy <- data[data$avg.bleach == 0,c('Long','Lat')]
xy <- SpatialPointsDataFrame(data=xy,coords=xy[c('Long','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy, cex=.7)

temp <- subset(data, avg.bleach > 0)
temp <- temp[with(temp, order(temp$avg.bleach)),]
xy <- temp[c('Long','Lat')]
xy <- SpatialPointsDataFrame(data=xy,coords=xy[c('Long','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy, cex=.7)
points(xy, cex=.9, pch=19, col=pal[temp$avg.bleach])
windowsFonts(A=windowsFont('Arial Unicode MS'))
text(-7868896,-2922012,'Indian Ocean',cex=1.0,family='A')
text(9438742,487176,'Pacific Ocean',cex=1.0,family='A')
north.arrow(x=(-16654136+111319.4*320), y=1615153*2, len=(111319.4*2), lab="N", cex=.7)

#legend
plotrix::color.legend(9684797.171+25e5,-28*111319.4666666667,15807371.62+25e5,-23.5*111319.4666666667,legend=c(1,100),rect.col=pal,cex=1)
points(9239521.171+25e5,-25.75*111319.4666666667,lwd=1.5)
text(9239521.171+25e5,-20.3*111319.4666666667,'0')
text(((15807371.62+25e5)-(9684797.171+25e5))/2+(9684797.171+25e5),-18*111319.4666666667,"Bleaching %")
dev.off()

# Figure 2 ----------------------------------------------------------------

data <- data[!is.na(data$Diversity),]

# standardize explanatory variables
X_raw <- data[c('Diversity','Latitude.Degrees','Year','Depth','ClimSST','Temperature_Kelvin','Temperature_Maximum','Temperature_Kelvin_Standard_Deviation','SSTA','SSTA_Minimum','SSTA_Maximum','SSTA_DHW','SSTA_Frequency','SSTA_Frequency_Standard_Deviation','TSA_Frequency','TSA_Frequency_Standard_Deviation','TSA_DHW_Standard_Deviation','rate_of_SST_change')]
X_raw$Latitude.Degrees <- abs(X_raw$Latitude.Degrees)

standardize_function<-function(x){
  x.standardized=(x-mean(na.omit(x)))/sd(na.omit(x))
  return(x.standardized)
}

#Diversity
X_raw$Diversity.standardized<-standardize_function(X_raw$Diversity)

#Latitude
#X_raw$Latitude<-X_raw$Latitude.Degrees
X_raw$Latitude<-abs(X_raw$Latitude.Degrees)
X_raw$Latitude.standardized<- standardize_function(X_raw$Latitude)
#ClimSST
X_raw$ClimSST.standardized<-standardize_function(X_raw$ClimSST)
#Temperature
X_raw$Temperature_Kelvin.standardized<- standardize_function(X_raw$Temperature_Kelvin)
#Temperature_Maximum
X_raw$Temperature_Maximum.standardized<- standardize_function(X_raw$Temperature_Maximum)
#Depth
X_raw$Depth.standardized<-standardize_function(X_raw$Depth)
#Year
X_raw$Year.standardized<-standardize_function(X_raw$Year)
#Latitude
X_raw$Latitude.standardized<-standardize_function(X_raw$Latitude.Degrees)
#SSTA
X_raw$SSTA.standardized<- standardize_function(X_raw$SSTA)
#SSTA_Minimum
X_raw$SSTA_Minimum.standardized<- standardize_function(X_raw$SSTA_Minimum)
#SSTA_Maximum
X_raw$SSTA_Maximum.standardized<- standardize_function(X_raw$SSTA_Maximum)
#SSTA_DHW
X_raw$SSTA_DHW.standardized<-standardize_function(X_raw$SSTA_DHW)
#SSTA_Frequency
X_raw$SSTA_Frequency.standardized<-standardize_function(X_raw$SSTA_Frequency)
#TSA_Frequency
X_raw$TSA_Frequency.standardized<- standardize_function(X_raw$TSA_Frequency)
#rate_of_SST_change
X_raw$rate_of_SST_change.standardized<- standardize_function(as.numeric(X_raw$rate_of_SST_change))


X <- model.matrix(~ Diversity.standardized + Latitude.standardized + Year.standardized + 
                    Depth.standardized + ClimSST.standardized + 
                    Temperature_Kelvin.standardized + 
                    Temperature_Maximum.standardized + 
                    Temperature_Kelvin_Standard_Deviation + 
                    SSTA.standardized + 
                    SSTA_Minimum.standardized + SSTA_Maximum.standardized + 
                    SSTA_DHW.standardized + SSTA_Frequency.standardized + 
                    SSTA_Frequency_Standard_Deviation + 
                    TSA_Frequency.standardized + 
                    TSA_Frequency_Standard_Deviation + 
                    TSA_DHW_Standard_Deviation +
                    rate_of_SST_change.standardized, 
                  data =X_raw)

K <- ncol(X)

### New part
Site <- as.factor(as.numeric(as.factor(as.numeric(data$Reef.ID))))
Nre<- length(unique(Site))

#do average bleaching instead of binary bleaching
win.data <- list(Y    = round(data$avg.bleach,0),
                 N    = nrow(data),
                 X    = X,
                 K    = K,
                 #LogN = log(nrow(Bleaching_Data_with_cortad_variables_ecoregion)),
                 # New part
                 Site=Site, #Random effects identification 
                 Nre=Nre # Number of random effects
                 #End of New part
)

# standardize_function<-function(x){
#   x.standardized=(x-mean(na.omit(x)))/sd(na.omit(x))
#   return(x.standardized)
# }
# 
# X_standardized <- X_raw; for(i in 1:ncol(X_raw)) X_standardized[,i] <- standardize_function(X_raw[,i])
# 
# X <- model.matrix(~ as.matrix(X_standardized))

# data inputs
# win.data <- list(Y    = round(data$avg.bleach,0),
#                  N    = nrow(data),
#                  X    = X,
#                  K    = ncol(X),
#                  Site = as.numeric(as.factor(as.character(data$Reef.ID))), 
#                  Nre = length(unique(as.numeric(as.factor(as.character(data$Reef.ID)))))
# )

# JAGS model code
sink("GLMM.txt")
cat("
    model{
    #1. Priors
    for (i in 1:K) { beta[i]  ~ dnorm(0, 0.0001) }
    
    #Priors for random intercepts
    for (i in 1:Nre) {a[i] ~ dnorm(0, tau)}
    
    #Prior for tau_Site
    sigma_Site ~ dgamma(0.001, 0.001)
    tau <- 1 / (sigma_Site * sigma_Site)
    
    #Priors for k (size)
    size ~ dunif(0, 20)
    
    #2. Likelihood
    for (i in 1:N) {
    Y[i] ~  dnegbin(p[i], size)
    p[i] <- size / (size + mu[i]) 
    log(mu[i]) <- eta[i]
    eta[i]     <- inprod(beta[], X[i,]) + a[Site[i]]
    
    
    #3. Model checks
    #Pearson residuals
    Exp[i] <- mu[i]
    Var[i] <- mu[i] + mu[i] * mu[i] / size
    E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])
    
    #Simulated data with mean/variance taken from the fitted model
    #See text under block B, below.
    YNew[i] ~  dnegbin(p[i], size)                     
    
    #Pearson residual for predicted data     
    ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i])
    
    #Squared residuals
    D[i]    <- pow(E[i], 2)
    DNew[i] <- pow(ENew[i], 2)
    }         
    
    #Sum of squared Pearson residuals:
    Fit     <- sum(D[1:N])
    #Sum of squared predicted Pearson residuals:  
    FitNew  <- sum(DNew[1:N]) 
    
    #Posterior predictive checks
    Y.mean <- mean(Y[])
    YNew.mean <- mean(YNew[])
    pval.mean <- step(YNew.mean-Y.mean) 
    
    #R-squared
    varF <- sd(Y)^2
    varE <- sd(Y - YNew)^2
    R2 <- varF/(varF+varE)  
    }
    ",fill = TRUE)
sink()

# run model
inits  <- function () {
  list(beta = rnorm(ncol(X), 0, 0.1),
       size = runif(1, 0, 20 ) )}
params <- c("beta", "E", "Fit", "FitNew", "YNew", "pval.mean", "R2")
params <- c("beta")

G1 <- jags(data       = win.data,
           inits      = inits,
           parameters = params,
           model      = "GLMM.txt",
           n.thin     = 1,
           n.chains   = 3,
           n.burnin   = 4000,
           n.iter     = 5000)

vars <- c("beta[1]", "beta[2]","beta[3]","beta[4]","beta[5]", "beta[6]", "beta[7]", "beta[8]", "beta[9]", "beta[10]", "beta[11]", "beta[12]", "beta[13]", "beta[14]", "beta[15]", "beta[16]", "beta[17]","beta[18]", "beta[19]")

OUT1 <- MyBUGSOutput(G1$BUGSoutput, vars)
print(OUT1, digits =3)
summary(G1)

# coef plot
#double check down vs up, the labels got flipped?
G2_df=data.frame(variableG2=c("Intercept", "Diversity", "Latitude", "Year", "Depth", "ClimSST", "SST", "SST_Max", "SST_stdev", "SSTA", "SSTA_Min", "SSTA_Max", "SSTA_DHW", "SSTA_Freq", "SSTA_Freq_stdev", "TSA_Freq", "TSA_Freq_stdev", "TSA_DHW_stdev", "Rate_of_SST_change"),MeanG2=OUT1[,1],Down=OUT1[,3],Up=OUT1[,4],Down_quarter=OUT1[,5], Up_quarter=OUT1[,6])
#write.csv(G2_df, file = file.path(home,'output','beta_est_20181208.csv'),row.names=F)
#G2_df<-read.csv(file.path(home,'output','beta_est_20181208.csv'))
G2_df$color <- ("white")
G2_df$color[(G2_df$MeanG2 > 0) & (G2_df$Down>0)] <- 'red'
G2_df$color[(G2_df$MeanG2 < 0) & (G2_df$Up<0)] <- 'blue'

tiff(file=file.path(home,'output','Figure_2_positive_intercept.tif'),height=2000,width=2700,res=300)
ggplot(G2_df,aes(x=reorder(variableG2, MeanG2), MeanG2)) +
  geom_errorbar(aes(ymax=G2_df$Up, ymin=G2_df$Down), width=0) +
  geom_errorbar(aes(ymax=G2_df$Up_quarter, ymin=G2_df$Down_quarter), width=0, size=1.3) +
  geom_point(pch=21, size=3, fill=G2_df$color, color="black") +
  coord_flip() +
  theme_grey(base_size=15) +
  guides(colour=FALSE)+
  geom_hline(yintercept=0, linetype="dashed", color="gray") +
  labs(y=expression(paste("Estimated ",gamma," coefficients")), x="")
dev.off()

# Model checks
#MyBUGSOutput(G1$BUGSoutput, 'R2')
#MyBUGSOutput(G1$BUGSoutput, 'pval.mean')


# Figure 3 ----------------------------------------------------------------
ecor.sg.fil <-subset(ecor.sg, ERG=='ERG001' | ERG=='ERG002' | ERG=='ERG003' |
                       ERG=='ERG005' | ERG=='ERG006' | ERG=='ERG007'
                     | ERG=='ERG009' | ERG=='ERG010' | ERG=='ERG012'
                     | ERG=='ERG013' | ERG=='ERG016' | ERG=='ERG017' 
                     | ERG=='ERG018' | ERG=='ERG020' | ERG=='ERG021' 
                     | ERG=='ERG025' | ERG=='ERG026' | ERG=='ERG028' 
                     | ERG=='ERG031' | ERG=='ERG032' | ERG=='ERG033' 
                     | ERG=='ERG034' | ERG=='ERG035' | ERG=='ERG036' 
                     | ERG=='ERG037' | ERG=='ERG038' | ERG=='ERG041' 
                     | ERG=='ERG042' | ERG=='ERG044' | ERG=='ERG045' 
                     | ERG=='ERG046' | ERG=='ERG047' | ERG=='ERG048' 
                     | ERG=='ERG049' | ERG=='ERG050' | ERG=='ERG051' 
                     | ERG=='ERG053' | ERG=='ERG055' | ERG=='ERG056' 
                     | ERG=='ERG057' | ERG=='ERG058' | ERG=='ERG061' 
                     | ERG=='ERG062' | ERG=='ERG065' | ERG=='ERG066' 
                     | ERG=='ERG069' | ERG=='ERG070' | ERG=='ERG071' 
                     | ERG=='ERG072' | ERG=='ERG073' | ERG=='ERG074' | ERG=='ERG102'
                     | ERG=='ERG077' | ERG=='ERG078' | ERG=='ERG080' 
                     | ERG=='ERG081' | ERG=='ERG091' | ERG=='ERG093' 
                     | ERG=='ERG106' | ERG=='ERG107' 
                     | ERG=='ERG108' | ERG=='ERG111' | ERG=='ERG112' 
                     | ERG=='ERG113' | ERG=='ERG114' | ERG=='ERG115'
                     | ERG=='ERG116' | ERG=='ERG123' | ERG=='ERG129' 
                     | ERG=='ERG136' | ERG=='ERG137' | ERG=='ERG138' 
                     | ERG=='ERG139' | ERG=='ERG140' | ERG=='ERG141' 
                     | ERG=='ERG144' | ERG=='ERG146')
cols<-c(rev(brewer.pal(5,"Blues")), brewer.pal(5,"Reds"))


cat("model{
    for (i in 1:N){  
    env[i] ~ dnorm(0,tau) 
    m[i] ~ dbin(p[i],n[i]) 
    p[i] <- 1-exp(-h[i]*dt) 
    log(h[i]) <- Beta1+ Beta2[i]*(temp[i]-mean(temp[])/sd(temp[])) + env[i] 
    Beta2[i]~dnorm(0.0, 1.0E-6) I(-2,2)
    }
    
    B_slp <- sum(Beta2[])/N # relationship; slope of Beta
    
    #Priors
    Beta1 ~ dnorm (0, 1.0E-6)
    tau <- 1/sd/sd
    sd ~ dunif(0, 0.5)
    
    }",file='hazard.txt')


ecor_data$sst_per_ecoregion_standardized<-standardize_function(ecor_data$sst_per_ecoregion)
ecor_data$diversity_standardized<-standardize_function(ecor_data$diversity)


### B - SST_stdev
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$tsa_freq_per_ecoregion, dt=1)
##hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$sst_per_ecoregion, dt=1)
##hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$sst_per_ecoregion_standardized, dt=1)
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$ssta_dhw_per_ecoregion, dt=1)
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$rate_of_sst_change_per_ecoregion, dt=1)
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$sst_stdev_per_ecoregion, dt=1)
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$ssta_freq_stdev_per_ecoregion, dt=1)
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$tsa_freq_stdev_per_ecoregion, dt=1)
#hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$tsa_dhw_stdev_per_ecoregion, dt=1)
hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$diversity, dt=1)
hdata <- list(N=c(77), m=round(ecor_data$average_bleaching_mean_per_ecoregion,0), n=ecor_data$n, temp=ecor_data$diversity_standardized, dt=1)


params<- c("Beta2")
hazard.jags <- jags(data=hdata,
                    parameters.to.save=params,
                    inits=list(c(Beta1=1), c(Beta1=.5), c(Beta1=.1)),
                    model.file="hazard.txt",
                    n.chains=3,
                    n.thin=1,
                    n.burnin=500,
                    n.iter=5000)
betas_mean <- hazard.jags$BUGSoutput$mean$Beta2
betas_median <- hazard.jags$BUGSoutput$median$Beta2

#G2_df=data.frame(variableG2=c(seq(1:77)),MeanG2=hazard.jags$BUGSoutput$summary[1:77,1],Up=hazard.jags$BUGSoutput$summary[1:77,7],Down=hazard.jags$BUGSoutput$summary[1:77,3])
G2_df=data.frame(variableG2=c(seq(1:77)),MedianG2=hazard.jags$BUGSoutput$summary[1:77,5],Up=hazard.jags$BUGSoutput$summary[1:77,7],Down=hazard.jags$BUGSoutput$summary[1:77,3])
G2_df_full=data.frame(hazard.jags$BUGSoutput$summary)

#brks<-seq(from=(-1*(max(abs(min(betas_mean)), abs(max(betas_mean))))), to=(max(abs(min(betas_mean)), abs(max(betas_mean)))), length.out=10)
#grps<-(cut(betas_mean, brks, include.lowest=TRUE))
brks<-seq(from=(-1*(max(abs(min(betas_median)), abs(max(betas_median))))), to=(max(abs(min(betas_median)), abs(max(betas_median)))), length.out=10)
grps<-(cut(betas_median, brks, include.lowest=TRUE))

#tiff(file=file.path(home,'output','Figure_3_JAGS_tsa_freq_per_ecoregion_median.tif'),height=800,width=3000,res=300)
##tiff(file=file.path(home,'output','Figure_3_JAGS_sst_per_ecoregion_median.tif'),height=800,width=3000,res=300)
##tiff(file=file.path(home,'output','Figure_3_JAGS_sst_per_ecoregion_standardized_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_ssta_dhw_per_ecoregion_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_sst_rate_of_change_per_ecoregion_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_sst_stdev_per_ecoregion_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_ssta_freq_stdev_per_ecoregion_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_tsa_freq_stdev_per_ecoregion_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_tsa_dhw_stdev_per_ecoregion_median.tif'),height=800,width=3000,res=300)
#tiff(file=file.path(home,'output','Figure_3_JAGS_diversity.tif'),height=800,width=3000,res=300)
tiff(file=file.path(home,'output','Figure_3_JAGS_diversity_standardized.tif'),height=800,width=3000,res=300)
par(mgp=c(0.5,0.6,0), mar=c(1,1,1,1))
plot(wlrd.p,ylim=c(-4400000,4400000),xlim=c(-2000000,2000000), col='grey90',border='grey70')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60°'),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
box()
plot(ecor.sg.fil,col=cols[unclass(grps)],add=T)
windowsFonts(A=windowsFont('Arial Unicode MS'))
text(-7868896,-2922012,'Indian Ocean',cex=1.0,family='A')
text(9438742,487176,'Pacific Ocean',cex=1.0,family='A')
plotrix::color.legend(9684797.171+25e5,-28*111319.4666666667,15807371.62+25e5,-23.5*111319.4666666667,legend=c(round(min(brks), digits=1),round(max(brks), digits=1)),rect.col=cols,cex=1)
north.arrow(x=(-16654136+111319.4*320), y=1615153*2, len=(111319.4*2), lab="N", cex=.7)
dev.off()


#params<- c("h", "p", "Beta2", "B")
params<- c("Beta2")
model.file<- system.file(package="R2OpenBUGS", "model", "Hazard Function.txt")
hazard.bugs<- bugs(n.burnin=5000,
                    data=hdata,
                    inits=list(Beta1=c(1), Beta1=c(.5), Beta1=c(.1)),
                    params,
                    model.file=file.path(home,'Hazard Function.txt'),
                    n.chains=3,
                    n.iter=50000, debug=FALSE,DIC=TRUE)
betas_mean <- hazard.bugs$summary[1:77,1]
betas_median <- hazard.bugs$summary[1:77,5]

# coef plot
#G2_df=data.frame(variableG2=c(seq(1:77)),MeanG2=hazard.bugs$summary[1:77,1],Up=hazard.bugs$summary[1:77,7],Down=hazard.bugs$summary[1:77,3])
G2_df=data.frame(variableG2=c(seq(1:77)),MedianG2=hazard.bugs$summary[1:77,5],Up=hazard.bugs$summary[1:77,7],Down=hazard.bugs$summary[1:77,3])
G2_df_full=data.frame(hazard.bugs$summary)
write.csv(G2_df, file = file.path(graphs_directory,'BUGS_tsa_dhw_stdev_per_ecoregion.csv'),row.names=F)
write.csv(G2_df_full, file = file.path(graphs_directory,'BUGS_tsa_dhw_stdev_stdev_per_ecoregion_all_values.csv'),row.names=F)


G2_df$color <- ("white")
G2_df$color[G2_df$Up < 0] <- 'blue'
G2_df$color[G2_df$Down > 0] <- 'red'

tiff(file=file.path(graphs_directory,'BUGS_tsa_freq_per_ecoregion.tif'),height=2000,width=2700,res=300)
ggplot(G2_df,
       aes(x=reorder(variableG2, MedianG2), MedianG2)) +
       #aes(x=reorder(variableG2, MeanG2), MeanG2)) +
  geom_errorbar(aes(ymax=G2_df$Up, ymin=G2_df$Down), width=0) +
  geom_point(size=3, color=G2_df$color) +
  coord_flip() +
  theme_grey(base_size=15) +
  guides(colour=FALSE)+
  geom_hline(yintercept=0, linetype="dashed", color="gray") +
  labs(y=expression(paste("Estimated ",gamma," coefficients")), x="")
dev.off()


#brks<-seq(from=(-1*(max(abs(min(betas_mean)), abs(max(betas_mean))))), to=(max(abs(min(betas_mean)), abs(max(betas_mean)))), length.out=10)
#grps<-(cut(betas_mean, brks, include.lowest=TRUE))

brks<-seq(from=(-1*(max(abs(min(betas_median)), abs(max(betas_median))))), to=(max(abs(min(betas_median)), abs(max(betas_median)))), length.out=10)
grps<-(cut(betas_median, brks, include.lowest=TRUE))


#tiff(file=file.path(home,'output','Figure_3_BUGS_tsa_freq_per_ecoregion_median.tif'),height=800,width=3000,res=300)
tiff(file=file.path(home,'output','Figure_3_BUGS_sst_per_ecoregion_median.tif'),height=800,width=3000,res=300)
par(mgp=c(0.5,0.6,0), mar=c(1,1,1,1))
plot(wlrd.p,ylim=c(-4400000,4400000),xlim=c(-2000000,2000000), col='grey90',border='grey70')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°E','180°','60°W'),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°N','0°','23°S'),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
box()
plot(ecor.sg.fil,col=cols[unclass(grps)],add=T)
#text(-18987920,37*111319.4666666667,'a',cex=1.3)
windowsFonts(A=windowsFont('Arial Unicode MS'))
text(-7868896,-2922012,'Indian Ocean',cex=0.7,family='A')
text(9438742,487176,'Pacific Ocean',cex=0.7,family='A')
plotrix::color.legend(9684797.171+25e5,-28*111319.4666666667,15807371.62+25e5,-23.5*111319.4666666667,legend=c(round(min(brks), digits=1),round(max(brks), digits=1)),rect.col=cols,cex=1)
north.arrow(x=(-16654136+111319.4*320), y=1615153*2, len=(111319.4*2), lab="N", cex=.7)
dev.off()


# Figure 4 ----------------------------------------------------------------
Bpresent <- subset(data, Average_bleaching > 0) 
Bpres1999 <- subset(Bpresent, Year > 1999) 
Bpres2002 <- subset(Bpresent, Year > 2001) 

plot(Average_bleaching ~ Year, data=Bpres2002)
plot(Average_bleaching ~ as.factor(Year), data=Bpres2002)

tglm <- glm.nb(Average_bleaching ~ as.factor(Year), data=Bpres2002)
newdat <- data.frame(Year=as.factor(unique(Bpres2002$Year)))
tglm.pred <- predict(tglm,newdata=newdat,se.fit=T)
plot(as.numeric(newdat$Year),tglm.pred$fit,type='p')
plotCI(as.numeric(newdat$Year),tglm.pred$fit,ui=tglm.pred$fit+1.96*tglm.pred$se.fit,li=tglm.pred$fit-1.96*tglm.pred$se.fit)

tglm <- glm.nb(Average_bleaching ~ as.numeric(Year), data=Bpres2002)
summary(tglm)
tglm.pred <- predict(tglm,newdata=newdat,se.fit=T)
plot(Average_bleaching ~ as.numeric(as.factor(Bpres2002$Year)), data=Bpres2002)
points(as.numeric(newdat$Year),tglm.pred$fit,type='l',col='red')

all2002 <- subset(data, Year > 2001)
all2002$Year.num <- as.numeric(as.factor(all2002$Year))
tglm <- glm(bleach.prez ~ Year.num,family='binomial', data=all2002)
summary(tglm)
newdat <- data.frame(Year.num=unique(all2002$Year.num),Year=unique(all2002$Year))
newdat <- arrange(newdat,Year)
tglm.pred <- predict(tglm,newdata=newdat[c('Year.num')],se.fit=T)
library(boot)
plot(newdat$Year,inv.logit(tglm.pred$fit),type='l',col='red')

tiff(file=file.path(home,'Figure4_alt.tif'),height=1800,width=2200,res=300)
par(mfrow=c(1,1),mgp=c(2.2,0.7,0),mar=c(3.5,3.5,2,3.5))
plot(Average_bleaching ~ as.factor(Year), data=Bpres2002, outline=FALSE, lty=1, staplewex=0, boxwex=0.8, boxlwd=1, medlwd=1,col='grey90',xlab='Year',ylab='% Bleaching')
par(new=T)
plot(newdat$Year.num,inv.logit(tglm.pred$fit),type='l',col='red',xaxt='n',yaxt='n',lwd=2,ylab='',xlab='',ylim=c(0.25,.7))
polygon(c(newdat$Year.num,rev(newdat$Year.num)),c(inv.logit(tglm.pred$fit+1.96*tglm.pred$se.fit),rev(inv.logit(tglm.pred$fit-1.96*tglm.pred$se.fit))),border=F,col=rgb(255,0,0,100,max=255))
axis(4)
mtext('Probability Bleaching',side=4,line=2)
dev.off()


# Figure 5 ----------------------------------------------------------------
tiff(file=file.path(home,'Figure5.tif'),height=1800,width=4000,res=300)
par(mfrow=c(1,2),mgp=c(2.2,0.7,0))
# A - temperature at bleaching over time

plot(Bpres1999$Year,Bpres1999$Temp, cex=1, cex.axis=1.5, cex.lab=1.5,ylab=expression("Bleaching temperature "~degree~C), xlab="Year",pch=21,bg='white',col='darkgrey',bty='l',ylim=c(18,34.5))

mod <- lm(Temp ~ Year, data = Bpres1999)
summary(mod)
newx <- data.frame(Year=seq(2000, 2017, length.out=100))
preds <- predict(mod, newdata = newx, 
                 interval = 'confidence')
polygon(c(rev(newx$Year), newx$Year), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
lines(newx$Year, preds[ ,1], lwd=2, col = 'black')
lines(newx$Year, preds[ ,3], lwd=2, lty = 'dashed', col = 'black')
lines(newx$Year, preds[ ,2], lwd=2,lty = 'dashed', col = 'black')

text(2000,34.5,'A',cex=1.5)

# B - histogram of temperature at bleaching before and after 2007
BpresEarly=subset(Bpresent, Year < 2007 & Temp > 24, select=c(Temp, Average_bleaching, Year)) 
BpresLate=subset(Bpresent, Year >= 2007 & Temp > 24, select=c(Temp, Average_bleaching, Year)) 

hist(BpresEarly$Temp,col=rgb(107,176,192,100,max=255), xlab=expression("Bleaching temperature"~degree~C), ylab= "Probablity density", main="", freq=F, ylim= c(0, 0.3), cex.lab=1.5, cex.axis=1.5)
hist(BpresLate$Temp, freq=F, col=rgb(242,26,0,100,max=255), main="", add=T)

dearly <- fitdistr(BpresEarly$Temp,'weibull')
x=seq(24,34,.01) 
y=dweibull(x,dearly$estimate['shape'],dearly$estimate['scale']) # Shape and scale here, 
lines(x,y, col=rgb(107,176,192,255,max=255), lwd=2)

dlate <- fitdistr(BpresLate$Temp,'weibull')
y=dweibull(x,dlate$estimate['shape'],dlate$estimate['scale']) # Shape and scale here, 
lines(x,y, col=rgb(242,26,0,255,max=255), lwd=2)

text(24,.3,'B',cex=1.5)

dev.off()

#### compare means based on normal distribution
data$period <- NA
data$period[data$Year < 2007 & data$Temp > 24] <- 'early'
data$period[data$Year >= 2007 & data$Temp > 24] <- 'late'

summary(lm(data$Temp~as.factor(data$period)))

hist(data$Temp[data$period=='early'],col=rgb(0,0,1,1, alpha=0.3), xlab="Degrees Celsius", ylab= "Probablity density", main="", freq=F, ylim= c(0, 0.3), cex.lab=1.5, cex.axis=1.5)
hist(data$Temp[data$period=='late'], freq=F, col=rgb(1,0,0,1, alpha=0.3), main="", add=T)
lines(x,dnorm(x,mean(data$Temp[data$period=='early'],na.rm=T),sd(data$Temp[data$period=='early'],na.rm=T)), col="blue", lwd=2)
lines(x,dnorm(x,mean(data$Temp[data$period=='late'],na.rm=T),sd(data$Temp[data$period=='late'],na.rm=T)), col="red", lwd=2)

#### compare difference in weibull distributions
poolw <- fitdistr(data$Temp[data$Temp > 24],'weibull')
logLik_sum <- dearly$loglik + dlate$loglik
logLik_pooled <- poolw$loglik
pchisq(2*(logLik_sum-logLik_pooled),df=2,lower.tail=FALSE) 


# SOM ---------------------------------------------------------------------


