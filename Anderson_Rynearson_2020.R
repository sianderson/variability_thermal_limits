# Ecological Simulation from Anderson and Rynearson, 2020
# Stephanie I. Anderson updated 02/02/21

# Load Libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(gridExtra)
library(lubridate)

##################################################################
# Creating a Temperature Function to Describe Narragansett Bay, RI

## Load Data
data <- read.csv("data/NB_SST.csv", header=TRUE)
data <-melt(data, id=c("Date", "Temp"))

## Reformat date to be readable
# Strips week and year from date format to create separate columns
Week <- format(as.POSIXct(strptime(data$Date,"%m/%d/%Y",tz="")) ,format = "%W")
Year <- format(as.POSIXct(strptime(data$Date,"%m/%d/%Y",tz="")) ,format = "%Y")
data$Week <- Week
data$Year <- Year
data$Year <- as.numeric(as.character(data$Year))
week.mean <- aggregate(Temp~Week, data, mean)
week.mean$Week <- as.numeric(as.character(week.mean$Week))
week.mean$Temp <- as.numeric(as.character(week.mean$Temp))

# Assess temperature from 2008 to 2015
data_sub <-subset(data, Year > 8 & Year < 15 )

# find an average weekly temperature
mean <- aggregate(Temp~Week, data_sub, mean)
mean$Week <- as.numeric(as.character(mean$Week))
mean$Temp <- as.numeric(as.character(mean$Temp))

# Fit a sinusoidal function to temperature data
SSTk <- lm(Temp ~ (sin(2*pi*Week/52)+cos(2*pi*Week/52)),data=mean)

# plot of results
# points indicate all temperature recordings for that week over the time period
# line indicates average temperature
plot(data_sub$Week, data_sub$Temp, col="gray50", xlab="", ylab="")
lines(mean$Week, SSTk$fitted.values, col=1, lwd=2)
title(ylab=expression(bold("Temperature (ºC)")), xlab=expression(bold("Week")), line=2)

# create a time vector based on a weekly scale
x=seq(from=1, to=365, by=0.1)

# Save fitted function (from SSTk) to new temp function
temp<-function(x){
  r<-(12.07626 + -6.21474  * sin(x/58.09155) + -8.62064 * cos(x/58.09155)) 
  r
}

##################################################################
# Fit curves to growth data from each species
data<-read.csv("data/Growth.csv")
Smar <- subset(data, Species=="S. marinoi")
Spse <- subset(data, Species=="S. pseudocostatum")
Sgre <- subset(data, Species=="S. grethae")
Smen <- subset(data, Species=="S. menzelii")
Sdoh <- subset(data, Species=="S. dohrnii")

# curve function comes from Norberg et al. 2004
nbcurve<-function(temp,opt,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-opt)/(w/2))^2)
  res
}

# non-linear least squares estimations
m1<-nls(Growth~nbcurve(Temp,o,w,a,b),start=list(w=5,o=10,a=1.5,b=0.02),data=Smar)
m1b<-nls(Growth~nbcurve(Temp,o,w,a,b),start=list(w=34.6,o=20.6,a=1.41,b=0.02),data=Spse)
m1c<-nls(Growth~nbcurve(Temp,o,w,a,b),start=list(w=10,o=15,a=0.04,b=0.11),data=Sgre)
m1d<-nls(Growth~nbcurve(Temp,o,w,a,b),start=list(w=37,o=20,a=1.5,b=0.02),data=Sdoh)
m1e<-nls(Growth~nbcurve(Temp,o,w,a,b),start=list(w=20,o=15,a=0.4,b=0.12),data=Smen)

# save coefficients
cfs<-as.list(coef(m1))
cfs2<-as.list(coef(m1b))
cfs3<-as.list(coef(m1c))
cfs4<-as.list(coef(m1d))
cfs5<-as.list(coef(m1e))

##################################################################
# Calculating growth rates throughout the year, based on temperature

t = 0:365
smen=nbcurve(temp(t),cfs5$o,cfs5$w,cfs5$a,cfs5$b)
smar=nbcurve(temp(t),cfs$o,cfs$w,cfs$a,cfs$b)
sgre=nbcurve(temp(t),cfs3$o,cfs3$w,cfs3$a,cfs3$b)
spse=nbcurve(temp(t),cfs2$o,cfs2$w,cfs2$a,cfs2$b)
sdoh=nbcurve(temp(t),cfs4$o,cfs4$w,cfs4$a,cfs4$b)

# visualization
plot(nbcurve(temp(0:365),cfs5$o,cfs5$w,cfs5$a,cfs5$b),add=T, col='darkorange1', lwd=3, ylim=c(0,2.0), type="l", xlab="", ylab="", xaxt='n')
curve(nbcurve(temp(x),cfs$o,cfs$w,cfs$a,cfs$b),0, 365,add=T,col='lightskyblue', lwd=3)
curve(nbcurve(temp(x),cfs3$o,cfs3$w,cfs3$a,cfs3$b), 0, 365,add=T,col='goldenrod1', lwd=3)
curve(nbcurve(temp(x),cfs2$o,cfs2$w,cfs2$a,cfs2$b),0, 365,add=T,col='firebrick3', lwd=3)
curve(nbcurve(temp(x),cfs4$o,cfs4$w,cfs4$a,cfs4$b),0, 365,add=T,col='royalblue3', lwd=3)
legend(-10, 2.1, c("S. marinoi", "S. dohrnii", "S. menzelii", "S. grethae","S. pseudocostatum"),
       col = c('lightskyblue','royalblue3','darkorange1','goldenrod1','firebrick3', 6), 
       pch = c(16, 17, 19, 15, 18),
       merge = FALSE, bg = "white",
       bty="n", 
       text.font = 3,
       cex=0.9, 
       y.intersp=0.8)
axis(side=1,at=c(0,91,182,274, 366),labels=c("Jan","Apr", "Jul", "Oct", "Jan"))
title(ylab=expression(bold("Specific Growth Rate (d"^"-1" *")")), xlab=expression(bold("Month")), line=2)

##################################################################
# Conduct Ecological Simulation

N0 = 1000  #initial population size
t = 0:1825 # five years
times= 1826

# create empty vectors to populate
N = vector(length = times)
Nm = vector(length = times)
Nd = vector(length = times)
Np = vector(length = times)
Ng = vector(length = times)

# Start all species with the same initial population size
N[1] = N0
Nm[1] = N0
Nd[1] = N0
Np[1] = N0
Ng[1] = N0

# include time step
dt=0.2
lambda=0.00001

# Population dependent on thermal growth rates and density dependent control
# Grazing based on Ivlev, 1945
# Refer to Anderson and Rynearson, 2020 for more information
for (t in 2:times) {
  N[t] = N[t - 1] * exp((nbcurve(temp(t),cfs$o,cfs$w,cfs$a,cfs$b)-1.8*(1-exp(-lambda*(N[t-1]))))*dt)
  Nm[t] = Nm[t - 1] * exp((nbcurve(temp(t),cfs5$o,cfs5$w,cfs5$a,cfs5$b)-1.8*(1-exp(-lambda*(Nm[t-1]))))*dt)
  Nd[t] = Nd[t - 1] * exp((nbcurve(temp(t),cfs4$o,cfs4$w,cfs4$a,cfs4$b)-1.8*(1-exp(-lambda*(Nd[t-1]))))*dt)
  Np[t] = Np[t - 1] * exp((nbcurve(temp(t),cfs2$o,cfs2$w,cfs2$a,cfs2$b)-1.8*(1-exp(-lambda*(Np[t-1]))))*dt)
  Ng[t] = Ng[t - 1] * exp((nbcurve(temp(t),cfs3$o,cfs3$w,cfs3$a,cfs3$b)-1.8*(1-exp(-lambda*(Ng[t-1]))))*dt)
}

##################################################################
# Plot results
# Abundance over time
plot(1:times, N, log="y",type="l", lwd=2, las = 1, col='lightskyblue', ylab="", xlab="", ylim=c(1e1,1e9))
points(1:times, Nm, log="y",type="l", lwd=2, las = 1, col='darkorange1')
points(1:times, Nd, log="y",type="l", lwd=2, las = 1, col='royalblue3')
points(1:times, Np, log="y",type="l", lwd=2, las = 1, col='firebrick3')
points(1:times, Ng, log="y",type="l", lwd=2, las = 1, col='goldenrod1')
par(mar= c(5, 5, 4, 2))
title(ylab=expression(bold("Concentration (Cells/L)")), xlab=expression(bold("Day")), line=3)
legend(-40, 1e9, c("S. marinoi", "S. dohrnii", "S. menzelli", "S. grethae","S. pseudocostatum"),
       col = c('lightskyblue','royalblue3','darkorange1','goldenrod1','firebrick3', 6), 
       pch = c(16, 17, 19, 15, 18),
       merge = FALSE, bg = "white",
       bty="n", 
       text.font = 3,
       cex=0.9)

# Plot as percent compositon over time
# combine vectors into single data frame
df=data.frame(N, Nm, Nd, Np, Ng, t)

# Calculate total population and percent compositon
df <- ddply(df,. (t), mutate, all=sum(N, Nm, Nd, Np, Ng))
df <- melt(df, id=c("t","all"))
df <- ddply(df,. (t, variable), mutate, percent = 100*(value/all))

# Plot results
ggplot(df, aes(x=t, y=percent/100, colour=variable))+
  xlim(1095, 1825)+ #Ignore first 3 years to allow for model stabilization
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.6)) +
  geom_line(size=1)+
  guides(color=FALSE)+
  theme_classic()+
  scale_colour_manual(values=c("Ng"="goldenrod1", 
                               "N"="lightskyblue",
                               "Nm"="darkorange1",
                               "Nd"="royalblue3",
                               "Np"="firebrick3")) + 
  labs(y="Species Composition", x="Day") 

