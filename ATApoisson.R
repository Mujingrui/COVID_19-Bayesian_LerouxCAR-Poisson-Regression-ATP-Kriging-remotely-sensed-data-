setwd("E:/ResearchWork2")
library(rgdal)
library(sp)
library(atakrig)


census_PHU_ad2 <- read.csv("F:/Research Work/COVID-2019/Canada/ON/data_preparation/census_PHU_adjusted.csv",header = TRUE,sep = ",")
covid_data2 <- read.csv("F:/Research Work/Daily Cases/covid_data_biweekly_rate.csv",sep=",",header = TRUE)
covid_census2 <- merge(covid_data2,census_PHU_ad,by.x = "HRUID20181.x",by.y = "GEO_CODE")
covid_census2_order <- covid_census2[order(covid_census2$Date),]

######
covid_census2_order$cases_new <- round(covid_census2_order$Age_Adjusted_Rate1*(covid_census2_order$Population2019/100000),0)
ATA_lst_data2 <- list(NULL)
length(ATA_lst_data2) <- 32
ATA_lst_model <- list(NULL)
length(ATA_lst_model) <- 32
ATA_lst_result <- list(NULL)
length(ATA_lst_result) <- 32
PHU_boundary <- rgdal::readOGR(dsn = "F:/Research Work/COVID-2019/Canada/ON/Ministry_of_Health_Public_Health_Unit_Boundary-shp",layer = "Ministry_of_Health_Public_Health_Unit_Boundary")

for(i in 1:32){
  ATA_lst_data2[[i]] <- merge(PHU_boundary,covid_census2_order[c(((i-1)*34+1):(34*i)),],by.x = "PHU_ID", by.y = "PHU_ID1")
  ATA_lst_data2[[i]] <- spTransform(ATA_lst_data2[[i]], CRS("+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
  ATA_lst_data2$response <- ATA_lst_data2[[i]]$cases_new/(ATA_lst_data2[[i]]$Population2019/100000)
  names(ATA_lst_data2[[i]]@data)[18] <- "X65_Years._over1"
}
ml12 <- glm(cases_new ~ offset(log(Population2019/100000))+Population_Density+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data2[[1]]@data)
ATA_lst_data2[[1]]$residuals <- ml12$residuals
ATA_lst_data2[[1]]$PHU_ID <- as.character(ATA_lst_data2[[1]]$PHU_ID)
obs.discrete <- discretizePolygon(ATA_lst_data2[[1]],cellsize=15000,id="PHU_ID",value="residuals")
pointsv <- deconvPointVgm(obs.discrete,model="Exp",ngroup=10,rd=0.75,fig=TRUE)
pred.cv <- ataKriging.cv(obs.discrete, nfold=length(ATA_lst_data2[[1]]), pointsv, showProgress = TRUE)
prediction_frame <- ATA_lst_data2[[2]]@data[,c(1,3,11:12,17:28)]
prediction_pre <- merge(prediction_frame,pred.cv,by.x = "PHU_ID",by.y = "areaId")
#pre_matrix <- as.matrix(prediction_pre[,9:19])
for(i in 1:34){
  prediction_pre$local_risk[i] <- exp(ml12$coefficients[1]*prediction_frame$Population_Density[i]+ml12$coefficients[2]*prediction_frame$X65_Years._over1[i]+ml12$coefficients[3]*prediction_frame$Apartment_Proportion[i]+ml12$coefficients[4]*prediction_frame$Low_income[i]+ml12$coefficients[5]*prediction_frame$X18_65Years[i]+ml12$coefficients[6]*prediction_frame$X65Years_over[i]+ml12$coefficients[7]*prediction_frame$Other_Proportion[i]+ml12$coefficients[8]*prediction_frame$Apartment2_Proportion[i]+ml12$coefficients[9]*prediction_frame$Average_Household_size[i]+ml12$coefficients[10]*prediction_frame$Self_employed_Proportion[i]+ml12$coefficients[11]*prediction_frame$Health_occupations_proportion[i]+prediction_pre$pred[i])
}
####
library(gstat)
plotDeconvVgm <- function(v, main=NULL, posx=NULL, posy=NULL, lwd=2, showRegVgm=FALSE) {
  plotDVgm <- function(v, main) {
    xlim <- c(0, 1.1 * max(v$experientialAreaVariogram$dist))
    xx <- seq(0, xlim[2], length=100)
    xx[1] <- 1.0e-3
    yy <- variogramLine(v$areaVariogram, covariance=FALSE, dist_vector=xx)$gamma
    yy2 <- variogramLine(v$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma
    
    ylim <- c(min(0, min(v$experientialAreaVariogram$gamma)),
              1.1 * max(c(v$experientialAreaVariogram$gamma, yy, yy2)))
    #plot(v$experientialAreaVariogram$dist, v$experientialAreaVariogram$gamma,
        # xaxs="i", yaxs="i", col="blue", xlim=xlim, ylim=ylim, ann = FALSE,
        # xlab="distance", ylab="semivariance", main=main)
    lines(xx, yy, col="blue", lwd=lwd)
    lines(xx, yy2, col="red", lwd=lwd)
    if(showRegVgm)
      lines(v$regularizedAreaVariogram$dist, v$regularizedAreaVariogram$gamma, col="black", type="l", lty=2, lwd=lwd)
    mtext(main, side = 3, line = 0.2, cex = 0.8)
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar=c(2.5, 3.5, 1, .5)) # reduce the margins around the figure
  par(mgp=c(1.5, .5, 0)) #  reduce the spacing between the figure plotting region and the axis labels
  par(oma=c(1, 1, 2.5, 0)) #  add an outer margin to the top of the graph
  
  if (hasName(v, "regularizedAreaVariogram")) {
    plotDVgm(v, main = "")
    if(is.null(posx)) posx <- "bottomright"
  } else {
    m0 <- (sqrt(1+8*length(v))-1)/2
    n <- sort(names(v)[1:m0])
    m1 <- m0*(m0-1)/2 + m0
    m <- matrix(0, nrow=m0, ncol=m0)
    m[lower.tri(m, diag = TRUE)] <- seq(m1)
    m[1,2:m0] <- m1 + 1
    layout(mat=m)
    
    for (i in 1:m0) {
      plotDVgm(v[[n[i]]], main = n[i])
      for (j in 1:m0) {
        if(j > i) {
          n2 <- crossName(n[i],n[j])
          plotDVgm(v[[n2]], main = n2)
        }
      }
    }
    plot.new()
    if(is.null(posx)) posx <- "left"
  }
  
  if(showRegVgm) {
    legend(posx, posy, bty="n",
           legend=c("Empirical area variogram", "Fitted area variogram",
                    "Deconvoluted point variogram", "Regularized area variogram"),
           pch=c(1,NA,NA,NA), lty=c(NA,1,1,2), lwd=c(NA,lwd,lwd,lwd),
           col=c("blue","blue","red","black"))
  } else {
    legend(posx, posy, bty="n",
           legend=c("Empirical area variogram", "Fitted area variogram",
                    "Deconvoluted point variogram"),
           pch=c(1,NA,NA), lty=c(NA,1,1), lwd=c(NA,lwd,lwd),
           col=c("blue","blue","red","black"))
  }
  
  if(is.null(main)) main <- "Deconvoluted variogram"
  mtext(main, outer = TRUE,  side = 3, cex = 1.2, line = 0)
  mtext("distance", outer = TRUE,  side = 1, cex = 1.05, line = 0)
  mtext("semivariance", outer = TRUE,  side = 2, cex = 1.05, line = -1)
}
plotDeconvVgm(pointsv)

plotDeconvVgm2 <- function(v, main=NULL, posx=NULL, posy=NULL, lwd=2, showRegVgm=FALSE) {
  plotDVgm <- function(v, main) {
    xlim <- c(0, 1.1 * max(v$experientialAreaVariogram$dist))
    xx <- seq(0, xlim[2], length=100)
    xx[1] <- 1.0e-3
    yy <- variogramLine(v$areaVariogram, covariance=FALSE, dist_vector=xx)$gamma
    yy2 <- variogramLine(v$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma
    
    ylim <- c(min(0, min(v$experientialAreaVariogram$gamma)),
              1.1 * max(c(v$experientialAreaVariogram$gamma, yy, yy2)))
    #plot(v$experientialAreaVariogram$dist, v$experientialAreaVariogram$gamma,
    # xaxs="i", yaxs="i", col="blue", xlim=xlim, ylim=ylim, ann = FALSE,
    # xlab="distance", ylab="semivariance", main=main)
    lines(xx, yy, col="blue", lwd=lwd)
    lines(xx, yy2, col="red", lwd=lwd)
    if(showRegVgm)
      lines(v$regularizedAreaVariogram$dist, v$regularizedAreaVariogram$gamma, col="black", type="l", lty=2, lwd=lwd)
    mtext(main, side = 3, line = 0.2, cex = 0.8)
  }}
plotDVgm2(pointsv)
#demo <- data.frame(xx, yy)
yy2 <- variogramLine(pointsv$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma
demo <- data.frame(xx, yy,yy2)
demo_long <- melt(demo, id="xx")  # convert to long format

library(ggplot2)
library(reshape)
ggplot(data=demo_long,
       aes(x=xx, y=value, colour=variable)) +
  geom_line()
####
ml22 <- glm(cases_new ~ offset(log(Population2019/100000))+Population_Density+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data2[[2]]@data)
ATA_lst_data2[[2]]$residuals <- ml22$residuals
ATA_lst_data2[[2]]$PHU_ID <- as.character(ATA_lst_data2[[2]]$PHU_ID)
obs.discrete22 <- discretizePolygon(ATA_lst_data2[[2]],cellsize=15000,id="PHU_ID",value="residuals")
pointsv22 <- deconvPointVgm(obs.discrete22,model="Exp",ngroup=35,rd=0.75,fig=TRUE)
predictionLocations2 <- ATA_lst_data[[2]]
pred.discrete2 <- discretizePolygon(predictionLocations2,cellsize=15000,id="PHU_ID")
pred2 <- ataKriging(obs.discrete2,pred.discrete2,pointsv2$pointVariogram)
pred.cv <- ataKriging.cv(obs.discrete2, nfold=length(ATA_lst_data2[[2]]), pointsv22, showProgress = TRUE)
View(ATA_lst_data[[2]])
ml22 <- glm(cases_new ~ offset(log(Population2019/100000))+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[2]]@data)
ATA_lst_data[[2]]$residuals2 <- ml22$residuals
#ATA_lst_data[[2]]$PHU_ID <- as.character(ATA_lst_data[[2]]$PHU_ID)
obs.discrete22 <- discretizePolygon(ATA_lst_data[[2]],cellsize=15000,id="PHU_ID",value="residuals2")
pointsv22 <- deconvPointVgm(obs.discrete22,model="Exp",ngroup=16,rd=0.75,fig=TRUE)
predictionLocations2 <- ATA_lst_data[[2]]
# pred.discrete2 <- discretizePolygon(predictionLocations2,cellsize=15000,id="PHU_ID")
# pred22 <- ataKriging(obs.discrete22,pred.discrete2,pointsv22$pointVariogram)
pred.cv2 <- ataKriging.cv(obs.discrete22, nfold=length(ATA_lst_data[[2]]), pointsv22, showProgress = TRUE)
prediction_frame <- ATA_lst_data[[2]]@data[,c(1,3,11:15,17,18:29,32)]
prediction_pre <- merge(prediction_frame,pred.cv2,by.x = "PHU_ID",by.y = "areaId")
pre_matrix <- as.matrix(prediction_pre[,9:19])
for(i in 1:34){
  prediction_pre$local_risk[i] <- exp(ml22$coefficients[1]*prediction_frame$X65_Years._over1[i]+ml22$coefficients[2]*prediction_frame$Apartment_Proportion[i]+ml22$coefficients[3]*prediction_frame$Low_income[i]+ml22$coefficients[4]*prediction_frame$X18_65Years[i]+ml22$coefficients[5]*prediction_frame$X65Years_over[i]+ml22$coefficients[6]*prediction_frame$Other_Proportion[i]+ml22$coefficients[7]*prediction_frame$Apartment2_Proportion[i]+ml22$coefficients[8]*prediction_frame$Average_Household_size[i]+ml22$coefficients[9]*prediction_frame$Self_employed_Proportion[i]+ml22$coefficients[10]*prediction_frame$Health_occupations_proportion[i]+prediction_pre$pred[i])
}
 
###
ml3 <- glm(cases_new ~ offset(log(Population2019/100000))+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[3]]@data)
ATA_lst_data[[3]]$residuals <- ml3$residuals
ATA_lst_data[[3]]$PHU_ID <- as.character(ATA_lst_data[[3]]$PHU_ID)
obs.discrete3 <- discretizePolygon(ATA_lst_data[[3]],cellsize=15000,id="PHU_ID",value="residuals")
pointsv3 <- deconvPointVgm(obs.discrete3,model="Exp",ngroup=10,rd=0.75,fig=TRUE)
library(reshape)
library(ggplot2)
plotDVgm2 <- function(v){
  xlim <- c(0, 1.1 * max(v$experientialAreaVariogram$dist))
  xx <- seq(0, xlim[2], length=100)
  xx[1] <- 1.0e-3
  yy <- variogramLine(v$areaVariogram, covariance=FALSE, dist_vector=xx)$gamma
  yy2 <- variogramLine(v$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma
  
  ylim <- c(min(0, min(v$experientialAreaVariogram$gamma)),
            1.1 * max(c(v$experientialAreaVariogram$gamma, yy, yy2)))
  
  demo <- data.frame(xx, yy,yy2)
  colnames(demo) <- c("Distance","Fitted Area","Devoncoluted Point")
  demo_long <- melt(demo, id="Distance")  # convert to long format
  colnames(demo_long)[2] <- c("Semivariogram")
  #library(ggplot2)
  ggplot(data=demo_long,
         aes(x=Distance, y=value, colour=Semivariogram)) +
         geom_line()+
         theme_bw()
         #+theme(legend.position = "bottom")
         #scale_scale_continuous(name = "Semivariograms",
                            # breaks = c("yy","yy2"),
                            # labels = c("Fitted Area","Deconvoluted Point")
    
}
plot3 <- plotDVgm2(pointsv3)


predictionLocations3 <- ATA_lst_data[[3]]
pred.discrete3 <- discretizePolygon(predictionLocations3,cellsize=15000,id="PHU_ID")
pred3 <- ataKriging(obs.discrete3,pred.discrete3,pointsv3$pointVariogram)
pred.cv3 <- ataKriging.cv(obs.discrete3, nfold=length(ATA_lst_data[[3]]), pointsv3, showProgress = TRUE)





pred33 <- atpKriging(obs.discrete3,centsid_census[,2:3],pointsv3$pointVariogram)
pred33 <- merge(pred33,rectangles_census@data,by.x = "areaId", by.y = "ID")
pred33$cases_estimated <- round(exp(m13$coefficients[1]+m13$coefficients[2]*pred33$Other_Proportion+m13$coefficients[3]*pred33$Apartment2_Proportion+
                                m13$coefficients[4]*pred33$Avr_H_S+m13$coefficients[5]*pred33$Self_employed_Proportion+m13$coefficients[6]*pred33$Health_occupations_proportion+pred33$pred),0)
rectangles_pred33 <- merge(rectangles_census,pred33,by.x = "ID",by.y = "areaId")
l_dww <- tm_shape(rectangles_pred33) + tm_fill("cases_estimated", breaks=c(0,10,30,50,100,300,500,1000,2000,4000,6000,7000), Palette = "Reds", title = "Estimated Infected Risk \n (cases per 100,000)") +
  tm_borders(alpha=.1) +
  tm_facets(by = "Date")+
  tm_layout(legend.text.size = 0.4,legend.title.size = 1, frame = FALSE, legend.outside = TRUE,
            legend.show = TRUE, legend.hist.height = 2)


ml4 <- glm(cases_new ~ offset(log(Population2019/100000))+Population_Density+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[4]]@data)
ml5 <- glm(cases_new ~ offset(log(Population2019/100000))+Population_Density+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[5]]@data)
for(i in 1:32){
  ATA_lst_model[[i]] <- glm(cases_new ~ offset(log(Population2019/100000))+Population_Density+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[i]]@data)
}
ATA_lst_model2 <- list(NULL)
length(ATA_lst_model2) <- 32
for(i in 1:32){
  ATA_lst_model2[[i]] <- glm(cases_new ~ offset(log(Population2019/100000))+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[i]]@data)
}

ATA_lst_model3 <- list(NULL)
length(ATA_lst_model3) <- 32
obs.discrete <- list(NULL)
length(obs.discrete) <- 32

for(i in 1:32){
  ATA_lst_model3[[i]] <- glm(cases_new ~ offset(log(Population2019/100000))+Population_Density+X65_Years._over1+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[i]]@data)
}
for(i in 1:32){
  ATA_lst_data[[i]]$residuals <- ATA_lst_model3[[i]]$residuals
  ATA_lst_data[[i]]$PHU_ID <- as.character(ATA_lst_data[[i]]$PHU_ID)
  obs.discrete[[i]] <- discretizePolygon(ATA_lst_data[[i]],cellsize=15000,id="PHU_ID",value="residuals")
}
pointsv1 <- deconvPointVgm(obs.discrete[[1]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv2 <- deconvPointVgm(obs.discrete[[2]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv3 <- deconvPointVgm(obs.discrete[[3]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv4 <- deconvPointVgm(obs.discrete[[4]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv5 <- deconvPointVgm(obs.discrete[[5]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv6 <- deconvPointVgm(obs.discrete[[6]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv7 <- deconvPointVgm(obs.discrete[[7]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv8 <- deconvPointVgm(obs.discrete[[8]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv9 <- deconvPointVgm(obs.discrete[[9]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv10 <- deconvPointVgm(obs.discrete[[10]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv11 <- deconvPointVgm(obs.discrete[[11]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv12 <- deconvPointVgm(obs.discrete[[12]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv13 <- deconvPointVgm(obs.discrete[[13]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv14 <- deconvPointVgm(obs.discrete[[14]],model="Sph",ngroup=14,rd=0.75,fig=TRUE)
pointsv15 <- deconvPointVgm(obs.discrete[[15]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv16 <- deconvPointVgm(obs.discrete[[16]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv17 <- deconvPointVgm(obs.discrete[[17]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv18 <- deconvPointVgm(obs.discrete[[18]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv19 <- deconvPointVgm(obs.discrete[[19]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv20 <- deconvPointVgm(obs.discrete[[20]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv21 <- deconvPointVgm(obs.discrete[[21]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv22 <- deconvPointVgm(obs.discrete[[22]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv23 <- deconvPointVgm(obs.discrete[[23]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv24 <- deconvPointVgm(obs.discrete[[24]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv25 <- deconvPointVgm(obs.discrete[[25]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv26 <- deconvPointVgm(obs.discrete[[26]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv27 <- deconvPointVgm(obs.discrete[[27]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv28 <- deconvPointVgm(obs.discrete[[28]],model="Exp",ngroup=14,rd=0.75,fig=TRUE)
pointsv29 <- deconvPointVgm(obs.discrete[[29]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv30 <- deconvPointVgm(obs.discrete[[30]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv31 <- deconvPointVgm(obs.discrete[[31]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)
pointsv32 <- deconvPointVgm(obs.discrete[[32]],model="Exp",ngroup=12,rd=0.75,fig=TRUE)

ATA_lst_point <- list(NULL)
length(ATA_lst_point) <- 32
ATA_lst_point[[1]] <- pointsv1
ATA_lst_point[[2]] <- pointsv2
ATA_lst_point[[3]] <- pointsv3
ATA_lst_point[[4]] <- pointsv4
ATA_lst_point[[5]] <- pointsv5
ATA_lst_point[[6]] <- pointsv6
ATA_lst_point[[7]] <- pointsv7
ATA_lst_point[[8]] <- pointsv8
ATA_lst_point[[9]] <- pointsv9
ATA_lst_point[[10]] <- pointsv10
ATA_lst_point[[11]] <- pointsv11
ATA_lst_point[[12]] <- pointsv12
ATA_lst_point[[13]] <- pointsv13
ATA_lst_point[[14]] <- pointsv14
ATA_lst_point[[15]] <- pointsv15
ATA_lst_point[[16]] <- pointsv16
ATA_lst_point[[17]] <- pointsv17
ATA_lst_point[[18]] <- pointsv18
ATA_lst_point[[19]] <- pointsv19
ATA_lst_point[[20]] <- pointsv20
ATA_lst_point[[21]] <- pointsv21
ATA_lst_point[[22]] <- pointsv22
ATA_lst_point[[23]] <- pointsv23
ATA_lst_point[[24]] <- pointsv24
ATA_lst_point[[25]] <- pointsv25
ATA_lst_point[[26]] <- pointsv26
ATA_lst_point[[27]] <- pointsv27
ATA_lst_point[[28]] <- pointsv28
ATA_lst_point[[29]] <- pointsv29
ATA_lst_point[[30]] <- pointsv30
ATA_lst_point[[31]] <- pointsv31
ATA_lst_point[[32]] <- pointsv32
saveRDS(ATA_lst_point,"E:/ResearchWork2/ATAPoisson_point/ATPPoisson_point.RDS")
saveRDS(ATA_lst_model3,"E:/ResearchWork2/ATAPoisson_point/ATPPoisson_glm.RDS")
saveRDS(ATA_lst_data,"E:/ResearchWork2/ATAPoisson_point/ATPPoisson_data.RDS")

ATA_lst_model3 <- readRDS("E:/ResearchWork2/ATAPoisson_point/ATPPoisson_glm.RDS")
ATA_lst_point <- readRDS("E:/ResearchWork2/ATAPoisson_point/ATPPoisson_point.RDS")


pointsv3 <- deconvPointVgm(obs.discrete3,model="Exp",ngroup=10,rd=0.75,fig=TRUE)
predictionLocations3 <- ATA_lst_data[[3]]
pred.discrete3 <- discretizePolygon(predictionLocations3,cellsize=15000,id="PHU_ID")
pred3 <- ataKriging(obs.discrete3,pred.discrete3,pointsv3$pointVariogram)

pred33 <- atpKriging(obs.discrete3,centsid_census[,2:3],pointsv3$pointVariogram)

##############
ATA_lst_data <- list(NULL)
length(ATA_lst_data) <- 31
ATA_lst_model <- list(NULL)
length(ATA_lst_model) <- 31
ATA_lst_result <- list(NULL)
length(ATA_lst_result) <- 31
covid_census_order <- covid_census[order(covid_census$Date),]
covid_census_order$cases_new <- round(covid_census_order$Age_Adjusted_Rate*(covid_census_order$POP2019/100000),0)

for(i in 1:31){
  ATA_lst_data[[i]] <- merge(PHU_boundary,covid_census_order[c(((i-1)*34+1):(34*i)),],by.x = "PHU_ID", by.y = "PHU_ID")
  ATA_lst_data[[i]] <- spTransform(ATA_lst_data[[i]], CRS("+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
  #ATA_lst_data[[i]]$response <- ATA_lst_data[[i]]$cases_new/(ATA_lst_data[[i]]$POP2019/100000)
  #names(ATA_lst_data[[i]]@data)[18] <- "X65_Years._over1"
}
ml1 <- glm(cases_new ~ offset(log(POP2019/100000))+X65_Years._over+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[1]]@data)
ATA_lst_data[[1]]$residuals <- ml1$residuals
ATA_lst_data[[1]]$PHU_ID <- as.character(ATA_lst_data[[1]]$PHU_ID)
obs.discrete1 <- discretizePolygon(ATA_lst_data[[1]],cellsize=15000,id="PHU_ID",value="residuals")
pointsv1 <- deconvPointVgm(obs.discrete1,model="Exp",ngroup=25,fixed.range = 450000,rd=0.75,fig=TRUE)

plotDVgm2 <- function(v){
  xlim <- c(0, 1.1 * max(v$experientialAreaVariogram$dist))
  xx <- seq(0, xlim[2], length=100)
  xx[1] <- 1.0e-3
  yy <- variogramLine(v$areaVariogram, covariance=FALSE, dist_vector=xx)$gamma
  yy2 <- variogramLine(v$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma
  
  ylim <- c(min(0, min(v$experientialAreaVariogram$gamma)),
            1.1 * max(c(v$experientialAreaVariogram$gamma, yy, yy2)))
  
  demo <- data.frame(xx, yy,yy2)
  colnames(demo) <- c("Distance","Fitted Area","Devoncoluted Point")
  demo_long <- melt(demo, id="Distance")  # convert to long format
  colnames(demo_long)[2] <- c("Semivariogram")
  #library(ggplot2)
  ggplot(data=demo_long,
         aes(x=Distance, y=value, colour=Semivariogram)) +
    geom_line()+
    theme_bw()
  #+theme(legend.position = "bottom")
  #scale_scale_continuous(name = "Semivariograms",
  # breaks = c("yy","yy2"),
  # labels = c("Fitted Area","Deconvoluted Point")
  
}
plot1 <- plotDVgm2(pointsv1)
pred.cv1 <- ataKriging.cv(obs.discrete1, nfold=length(ATA_lst_data[[1]]), pointsv1, showProgress = TRUE)
prediction_frame1 <- ATA_lst_data[[1]]@data[,c(1,12,13,31:41)]
prediction_pre1 <- merge(prediction_frame1,pred.cv1,by.x = "PHU_ID",by.y = "areaId")
local_risk <- rep(0,34)
for(i in 1:34){
  prediction_pre1$local_risk[i] <- exp(ml1$coefficients[1]*prediction_frame1$X65_Years._over[i]+ml1$coefficients[2]*prediction_frame1$Apartment_Proportion[i]+ml1$coefficients[3]*prediction_frame1$Low_income[i]+ml1$coefficients[4]*prediction_frame1$X18_65Years[i]+ml1$coefficients[5]*prediction_frame1$X65Years_over[i]+ml1$coefficients[6]*prediction_frame1$Other_Proportion[i]+ml1$coefficients[7]*prediction_frame1$Apartment2_Proportion[i]+ml1$coefficients[8]*prediction_frame1$Average_Household_size[i]+ml1$coefficients[9]*prediction_frame1$Self_employed_Proportion[i]+ml1$coefficients[10]*prediction_frame1$Health_occupations_proportion[i]+prediction_pre1$pred[i])
}



ml2 <- glm(cases_new ~ offset(log(POP2019/100000))+X65_Years._over+Apartment_Proportion+Low_income+X18_65Years+X65Years_over+Other_Proportion+Apartment2_Proportion+Average_Household_size+Self_employed_Proportion+Health_occupations_proportion-1, family=poisson(link = "log"), data=ATA_lst_data[[2]]@data)
ATA_lst_data[[2]]$residuals <- ml2$residuals
ATA_lst_data[[2]]$PHU_ID <- as.character(ATA_lst_data[[2]]$PHU_ID)
obs.discrete2 <- discretizePolygon(ATA_lst_data[[2]],cellsize=15000,id="PHU_ID",value="residuals")
pointsv2 <- deconvPointVgm(obs.discrete1,model="Exp",ngroup=25,rd=0.75,fixed.range = 150000, fig=TRUE)
plot2 <- plotDVgm2(pointsv2)
pred.cv2 <- ataKriging.cv(obs.discrete2, nfold=length(ATA_lst_data[[2]]), pointsv2, showProgress = TRUE)
prediction_frame2 <- ATA_lst_data[[2]]@data[,c(1,12,13,31:41)]
prediction_pre2 <- merge(prediction_frame2,pred.cv2,by.x = "PHU_ID",by.y = "areaId")
local_risk <- rep(0,34)
for(i in 1:34){
  prediction_pre2$local_risk[i] <- exp(ml1$coefficients[1]*prediction_frame2$X65_Years._over[i]+ml1$coefficients[2]*prediction_frame2$Apartment_Proportion[i]+ml1$coefficients[3]*prediction_frame2$Low_income[i]+ml1$coefficients[4]*prediction_frame2$X18_65Years[i]+ml1$coefficients[5]*prediction_frame1$X65Years_over[i]+ml1$coefficients[6]*prediction_frame2$Other_Proportion[i]+ml1$coefficients[7]*prediction_frame2$Apartment2_Proportion[i]+ml1$coefficients[8]*prediction_frame2$Average_Household_size[i]+ml1$coefficients[9]*prediction_frame2$Self_employed_Proportion[i]+ml1$coefficients[10]*prediction_frame2$Health_occupations_proportion[i]+prediction_pre2$pred[i])
}