setwd("E:/ResearchWork2")
covid_biweekly_census <- read.csv("E:/ResearchWork2/data_preparation/covid_biweekly_census.csv",sep=",",header=TRUE)
covid_biweekly_census <- covid_biweekly_census[order(covid_biweekly_census$Date),]
covid_biweekly_census <- covid_biweekly_census[order(covid_biweekly_census$Geography1.x),]

library(rgdal)
library(sp)
library(sf)
library(tigris)
library(INLA)
library(spdep)

PHU_boundary <- rgdal::readOGR(dsn = "F:/Research Work/COVID-2019/Canada/ON/Ministry_of_Health_Public_Health_Unit_Boundary-shp",layer = "Ministry_of_Health_Public_Health_Unit_Boundary")
covid_data2_1 <- subset(covid_biweekly_census,PHU_ID1 == 2226)
ID1 <- rep(1,nrow(covid_data2_1))
for(i in 2:34){
  ID1 <- c(ID1,rep(i,32))
}

date_num_str_int <- rep(c(1:32),34)
date_num_unstr_int <- rep(c(1:32),34)

date_num_str <- rep(c(1:32),34)
date_num_unstr <- rep(c(1:32),34)

region_date_num <- c(1:1088)

covid_biweekly_census$ID1 <- ID1
covid_biweekly_census$date_num_str_int <- date_num_str_int
covid_biweekly_census$date_num_unstr_int <- date_num_unstr_int
covid_biweekly_census$date_num_str <- date_num_str
covid_biweekly_census$date_num_unstr <- date_num_unstr
covid_biweekly_census$region_date_num <- region_date_num
covid_biweekly_census$region <- ID1
covid_biweekly_census$cases_new <- round(covid_biweekly_census$Age_Adjusted_Rate1*(covid_biweekly_census$Population2019/100000),0)
covid_censuss <- subset(covid_biweekly_census,Date == "2020-03-31")
PHU_covid_census <- merge(PHU_boundary,covid_censuss,by.x = "PHU_ID",by.y = "PHU_ID1")

nb_ON1 <- poly2nb(PHU_covid_census)
head(nb_ON1)
nb2INLA("map.adj1", nb_ON1)
g_ON1 <- inla.read.graph(filename = "map.adj1")

######
S=34
T=32
temp <- poly2nb(PHU_covid_census)
nb2INLA("map.adj1", temp)
#kenya.adj <- paste(getwd(),"/kenya.graph",sep="")
H <- inla.read.graph(filename="map.adj1")
####
D1 <- diff(diag(T),differences=1)
Q.gammaRW1 <- t(D1)%*%D1
D2 <- diff(diag(T),differences=2)
Q.gammaRW2 <- t(D2)%*%D2
####
Q.xi <- matrix(0, H$n, H$n)
for (i in 1:H$n){
  Q.xi[i,i]=H$nnbs[[i]]
  Q.xi[i,H$nbs[[i]]]=-1
}
Q.Leroux <- diag(S)-Q.xi
R1 <- kronecker(Q.gammaRW1,diag(S))
############
#r.def <- S
#########
source("F:/Research Work/Readings/Bayesian/Ugarte_et_al_(2)/Ugarte_et_al_(2014)/kronecker_nullspace.R")

Type <- "TypeIV"	# or Type <- "TypeII"/"TypeIII"/"TypeIV"

if(Type=="TypeII") {
  null.space <- kronecker.null.space(diag(S),Q.gammaRW1)
  R <- null.space[[1]]
  R <- as.matrix(R)
  A_delta <- as.matrix(null.space[[2]])
}
if(Type=="TypeIII") {
  null.space <- kronecker.null.space(Q.xi,diag(T))
  R <- null.space[[1]]
  R <- as.matrix(R)
  A_delta <- as.matrix(null.space[[2]])
}
if(Type=="TypeIV") {
  null.space <- kronecker.null.space(Q.xi,Q.gammaRW1)
  R <- null.space[[1]]
  R <- as.matrix(R)
  A_delta <- as.matrix(null.space[[2]])
}

rank.def <- dim(A_delta)[1]

#########
ddd <- Q.xi%x%Q.gammaRW1
ddd_eigen <- eigen(ddd)
ddd_R <- ddd_eigen$vectors[1024:1088,]

#######
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*beta+(b-1)*(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"

#e <- kronecker(matrix(1,1,T),diag(S))
e <- eigen(R1)
A.constr <- e$vectors[1055:1088,]



formula4 <- cases_new~Apartment_Proportion+Other_Proportion+Apartment2_Proportion+Average_Household_size+Low_income+X0_17Years+X18_65Years+X65Years_over+Self_employed_Proportion+Population_Density+Health_occupations_proportion+X65_Years._over+
               f(region, model="generic1", Cmatrix= Q.Leroux, constr=TRUE,
               hyper=list(prec=list(prior=sdunif),
               beta=list(prior=lunif)))+
               f(date_num_str, model="rw1")+
               f(date_num_unstr,model = "iid",constr=TRUE,
               hyper=list(prec=list(prior=sdunif)))+
               f(region_date_num,model="generic0", Cmatrix=ddd, rankdef=dim(ddd_R)[1], constr=TRUE,
               hyper=list(prec=list(prior=sdunif)),
               extraconstr=list(A=ddd_R, e=rep(0,dim(ddd_R)[1])))

model10a4<-inla(formula4, family="poisson", data=covid_biweekly_census, E = Population2019/100000,
               
               control.predictor=list(compute=TRUE),
               
               control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
               
               control.inla=list(strategy="simplified.laplace"),
               #verbose=TRUE
)
####

covid_biweekly_census_pre2 <- subset(covid_biweekly_census,Date>="2021-05-31"&Date<="2021-06-30")
covid_biweekly_census_pre2$cases_new <- NA
covid_biweekly_census_pre1 <- subset(covid_biweekly_census,Date < "2021-05-31")
covid_biweekly_census_pre <- rbind(covid_biweekly_census_pre1,covid_biweekly_census_pre2)
covid_biweekly_census_pre <- covid_biweekly_census_pre[order(covid_biweekly_census_pre$Date),]
covid_biweekly_census_pre <- covid_biweekly_census_pre[order(covid_biweekly_census_pre$Geography1.x),]
######
link = rep(NA, 1088)
link[which(is.na(covid_biweekly_census_pre$cases_new))] = 1
#####

model10a4_pre<-inla(formula4, family="poisson", data=covid_biweekly_census_pre, E = Population2019/100000,
                
                control.predictor=list(link = link),
                
                #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                
                control.inla=list(strategy="simplified.laplace"),
                #verbose=TRUE
)
########
covid_biweekly_census_pre22 <- subset(covid_biweekly_census,Date>="2021-06-14"&Date<="2021-06-30")
covid_biweekly_census_pre22$cases_new <- NA
covid_biweekly_census_pre11 <- subset(covid_biweekly_census,Date < "2021-06-14")
covid_biweekly_census_pre3 <- rbind(covid_biweekly_census_pre11,covid_biweekly_census_pre22)
covid_biweekly_census_pre3 <- covid_biweekly_census_pre3[order(covid_biweekly_census_pre3$Date),]
covid_biweekly_census_pre3 <- covid_biweekly_census_pre3[order(covid_biweekly_census_pre3$Geography1.x),]
######
link1 = rep(NA, 1088)
link1[which(is.na(covid_biweekly_census_pre3$cases_new))] = 1
###
model10a4_pre2<-inla(formula4, family="poisson", data=covid_biweekly_census_pre3, E = Population2019/100000,
                    
                    control.predictor=list(link = link1),
                    
                    #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                    
                    control.inla=list(strategy="simplified.laplace"),
                    #verbose=TRUE
)
#####
covid_biweekly_census_pre222 <- subset(covid_biweekly_census,Date=="2021-06-30")
covid_biweekly_census_pre222$cases_new <- NA
covid_biweekly_census_pre111 <- subset(covid_biweekly_census,Date < "2021-06-30")
covid_biweekly_census_pre4 <- rbind(covid_biweekly_census_pre111,covid_biweekly_census_pre222)
covid_biweekly_census_pre4 <- covid_biweekly_census_pre4[order(covid_biweekly_census_pre4$Date),]
covid_biweekly_census_pre4 <- covid_biweekly_census_pre4[order(covid_biweekly_census_pre4$Geography1.x),]
######
link2 = rep(NA, 1088)
link2[which(is.na(covid_biweekly_census_pre4$cases_new))] = 1
####
model10a4_pre3<-inla(formula4, family="poisson", data=covid_biweekly_census_pre4, E = Population2019/100000,
                     
                     control.predictor=list(link = link2),
                     
                     #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                     
                     control.inla=list(strategy="simplified.laplace"),
                     #verbose=TRUE
)
############

ddd_2 <- diag(S)%x%Q.gammaRW1
ddd2_eigen <- eigen(ddd_2)
ddd2_R <- ddd2_eigen$vectors[1055:1088,]

formula2 <- cases_new~Apartment_Proportion+Other_Proportion+Apartment2_Proportion+Average_Household_size+Low_income+X0_17Years+X18_65Years+X65Years_over+Self_employed_Proportion+Population_Density+Health_occupations_proportion+X65_Years._over+
  f(region, model="generic1", Cmatrix= Q.Leroux, constr=TRUE,
    hyper=list(prec=list(prior=sdunif),
               beta=list(prior=lunif)))+
  f(date_num_str, model="rw1")+
  f(date_num_unstr,model = "iid",constr=TRUE,
    hyper=list(prec=list(prior=sdunif)))+
  f(region_date_num,model="generic0", Cmatrix=ddd_2, rankdef=dim(ddd2_R)[1], constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=ddd2_R, e=rep(0,dim(ddd2_R)[1])))

model10a2<-inla(formula2, family="poisson", data=covid_biweekly_census, E = Population2019/100000,
                
                control.predictor=list(compute=TRUE),
                
                control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                
                control.inla=list(strategy="simplified.laplace"),
                #verbose=TRUE
)

model10a2_pre <-inla(formula2, family="poisson", data=covid_biweekly_census_pre, E = Population2019/100000,
                
                control.predictor=list(link = link),
                
                #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                
                control.inla=list(strategy="simplified.laplace"),
                #verbose=TRUE
)

model10a2_pre1 <-inla(formula2, family="poisson", data=covid_biweekly_census_pre3, E = Population2019/100000,
                     
                     control.predictor=list(link = link1),
                     
                     #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                     
                     control.inla=list(strategy="simplified.laplace"),
                     #verbose=TRUE
)

model10a2_pre2<-inla(formula2, family="poisson", data=covid_biweekly_census_pre4, E = Population2019/100000,
                     
                     control.predictor=list(link = link2),
                     
                     #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                     
                     control.inla=list(strategy="simplified.laplace"),
                     #verbose=TRUE
)
#######
ddd_3 <- Q.xi%x%diag(T)
ddd3_eigen <- eigen(ddd_3)
ddd3_R <- ddd3_eigen$vectors[1057:1088,]
formula3 <- cases_new~Apartment_Proportion+Other_Proportion+Apartment2_Proportion+Average_Household_size+Low_income+X0_17Years+X18_65Years+X65Years_over+Self_employed_Proportion+Population_Density+Health_occupations_proportion+X65_Years._over+
  f(region, model="generic1", Cmatrix= Q.Leroux, constr=TRUE,
    hyper=list(prec=list(prior=sdunif),
               beta=list(prior=lunif)))+
  f(date_num_str, model="rw1")+
  f(date_num_unstr,model = "iid",constr=TRUE,
    hyper=list(prec=list(prior=sdunif)))+
  f(region_date_num,model="generic0", Cmatrix=ddd_3, rankdef=dim(ddd3_R)[1], constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=ddd3_R, e=rep(0,dim(ddd3_R)[1])))

model10a3<-inla(formula3, family="poisson", data=covid_biweekly_census, E = Population2019/100000,
                
                control.predictor=list(compute=TRUE),
                
                control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                
                control.inla=list(strategy="simplified.laplace"),
                #verbose=TRUE
)

model10a3_pre <-inla(formula3, family="poisson", data=covid_biweekly_census_pre, E = Population2019/100000,
                     
                     control.predictor=list(link = link),
                     
                     #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                     
                     control.inla=list(strategy="simplified.laplace"),
                     #verbose=TRUE
)

model10a3_pre1 <-inla(formula3, family="poisson", data=covid_biweekly_census_pre3, E = Population2019/100000,
                      
                      control.predictor=list(link = link1),
                      
                      #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                      
                      control.inla=list(strategy="simplified.laplace"),
                      #verbose=TRUE
)

model10a3_pre2<-inla(formula3, family="poisson", data=covid_biweekly_census_pre4, E = Population2019/100000,
                     
                     control.predictor=list(link = link2),
                     
                     #control.compute=list(dic=TRUE,cpo = TRUE,waic = TRUE),
                     
                     control.inla=list(strategy="simplified.laplace"),
                     #verbose=TRUE
)

# pre_1 <- list(NULL)
# pre_1 <- length(3)
# pre_1[1] <- model10a2_pre
# pre_1[2] <- model10a3_pre
# pre_1[3] <- model10a4_pre
pre_type2_1 <- cbind(covid_biweekly_census,model10a2_pre2$summary.fitted.values)
pre_type3_1 <- cbind(covid_biweekly_census,model10a3_pre2$summary.fitted.values)
pre_type4_1 <- cbind(covid_biweekly_census,model10a4_pre3$summary.fitted.values)

write.csv(pre_type2_1,"E:/ResearchWork2/prediction/prediction1_type2.csv")
write.csv(pre_type3_1,"E:/ResearchWork2/prediction/prediction1_type3.csv")
write.csv(pre_type4_1,"E:/ResearchWork2/prediction/prediction1_type4.csv")

pre_type2_2 <- cbind(covid_biweekly_census,model10a2_pre1$summary.fitted.values)
pre_type3_2 <- cbind(covid_biweekly_census,model10a3_pre1$summary.fitted.values)
pre_type4_2 <- cbind(covid_biweekly_census,model10a4_pre2$summary.fitted.values)

write.csv(pre_type2_2,"E:/ResearchWork2/prediction/prediction2_type2.csv")
write.csv(pre_type3_2,"E:/ResearchWork2/prediction/prediction2_type3.csv")
write.csv(pre_type4_2,"E:/ResearchWork2/prediction/prediction2_type4.csv")


pre_type2_3 <- cbind(covid_biweekly_census,model10a2_pre$summary.fitted.values)
pre_type3_3 <- cbind(covid_biweekly_census,model10a3_pre$summary.fitted.values)
pre_type4_3 <- cbind(covid_biweekly_census,model10a4_pre$summary.fitted.values)

write.csv(pre_type2_3,"E:/ResearchWork2/prediction/prediction3_type2.csv")
write.csv(pre_type3_3,"E:/ResearchWork2/prediction/prediction3_type3.csv")
write.csv(pre_type4_3,"E:/ResearchWork2/prediction/prediction3_type4.csv")


