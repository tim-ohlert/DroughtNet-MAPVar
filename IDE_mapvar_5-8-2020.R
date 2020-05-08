library(tidyverse)
library(plyr)
library(visreg)
library(cowplot)


d.fun <- function (var, k=0){
  var<-var[!is.na(var)]
  if(length(var)<3){return(NA)}
  if(sd(var)==0){return(0)}
  if (min(var)<0){var <- var + abs(min(var)) + 0.01*(max(var)-min(var))}
  if (length(var)<2){return(NA)}else{
    var <- var + k
    aux <- numeric(length(var)-1)
    for (i in 1:(length(var)-1)){
      aux [i] <- abs(log(var[i+1]/var[i]))
    }
  }
  return(mean(aux, na.rm=T))}

#DROUGHTNET

TPAMerge50yrs <- read.csv("C:/Users/ohler/Dropbox/IDE Meeting_Oct2019/data/precip/Yearly Precip_TPA/TPAMerge50yrs.csv")

climate_frame <- ddply(TPAMerge50yrs, .(site_code),
                       function(x)data.frame(
                         d.var = d.fun(x$totalPRE),
                         MAP = mean(x$totalPRE),
                         cv = sd(x$totalPRE)/mean(x$totalPRE),
                         sd = sd(x$totalPRE)
                       ))


reduced_npp <- read.csv("C:/Users/ohler/Dropbox/IDE Meeting_Oct2019/Full Biomass/anpp_clean_trt_ppt_02-26-2020.csv")

#only use plots that are not manipulated
control_biomass_DN <- reduced_npp %>% 
  subset( n_treat_days < 30 | trt == "Control")%>%
  subset(ANPP == "Yes" )%>%
  ddply( c("site_code", "year", "plot", "subplot"),
         function(x)data.frame(
           biomass = sum(x$mass)
         )) %>%
  ddply( c("site_code", "year", "plot"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  ddply( c("site_code", "year"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  ddply( c("site_code"),
         function(x)data.frame(
           biomass = mean(x$biomass),
           n.year = length(x$year)
         )) %>%
  subset(n.year > 2)%>%
  merge(climate_frame, by = "site_code",all=TRUE)%>%
  subset(MAP < 1200)

control_biomass_DN <- control_biomass_DN[!is.na(control_biomass_DN$n.year),]
control_biomass_DN$network <- "DroughtNet"

#NUTNET


CRU <- read.csv("C:/Users/ohler/Dropbox/NutNet data/climate/CRU/CRU-annual_2018-07-06.csv")
CRU <- subset(CRU, year>=1965 & year<=2014)

climate_frame <- ddply(CRU, .(site_code),
                       function(x)data.frame(
                         d.var = d.fun(x$precip),
                         MAP = mean(x$precip),
                         cv = sd(x$precip)/mean(x$precip),
                         sd = sd(x$precip)
                       ))



full.biomass <- read.csv("C:/Users/ohler/Dropbox/NutNet data/full-biomass-30-April-2020.csv")


control_biomass_NN <- full.biomass %>%
  subset(live != 0)%>%
  subset(year_trt == 0 | trt=="Control")%>%
  #block, plot, subplot
  ddply( c("site_code", "year", "plot", "subplot"),
         function(x)data.frame(
           biomass = sum(x$mass)
         )) %>%
  ddply( c("site_code", "year", "plot"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  ddply( c("site_code", "year"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  ddply( c("site_code"),
         function(x)data.frame(
           biomass = mean(x$biomass),
           n.year = length(x$year)
         )) %>%
  subset(n.year > 2)%>%
  merge(climate_frame, by = "site_code",all.x=TRUE)%>%
  subset(MAP < 1500)

control_biomass_NN$network <- "NutNet"

dn_nn <- rbind(control_biomass_DN,control_biomass_NN)


##
#MAP models
map <- lm(biomass~MAP, data=dn_nn)
map.2 <- lm(biomass~MAP+I(MAP^2), data=dn_nn)
d <- lm(biomass~d.var, data=dn_nn)
map.d <- lm(biomass~MAP * d.var, data=dn_nn)

AIC(map,map.2,d,map.d)
summary(map)
summary(map.2)
summary(d)
summary(map.d)
visreg(map.2)
ggplot(dn_nn, aes(x=MAP, y=biomass))+
  geom_smooth(method="lm")+
  geom_point(aes(size=d.var),shape=19)+
  ylab("ANPP (g/m2)")+
  theme_cowplot()


map.resid <- resid(map)
plot(dn_nn$d.var, map.resid)
abline(0,0)

##MAP variability correlation
d_map<- lm(d.var~MAP, data=dn_nn)
d_map.2<- lm(d.var~MAP+I(MAP^2), data=dn_nn)
d_map.3<- lm(d.var~MAP+I(MAP^2)+I(MAP^3), data=dn_nn)
AIC(d_map,d_map.2,d_map.3)
summary(d_map)

visreg(d_map)


##CV - D check
cv_d <- lm(d.var~cv, data = dn_nn)
summary(cv_d)
visreg(cv_d)

###################
##Annual precipitation


#DroughtNet
annual_control_DN <- reduced_npp %>% 
  subset( n_treat_days < 30 | trt == "Control")%>%
  subset(ANPP == "Yes")%>%
  ddply( c("site_code", "year", "plot", "subplot", "ppt"),
         function(x)data.frame(
           biomass = sum(x$mass)
         )) %>%
  ddply( c("site_code", "year", "plot", "ppt"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  ddply( c("site_code", "year", "ppt"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  merge(climate_frame, by = "site_code",all=TRUE)%>%
  subset(MAP < 1200)

annual_control_DN <- annual_control_DN[!is.na(annual_control_DN$ppt),]




#NutNet
weather_annual <- read.csv("C:/Users/ohler/Dropbox/NutNet data/climate/WeatherStation/Weather_annual_20191110.csv")


annual_control_NN <- full.biomass %>%
  subset(live != 0)%>%
  subset(year_trt == 0 | trt=="Control")%>%
  #block, plot, subplot
  ddply( c("site_code", "year", "plot", "subplot"),
         function(x)data.frame(
           biomass = sum(x$mass)
         )) %>%
  ddply( c("site_code", "year", "plot"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         )) %>%
  ddply( c("site_code", "year"),
         function(x)data.frame(
           biomass = mean(x$biomass)
         ))%>%
  merge(climate_frame, by = "site_code",all.x=TRUE)%>%
  merge(weather_annual, by = c("site_code","year"))%>%
  subset(MAP < 1500)

annual_control_NN <- annual_control_NN[!is.na(annual_control_NN$d.var),]




dn_nn_annual <- rbind(annual_control_DN,annual_control_NN)



ppt <- lmer(biomass~ppt+(1|site_code), data=dn_nn_annual)
ppt.2 <- lmer(biomass~ppt+I(ppt^2)+(1|site_code),data=dn_nn_annual)
d <- lmer(biomass~d.var+(1|site_code), data=dn_nn_annual)
ppt.d <- lmer(biomass~ppt * d.var+(1|site_code), data=dn_nn_annual)


AIC(ppt,ppt.2,d,ppt.d)
summary(ppt)
summary(d)
summary(ppt.d)
visreg(ppt, ylab = "NPP")
visreg(d)















