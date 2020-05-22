### EXAMPLE OF USE - RAVENNA
source("clisagri.r")

#read meteo data
meteo = read.csv("MeteoRavenna.csv", sep=",", header = TRUE)
# convert dates in meteo to class Date
meteo$DAY = as.Date(meteo$DAY, format="%d/%m/%Y")

######### calculate all indicators for fixed period s#########
ravenna = clisagri(meteo, lat = 41.5, drought.coef=NULL)

######### calculate indicators based on development stages DVS #########
source("phenology.r")
parameters = read.csv("ParametersRavenna.csv", header=TRUE, sep=",")
sowing = read.csv("SowingRavenna.csv", header=TRUE, sep=",")
sowing$DAY = as.Date(sowing$DAY, format="%d/%m/%Y")
meteo = phenology(meteo, parameters, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS)
ravenna.dvs = agro.clim(meteo,41.5,sowing=sowing$DAY)

######## phenological simulation tailored for breeders #########
source("phenology.breeder.r")
parameters = read.csv("ParametersRavenna.breeder.csv", header=TRUE, sep=",")
meteo = phenology.breeder(meteo, parameters, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS, dvs.end = 6)
ravenna.dvs = clisagri(meteo,41.5,sowing=sowing$DAY)


##########  Figure 2 in manuscript ####### 
#### development stage simulation using standard phenology model ######
source("hazard.map.r")

source("phenology.r")
parameters = read.csv("ParametersRavenna.csv", header=TRUE, sep=",")

dvs.ravenna = dvs.calc(meteo=meteo.ravenna, sowing="2002-10-30", parameters=parameters, latitude=44.42, r.tsum1=seq(600,1000,by=20), r.tsum2=seq(600, 1000, by=20), parallel=TRUE, ncores=10)


