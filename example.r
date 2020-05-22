### EXAMPLE OF USE - RAVENNA
source("R/clisagri.r")

#read meteo data
meteo = read.csv("data/MeteoRavenna.csv", sep=",", header = TRUE)
# convert dates in meteo to class Date
meteo$DAY = as.Date(meteo$DAY, format="%d/%m/%Y")

######### calculate all indicators for fixed period s#########
ravenna = clisagri(meteo, lat = 41.5, drought.coef=NULL)

######### calculate indicators based on development stages DVS #########
source("R/phenology.r")
parameters = read.csv("data/ParametersRavenna.csv", header=TRUE, sep=",")
sowing = read.csv("data/SowingRavenna.csv", header=TRUE, sep=",")
sowing$DAY = as.Date(sowing$DAY, format="%d/%m/%Y")
meteo = phenology(meteo, parameters, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS)
ravenna.dvs = clisagri(meteo,41.5,sowing=sowing$DAY)

######## phenological simulation tailored for breeders #########
source("R/phenology.breeder.r")
parameters = read.csv("data/ParametersRavenna.breeder.csv", header=TRUE, sep=",")
meteo = phenology.breeder(meteo, parameters, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS, dvs.end = 6)
ravenna.dvs = clisagri(meteo,41.5,sowing=sowing$DAY)


######### reproduction of figures from manuscript Ceglar et al. (Clisagri: An R package for agro-climate services)
######### TBD
