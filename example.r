### EXAMPLE OF USE - RAVENNA
sapply(list.files("R/", full.names = TRUE), source)

########read meteo data
meteo = read.csv("data/MeteoRavenna.csv", sep=",", header = TRUE)
######## convert dates in meteo to class Date
meteo$DAY = as.Date(meteo$DAY, format="%d/%m/%Y")

######### calculate all   indicators for fixed periods #########
ravenna = clisagri(meteo, lat = 41.5, drought.coef=NULL, types = c(1:16))

######### calculate indicators based on development stages DVS #########
parameters = read.csv("data/ParametersRavenna.csv", header=TRUE, sep=",")
sowing = read.csv("data/SowingRavenna.csv", header=TRUE, sep=",")
sowing$DAY = as.Date(sowing$DAY, format="%d/%m/%Y")
meteo = phenology(meteo, parameters, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS)
ravenna.dvs = clisagri(meteo,41.5,sowing=sowing$DAY)

######## value of indicator 16 based on different combinations of TSUM1/TSUM2
type.6.ravenna = hazard.calc(meteo = meteo, 
                              sowing = sowing, 
                              parameters = parameters, 
                              latitude = 41, 
                              r.tsum1 = seq(600,1000,by=50), 
                              r.tsum2 = seq(600,1000,by=50), 
                              type = 6, parallel = TRUE, 
                              ncores = 2)
hazard.map(type.6.ravenna, type = 6, year = 2012)
tser = clim.plot(type.6.ravenna, type = 6)

######## phenological simulation tailored for multiphase phenological model #########
parameters.breed = read.csv("data/ParametersRavenna.breeder.csv", header=TRUE, sep=",")
meteo = read.csv("data/MeteoRavenna.csv", sep=",", header = TRUE)
meteo$DAY = as.Date(meteo$DAY, format="%d/%m/%Y")
meteo = phenology.breeder(meteo, parameters.breed, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS, dvs.end = 6)
ravenna.dvs.breed = clisagri(meteo,41.5,sowing=sowing$DAY)
