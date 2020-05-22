% File src/library/base/man/clisagri.Rd
\name{clisagri}
\alias{clisagri}
\title{Calculate a set of agro-climatic indicators quantifying (un)favourable climatic conditions during sensitive stages of durum wheat growth}
\description{
  Function \code{clisagri} provides a set of agro-climatic indicators characterizing weather conditions during sensitive stages of durum wheat growth. Three groups of indicators are distinguished: hydrological balance, excessive wetness and temperature indicators (heat stress and cold shocks). The timing of sensitive phases can be determined based on static dates or dynamic dates based on phenological development stage. Phenological stages can be calculated using \code{phenology}, for generic wheat variety or user specified variety (see \code{phenology} and Ceglar et al. (2019) for more details).\n
  Input meteorological data is provided in \code{meteo} data frame, where the date (DAY, YYYY-MM-DD), maximum daily temperature (TEMPERATURE_MAX, deg. C), minimum daily temperature (TEMPERATURE_MIN, deg. C), average daily temperature (TEMPERATURE_AVG, deg. C), daily precipitation (PRECIPITATION, mm) and optionaly crop development stage (DVS, unitless) are provided. 
  Static timing of sensitive development stages is pre-defined based on average inter-annual timing of occurrence, and is generally limited to specific months of pre-sowing and growing period of winter wheat. Current version has a set of default values for specified months. Dynamic phenological timing can be provided either by user or calculated using \code{phenology}. In the first case, user needs to add column DVS in \code{meteo} data frame, corresponding to development BBCH development stage of wheat for each day of the growing season. In the latter case, user needs to specify a minimum set of variety specific parameters (see \code{phenology} for more details). \n
  Indicators for hydrological balance are based on calculated Standardized precipitation evapotranspiration index (SPEI). These are calculated for static periods, if no DVS is provided, and belong to types 1-6: 1 - hydrological balance during pre-sowing period (September and October), 2 - hydrological around sowing and emergence (October and November), 3 - hydrological balance during the tillering period (December-March), 4 - hydrological balance between stem elongation and booting period (April-May), 5 - hydrological balance between heading and physiological maturity (May-June), and 6 - hydrological balance between sowing and physiological maturity (November-June). In case DVS is provided, the same hydrological balance indicators are calculated during specified developent stages: 2 - between BBCH=0 and BBCH=9, 3 - between BBCH=20 and BBCH=29, 4 - between BBCH=30 and BBCH=49, 5 - between BBCH=50 and BBCH=89, and 6 - between BBCH=0 and BBCH=89. \n
  Indicators for excessive rainfall belong to types between 7 and 13. Type 7 can only be calculated for static period as it characterizes pre-sowing conditions (rainfall amount between September and October). Other types are calculated for static periods, if no DVS is provided: 8 - number of days with rain above 10 mm during tillering period (December - March), 9 - number of days with rain above 40 mm during tillering (December - March), 10 - number of days with rain above 5 mm between heading and physiological maturity (April-June), 11 - number of days with rain above 40 mm between heading and physiological maturity, 12 - maximum number of consecutive days with rain above 5 mm between heading and end of flowering (April - May), and 13 - maximum number of consecutive days with rain above 5 mm between end of flowering and physiological maturity (May - June). In case DVS is provided, the same excessive rainfall indicators are calculated during specified developent stages: types 8 and 9 - between BBCH=20 and BBCH=29, types 10 and 11 - between BBCH=51 and BBCH=89, type 12 bwetween BBCH=51 and BBCH=69, and type 13 betwen BBCH=69 and BBCH=89. \n
  Cold stress indicator belongs to type 14, and is based on number of days with minimum daily temperature below 2 deg. C between booting and flowering. In case no DVS is provided, the period between April and May is considered, while period between BBCH=41 and BBCH=69 is considered when DVS is provided. \n
  Heat stress indicators belong to types 15 and 16: 15 - number of hot days with maximum daily temperature above 28 deg. C between end of stem elongation and end of flowering (March - May), 16 - number of hot days with maximum daily temperature above 28 deg. C during grain filling period (May - June). If DVS is provided, these indicators are calculated for the following periods: 15 - between BBCH=39 and BBCH=69, 16 - between BBCH=69 and BBCH=89. 
}
\usage{
  clisagri(meteo, lat, types, sowing)
}
\arguments{
  \item{meteo}{a data frame containing the following columns: DAY (YYYY-MM-DD), TEMPERATURE_MAX, TEMPERATURE_MIN, TEMPERATURE_AVG, PRECIPITATINO, DVS (optionaly)}
  \item{lat}{location latitude}
  \item{types}{specified types of indicators to be calculated; by default all indicators (1-16) are calculated }
  \item{sowing}{sowing dates for winter wheat}
}
\value{data frame with the following columns: type (type of indicator), year (always harvest year), prec (precipitation cumulates, mm), etp (potential evapotranspiration, mm), value (value of indicator, units depend on type of indicator)}
\seealso{
  \code{\link{phenology}}.
}
\examples{
#read meteo data
meteo = read.csv("MeteoFoggia.csv", sep=",", header = TRUE)

# convert dates in meteo to class Date
meteo$DAY = as.Date(meteo$DAY, format="%d/%m/%Y")

######### calculate all indicators for fixed period s#########
foggia = clisagri(meteo, lat = 41.5, drought.coef=NULL)

######### calculate indicators based on development stages DVS #########
parameters = read.csv("ParametersFoggia.csv", header=TRUE, sep=",")
sowing = read.csv("SowingFoggia.csv", header=TRUE, sep=",")
sowing$DAY = as.Date(sowing$DAY, format="%d/%m/%Y")

meteo = phenology(meteo, parameters, sowing$DAY, 41.5)
meteo$DVS = DVS2BBCH(meteo$DVS)
foggia.dvs = clisagri(meteo,41.5,sowing=sowing$DAY)

}


\name{phenology}
\alias{phenology}
\title{Calculate BBCH development stage of winter wheat based on observed meteorological conditions}
\description{
  Function \code{phenology} calculates development stage (BBCH) based on input daily temperature data. Phenological development of winter wheat is simulated based on the temperature accumulated above a base temperature (for winter wheat considered to be 0 deg.C). Development rate remains constant above a certain maximum effective temperature Tmax (for winter wheat considered to be 30 deg.C) and is corrected for the effects of vernalization and photoperiod (Ceglar et al., 2019). \n
  Parameters needed to calculate development stage need to be provided in f.variety data frame, which contains columns: PARAMETER_CODE, PARAMETER_XVALUE and PARAMETER_YVALUE. A set of 29 parameters need to be provided (see example parameter file). The description of parameters can be obtained by calling \code{phenology.parameters()} function. Generally, the most important parameters identifying different varieties are TSUM1 and TSUM2 (thermal requirements for vegetative and reproductive stages, respectivelly), while other parameters can maintain default values provided in parameter file. In case phenological observations are available, one can perform sensitivity analysis and calibration of other parameters as well (functionality not yet implemented in current version of package). 
}
\usage{
  phenology(meteo, f.variety, sowing, lat)
}
\arguments{
  \item{meteo}{a data frame containing the following columns: DAY (YYYY-MM-DD), TEMPERATURE_MAX, TEMPERATURE_MIN, TEMPERATURE_AVG, PRECIPITATINO, DVS (optionaly)}
  \item{f.variety}{a data frame containing the parameter values (depending on wheat variety) for phenological model. See \code{phenology.parameters() for necessary parameter table and values. }
  \item{sowing}{sowing dates for winter wheat}
  \item{lat}{location latitude}
}
\value{data frame with the same columns as meteo input, with additional column DVS specifying calculated development stage. DVS values range from 0 (sowing) - 2 (maturity). }
\seealso{
  \code{\link{phenology.parameters}}.
}
\examples{
#read meteo data
meteo = read.csv("MeteoFoggia.csv", sep=",", header = TRUE)

# convert dates in meteo to class Date
meteo$DAY = as.Date(meteo$DAY, format="%d/%m/%Y")


######### calculate indicators based on development stages DVS #########
parameters = read.csv("ParametersFoggia.csv", header=TRUE, sep=",")
sowing = read.csv("SowingFoggia.csv", header=TRUE, sep=",")
sowing$DAY = as.Date(sowing$DAY, format="%d/%m/%Y")

meteo = phenology(meteo, parameters, sowing$DAY, 41.5)

# convert simulated DVS to BBCH scale
meteo$DVS = DVS2BBCH(meteo$DVS) 

}



\name{DVS2BBCH}
\alias{DVS2BBCH}
\title{Convert between \code{phenology} development stage codes to BBCH development stage}
\description{
  Converter between WOFOST development stage (ranging from 0 to 2) to BBCH scale. 
}
\usage{
  DVS2BBCH(dvs)
}
\arguments{
  \item{dvs}{development stage as simulated by \code{phenology}}
}
\value{development stage in BBCH scale}
\seealso{
  \code{\link{phenology}}.
}
\examples{
}

\name{phenology.parameters}
\alias{phenology.parameters}
\title{Description of parameters used by phenological model \code{phenology} implemented in this package. }
\description{
  Display description of parameters used by phenological model \code{phenology} implemented in this package. See WOFOST manual (Boogard et al., 2014) for more detailed description. 
}
\usage{
  phenology.parameters()
}
\arguments{
}
\value{}
\seealso{
  \code{\link{phenology}}.
}
\examples{
}

