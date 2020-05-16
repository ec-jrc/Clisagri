#
# Copyright 2019 European Union
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission –
# subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at: https://joinup.ec.europa.eu/software/page/eupl
#
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is
# istributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
# express or implied. See the Licence for the specific language governing permissions and limitations under
# the Licence.
#

#*********************************************************
#             RAGROCLIM
#*********************************************************

#*********************************************************
# Climate calculation ???
#
# Parameters :
#  - meteo        [INPUT] 	- daily meteorological data
#  - lat          [INPUT] 	- latitude of location (degrees)
#  - types        [INPUT]   - vector with 16 types of indicators
#  - sowing       [INPUT]   - vector of sowing dates, to be provided if calculation is based on DVS)
#  - drought.coef [INPUT]   - not yet implemented
# RETURN :
# [DATA.FRAME]    - with following structure :
#                     type of indicator
#                     year
#                     cumulate of precipitation in relevant period
#                     cumulate of reference evapotranspiration in relevant period
#                     value of indicator
#*********************************************************
agro.clim = function(meteo, lat, types=c(1:16), sowing=NULL,  drought.coef=NULL)
{
  require(lubridate)
  
  if(!is.null(drought.coef))
    load(drought.coef, envir=environment(agro.clim))

  if(!is.null(meteo$DVS) && is.null(sowing))
    return("SOWING DATE needs to be present if calculation is based on DVS")

  meteo$YEAR=as.integer(format(meteo$DAY, "%Y"))
  meteo$MONTH=as.integer(format(meteo$DAY, "%m"))

  years = unique(meteo$YEAR)

  indexes = c()
  for(type in types)
  {
    if(type >=1 && type <=6)
      meteo.ind = agro.clim.type_1_6(meteo, lat, type, years=years, sowing=sowing)
    if(type == 7)
      meteo.ind = agro.clim.type_7(meteo, lat, type, years=years,  sowing=sowing)
    if(type == 8)
      meteo.ind = agro.clim.type_8(meteo, lat, type, years=years,  sowing=sowing)
    if(type == 9)
      meteo.ind = agro.clim.type_9(meteo, lat, type, years=years,  sowing=sowing)
    if(type == 10)
      meteo.ind = agro.clim.type_10(meteo, lat, type, years=years,  sowing=sowing)
    if(type == 11)
      meteo.ind = agro.clim.type_11(meteo, lat, type, years=years,  sowing=sowing)
    if(type >=12 && type <=13)
      meteo.ind = agro.clim.type_12_13(meteo, lat, type, years=years,  sowing=sowing)
    if(type == 14)
      meteo.ind = agro.clim.type_14(meteo, lat, type, years=years, sowing=sowing)
    if(type >=15 && type <=16)
      meteo.ind = agro.clim.type_15_16(meteo, lat, type, years=years,  sowing=sowing)

    indexes = rbind(indexes, meteo.ind)
  }
  indexes = indexes[order(indexes$year, indexes$type),]
  return(indexes)
  
}

  

agro.clim.type_1_6 = function(meteo, lat, type, years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type,
                                                year=year,
                                                prec=sum(meteo$PRECIPITATION[which_],
                                                         na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat),
                                                             tmean=meteo$TEMPERATURE_AVG[which_],
                                                             trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],
                                                             rr=meteo$PRECIPITATION[which_]))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA))
    }
    
    meteo.ind = cbind(meteo.ind, value=spei.period(meteo.ind, drought.coef=get0("d.m", ifnotfound = NULL), type=type))
  }
  
  if(!is.null(meteo$DVS) && type > 1)
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type,
                                                year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_],
                                                           na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_],lat),
                                                             tmean=meteo$TEMPERATURE_AVG[which_],
                                                             trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],
                                                             rr=meteo$PRECIPITATION[which_]))))
      else
        meteo.ind = rbind(meteo.ind,
                          data.frame(type=type,
                                     year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                     prec=NA,
                                     etp=NA))
    }
    
    meteo.ind = cbind(meteo.ind, value=spei.period(meteo.ind, drought.coef=get0("d.dvs", ifnotfound = NULL), type=type))
  }
  return(meteo.ind)
}
  
agro.clim.type_7 = function(meteo, lat, type, years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type,
                                                year=year,
                                                prec=sum(meteo$PRECIPITATION[which_],
                                                         na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat),
                                                             tmean=meteo$TEMPERATURE_AVG[which_],
                                                             trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],
                                                             rr=meteo$PRECIPITATION[which_])),
                                                value = sum(meteo$PRECIPITATION[which_], na.rm=TRUE)))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA,value=NA))
    }
  }
  return(meteo.ind)
}

agro.clim.type_8 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=10)) ))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA,value=NA))
    }
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat),
                                                             tmean=meteo$TEMPERATURE_AVG[which_],
                                                             trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],
                                                             rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=10))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type,
                                                year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec=NA,
                                                etp=NA,
                                                value=NA))
      
    }
  }
  return(meteo.ind)
}

agro.clim.type_9 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=40)) ))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA,value=NA))
    }
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=40))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])), 
                                                prec=NA, etp=NA, value=NA))
      
    }
  }
  return(meteo.ind)
}

agro.clim.type_10 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=5)) ))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA,value=NA))
    }
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=5))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])), prec=NA, etp=NA, value=NA))
      
    }
  }
  return(meteo.ind)
}

agro.clim.type_11 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=40)) ))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA,value=NA))
    }
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = length(which(meteo$PRECIPITATION[which_]>=40))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec=NA, etp=NA, value=NA))
      
    }
  }
  return(meteo.ind)
}

agro.clim.type_12_13 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
      {
        ml = rle(as.vector(meteo$PRECIPITATION[which_]>=5))
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = ifelse(length(which(ml$values==TRUE)) > 0, max(ml$lengths[which(ml$values==TRUE)]), 0)))
      } else {
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA,value=NA))
      }
    }
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
      {
        ml = rle(as.vector(meteo$PRECIPITATION[which_]>=5))
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value = ifelse(length(which(ml$values==TRUE)) > 0, max(ml$lengths[which(ml$values==TRUE)]), 0)))
      } else {
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec=NA, etp=NA, value=NA))
      }
    }
  }
  return(meteo.ind)
}

agro.clim.type_14 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value=length(which(meteo$TEMPERATURE_MIN[which_]<=2))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA, value=NA))
    }
    
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value=length(which(meteo$TEMPERATURE_MIN[which_]<=2))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec=NA, etp=NA, value=NA))
      
    }
  }
  return(meteo.ind)
}

agro.clim.type_15_16 = function(meteo, lat, type,  years, sowing=NULL)
{
  a = assign.type(type)
  meteo.ind = c()
  if(is.null(meteo$DVS))
  {
    for(year in years) {
      which_ = check.completeness.m(a, meteo, year)
      
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value=length(which(meteo$TEMPERATURE_MAX[which_]>=hs(type)))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=year, prec=NA, etp=NA, value=NA))
    }
    
  }
  
  if(!is.null(meteo$DVS))
  {
    b = assign.dvs(type)
    for(i in 1:length(sowing))
    {
      which_ = check.completeness.dvs(b,meteo, sowing[i])
      if(length(which_)>0)
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year = ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec = sum(meteo$PRECIPITATION[which_], na.rm=TRUE),
                                                etp=sum(evap(ra=raest(meteo$MONTH[which_], lat), tmean=meteo$TEMPERATURE_AVG[which_], trange=meteo$TEMPERATURE_MAX[which_]-meteo$TEMPERATURE_MIN[which_],rr=meteo$PRECIPITATION[which_])),
                                                value=length(which(meteo$TEMPERATURE_MAX[which_]>=hs(type)))))
      else
        meteo.ind = rbind(meteo.ind, data.frame(type=type, year=ifelse(yday(sowing[i]) > 180, year(sowing[i])+1, year(sowing[i])),
                                                prec=NA, etp=NA, value=NA))
    }
  }
  return(meteo.ind)
}

#*********************************************************
# Check completness
#
# Parameters :
#  - a            [INPUT]   - start and end month of period for calculating indicator value
#  - meteo        [INPUT] 	- daily meteorological data
#  - year         [INPUT] 	- year
# RETURN : index vector for reading values from meteo data frame
#*********************************************************
check.completeness.m = function(a, meteo, year)
{
  if(a$start.month <= a$end.month && a$start.month >8)
  {
    which_ = which((meteo$YEAR==year-1 & meteo$MONTH >= a$start.month & meteo$MONTH <= a$end.month))
    months = unique(meteo$MONTH[which_])
    if(length(which_) > 0)
      if(length(months) != (-a$start.month+a$end.month+1))
        which_ = c()
  }

  if(a$start.month >= a$end.month)
  {
    which_ = which((meteo$YEAR==year-1 & meteo$MONTH >= a$start.month) | (meteo$YEAR==year & meteo$MONTH <= a$end.month))
    months = unique(meteo$MONTH[which_])
    if(length(which_) > 0)
      if(length(months) != (13-a$start.month+a$end.month))
        which_ = c()
  }

  if(a$start.month <= a$end.month && a$start.month <=8)
  {
    which_ = which(meteo$YEAR==year & meteo$MONTH >= a$start.month & meteo$MONTH <= a$end.month)
    months = unique(meteo$MONTH[which_])
    if(length(which_) > 0)
      if(length(months) != (a$end.month-a$start.month + 1))
        which_ = c()
  }
  return(which_)
}

#*********************************************************
# HEat stress critical temperatures
#*********************************************************
hs = function(type)
{
  hs = vector(length=16)
  hs[1:14] = NA
  hs[15] = 31
  hs[16] = 35
  return(hs[type])
}

#*********************************************************
# Check completness
#
# Parameters :
#  - b            [INPUT]   - start and end BBCH for calculating indicator value
#  - meteo        [INPUT] 	- daily meteorological data
#  - year         [INPUT] 	- year
# RETURN : index vector for reading values from meteo data frame
#*********************************************************
check.completeness.dvs = function(b, meteo, sow)
{

  which_ = which((meteo$DAY >= sow & meteo$DAY <= sow+360 & meteo$DVS >= b$start.dvs & meteo$DVS <= b$end.dvs))
#    if(max(meteo$DVS[which_],na.rm=TRUE) < (b$end.dvs-0.1))
#      which_ = c()
  if(length(which_)>0 && min(meteo$DVS[which_],na.rm=TRUE) < b$start.dvs)
    which_ = c()

  return(which_)
}


#*********************************************************
# Assign type
#
# Parameters :
#  - type            [INPUT]   - type of parameter 
# RETURN : list containing the start and end months charancterizing the period for indicator calculation
#*********************************************************
assign.type = function(type)
{
  period=matrix(data=c(9,10,12,4,5,11,9,12,12,4,4,4,5,4,3,5,10,11,3,4,6,6,10,3,3,6,6,5,6,5,5,6), ncol=2)
  return(list(start.month=period[type,1],end.month=period[type,2]))
}

#*********************************************************
# Assign DVS
#
# Parameters :
#  - type            [INPUT]   - type of parameter 
# RETURN : list containing the start and end BBCH values, charancterizing the development stage for indicator calculation
#*********************************************************
assign.dvs = function(type)
{
  period=matrix(data=c(NA,0,20,30,50,0,NA,20,20,51,51,51,69,41,39,69,NA,9,29,49,89,89,NA, 29,29,89,89,69,89,69,69,89), ncol=2)
  return(list(start.dvs=period[type,1],end.dvs=period[type,2]))
}

#*********************************************************
# HDD Period
#
# Parameters :
#  - meteo          [INPUT]   - meteo data frame
#  - thr            [INPUT]   - temperature threshold for calculating heat degree days
# RETURN : heat degree days 
#*********************************************************
hdd.period = function(meteo,thr)
{
  meteo$TEMPERATURE_MAX[meteo$TEMPERATURE_MAX<=thr] = 0
  return(sum(meteo$TEMPERATURE_MAX,na.rm=TRUE))
}

#*********************************************************
# SPEI Period
#
# Parameters :
#  - meteo.spei    [INPUT]   - input meteorological data fro calculation of SPEI
#  - type          [INPUT]   - type of drought indicator
#  - drought.coef  [INPUT]   - not yet implemented
# RETURN : SPEI values for specified type of drought indicator
#*********************************************************
spei.period = function(meteo.spei, type, drought.coef=NULL)
{
  meteo.spei$spei = meteo.spei$prec - meteo.spei$etp
  meteo.spei$spei.d = NA
  if(length(which(is.na(meteo.spei$spei))) < 0.5*length(meteo.spei$spei))
  {
    if(!is.null(drought.coef))
    {
      mycdf.clino = drought.coef[[type]][[1]]
      pze = drought.coef[[type]][[2]]
      max.obs = drought.coef[[type]][[3]]
      min.obs = drought.coef[[type]][[4]]
    }
    else
    {
      mycdf.clino = ecdf(c(meteo.spei$spei))
      pze = length(which(meteo.spei$spei==0))/length(meteo.spei$spei)
      min.obs = min(meteo.spei$spei, na.rm=TRUE)
      max.obs = max(meteo.spei$spei, na.rm=TRUE)
    }

  meteo.spei$spei[meteo.spei$spei<= min.obs] = min.obs + 0.01
  meteo.spei$spei[meteo.spei$spei>= max.obs] = max.obs - 0.01

  for(year in meteo.spei$year)
    meteo.spei$spei.d[meteo.spei$year==year]=qnorm(pze+(1-pze)*pnorm(qnorm(mycdf.clino(meteo.spei$spei[meteo.spei$year==year]))))

  my_inf=which(is.infinite(meteo.spei$spei.d))
  meteo.spei$spei.d[my_inf]=NA
  }

  return(meteo.spei$spei.d)
}

#*********************************************************
# SPI Period
#
# Parameters :
#  - meteo.spei    [INPUT]   - input meteorological data fro calculation of SPI
#  - type          [INPUT]   - type of drought indicator
#  - drought.coef  [INPUT]   - not yet implemented
# RETURN : SPI values for specified type of drought indicator
#*********************************************************
spi.period = function(meteo.spei, type, drought.coef=NULL)
{

  meteo.spei$spei = meteo.spei$prec
  meteo.spei$spei.d = NA

  if(!is.null(drought.coef))
  {
    mycdf.clino = drought.coef[[type]][[1]]
    pze = drought.coef[[type]][[2]]
    max.obs = drought.coef[[type]][[3]]
    min.obs = drought.coef[[type]][[4]]
  }
  else
  {
    mycdf.clino = ecdf(c(meteo.spei$spei))
    pze = length(which(meteo.spei$spei==0))/length(meteo.spei$spei)
    min.obs = min(meteo.spei$spei, na.rm=TRUE)
    max.obs = max(meteo.spei$spei, na.rm=TRUE)
  }


  meteo.spei$spei[meteo.spei$spei<= min.obs] = min.obs + 0.01
  meteo.spei$spei[meteo.spei$spei>= max.obs] = max.obs - 0.01

  for(year in meteo.spei$year)
    meteo.spei$spei.d[meteo.spei$year==year]=qnorm(pze+(1-pze)*pnorm(qnorm(mycdf.clino(meteo.spei$spei[meteo.spei$year==year]))))

  my_inf=which(is.infinite(meteo.spei$spei.d))
  meteo.spei$spei.d[my_inf]=NA

  return(meteo.spei$spei.d)
}


#*********************************************************
# Global radiation
#
# Parameters :
#  - c      [INPUT]   - vector of months
#  - lat    [INPUT]   - latitude
# RETURN :
#   value of extraterrestrial solar radiation
#*********************************************************
raest=function(c,lat){
  J <- as.integer(30.5 * c - 14.6)
  delta <- 0.409 * sin(0.0172 * J - 1.39)
  dr <- 1 + 0.033 * cos(0.0172 * J)
  latr <- lat/57.2957795
  sset <- -tan(latr) * tan(delta)
  omegas <- sset * 0
  omegas[sset >= {-1} & sset <= 1] <- acos(sset[sset >= { -1 } & sset <= 1])
  omegas[sset < { -1 }] <- max(omegas)
  Ra <- 37.6 * dr * (omegas * sin(latr) * sin(delta) +  cos(latr) * cos(delta) * sin(omegas))
  return(Ra)
}

#*********************************************************
# Evapotraspiration
#
# Parameters :
#  - ra      [INPUT]   - global radiation
#  - tmean   [INPUT]   - temperature mean
#  - trange  [INPUT]   - temperature range
# RETURN :
#   value of globalevapotraspiration
#*********************************************************
#evapotranspiration
evap=function(ra,tmean,trange,rr){
  trange[trange<=0] = 0.1
  ev = vector(length=length(tmean))
  #ev = 0.0013*0.408*ra*(tmean+17)*(trange-0.01239)^0.76
  #ev = 0.0013*0.408*ra*(tmean+17)*(trange-0.0123*rr)^0.76
  w_ = which(trange>=0.0874*rr)
  ev[w_] = 0.0019*0.408*ra[w_]*(tmean[w_]+21.05849)*(trange[w_]-0.0874*rr[w_])^0.6278
  w_1 =  which(trange<0.0874*rr)
  if(length(w_1)>0)
    ev[w_1] = 0
  return(ev)
}



#*********************************************************
# Convert DVS value from phenological model to BBCH scale
#
# Parameters :
#  - dvs      [INPUT]   - simualted development stage using phenology model
# RETURN : converted BBCH development stages
#
#*********************************************************
DVS2BBCH = function(dvs, dvs.end=2)
{
    if(dvs.end>2) {
      BBCH = matrix(data=c(0,1,2,3,4,5,6,0,10,20,30,40,60,89), ncol=2)
    }
    else {
      BBCH = matrix(data=c(0,0.01,0.10,0.34,0.70,0.81,0.90,1.11,2, 0,9,20,30,40,50,60,70,89), ncol=2)
    }
    BBCH_ = seq(0,dvs.end,by=0.005)
    DVS_ = c()
    for(i in 1:length(BBCH_))
    {
        which_0 = tail(which(BBCH[,1] <= BBCH_[i]), n=1)
        if(BBCH_[i] < dvs.end)
            DVS_ = c(DVS_, ((BBCH[which_0+1,2]-BBCH[which_0,2]) / (BBCH[which_0+1,1]-BBCH[which_0,1])) * (BBCH_[i] - BBCH[which_0,1]) + BBCH[which_0,2])
        else
            DVS_ = c(DVS_, 89)
    }

    dvs.0 = unlist(apply(X=t(dvs), FUN=function(X,DVS,BBCH) ifelse(is.na(X), NA, DVS[which.min(abs(X-BBCH_))]), MARGIN=2, DVS=DVS_,BBCH=BBCH_))
    return(dvs.0)
}

# get observed date of specified DVS
dvs.date = function(meteo, dvs=69)
{
  require(lubridate)
  meteo$YEAR=as.integer(format(meteo$DAY, "%Y"))
  meteo$MONTH=as.integer(format(meteo$DAY, "%m"))
  
  years = unique(meteo$YEAR)
 
  dvs_ = data.frame()
  for(year in years)
  {
    meteo_ = meteo[meteo$YEAR==year,]
    which_ = which.min(abs(meteo_$DVS-dvs))
    diff = abs(meteo_$DVS-dvs)[which_]
    if(length(which_) > 0 && diff < 2)
      dvs_ = rbind(dvs_, data.frame(year=year,dvs=dvs,doy=yday(meteo_$DAY[which.min(abs(meteo_$DVS-dvs))])))
    else
      dvs_ = rbind(dvs_, data.frame(year=year,dvs=dvs,doy=NA))
  }
  return(dvs_)
}

