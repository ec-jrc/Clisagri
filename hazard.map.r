hazard.map = function(data, type, year=NULL)
{
  require(ggplot2)
  require(scales)
  require(RColorBrewer)
  require(metR)
  
  if(type >= 2 && type <= 6){
    colors = brewer.pal(n=7,name="RdBu")
    breaks = c(-2,-1.5,-1,1,1.5,2)
    labels=c(-2,"",-1,1,"",2)
    limits=c(-max(abs(data$value),2, na.rm=TRUE), max(abs(data$value),2, na.rm=TRUE))
    levels = c(-1,-1.5,-2)
  }
  if(type >= 8 && type <= 13){
    colors = brewer.pal(n=9,name="Blues")
    breaks = seq(0,9,by=1)
    labels=c("", 1,"",3, "", 5, "", 7, "", "9>")
    data$value[data$value>9]=9
    limits=c(0, max(abs(data$value), 9, na.rm=TRUE))
  }
  if(type >= 14 && type <= 16){
    colors = brewer.pal(n=9,name="Reds")
    mx = ceiling(max(data$value, na.rm=TRUE) / 9)
    breaks = seq(0,max(data$value, na.rm=TRUE),by=mx)
    labels=seq(0, max(data$value, na.rm=TRUE), by=mx)
    limits=c(0, max(abs(data$value), na.rm=TRUE))
  }
  
  if(is.null(year)) 
    year = unique(data$year)[1]

    data = data[data$year == year,]
    data$maturity.s = loess(maturity~tsum1+tsum2, degree = 2, data = data)$fitted
    data$flow.s = loess(flower~tsum1+tsum2, degree = 2, data = data)$fitted
    
    g1 = ggplot(data, aes(x=tsum1, y=tsum2)) + 
      geom_tile(aes(fill=value)) + 
      stat_contour(aes(z=maturity.s) ,geom="contour", colour="black")+
      stat_contour(aes(z=flow.s) ,geom="contour", colour="grey")+
      scale_fill_gradientn(breaks=breaks, colours=colors, labels=labels, limits=limits)+
      #geom_text_contour(aes(z=flow)) +
      geom_text_contour(aes(z=maturity.s), stroke=0.2,size=4)+
      geom_text_contour(aes(z=flow.s), stroke=0.2,size=4,color="darkgrey", check_overlap=TRUE)+
      coord_fixed()+
      coord_cartesian(xlim = c(min(data$tsum1), max(data$tsum1)), ylim = c(min(data$tsum2), max(data$tsum2))) +
      theme_bw()+
      xlab("TSUM1 [GDD]")+
      ylab("TSUM2 [GDD]")
    
    print(g1)
    
  return(g1)
}



hazard.calc = function(meteo, sowing, parameters, latitude, r.tsum1=seq(300,700,by=10), r.tsum2=seq(450,1050,by=10), type=2, parallel=TRUE, ncores=10, mode="clim")
{
  require(parallel)
  
  # create a combination matrix
  combin = matrix(nrow=length(r.tsum1)*length(r.tsum2), ncol=2)
  l=1
  for(j in 1:length(r.tsum1))
    for(k in 1:length(r.tsum2))
    {
      combin[l,] = c(r.tsum1[j], r.tsum2[k])
      l=l+1
    }
  
  # calculate indicator of required type
  results = mclapply(X=1:dim(combin)[1], FUN=par.fun, mc.cores = ncores, meteo=meteo, type=type, combin=combin, sowing=sowing, parameters=parameters, lat=latitude)
  
  combin.res = c()
  for(i in 1:length(results))
    combin.res = rbind(combin.res, data.frame(tsum1=combin[i,1], tsum2=combin[i,2], results[[i]][,c(2,5,6,7)]))

  return(combin.res)
}

stage = function(meteo, dvs)
{
  days = c()
  meteo$YEAR=as.integer(format(meteo$DAY, "%Y"))
  for(year in unique(meteo$YEAR))
  {
    days = c(days, which.min(abs(meteo$DVS[meteo$YEAR==year]-dvs)))
  }
  days[1]=NA
  return(days[-length(days)])
  
}

dvs.calc = function(meteo, sowing, parameters, latitude, r.tsum1=seq(300,700,by=10), r.tsum2=seq(450,1050,by=10), parallel=TRUE, ncores=10, mode="clim") 
{
  require(parallel)
  if(class(sowing) != "Date")
   sowing = as.Date(sowing)
  
  # create a combination matrix
  combin = matrix(nrow=length(r.tsum1)*length(r.tsum2), ncol=2)
  l=1
  for(j in 1:length(r.tsum1))
    for(k in 1:length(r.tsum2))
    {
      combin[l,] = c(r.tsum1[j], r.tsum2[k])
      l=l+1
    }
  
  # calculate indicator of required type
  results = mclapply(X=1:dim(combin)[1], FUN=par.fun.dvs, mc.cores = ncores, meteo=meteo[meteo$DAY>=sowing & meteo$DAY <= (sowing+365),], combin=combin, sowing=sowing, parameters=parameters, lat=latitude)
  
  combin.res = c()
  for(i in 1:length(results))
    combin.res = rbind(combin.res, data.frame(tsum1=combin[i,1], tsum2=combin[i,2], results[[i]][,c(1,6)]))
  
  return(combin.res)
}

par.fun = function(X,meteo,lat,parameters,combin,sowing,type=1)
{
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE=="TSUM1"] = combin[X,1]
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE=="TSUM2"] = combin[X,2]
  meteo = phenology(meteo, parameters, sowing$DAY, latitude)

  meteo$DVS = DVS2BBCH(meteo$DVS, dvs.end = 2)
  foggia.dvs = agro.clim(meteo = meteo, types = type, lat = latitude, sowing=sowing$DAY)
  foggia.dvs$flower = stage(meteo,65)
  foggia.dvs$maturity = stage(meteo,89)
  return(foggia.dvs) 

}

par.fun.dvs = function(X,meteo,lat,parameters,combin,sowing)
{
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE=="TSUM1"] = combin[X,1]
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE=="TSUM2"] = combin[X,2]
  meteo = phenology(meteo, parameters, sowing, latitude)

  return(meteo)
}
