breeder.calc = function(x, meteo, sowing, parameters, types, latitude, year=NA, flowering=NA, maturity=NA)
{
  # calculate indicator of required type
  which_ = grep('TSUM_', parameters$PARAMETER_CODE)
  parameters$PARAMETER_YVALUE[which_] = x

  meteo = phenology.breeder(meteo, parameters, sowing$DAY, latitude)
  meteo$DVS = DVS2BBCH(meteo$DVS, dvs.end=6)
  
  loc.dvs = clisagri(meteo = meteo, types = types, lat = latitude, sowing=sowing$DAY)
  loc.dvs$flower = rep(stage(meteo,65), each=length(types))
  loc.dvs$maturity = rep(stage(meteo,89), each=length(types))
  
  if(!is.na(flowering))
  {
    loc.dvs$value[loc.dvs$flower<=flowering] = NA
  }
  
  if(!is.na(maturity))
    loc.dvs$value[loc.dvs$maturity<=maturity] = NA
  
  if(!is.na(year) & length(which(is.na(loc.dvs$value))) <= 0.25 * length(loc.dvs$value))
  {
    fit.hb = 0
    if(length(which(types>=2 & types <=6)) >= 1)
      for(tp in types[types>=2 & types <=6])
      {
        w_hb = which(loc.dvs$type==tp)
        w_hb_year = which(loc.dvs$type==tp & loc.dvs$year==year)
        fit.hb = fit.hb - (abs(loc.dvs$value[w_hb_year]) / max(abs(loc.dvs$value[w_hb]), na.rm=TRUE))^4
      }
    
    fit.nd = 0
    if(length(which(types>=8 & types <=16)) >= 1)
      for(tp in types[types>=8 & types <=16])
      {
        w_hb = which(loc.dvs$type==tp & !is.na(loc.dvs$value))
        w_hb_year = which(loc.dvs$type==tp & !is.na(loc.dvs$value) & loc.dvs$year==year)
        fit.nd = fit.nd - sqrt(1 / length.stage(meteo = meteo, ind = tp, sowing = sowing$DAY[loc.dvs$year[w_hb]==year]) * (loc.dvs$value[w_hb_year]))
      }
     fitness = (fit.nd + fit.hb)
  } else {
    fitness = -length(types)
  }
  
  return(fitness)
}

breeder.plot = function(data, year, types, parameters, x, sowing, latitude)
{
  require(ggplot2)
  # periods of different indicators
  period=matrix(data=c(NA,0,20,30,50,0,NA,20,20,51,51,51,69,41,39,69, NA,9,29,49,89,89,NA, 29,29,89,89,69,89,69,69,89), ncol=2)
  
  # random realization
  which_ = grep('TSUM_', parameters$PARAMETER_CODE)
  parameters$PARAMETER_YVALUE[which_] = x
  
  meteo = phenology.breeder(meteo, parameters, sowing$DAY, latitude)
  meteo$DVS.r = DVS2BBCH(meteo$DVS, dvs.end=6)
  
  meteo$DVS = meteo$DVS.r
  loc.dvs.r = clisagri(meteo = meteo, types = types, lat = latitude, sowing=sowing$DAY)
  
  # optimal solution
  parameters$PARAMETER_YVALUE[which_] = colMeans(data@solution)
  meteo = phenology.breeder(meteo, parameters, sowing$DAY, latitude)
  meteo$DVS.opt = DVS2BBCH(meteo$DVS, dvs.end=6)
  
  meteo$DVS = meteo$DVS.opt
  loc.dvs.opt = clisagri(meteo = meteo, types = types, lat = latitude, sowing=sowing$DAY)
  
  # dates of critical phenophases
  labels = data.frame(stage=c("SOW", "EM", "TIL", "STEM", "BOOT", "ANTH", "MAT"), bbch=c(0, 11, 21, 31, 41, 69, 89))
  
  which_ = which(meteo$DAY>=sowing$DAY[which(unique(loc.dvs.opt$year)==year)] & meteo$DAY <=sowing$DAY[which(unique(loc.dvs.opt$year)==year)]+360)
  plot.data = data.frame(day=meteo$DAY[which_], value=meteo$DVS.r[which_], type="R", ind="DVS")
  plot.data = rbind(plot.data, data.frame(day=meteo$DAY[which_], value=meteo$DVS.opt[which_], type="OPT", ind="DVS"))
  plot.data = plot.data[!is.na(plot.data$value),]
  
  plot.data$stage = NA
  
  for(i in 1:(dim(labels)[1]-1))
    plot.data$stage[plot.data$value >= labels[i,2] & plot.data$value < labels[i+1,2]] = paste(labels[i,1],"-",labels[i+1,1],sep="")

  plot.ind = c()
  for(type in types)
  {
    a = data.frame(day=meteo$DAY[which_], dvs=meteo$DVS.r[which_], value=loc.dvs.r$value[loc.dvs.r$year==year & loc.dvs.r$type==type])
    a$value[a$dvs<period[type,1] | a$dvs>period[type,2] | is.na(a$dvs)] = NA
    plot.ind = rbind(plot.ind, data.frame(day=meteo$DAY[which_], value=a$value, option="R", type=type, ind=ifelse(type<=6, "HYDRO. BALANCE", "NUMBER OF DAYS")))
    
    a = data.frame(day=meteo$DAY[which_], dvs=meteo$DVS.opt[which_], value=loc.dvs.opt$value[loc.dvs.opt$year==year & loc.dvs.opt$type==type])
    a$value[a$dvs<period[type,1] | a$dvs>period[type,2] | is.na(a$dvs)] = NA
    plot.ind = rbind(plot.ind, data.frame(day=meteo$DAY[which_], value=a$value, option="OPT", type=type, ind=ifelse(type<=6, "HYDRO. BALANCE", "NUMBER OF DAYS")))    
  }
  
  plot.data$stage <- factor(plot.data$stage , levels = unique(plot.data$stage))
  plot.data$type <- factor(plot.data$type , levels = c("OPT", "R"))
  xlim = c(min(plot.data$day), max(plot.data$day))
  g1 = ggplot(plot.data, aes(x=day, y=type, fill=stage)) + 
    geom_tile(aes(height=0.5)) + 
    #guides(fill = FALSE) +
    scale_x_date(date_breaks = "1 month", labels = date_format("%Y-%m"), limits=xlim)
  
  colors = brewer.pal(n=7,name="RdBu")
  breaks = c(-2,-1.5,-1,1,1.5,2)
  labels=c(-2,"",-1,1,"",2)
  levels = c(-1,-1.5,-2)
  g2 = ggplot(data=plot.ind[plot.ind$ind=="HYDRO. BALANCE",], aes(x=day, y=option, fill=value)) + 
    geom_tile(aes(height=0.5)) + 
    scale_fill_gradientn(breaks=breaks, colours=colors, labels=labels, limits=c(-2,2))+
    theme_bw()+
    xlab("Date")+
    ylab("")+
    facet_wrap(~type, ncol = 1)+
    scale_x_date(date_breaks = "1 month", labels = date_format("%Y-%m"), limits=xlim)
  
  colors = brewer.pal(n=9,name="Reds")
  mx = ceiling(max(plot.ind$value[plot.ind$type>=8], na.rm=TRUE) / 9)
  breaks = seq(0,max(plot.ind$value[plot.ind$type>=8], na.rm=TRUE),by=mx)
  labels=seq(0, max(plot.ind$value[plot.ind$type>=8], na.rm=TRUE), by=mx)
  limits=c(0, max(abs(plot.ind$value[plot.ind$type>=8]), na.rm=TRUE))
  g3 = ggplot(data=plot.ind[plot.ind$ind=="NUMBER OF DAYS",], aes(x=day, y=option, fill=value)) + 
    geom_tile(aes(height=0.5)) + 
    scale_fill_gradientn(breaks=breaks, colours=colors, labels=labels, limits=limits)+
    theme_bw()+
    xlab("Date")+
    ylab("")+
    facet_wrap(~type, ncol = 1)+
    scale_x_date(date_breaks = "1 month", labels = date_format("%Y-%m"), limits=xlim)
    
  
  plot<-plot_grid(g1, g2, g3, align='v', vjust=1, scale = 1, ncol = 1, axis="rlbt")
  
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

length.stage = function(meteo, ind, sowing)
{
  dvs_ = assign.dvs(ind)
  which_ = which(meteo$DAY==sowing)
  meteo_ = meteo[which_:(which_+365),]
  return(length(which(meteo_$DVS >= dvs_$start.dvs & meteo_$DVS <= dvs_$end.dvs)))
}

dvs.calc = function(meteo, sowing, parameters, latitude, r.tsum1=seq(300,700,by=10), r.tsum2=seq(450,1050,by=10), parallel=TRUE, ncores=10, mode="clim") 
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
  results = mclapply(X=1:dim(combin)[1], FUN=par.fun.dvs, mc.cores = ncores, meteo=meteo[meteo$DAY>=sowing & meteo$DAY <= (sowing+365),], type=type, combin=combin, sowing=sowing, parameters=parameters, lat=latitude)
  
  combin.res = c()
  for(i in 1:length(results))
    combin.res = rbind(combin.res, data.frame(tsum1=combin[i,1], tsum2=combin[i,2], results[[i]][,c(1,6)]))
  
  return(combin.res)
}

par.fun = function(X,meteo,lat,parameters,tsums,sowing,type)
{
  which_ = grep('TSUM_', parameters$PARAMETER_CODE)
  parameters$PARAMETER_YVALUE[which_] = tsums[,X]
  
  meteo = phenology.breeder(meteo, parameters, sowing$DAY, latitude)
  meteo$DVS = DVS2BBCH(meteo$DVS, dvs.end=6)
  
  foggia.dvs = clisagri(meteo = meteo, types = type, lat = latitude, sowing=sowing$DAY)
  foggia.dvs$flower = rep(stage(meteo,65), each=length(type))
  foggia.dvs$maturity = rep(stage(meteo,89), each=length(type))
  return(foggia.dvs) 
  
}

par.fun.dvs = function(X,meteo,lat,parameters,combin,sowing,type=1)
{
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE=="TSUM1"] = combin[X,1]
  parameters$PARAMETER_XVALUE[parameters$PARAMETER_CODE=="TSUM2"] = combin[X,2]
  meteo = phenology(meteo, parameters, sowing, latitude)
  
  return(meteo)
}
