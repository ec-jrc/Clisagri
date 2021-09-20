clim.plot = function(data, type=NA)
{
  require(ggplot2)
  require(scales)

  tsum1 = unique(data$tsum1)
  tsum2 = unique(data$tsum2)
  
  z <- outer(seq(-1,1,length=length(tsum1)), seq(-1,1,length=length(tsum2)), FUN = fun_xy)
  
  Z = cbind(expand.grid(tsum1, tsum2), as.vector(z[,]))
  names(Z) = c("tsum1", "tsum2", "color")

  gz=ggplot(Z, aes(x=tsum1,y=tsum2)) +
    geom_tile(fill=Z$color)+
    xlab("TSUM1 [GDD]")+
    ylab("TSUM2 [GDD]")+
    theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
    scale_x_continuous(limits = c(min(tsum1)-10, max(tsum1)+10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(tsum2)-10,max(tsum2)+10), expand = c(0, 0)) 
  
  data$color = NA
  for(i in 1:length(tsum1))
    for(j in 1:length(tsum2))
      data$color[which(data$tsum1==tsum1[i]&data$tsum2==tsum2[j])] = z[i,j]
  
  g2 = ggplot(data, aes(x=year, y=value)) + 
    geom_jitter(colour=data$color, size=length(unique(data$year)) * 2 / (length(tsum1)*length(tsum2)), shape=1, alpha=0.8) +
    theme(legend.position = "none") +
    theme_bw() +
    ggtitle(paste("Indicator type: ", type))+
    ylab("Value")+
    xlab("Year")
  
  print(g2) 
  dev.new() 
  print(gz)
  
  return(list(codes=gz, timeseries=g2))
}

dvs.plot = function(data)
{
  
  require(ggplot2)
  require(scales)

  tsum1 = unique(data$tsum1)
  tsum2 = unique(data$tsum2)
  
  z <- outer(seq(-1,1,length=length(tsum1)), seq(-1,1,length=length(tsum2)), FUN = fun_xy)
  
  Z = cbind(expand.grid(tsum1, tsum2), as.vector(z))
  names(Z) = c("tsum1", "tsum2", "color")
  
  gz=ggplot(Z, aes(x=tsum1,y=tsum2)) +
    geom_tile(fill=Z$color)+
    xlab("TSUM1 [GDD]")+
    ylab("TSUM2 [GDD]")+
    theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
    scale_x_continuous(limits = c(min(tsum1)-10, max(tsum1)+10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(tsum2)-10,max(tsum2)+10), expand = c(0, 0))
  
  data$color = NA
  for(i in 1:length(tsum1))
    for(j in 1:length(tsum2))
      data$color[which(data$tsum1==tsum1[i]&data$tsum2==tsum2[j])] = z[i,j]
  
  dvs.stages = c(0.01, 0.1, 0.34, 0.72, 0.82, 0.92, 1.16, 2.0)
  dvs.stages = data.frame(DVS=dvs.stages, NAME=c("emergence", "tillering", "stem elongation", "booting", "heading", "flowering", "grain filling", "maturity"))
  
  new_line <- element_line(color = "darkgrey", 
                           size = 0.5,
                           linetype = 1,
                           lineend = "round")
  
  
  g1 = ggplot(data, aes(x=DAY, y=DVS, group=color)) +
    geom_path(color=data$color, size=0.2, alpha=0.5) +
    theme(legend.position = "none") +
    theme_bw() +
    ylab("Development stage")+
    xlab("DAY") +
    scale_x_date(breaks=pretty_breaks(), date_minor_breaks="1 month") +
    scale_y_continuous(name="Development stage", breaks=dvs.stages[,1], labels=dvs.stages[,2]) + 
    theme(panel.grid = new_line)+
    xlim(min(data$DAY), max(data$DAY[!is.na(data$DVS)]))+
    coord_cartesian(ylim = c(0,2)) 
  
  print(gz)
  dev.new()
  print(g1)
  
  return(list(codes=gz, dvs.plot=g1))
}

fun_xy <- function(x, y){
    
  R <- (x+1)/2
  G <- (1-y)/2
  B <- (1-x)/2
  A <- 1- 0.5*exp(-(x^2+y^2)/0.2)
  
  rgb(R, G, B, A)
    
}
