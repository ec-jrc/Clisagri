phenology.breeder = function(meteo, f.variety, sowing, lat)
{

  # definition of parameters
  Parameters = list()
  Parameters$TEFFMX = -99  # Max eff temperature for emergence
#  Parameters$TSUM1  = -99  # Temperature sum emergence to anthesis
#  Parameters$TSUM2  = -99  # Temperature sum anthesis to maturity
  Parameters$IDSL   = -99  # Switch for photoperiod (1) and vernalisation (2)
  Parameters$DLO    = -99  # Optimal day length for phenol. development
  Parameters$DLC    = -99  # Critical day length for phenol. development
  Parameters$DVSI   = -99  # Initial development stage
  Parameters$DVSEND = -99  # Final development stage
  Parameters$VERNSAT = -99
  Parameters$VERNBASE = -99
  Parameters$VERNDVS = -99

  #Parameters$DTSMTB = AfgenTrait() # Temperature response function for phenol.
  Parameters$crop_start_type = NA
  Parameters$crop_end_type=NA
  
  #TSUM dependent on the development stage
  Parameters$TSUM = matrix(nrow=10,ncol=2)
  Parameters$TBASE = matrix(nrow=10,ncol=2)
  Parameters$DVS_PERIOD = matrix(nrow=10,ncol=2)

  # definition of rate variables
  RateVariables = list()
  RateVariables$DTSUME = -99  # increase in temperature sum for emergence
  RateVariables$DTSUM  = -99  # increase in temperature sum
  RateVariables$DVR    = -99  # development rate
  RateVariables$VERNR = -99
  RateVariables$VERNFAC = -99
  RateVariables$DVRED = -99

  # definition of state variables
  StateVariables = list()
  StateVariables$DVS   = -99  # Development stage
  StateVariables$TSUM  = -99  # Temperature sum state

  # States which register phenological events
  StateVariables$DOS = -99 # Day of sowing
  StateVariables$DOE = -99 # Day of emergence
  StateVariables$DOA = -99 # Day of anthesis
  StateVariables$DOM = -99 # Day of maturity
  StateVariables$DOH = -99 # Day of harvest
  StateVariables$VERN = -99

  StateVariables$DOV = NA
  StateVariables$ISVERNALISED = TRUE

  StateVariables$STAGE = NA

  #########
  self = list(params=Parameters, rates=RateVariables, states=StateVariables, force_vernalisation=FALSE)

  # initialize DVS
  meteo$DVS = NA
  for(i in 1:length(sowing))
  {
    which_ = which(meteo$DAY==sowing[i])
    if(length(which_)>0)
    {
      meteo_ = meteo[which_:ifelse((which_+365) <= dim(meteo)[1],which_+365,dim(meteo)[1]),]
      self = initialize(self, sowing[i], f.variety, start.type="sowing", end.type="maturity")  
  
      days = 1
      drv = c()
      drv$LAT = lat
      while(self$states$DVS <= self$params$DVSEND && days < dim(meteo_)[1])
      {
        drv$TEMP = meteo_$TEMPERATURE_AVG[days]#(meteo_$TEMPERATURE_MAX[days] + meteo_$TEMPERATURE_MIN[days]) / 2
        self = calc_rates(self, meteo_$DAY[days], drv)
        self = integrate.dvs(self, meteo_$DAY[days])
        meteo$DVS[meteo$DAY==meteo_$DAY[days]] = self$states$DVS
        days = days + 1
      }
    }
  }

  return(meteo)
}

# crop start type:   "sowing", "emergence"
# crop_end_type :    "maturity", "harvest", "earliest"
initialize = function(self, day, f.variety, start.type="sowing", end.type="maturity")
{

  self$params$TEFFMX = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="TEFFMX"]
  self$params$DVSI = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="DVSI"]
  self$params$IDSL = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="IDSL"]
  self$params$DLO = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="DLO"]
  self$params$DLC = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="DLC"]
  self$params$DVSEND = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="DVSEND"]
  self$params$crop_start_type = start.type
  self$params$crop_end_type = end.type
  
  which_ = grep('TSUM_', f.variety$PARAMETER_CODE)
  self$params$TSUM[1:length(which_),1] = f.variety$PARAMETER_XVALUE[which_]
  self$params$TSUM[1:length(which_),2] = f.variety$PARAMETER_YVALUE[which_]
  self$params$TSUM =  self$params$TSUM[1:length(which_),]
  
  which_ = grep('TBASE_', f.variety$PARAMETER_CODE)
  self$params$TBASE[1:length(which_),1] = f.variety$PARAMETER_XVALUE[which_]
  self$params$TBASE[1:length(which_),2] = f.variety$PARAMETER_YVALUE[which_]
  self$params$TBASE =  self$params$TBASE[1:length(which_),]
  
  which_ = grep('DVS_', f.variety$PARAMETER_CODE)
  self$params$DVS_PERIOD[1:length(which_),1] = f.variety$PARAMETER_XVALUE[which_]
  self$params$DVS_PERIOD[1:length(which_),2] = f.variety$PARAMETER_YVALUE[which_]
  self$params$DVS_PERIOD = self$params$DVS_PERIOD[1:length(which_),]
  
  # Define initial states
  states.curr = get_initial_stage(self,day)
  
  self$states$DOS = states.curr$DOS # Day of sowing
  self$states$DOE = states.curr$DOE # Day of emergence
  self$states$DVS = 0
  self$states$TSUM = 0
  self$states$STAGE = states.curr$STAGE
  
  self$rates$DTSUM = 0
  self$rates$DTSUME = 0
  self$rates$DVR = 0
  self$rates$DVRED = 1

  
  if(self$params$IDSL >= 2) {
    self$params$VERNDVS = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="VERNDVS"]
    self$params$VERNBASE = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="VERNBASE"]
    self$params$VERNSAT = f.variety$PARAMETER_XVALUE[f.variety$PARAMETER_CODE=="VERNSAT"]
    
    self$rates$VERNR = 0
    self$rates$VERNFAC = 0
    self$states$VERN = 0
    self$states$DOV = NA
    self$states$ISVERNALISED = FALSE
  }

  return(self)   
}

calc_rates = function(self, day, drv)
{
  # Day length sensitivity
  DVRED = 1
  if(self$params$IDSL >= 1) {
    DAYLP = as.numeric(daylength.1(day, drv$LAT))
    DVRED = (DAYLP - self$params$DLC)/(self$params$DLO - self$params$DLC)
    if(DVRED < 0)
      DVRED = 0
    if(DVRED > 1)
      DVRED = 1
  }
  
  # Vernalisation
  VERNFAC = 1.
  if(self$params$IDSL >= 2) {
    if(self$states$STAGE == 1) {
      self = calc.rates.vernalisation(self, day, drv)
      VERNFAC = self$rates$VERNFAC
    }
  }
  # Development rates
  if(self$states$STAGE == 0) {
    self$rates$DTSUM = dtsmtb(drv$TEMP, self)
    self$rates$DVR = self$rates$DTSUM/self$params$TSUM[self$params$TSUM[,1]==floor(self$states$DVS),2]
  }
  if(self$states$STAGE == 1) {
    self$rates$DTSUM = dtsmtb(drv$TEMP, self) * VERNFAC * DVRED
    self$rates$DVR = self$rates$DTSUM/self$params$TSUM[self$params$TSUM[,1]==floor(self$states$DVS),2]
  }
  if(self$states$STAGE == 2) {
    self$rates$DTSUM = dtsmtb(drv$TEMP, self)
    self$rates$DVR = self$rates$DTSUM/self$params$TSUM[self$params$TSUM[,1]==floor(self$states$DVS),2]
  }
  
  if(is.na(self$states$STAGE))
    error = 1
  return(self)
}

integrate.dvs = function(self, day, delt=1.0)
{
  error = 0
  #Updates the state variable and checks for phenologic stages
  # Integrate vernalisation module
  
  if(self$params$IDSL >= 2) {
    if(self$states$STAGE == 1) {
      self = integrate.vernalisation(self, day, delt)
    }
  }

  # Integrate phenologic states
  self$states$DVS = self$states$DVS + self$rates$DVR
  self$states$TSUM = self$states$TSUM + self$rates$DTSUM

  # Check if a new stage is reached
  self = next_stage(self)
   
  if(is.na(self$states$STAGE)) {
    msg = c("No STAGE defined in phenology submodule")
    error = 1
  }
  return(self)
}
  
get_initial_stage = function(self, day)
{
  error = 0
  
  if(self$params$crop_start_type == "emergence") {
    STAGE = 1
  }
  if(self$params$crop_start_type == "sowing") {
    STAGE = 0
  }
  if(self$params$crop_start_type != "emergence" || self$params$crop_start_type != "sowing")
    error = 1
  
  return(list(STAGE = STAGE))
}

next_stage = function(self) 
{
  curr.STAGE = self$states$STAGE

  if(curr.STAGE < 2) {
    self$states$STAGE = self$params$DVS_PERIOD[self$params$DVS_PERIOD[,1]==floor(self$states$DVS),2]
  }
  if(curr.STAGE == 2) {
    msg = "Cannot move to next phenology stage: maturity already reached!"
    error = 1
  }
  
  return(self)
}

calc.rates.vernalisation = function(self, day, drv)
{
  DVS = self$states$DVS
  
  if(!self$states$ISVERNALISED) {
    if(DVS < self$params$VERNDVS) {
      self$rates$VERNR = vernrtb(drv$TEMP)
      r = (self$states$VERN - self$params$VERNBASE)/(self$params$VERNSAT-self$params$VERNBASE)
      if(r < 0)
        r = 0
      if(r > 1)
        r = 1
      self$rates$VERNFAC = r
    } else {
      self$rates$VERNR = 0.
      self$rates$VERNFAC = 1.0
      self$force_vernalisation = TRUE
    }
  } else {
    self$rates$VERNR = 0.
    self$rates$VERNFAC = 1.0
  }
  return(self)
}

integrate.vernalisation = function(self, day, delt=1.0)
{

  self$states$VERN = self$states$VERN + self$rates$VERNR
  
  if(self$states$VERN >= self$params$VERNSAT) {
     # Vernalisation requirements reached
     self$states$ISVERNALISED = TRUE
     if(is.na(self$states$DOV)) {
       self$states$DOV = day
     }
  }
  if(self$force_vernalisation && self$states$VERN < self$params$VERNSAT) {  # Critical DVS for vernalisation reached
    # Force vernalisation, but do not set DOV
    self$states$ISVERNALISED = TRUE
    # Write log message to warn about forced vernalisation
  }
  if(self$states$VERN < self$params$VERNSAT && !self$force_vernalisation)  # Reduction factor for phenologic development
    self$states$ISVERNALISED = FALSE
  
  return(self)
}

vernrtb = function(temp)
{
  params.v = matrix(data=c(-30,-8,-4,3, 10,17,20, 50,0,0,0,1,1,0,0,0), nrow=8, ncol=2)
  return(approx(params.v[,1], params.v[,2], temp)$y)
}

dtsmtb = function(temp, self) {
  tbase = ftbase(self)
  params.v = matrix(data=c(-50,tbase,30,45, 50, 0,0,30-tbase,30-tbase,30-tbase), nrow=5,ncol=2)
  return(approx(params.v[,1], params.v[,2], temp)$y)
}

ftsum = function(self) {
  return(self$params$TSUM[self$params$TSUM[,1]==floor(self$states$DVS),2])
}

ftbase = function(self) {
  return(self$params$TBASE[self$params$TBASE[,1]==floor(self$states$DVS),2])
}

daylength.1 = function(day, latitude, angle=-4)
{

  # Calculate day-of-year from date object day
  IDAY = as.numeric(format(day, "%j"))
    
  # constants
  RAD = 0.0174533
  PI = 3.1415926
    

    
   # calculate daylength
    ANGLE = angle
    LAT = latitude
    DEC = -asin(sin(23.45*RAD)*cos(2.*PI*((IDAY)+10.)/365.))
    SINLD = sin(RAD*LAT)*sin(DEC)
    COSLD = cos(RAD*LAT)*cos(DEC)
    AOB   = (-sin(ANGLE*RAD)+SINLD)/COSLD
    
    # daylength
    if(abs(AOB) <= 1.0) 
      DAYLP = 12.0*(1.+2.*asin((-sin(ANGLE*RAD)+SINLD)/COSLD)/PI)
    if(AOB > 1.0)
      DAYLP = 24.0
    if(abs(AOB) > 1.0 && AOB < 1.0)
      DAYLP =  0.0
    
    # store results in cache
    
    return(DAYLP)
    
}

