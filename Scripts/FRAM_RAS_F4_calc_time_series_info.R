# function to extract time series components (seasonality and trend) using FFT
# creates descriptive model data based on the extracted info


# input: vector containing abundance information of an organism over time 
calcTimeSeriesInfo=function(x){
  # calculate phase using complex argument information (angle between Re and Im)
  # calculate the difference to highest value in ref sinus (pi*0.5)
  # calculate the phase shift as the proportion of pi
  
  getPhase=function(z3,f=0){
    
    z3=fft(approx(z3,n=1e+5)$y)
    z3=z3[2:round(length(z3)/2)]
    if(f==0){
      f=which.max(abs(z3))
    }
    f=f
    #phase_shift x*pi
    ps=NaN
    ps=((Arg(z3)[f]+(pi/2))/pi)
    # if shift is higher than 0.5*pi than proportion is -2
    if((Arg(z3)[f]+(pi/2))/pi>1){
      ps=((Arg(z3)[f]+(pi/2))/pi)-2
    }
    return(c(phase=ps,freq=f))
  }
  
  
  #time set -> 2*pi == one year == 365 days
  # consider difference to full year
  #rounded number of years
  nyear=round(length(x)/365)
  #difference to it
  dtfy=1-(round(length(x)/365)-(length(x)/365))
  
  lastTP=(2*pi)
  ts=seq(0,lastTP,length.out = length(x))
  ts=ts[1:length(x)]
  #amplitude
  Amp=(max(x)-min(x))/2
  #yshift
  ys=max(x)-Amp
  phase_freq=getPhase(x)
  #### identify trend using linear regression
  fit=lm(x~ts)
  
  z=Amp*sin(phase_freq[2]*ts+phase_freq[1]*pi)+ys
  ## add trend
  zt=z+fit$fitted.values-fit$coefficients[1]

  
  ### extract max (hp) / min (tp) points 
  hp=((pi*0.5)/phase_freq[2]) - (phase_freq[1]*pi/phase_freq[2])
  tp=hp-pi/phase_freq[2]
  all_hp=c(hp,(hp+seq(phase_freq[2]+1)*(2*pi/phase_freq[2])))
  all_hp=all_hp[which(all_hp >=0 & all_hp <= 2*pi)]
  
  all_tp=c(tp,(tp+seq(phase_freq[2]+1)*(2*pi/phase_freq[2])))
  #all_tp=all_tp[which(all_tp >=0 & all_tp <= 2*pi)]
  all_tp=all_tp[which(all_tp <= 2*pi)]
  
  ### convert to x indices
  xtp=round(all_tp/(2*pi)*length(x))
  xhp=round(all_hp/(2*pi)*length(x))
  
  hp_tp=list(coef_hp=all_hp,coef_tp=all_tp,xhp=xhp,xtp=xtp)
  
  ## store all time series components in a list
  ts_comp=list(Amp=Amp, phase_freq=phase_freq, ys=ys,fit=fit,modData_sin=z,mod_dataTS=zt,hp_tp=hp_tp)
  
  return(ts_comp)  
}
