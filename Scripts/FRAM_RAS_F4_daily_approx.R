## daily interpolation function

# input: date vector and abundance matrix
# output: date vector with daily time interval and approximated abundance matrix + season info
# function contain a fix=T parameter which fixed the data to exactly 4 years (4*365)
# 29th of February will be removed 
daily_approx=function(dt,ee,fix=T){
  
  # daily approximation
  dtw=approx(x=dt,y=ee[1,],xout=seq(dt[1],dt[length(dt)],"day"))$x
  
  ee=apply(ee,1,function(x){
    yy=approx(x=dt,y=x,xout=seq(dt[1],dt[length(dt)],"day"))
    return(yy$y)
  })
  ee=t(ee)
  #remove february 29
  f29=grep("02-29",dtw)
  dtw=dtw[-f29]
  ee=ee[,-f29]
  #start
  #i=which(dtw==sort(dt)[1])
  #end
  #j=which(dtw==sort(dt)[length(dt)])
  #ee=ee[,i:j]
  #dtw=dtw[i:j]
  if(fix==TRUE){ 
    # fix to exactly 4 year from first day
    lastDay=dtw[1]+(4*365)
    li=which(dtw==lastDay)
    dtw=dtw[1:li]
    ee=ee[,1:li]
}
  
  # season
  # spring 01.03. bis 31.05.
  # summer 01.06. bis 31.08.
  # fall 01.09. bis 30.11.
  # winter 01.12. bis 28./29.02.
  mn=as.character(dtw,format="%m")
  season=rep("winter",length(mn))
  season[grep('03|04|05',mn)]="spring"
  season[grep('06|07|08',mn)]="summer"
  season[grep('09|10|11',mn)]="fall"
  
  return(list(dtw=dtw,ee=ee,season=season))
}
