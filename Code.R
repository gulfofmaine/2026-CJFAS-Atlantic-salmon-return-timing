
# loading packages
library(tidyverse)
library(lubridate)
library(mgcv)
library(itsadug)
library(MuMIn)
library(here)
library(gridExtra)
library(patchwork)
library(cowplot)
library(stringi)
library(knitr)
library(sf)
library(trend)
library(broom)
library(scales)
library(segmented)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggpattern)
library(zoo)
library(ggpubr)
library(conflicted)
conflict_prefer("filter","dplyr")
conflict_prefer("select","dplyr")
conflict_prefer("here","here")
theme_set(theme_bw())
options(na.action="na.fail")
options(knitr.kable.NA='')

# OISST data

# pre-processed 15-day smoothed OISST NW Atlantic lat 40,70 lon -75,-35
OISST<-readRDS("NWAtl19822020_smooth.Rds")
OISST$temp<-as.numeric(format(OISST$temp,digits=3))
colnames(OISST)[4]<-"Date"
colnames(OISST)[5]<-"Temp"

# winter phenology metrics
# seasonal year (Aug01 year i to July31 year i+1)
wdat<-OISST
wdat$seasonyear<-ifelse(month(wdat$Date)%in%c(1:7),year(wdat$Date),year(wdat$Date)+1)

# not using winter 1981-1982 or 2019-2020
wdat%>%filter(seasonyear%in%c(1983:2019))->wdat

# for every cell, identify the winter min temperature for every winter,
# then find the highest of the winter minima and add 0.5 = threshold
wdat%>%
  group_by(cell,seasonyear)%>%
  summarise(mntemp=min(Temp,na.rm=TRUE))%>%
  group_by(cell)%>%
  summarise(highwinmin=max(mntemp))%>%
  mutate(winthrsh=highwinmin+0.5)%>%
  select(-highwinmin)->winmin

# join with all data
left_join(wdat,winmin,by="cell")%>%
  mutate(isBelow=Temp<winthrsh)%>%
  mutate(year=year(Date),month=month(Date),yearday=yday(Date))->win

# filter out days that are below threshold and by yearday, 
# identify first, last, and total days
# total days not consecutive days but accumulated days
win%>%
  filter(isBelow==TRUE)%>%
  group_by(cell,seasonyear)%>%
  summarise(winsrt=first(yearday),winend=last(yearday),winlen=n())->wphenol

# if winter starts in prior calendar year, adjust start dates as negative
# days relative to jan 1st, so dec 31st would be day 0, dec 30th day -1 etc...
wphenol$isprior<-ifelse(wphenol$winsrt>=213,"prior",NA)

adjwinsrt<-c()
for(i in 1:nrow(wphenol)){
  z<-ifelse(!is.na(wphenol$isprior[i]),wphenol$winsrt[i]-366,wphenol$winsrt[i])
  adjwinsrt<-c(adjwinsrt,z)
}
wphenol$adjwinsrt<-adjwinsrt

# winter end
# if winter end in prior calendar year, adjust end dates as negative
wphenol$ispriorend<-ifelse(wphenol$winend>=213,"prior",NA)

# winter end, adjust end dates as negative
adjwinend<-c()
for(i in 1:nrow(wphenol)){
  z<-ifelse(!is.na(wphenol$ispriorend[i]),wphenol$winend[i]-366,wphenol$winend[i])
  adjwinend<-c(adjwinend,z)
}
wphenol$adjwinend<-adjwinend

# winter season length difference method
# difference method
wphenol%>%
  mutate(adjwinlen=adjwinend-adjwinsrt)->wphenol

wphenol$adjwinlen<-ifelse(wphenol$adjwinlen<=0,wphenol$winlen,wphenol$adjwinlen)

#saveRDS(wphenol,"wphenol_diffmeth.RDS")

# summer metrics
# for every cell, identify the summer max temperature for every summer,
# then find the lowest of the summer maxima and subtract 0.5 = threshold
OISST%>%
  group_by(cell,year=year(Date))%>%
  summarise(mxtemp=max(Temp,na.rm=TRUE))%>%
  group_by(cell)%>%
  summarise(mntemp=min(mxtemp))%>%
  mutate(summthrsh=mntemp-0.5)%>%
  select(-mntemp)->summax

# join with all data
left_join(OISST,summax,by="cell")%>%
  mutate(isAbove=Temp>summthrsh)%>%
  mutate(year=year(Date),month=month(Date),yearday=yday(Date))->summ

# filter out days that are above threshold and after day 100, 
# identify first, last, and total days
# total days not consecutive days but accumulated days
summ%>%
  filter(isAbove==TRUE)%>%
  filter(yearday%in%c(100:360))%>% # look for summer after day 100 and prior to day 360
  group_by(cell,year)%>%
  summarise(sumsrt=first(yearday),sumend=last(yearday),sumlen=n())->sphenol

# summer season length difference method
# difference method
sphenol%>%
  mutate(adjsumlen=sumend-sumsrt)->sphenol

#saveRDS(sphenol,"sphenol_diffmeth.RDS")

# relate cells to lon lat

# read in lon lat
lonlat<-OISST%>%select(lon,lat,cell)%>%distinct()

# relate to cell number
geomdat<-left_join(wphenol,lonlat,by="cell")
sgeomdat<-left_join(sphenol,lonlat,by="cell")

# write to csv
#write_csv(geomdat,"winterphenologyNWAtlantic19832019_diffmeth.csv")
#write_csv(sgeomdat,"summerphenologyNWAtlantic19822019_diffmeth.csv")

sphenol<-read_csv("summerphenologyNWAtlantic19822019_diffmeth.csv")
wphenol<-read_csv("winterphenologyNWAtlantic19832019_diffmeth.csv")

colnames(wphenol)[2]<-"year"
phenol<-left_join(wphenol,sphenol)

phenol<-st_as_sf(phenol,coords=c("lon","lat"),remove=FALSE)
st_crs(phenol)<-4326

# read in EPU shapefiles
reg2<-st_read("EPUs.shp")
reg2<-st_transform(reg2,crs=st_crs(4326))

# combined southern regions
reg2%>%filter(Name=="Scotian Shelf")->SS
SS<-st_intersection(phenol,SS)

# combined northern regions
reg2%>%filter(Name=="Northeast Newfoundland Shelf")->NN
NN<-st_intersection(phenol,NN)

cuezones<-bind_rows(NN,SS)

# two regions
# Southern: SNf, SS, GoM
# Northern: GB, FC, NNfS
regions<-tibble(region=c("southern","northern"),
                Name=c("Scotian Shelf","Northeast Newfoundland Shelf"))
left_join(cuezones,regions,by="Name")->cuezones

# calc averages
cuezones%>%
  group_by(region,year)%>%
  summarize(mnwinlen=mean(winlen),
            mnwinsrt=mean(adjwinsrt),
            mnwinend=mean(adjwinend),
            mnsumlen=mean(sumlen),
            mnsumsrt=mean(sumsrt),
            mnsumend=mean(sumend),
            mnadjwinlen=mean(adjwinlen),
            mnadjsumlen=mean(adjsumlen))->avg



# River data
# USGS 01034500 Penobscot River at West Enfield, Maine
Q<-read.table(
  "https://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=01034500&parameterCd=00060&startDT=1902-10-01&endDT=2019-12-31",
  sep="\t",
  header=TRUE,
  comment.char="#",
  skip=1,
  fill=TRUE,
  stringsAsFactors=FALSE)

Q<-Q[-1,]
Q$datetime<-as.Date(Q$datetime)
Q<-Q%>%select(datetime,X63753_00060_00003)
colnames(Q)<-c("date","flow")
Q$flow<-as.numeric(Q$flow)
Q%>%mutate(year=year(date),month=month(date),yearday=yday(date))->Q

# apply 15-day smooth
Q$flow<-rollmean(x=Q$flow,k=15,fill="extend")

# years with run timing data
Q%>%filter(year%in%c(1983:2019))->Q8319

# Winter-spring center of volume
# winter spring (january through may) center of volume
Q8319%>%
  filter(month%in%c(1:5))->wscv

# metrics
wscv%>%
  group_by(year)%>%
  mutate(cumulativeflow=cumsum(flow),
         maxflow=max(cumulativeflow),
         scaledflow=cumulativeflow/maxflow,
         meanyearflow=mean(flow))->wscv

# calc props.
wscv%>%
  group_by(year)%>%
  summarise(Q05=yearday[scaledflow>=0.05][1],
            Q10=yearday[scaledflow>=0.10][1],
            Q15=yearday[scaledflow>=0.15][1],
            Q20=yearday[scaledflow>=0.20][1],
            Q25=yearday[scaledflow>=0.25][1],
            Q30=yearday[scaledflow>=0.30][1],
            Q35=yearday[scaledflow>=0.35][1],
            Q40=yearday[scaledflow>=0.40][1],
            Q45=yearday[scaledflow>=0.45][1],
            Q50=yearday[scaledflow>=0.50][1],
            Q55=yearday[scaledflow>=0.55][1],
            Q60=yearday[scaledflow>=0.60][1],
            Q65=yearday[scaledflow>=0.65][1],
            Q70=yearday[scaledflow>=0.70][1],
            Q75=yearday[scaledflow>=0.75][1],
            Q80=yearday[scaledflow>=0.80][1],
            Q85=yearday[scaledflow>=0.85][1],
            Q90=yearday[scaledflow>=0.90][1],
            Q95=yearday[scaledflow>=0.95][1])->wscv_dat

# WSCV
wscv_dat%>%
  pivot_longer(cols=2:20,names_to="Qcuvol",values_to="yearday")%>%
  filter(Qcuvol=="Q50")%>%
  select(yearday)%>%
  pull()->wscv_var

# Annual center of volume
Q8319%>%
  group_by(year)%>%
  mutate(cumulativeflow=cumsum(flow),
         maxflow=max(cumulativeflow),
         scaledflow=cumulativeflow/maxflow)->scaletotalflow

# cumulative flow
scaletotalflow%>%
  group_by(year)%>%
  summarise(Q05=yearday[scaledflow>=0.05][1],
            Q10=yearday[scaledflow>=0.10][1],
            Q15=yearday[scaledflow>=0.15][1],
            Q20=yearday[scaledflow>=0.20][1],
            Q25=yearday[scaledflow>=0.25][1],
            Q30=yearday[scaledflow>=0.30][1],
            Q35=yearday[scaledflow>=0.35][1],
            Q40=yearday[scaledflow>=0.40][1],
            Q45=yearday[scaledflow>=0.45][1],
            Q50=yearday[scaledflow>=0.50][1],
            Q55=yearday[scaledflow>=0.55][1],
            Q60=yearday[scaledflow>=0.60][1],
            Q65=yearday[scaledflow>=0.65][1],
            Q70=yearday[scaledflow>=0.70][1],
            Q75=yearday[scaledflow>=0.75][1],
            Q80=yearday[scaledflow>=0.80][1],
            Q85=yearday[scaledflow>=0.85][1],
            Q90=yearday[scaledflow>=0.90][1],
            Q95=yearday[scaledflow>=0.95][1],)->qcuflow

# final var
Q50_var<-qcuflow$Q50

# River summer start phenology metric
# Run Timing data from Milford
dat<-read_csv("tblTendData2.csv")
dat<-dat%>%select(c(TendDate,SiteCode,WaterTemp))
colnames(dat)<-c("date","site","temperature")
dat$date<-as.Date(dat$date,format="%m/%d/%Y")
dat<-arrange(dat,date)
# filter out 1STILLW0.20
dat%>%filter(site!="1STILLW0.20")->dat
# daily mean
dat<-dat%>%
  group_by(date)%>%
  summarise(temperature=mean(temperature,na.rm=TRUE))

# NA in data
dat$temperature[which(dat$temperature=="NaN")]<-NA
dat<-dat%>%filter(!is.na(temperature))

# USGS temperature data
temp<-read_csv("PenobscotEddington.csv")
colnames(temp)<-c("date","temperature")

# date sequence
date<-tibble(date=seq.Date(as.Date("1978-01-01"),as.Date("2019-12-31"),"day"))

# merge run timing temperature to date sequence
riversites<-left_join(date,dat,by="date")
riversites$type<-NA
riversites$type[which(!is.na(riversites$temperature))]<-"Trap"

# format
temp$date<-as.Date(temp$date)

# merge Eddington temperature to riversites
riversites<-left_join(riversites,temp,by="date")
riversites%>%mutate(tmp=coalesce(temperature.x,temperature.y))->riversites
riversites$type[which(is.na(riversites$type)&!is.na(riversites$tmp))]<-"Gage"

# year, month, yearday
riversites<-riversites%>%mutate(year=year(date),month=month(date),yearday=yday(date))

# seasonal groupings
riversites$season<-ifelse(riversites$month%in%c(12,1,2),"djf",
                          ifelse(riversites$month%in%c(3:5),"mam",
                                 ifelse(riversites$month%in%c(6:8),"jja","son")))
riversites$season<-factor(riversites$season,levels=c("djf","mam","jja","son"))

# fill in remaining gaps with air temp model
airmod<-read_csv("Eddington_daily_obs.csv")

# kelvin to celsius, if negative make 0
airmod$daily_avg<-airmod$daily_avg-273.15
airmod$daily_avg<-ifelse(airmod$daily_avg<0,0,airmod$daily_avg)
colnames(airmod)<-c("daily_avg","date")

# merge airmod temperature to riversites
riversites<-left_join(riversites,airmod,by="date")
riversites%>%mutate(tmp=coalesce(tmp,daily_avg))->riversites
riversites$type[which(is.na(riversites$type)&!is.na(riversites$tmp))]<-"Air model"

# essential columns
riversites%>%select(date,tmp,year,yearday)->riverTdaily

# apply rolling 15-day smooth
riverTdaily$tmp15<-rollmean(x=riverTdaily$tmp,k=15,fill=NA)

# filter to 1983-2019 and find max temp by year
riverTdaily%>%
  filter(year%in%c(1983:2019))%>%
  group_by(year)%>%
  summarise(mxtmp=max(tmp15,na.rm=TRUE))->mxtmp

# lowest summer max
mxtmp%>%filter(mxtmp==min(mxtmp))->minmax # 1996 21.8

# threshold is minmax-0.5
threshold<-minmax$mxtmp-0.5

# calculate start, end, and length
riverTdaily%>%
  filter(year%in%c(1983:2019))%>%
  filter(tmp15>=threshold)%>%
  group_by(year)%>%
  summarise(rivsumsrt=first(yearday),
            rivsumend=last(yearday),
            rivsumlen=n())->riverphenol

# final var
rivsumsrt_var<-riverphenol$rivsumsrt

# Spring inflection point
# breakpoint est loop years
brkest<-c()
for(i in 1:37){
  # year
  yr<-c(1983:2019)[i]
  # jan through may
  rng<-c(1:151)
  riverTdaily%>%
    filter(year==yr&yearday%in%rng)%>%
    pull(tmp15)->ts
  mod<-lm(ts~rng)
  o<-segmented(mod)
  myf<-fitted(o)
  # breakpoint yearday
  brkest<-c(brkest,round(o$psi[2]))
}

# final var
sprinfl<-brkest

# Run timing data
runtiming<-read_csv("tblSalmonData.csv")

# use ID to extract year-month-day
dt<-c()
drcd<-c()
loc<-c()
for(i in 1:nrow(runtiming)){
  # split apart ID by "-" symbol
  q<-unlist(strsplit(runtiming$JoinID_tag[i],"-"))
  # store drainagecode (first element in q vector)
  drcd<-c(drcd,q[1])
  # store location (second element in q vector)
  loc<-c(loc,q[2])
  # find the 8 character date, insert "-" between yyyy-mm-dd
  if(nchar(q[5])==8){
    z<-q[5]
    stri_sub(z,5,4)<-"-"
    stri_sub(z,8,7)<-"-"
    dt<-c(dt,z)
    next
  }
  # find the 8 character date, insert "-" between yyyy-mm-dd
  if(length(q)>5){
    if(nchar(q[6])==8){
      z<-q[6]
      stri_sub(z,5,4)<-"-"
      stri_sub(z,8,7)<-"-"
      dt<-c(dt,z)
      next
    }
  }
  # find the 8 character date, insert "-" between yyyy-mm-dd
  if(nchar(q[4])==8){
    z<-q[4]
    stri_sub(z,5,4)<-"-"
    stri_sub(z,8,7)<-"-"
    dt<-c(dt,z)
    next
  }
  # if already separated by "-", find year alone
  q<-paste(q[which(q%in%as.character(1978:2017)):(length(q)-1)],collapse="-")
  # remove non-date info at end of ID
  if(length(unlist(strsplit(q,"-")))==4){q<-paste(unlist(strsplit(q,"-"))[-4],collapse="-")}
  # store as character vector
  dt<-c(dt,q)
}

# add DrainageCode column before ID column
runtiming<-add_column(runtiming,DrainageCode=drcd,.before="JoinID_tag")

# add TrapLocation column before ID column
runtiming<-add_column(runtiming,TrapLocation=loc,.before="JoinID_tag")

# add Date column after ID column
runtiming<-add_column(runtiming,Date=dt,.after="JoinID_tag")

# convert to Date format
dt<-as.Date(dt)

# add year column after Date column
runtiming<-add_column(runtiming,year=year(dt),.after="Date")

# add month column after year column
runtiming<-add_column(runtiming,month=month(dt),.after="year")

# add day column after month column
runtiming<-add_column(runtiming,day=day(dt),.after="month")

# add yearday column after day column
runtiming<-add_column(runtiming,yearday=yday(dt),.after="day")

# Filter only DrainageCode == "PN"
PNruntiming<-runtiming %>% filter(DrainageCode=="PN")

# Filter only IsRecap == 0
PNruntiming<-PNruntiming %>% filter(IsRecap=="No")

# Filter TrapLocation == "1MAINST47.50 or "1MAINST62.28"
PNruntiming<-PNruntiming %>% filter(TrapLocation%in%c("1MAINST47.50","1MAINST62.28"))

# Look at Origin types, Keep H only
PNruntiming<-PNruntiming %>% filter(Origin=="H")

# Look at SeaAge, simplify to 1,2,3 as new lowercase seaage column
# everything that isn't SeaAge==1 rename to multi sea winter (MSW)
PNruntiming<-add_column(PNruntiming,seaage=ifelse(PNruntiming$SeaAge==1,"1SW","MSW"))
PNruntiming<-PNruntiming %>% filter(seaage%in%c("1SW","MSW"))

# make Date object
PNruntiming$Date<-as.Date(PNruntiming$Date)

# daily counts
PNruntiming%>%
  filter(seaage=="MSW")%>%
  group_by(Date)%>%
  count()->x

# full date time series
Date<-seq.Date(from=as.Date("1978-01-01"),to=as.Date("2019-12-31"),by="day")
fulldate<-tibble(Date)
fulldate<-left_join(fulldate,x)
fulldate$n<-ifelse(is.na(fulldate$n),0,fulldate$n)

# 5-day smooth
x5day<-rollapply(data=fulldate$n,width=5,FUN=mean,partial=TRUE,align="center")
fulldate$MSW5day<-x5day

# 10-day smooth
x10day<-rollapply(data=fulldate$n,width=10,FUN=mean,partial=TRUE,align="center")
fulldate$MSW10day<-x10day

# calculate run metrics
colnames(fulldate)[2]<-"MSWraw"

# cumulative returns columns, 5%, 25%, 50%, 75%, and 95% total run, first days above 5, 25, 50, 75, and 95
fulldate%>%
  group_by(year=year(Date))%>%
  mutate(runtotraw=cumsum(MSWraw),runtot5day=cumsum(MSW5day),runtot10day=cumsum(MSW10day))%>%
  mutate(run05raw=max(runtotraw)*0.05,run055day=max(runtot5day)*0.05,run0510day=max(runtot10day)*0.05,
         run25raw=max(runtotraw)*0.25,run255day=max(runtot5day)*0.25,run2510day=max(runtot10day)*0.25,
         run50raw=max(runtotraw)*0.50,run505day=max(runtot5day)*0.50,run5010day=max(runtot10day)*0.50,
         run75raw=max(runtotraw)*0.75,run755day=max(runtot5day)*0.75,run7510day=max(runtot10day)*0.75,
         run95raw=max(runtotraw)*0.95,run955day=max(runtot5day)*0.95,run9510day=max(runtot10day)*0.95)%>%
  summarise(Q05raw=yday(Date[runtotraw>=run05raw][1]),
            Q055day=yday(Date[runtot5day>=run055day][1]),
            Q0510day=yday(Date[runtot10day>=run0510day][1]),
            Q25raw=yday(Date[runtotraw>=run25raw][1]),
            Q255day=yday(Date[runtot5day>=run255day][1]),
            Q2510day=yday(Date[runtot10day>=run2510day][1]),      
            Q50raw=yday(Date[runtotraw>=run50raw][1]),
            Q505day=yday(Date[runtot5day>=run505day][1]),
            Q5010day=yday(Date[runtot10day>=run5010day][1]),
            Q75raw=yday(Date[runtotraw>=run75raw][1]),
            Q755day=yday(Date[runtot5day>=run755day][1]),
            Q7510day=yday(Date[runtot10day>=run7510day][1]),
            Q95raw=yday(Date[runtotraw>=run95raw][1]),
            Q955day=yday(Date[runtot5day>=run955day][1]),
            Q9510day=yday(Date[runtot10day>=run9510day][1]))->MSWrunmetrics

# final
MSWrunmetrics$seaage<-"MSW"

# repeat for 1SW
# daily counts
PNruntiming%>%
  filter(seaage=="1SW")%>%
  group_by(Date)%>%
  count()->x

# full date time series
Date<-seq.Date(from=as.Date("1978-01-01"),to=as.Date("2019-12-31"),by="day")
fulldate<-tibble(Date)
fulldate<-left_join(fulldate,x)
fulldate$n<-ifelse(is.na(fulldate$n),0,fulldate$n)

# 5-day smooth
x5day<-rollapply(data=fulldate$n,width=5,FUN=mean,partial=TRUE,align="center")
fulldate$oneSW5day<-x5day

# 10-day smooth
x10day<-rollapply(data=fulldate$n,width=10,FUN=mean,partial=TRUE,align="center")
fulldate$oneSW10day<-x10day

# calculate run metrics
colnames(fulldate)[2]<-"oneSWraw"

# cumulative returns columns, 5%, 25%, 50%, 75%, and 95% total run, first days above 5, 25, 50, 75, and 95
fulldate%>%
  group_by(year=year(Date))%>%
  mutate(runtotraw=cumsum(oneSWraw),runtot5day=cumsum(oneSW5day),runtot10day=cumsum(oneSW10day))%>%
  mutate(run05raw=max(runtotraw)*0.05,run055day=max(runtot5day)*0.05,run0510day=max(runtot10day)*0.05,
         run25raw=max(runtotraw)*0.25,run255day=max(runtot5day)*0.25,run2510day=max(runtot10day)*0.25,
         run50raw=max(runtotraw)*0.50,run505day=max(runtot5day)*0.50,run5010day=max(runtot10day)*0.50,
         run75raw=max(runtotraw)*0.75,run755day=max(runtot5day)*0.75,run7510day=max(runtot10day)*0.75,
         run95raw=max(runtotraw)*0.95,run955day=max(runtot5day)*0.95,run9510day=max(runtot10day)*0.95)%>%
  summarise(Q05raw=yday(Date[runtotraw>=run05raw][1]),
            Q055day=yday(Date[runtot5day>=run055day][1]),
            Q0510day=yday(Date[runtot10day>=run0510day][1]),
            Q25raw=yday(Date[runtotraw>=run25raw][1]),
            Q255day=yday(Date[runtot5day>=run255day][1]),
            Q2510day=yday(Date[runtot10day>=run2510day][1]),      
            Q50raw=yday(Date[runtotraw>=run50raw][1]),
            Q505day=yday(Date[runtot5day>=run505day][1]),
            Q5010day=yday(Date[runtot10day>=run5010day][1]),
            Q75raw=yday(Date[runtotraw>=run75raw][1]),
            Q755day=yday(Date[runtot5day>=run755day][1]),
            Q7510day=yday(Date[runtot10day>=run7510day][1]),
            Q95raw=yday(Date[runtotraw>=run95raw][1]),
            Q955day=yday(Date[runtot5day>=run955day][1]),
            Q9510day=yday(Date[runtot10day>=run9510day][1]))->oneSWrunmetrics

# final
oneSWrunmetrics$seaage<-"1SW"

# bind seaages together
runmetrics<-bind_rows(MSWrunmetrics,oneSWrunmetrics)

# Combined OISST, river, and run timing data
alldat<-read_csv("Data_for_models.csv")

# use for models
subdat<-data.frame(alldat)

# variable names
codez<-tibble(cleanname=c("N.winlen","N.winsrt","N.winend",
                          "S.winlen","S.winsrt","S.winend",
                          "N.sumlen","N.sumsrt","N.sumend",
                          "S.sumlen","S.sumsrt","S.sumend",
                          "N.sumlen_lag",
                          "wscv","sprinfl","rivsumsrt","rivsumend","Q50"),
              modelname=c("s(northern_mnwinlen, k = 4)","s(northern_mnwinsrt, k = 4)","s(northern_mnwinend, k = 4)",
                          "s(southern_mnwinlen, k = 4)","s(southern_mnwinsrt, k = 4)","s(southern_mnwinend, k = 4)",
                          "s(northern_mnsumlen, k = 4)","s(northern_mnsumsrt, k = 4)","s(northern_mnsumend, k = 4)",
                          "s(southern_mnsumlen, k = 4)","s(southern_mnsumsrt, k = 4)","s(southern_mnsumend, k = 4)",
                          "s(northern_mnsumlen_lag, k = 4)",
                          "s(wscv, k = 4)","s(sprinfl, k = 4)","s(rivsumsrt, k = 4)","s(rivsumend, k = 4)","s(Q50, k = 4)"))

# Figure 1

# shapefiles
reg2<-st_read("EPUs.shp")
reg2<-st_transform(reg2,crs=st_crs(26919))
world<-ne_countries(scale="large",country=c("united states of america","canada","greenland"),returnclass="sf")
world<-st_transform(world,crs=st_crs(26919))

# plot
ggplot()+
  geom_sf(data=world,color="gray30",linewidth=0.7)+
  geom_sf(data=reg2,fill=c("#1b9e77","#e7298a"),color="transparent",alpha=0.1)+
  geom_sf(data=reg2,fill="transparent",color=c("#1b9e77","#e7298a"),linewidth=0.5)+
  geom_sf(data=world,color="gray30",linewidth=0.8)+
  coord_sf(xlim=c(13,250)*10000,ylim=c(45.3,78.5)*100000,expand=FALSE)+
  theme(text=element_text(size=17))->f1

f1

ggsave("Figure 1 basemap.pdf",f1,width=8.5,height=11,units="in",dpi=300)

# Figure 2

# percentile labels
perc.labs<-c("5th","25th","50th","75th","95th")
names(perc.labs)<-c("_Q05raw","_Q25raw","_Q50raw","_Q75raw","_Q95raw")

# plot
alldat%>%
  select(year:`1SW_Q95raw`)%>%
  pivot_longer(cols=MSW_Q05raw:`1SW_Q95raw`,names_to=c("seaage","percentile"),names_sep=3,values_to="yearday")%>%
  ggplot(aes(x=year,y=yearday,group=seaage))+
  stat_summary(fun.min=min,fun.max=max,fun=max,shape="",size=1,alpha=0.5)+
  geom_point(aes(fill=percentile),color="black",shape=21,size=2.5,alpha=1)+
  scale_x_continuous(breaks=seq(1980,2020,5))+
  labs(x="Year",y="Day of Year",color="Return Percentile",fill="Return Percentile")+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'),labels=perc.labs)+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'),labels=perc.labs)+
  facet_wrap(.~seaage,nrow=2)+
  theme(legend.position="bottom",strip.background=element_blank(),panel.grid.minor.x=element_blank(),
        text=element_text(size=15))+
  scale_y_continuous(breaks=c(91,105,121,135,152,166,182,196,213,227,244,258,274,288,305),
                     labels=c("Apr 1 (91)","Apr 15 (105)","May 1 (121)","May 15 (135)",
                              "Jun 1 (152)","Jun 15 (166)","Jul 1 (182)","Jul 15 (196)",
                              "Aug 1 (213)","Aug 15 (227)","Sep 1 (244)","Sep 15 (258)",
                              "Oct 1 (274)","Oct 15 (288)","Nov 1 (305)"))+
  coord_cartesian(ylim=c(121,288))->f2

f2

ggsave("Figure 2.pdf",f2,width=7,height=8,units="in",dpi=300)

# Figure 3

# read in data
wphenol<-readRDS("wphenol_diffmeth.Rds")
colnames(wphenol)[2]<-"year"
sphenol<-readRDS("sphenol_diffmeth.Rds")
phenol<-left_join(wphenol,sphenol)

# clean 
dat<-phenol%>%select(cell,lon,lat,year,adjwinsrt,adjwinend,adjwinlen,sumsrt,sumend,adjsumlen)
rm(phenol,sphenol,wphenol)

# overwintering area average by year
dat%>%
  filter(lon>=c(-59)&lon<=c(-45)&lat>=43&lat<=55)%>%
  group_by(year)%>%
  summarise(winsrt=mean(adjwinsrt),
            winend=mean(adjwinend),
            winlen=mean(adjwinlen),
            sumsrt=mean(sumsrt),
            sumend=mean(sumend),
            sumlen=mean(adjsumlen))->mndat

# Sen's slopes and p-values (take a bit to run)
dat[complete.cases(dat),]->dat
dat%>%
  group_by(cell)%>%
  summarise(winsrt_sens=sens.slope(adjwinsrt)[[1]],winsrt_pval=sens.slope(adjwinsrt)[[3]],
            winend_sens=sens.slope(adjwinend)[[1]],winend_pval=sens.slope(adjwinend)[[3]],
            winlen_sens=sens.slope(adjwinlen)[[1]],winlen_pval=sens.slope(adjwinlen)[[3]],
            sumsrt_sens=sens.slope(sumsrt)[[1]],sumsrt_pval=sens.slope(sumsrt)[[3]],
            sumend_sens=sens.slope(sumend)[[1]],sumend_pval=sens.slope(sumend)[[3]],
            sumlen_sens=sens.slope(adjsumlen)[[1]],sumlen_pval=sens.slope(adjsumlen)[[3]])->sens
sens<-left_join(sens,dat%>%select(cell,lon,lat)%>%distinct(),by="cell")

# rearrange data structure
sens%>%
  pivot_longer(cols=winsrt_sens:sumlen_pval,
               names_to=c("metric","sens"),
               names_pattern="(.*)_(.*)",
               values_to="value")%>%
  pivot_wider(names_from="sens",values_from="value")->sens

# sf projection
world<-ne_countries(scale="medium",returnclass="sf")
st_crs(world)<-4326

# rename metrics
sens$metric[which(sens$metric=="winsrt")]<-"Winter Start Sen's Slope"
sens$metric[which(sens$metric=="winend")]<-"Winter End Sen's Slope"
sens$metric[which(sens$metric=="winlen")]<-"Winter Length Sen's Slope"
sens$metric[which(sens$metric=="sumsrt")]<-"Summer Start Sen's Slope"
sens$metric[which(sens$metric=="sumend")]<-"Summer End Sen's Slope"
sens$metric[which(sens$metric=="sumlen")]<-"Summer Length Sen's Slope"
sens$metric<-factor(sens$metric,levels=c("Winter Start Sen's Slope",
                                         "Winter End Sen's Slope",
                                         "Winter Length Sen's Slope",
                                         "Summer Start Sen's Slope",
                                         "Summer End Sen's Slope",
                                         "Summer Length Sen's Slope"))

# northern and southern region shapefiles
reg2<-st_read("EPUs.shp")
reg2<-st_transform(reg2,crs=st_crs(4326))

sensletz<-tibble(letz=c("A","B","C","D","E","F"),metric=c("Winter Start Sen's Slope",
                                                          "Winter End Sen's Slope",
                                                          "Winter Length Sen's Slope",
                                                          "Summer Start Sen's Slope",
                                                          "Summer End Sen's Slope",
                                                          "Summer Length Sen's Slope"))

sens%>%left_join(sensletz,by="metric")->senslab

senslab$metric<-factor(senslab$metric,levels=c("Winter Start Sen's Slope",
                                               "Winter End Sen's Slope",
                                               "Winter Length Sen's Slope",
                                               "Summer Start Sen's Slope",
                                               "Summer End Sen's Slope",
                                               "Summer Length Sen's Slope"))

# Sens slope plots
ggplot()+
  geom_raster(data=senslab,aes(x=lon,y=lat,fill=sens))+
  scale_fill_gradient2(low="#7570b3",mid="#f7f7f7",high="#d95f02",midpoint=0,breaks=c(-5,0,3))+
  geom_sf(data=reg2,color=rep(c("#1b9e77","#F181B9"),6),fill="transparent",linewidth=0.8)+
  geom_text(data=senslab%>%filter(pval<0.05),aes(x=lon,y=lat,label="*"),
            size=1,color="gray60",nudge_x=0.05,nudge_y=-0.04)+
  geom_sf(data=world,color="gray30",size=0.2)+
  coord_sf(xlim=c(-74.9,-35.1),ylim=c(40.1,69.1),expand=FALSE)+
  scale_x_continuous(breaks=c(-70,-60,-50,-40))+
  facet_wrap(.~metric,ncol=3,nrow=2)+
  theme(legend.position="bottom",text=element_text(size=15),legend.text=element_text(size=15),
        strip.background=element_blank(),axis.title=element_blank())+
  labs(fill="Sen's Slope")->sensplt

# Time series
sensplt+geom_text(data=senslab,aes(x=-Inf,y=Inf,label=letz),hjust=-2.5,vjust=1.5,inherit.aes=FALSE)->sensplt

metric.labs<-c("Winter Length","Summer Length","Winter Start","Summer Start","Winter End","Summer End")
names(metric.labs)<-c("mnwinlen","mnsumlen","mnwinsrt","mnsumsrt","mnwinend","mnsumend")

czletz<-tibble(letz=c("G","H","I","J","K","L"),metric=c("mnwinsrt","mnwinend","mnwinlen","mnsumsrt","mnsumend","mnsumlen"))

subdat%>%
  select(c(year,northern_mnwinlen:southern_mnsumend))%>%
  pivot_longer(-year,names_to=c("region","metric"),names_sep="_")%>%
  left_join(czletz,by="metric")%>%
  mutate(metric=factor(metric,levels=c("mnwinsrt","mnwinend","mnwinlen","mnsumsrt","mnsumend","mnsumlen")))%>%
  mutate(region=str_replace(region,"^\\w{1}",toupper))->cuezonephenolletz

# G, H, J, K
cuezonephenolletz%>%
  filter(letz%in%c("G","H","J","K"))%>%
  ggplot(aes(x=year,y=value,group=region,color=region))+
  geom_line()+
  stat_smooth(geom="line",method="gam",se=FALSE,linewidth=1.5,alpha=0.8)+
  geom_point(size=1.5,alpha=0.8)+
  scale_x_continuous(breaks=c(1990,2000,2010,2020),minor_breaks=c(1985,1995,2005,2015))+
  guides(x=guide_axis(minor.ticks = TRUE))+
  scale_color_manual(values=rep(c("#F181B9","#1b9e77"),4))+
  labs(x="Year",y="Yearday",color="Region")+
  facet_wrap(.~metric,nrow=2,scales="free_y",labeller=labeller(metric=metric.labs))+
  theme(strip.background=element_blank(),legend.position="none",text=element_text(size=15),panel.grid.minor.x=element_blank(),legend.text=element_text(size=15),
        ,plot.margin=unit(c(5.5,12,5.5,5.5),"points"),
        axis.title.x=element_text(hjust=0.81))+
  geom_text(data=cuezonephenolletz%>%filter(letz%in%c("G","H","J","K")),aes(x=-Inf,y=Inf,label=letz),hjust=-1,vjust=1.5,inherit.aes=FALSE)->ghjk

# I, L
cuezonephenolletz%>%
  filter(letz%in%c("I","L"))%>%
  ggplot(aes(x=year,y=value,group=region,color=region))+
  geom_line()+
  stat_smooth(geom="line",method="gam",se=FALSE,linewidth=1.5,alpha=0.8)+
  geom_point(size=1.5,alpha=0.8)+
  scale_x_continuous(breaks=c(1990,2000,2010,2020),minor_breaks=c(1985,1995,2005,2015))+
  guides(x=guide_axis(minor.ticks = TRUE))+
  scale_color_manual(values=rep(c("#F181B9","#1b9e77"),2))+
  labs(x="",y="Number of Days",color="Region")+
  facet_wrap(.~metric,nrow=2,scales="free_y",labeller=labeller(metric=metric.labs))+
  theme(strip.background=element_blank(),legend.position="none",text=element_text(size=15),panel.grid.minor.x=element_blank(),legend.text=element_text(size=15),
        ,plot.margin=unit(c(5.5,12,5.5,5.5),"points"))+
  geom_text(data=cuezonephenolletz%>%filter(letz%in%c("I","L")),aes(x=-Inf,y=Inf,label=letz),hjust=-1,vjust=1.5,inherit.aes=FALSE)->il

ghjk+il+plot_layout(widths=c(2.15,1))->panels

cuezonephenolletz%>%
  filter(letz=="I")%>%
  ggplot(aes(x=year,y=value,group=region,color=region))+
  stat_smooth(geom="line",method="gam",se=FALSE,linewidth=1.5)+
  scale_color_manual(values=c("#F181B9","#1b9e77"),name="Region")+
  labs(x="Year",y="Number of Days",color="Region")+
  theme(legend.position="bottom",text=element_text(size=15),legend.text=element_text(size=15))->legplt

cowplot::get_plot_component(legplt,"guide-box",return_all=TRUE)[[2]]->leg

leg<-cowplot::get_legend(legplt)

panels/leg+plot_layout(heights=c(1,0.01))->p2

sensplt/p2+
  plot_layout(heights=c(2,1))->allplt

allplt

ggsave("Figure 3.pdf",allplt,width=10.5,height=15,units="in",dpi=300)

# 1SW 5th
# GAM
m1<-gam(X1SW_Q05raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinend,k=4)+
          s(southern_mnwinend,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnwinend, k = 4)`&&`s(northern_mnwinlen, k = 4)`)&&
              !(`s(southern_mnwinend, k = 4)`&&`s(southern_mnwinlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2","adjR^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat)

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$X1SW_Q05raw,
       Q05raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q05raw_preds-se.fit,ymax=Q05raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q05raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q05raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="1SW 5th Percentile Run Timing")+
  theme(text=element_text(size=14))->p1sw5th

p1sw5th

ggsave(filename="p1sw5th.png",plot=p1sw5th,dpi=300,width=7,height=3.5,units="in")

summary(gam(X1SW_Q05raw~
              s(northern_mnsumlen_lag,k=4)+
              s(northern_mnwinend,k=4)+
              s(southern_mnwinend,k=4)+
              s(northern_mnwinlen,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4)+
              s(sprinfl,k=4)+
              s(Q50,k=4)+
              s(rivsumsrt,k=4),
            data=subdat,method="REML"))

# 1
summary(gam(X1SW_Q05raw~
              s(northern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(X1SW_Q05raw~
              s(northern_mnwinlen,k=4)+
              s(southern_mnwinend,k=4),
            data=subdat,method="REML",select=TRUE))

# 3
summary(gam(X1SW_Q05raw~
              s(northern_mnwinlen,k=4)+
              s(rivsumsrt,k=4),
            data=subdat,method="REML",select=TRUE))

# 4
summary(gam(X1SW_Q05raw~
              s(southern_mnwinlen,k=4)+
              s(sprinfl,k=4),
            data=subdat,method="REML",select=FALSE))

# 5
summary(gam(X1SW_Q05raw~
              s(northern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 1SW 25th
# GAM
m1<-gam(X1SW_Q25raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinend,k=4)+
          s(southern_mnwinend,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnwinend, k = 4)`&&`s(northern_mnwinlen, k = 4)`)&&
              !(`s(southern_mnwinend, k = 4)`&&`s(southern_mnwinlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$X1SW_Q25raw,
       Q25raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q25raw_preds-se.fit,ymax=Q25raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q25raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q25raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="1SW 25th Percentile Run Timing")+
  theme(text=element_text(size=14))->p1sw25th

p1sw25th

ggsave(filename="p1sw25th.png",plot=p1sw25th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(X1SW_Q25raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="1SW 25th")+
      coord_cartesian(ylim=c(160,210))
    
  })
}

# 1
summary(gam(X1SW_Q25raw~
              s(northern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(X1SW_Q25raw~
              s(northern_mnwinlen,k=4)+
              s(southern_mnwinend,k=4),
            data=subdat,method="REML",select=TRUE))

# 3
summary(gam(X1SW_Q25raw~
              s(northern_mnwinlen,k=4)+
              s(Q50,k=4),
            data=subdat,method="REML",select=TRUE))

# 1SW 50th
# GAM
m1<-gam(X1SW_Q50raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinend,k=4)+
          s(southern_mnwinend,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnwinend, k = 4)`&&`s(northern_mnwinlen, k = 4)`)&&
              !(`s(southern_mnwinend, k = 4)`&&`s(southern_mnwinlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$X1SW_Q50raw,
       Q50raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q50raw_preds-se.fit,ymax=Q50raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q50raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q50raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="1SW 50th Percentile Run Timing")+
  theme(text=element_text(size=14))->p1sw50th

p1sw50th

ggsave(filename="p1sw50th.png",plot=p1sw50th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(X1SW_Q50raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="1SW 50th")+
      coord_cartesian(ylim=c(160,240))
    
  })
}

# 1
summary(gam(X1SW_Q50raw~
              s(northern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(X1SW_Q50raw~
              s(northern_mnwinlen,k=4)+
              s(southern_mnwinend,k=4),
            data=subdat,method="REML",select=TRUE))

# 1SW 75th
# GAM
m1<-gam(X1SW_Q75raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(northern_mnsumsrt,k=4)+
          s(southern_mnsumsrt,k=4)+
          s(northern_mnsumlen,k=4)+
          s(southern_mnsumlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnsumsrt, k = 4)`&&`s(northern_mnsumlen, k = 4)`)&&
              !(`s(southern_mnsumsrt, k = 4)`&&`s(southern_mnsumlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$X1SW_Q75raw,
       Q75raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q75raw_preds-se.fit,ymax=Q75raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q75raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q75raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="1SW 75th Percentile Run Timing")+
  theme(text=element_text(size=14))->p1sw75th

p1sw75th

ggsave(filename="p1sw75th.png",plot=p1sw75th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(X1SW_Q75raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="1SW 75th")+
      coord_cartesian(ylim=c(160,260))
    
  })
}

# 1
summary(gam(X1SW_Q75raw~
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen,k=4)+
              s(southern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(X1SW_Q75raw~
              s(northern_mnwinlen,k=4)+
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 3
summary(gam(X1SW_Q75raw~
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 4
summary(gam(X1SW_Q75raw~
              s(northern_mnwinlen,k=4)+
              s(Q50,k=4)+
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 5
summary(gam(X1SW_Q75raw~
              s(northern_mnwinlen,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 6
summary(gam(X1SW_Q75raw~
              s(northern_mnwinlen,k=4)+
              s(Q50,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 1SW 95th
# GAM
m1<-gam(X1SW_Q95raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(northern_mnsumsrt,k=4)+
          s(southern_mnsumsrt,k=4)+
          s(northern_mnsumlen,k=4)+
          s(southern_mnsumlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnsumsrt, k = 4)`&&`s(northern_mnsumlen, k = 4)`)&&
              !(`s(southern_mnsumsrt, k = 4)`&&`s(southern_mnsumlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)[[1]]

# predictions from each model
predz<-predict(pred.parms,newdata=subdat,se.fit=TRUE)

fit<-predz$fit
se.fit<-predz$se.fit

tibble(year=alldat$year,
       raw_obs=subdat$X1SW_Q95raw,
       Q95raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q95raw_preds-se.fit,ymax=Q95raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q95raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q95raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="1SW 95th Percentile Run Timing")+
  theme(text=element_text(size=14))->p1sw95th

p1sw95th

ggsave(filename="p1sw95th_1model.png",plot=p1sw95th,dpi=300,width=7,height=3.5,units="in")

# 1
summary(gam(X1SW_Q95raw~
              s(southern_mnwinlen,k=4)+
              s(rivsumsrt,k=4)+
              s(southern_mnsumsrt,k=4)+
              s(sprinfl,k=4),
            data=subdat,method="REML",select=TRUE))

# MSW 5th
# GAM
m1<-gam(MSW_Q05raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinend,k=4)+
          s(southern_mnwinend,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnwinend, k = 4)`&&`s(northern_mnwinlen, k = 4)`)&&
              !(`s(southern_mnwinend, k = 4)`&&`s(southern_mnwinlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$MSW_Q05raw,
       Q05raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q05raw_preds-se.fit,ymax=Q05raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q05raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q05raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="MSW 5th Percentile Run Timing")+
  theme(text=element_text(size=14))->pmsw5th

pmsw5th

ggsave(filename="pmsw5th.png",plot=pmsw5th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(MSW_Q05raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="MSW 5th")+
      coord_cartesian(ylim=c(120,180))
    
  })
}

# 1
summary(gam(MSW_Q05raw~
              s(northern_mnsumlen_lag,k=4)+
              s(Q50,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(MSW_Q05raw~
              s(northern_mnsumlen_lag,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# MSW 25th
# GAM
m1<-gam(MSW_Q25raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinend,k=4)+
          s(southern_mnwinend,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnwinend, k = 4)`&&`s(northern_mnwinlen, k = 4)`)&&
              !(`s(southern_mnwinend, k = 4)`&&`s(southern_mnwinlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$MSW_Q25raw,
       Q25raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q25raw_preds-se.fit,ymax=Q25raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q25raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q25raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="MSW 25th Percentile Run Timing")+
  theme(text=element_text(size=14))->pmsw25th

ggsave(filename="pmsw25th.png",plot=pmsw25th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(MSW_Q25raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="MSW 25th")+
      coord_cartesian(ylim=c(120,180))
    
  })
}

# 1
summary(gam(MSW_Q25raw~
              s(northern_mnsumlen_lag,k=4)+
              s(Q50,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(MSW_Q25raw~
              s(Q50,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 3
summary(gam(MSW_Q25raw~
              s(northern_mnsumlen_lag,k=4)+
              s(southern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 4
summary(gam(MSW_Q25raw~
              s(wscv,k=4)+
              s(southern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 5
summary(gam(MSW_Q25raw~
              s(wscv,k=4)+
              s(northern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 6
summary(gam(MSW_Q25raw~
              s(Q50,k=4)+
              s(southern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# MSW 50th
# GAM
m1<-gam(MSW_Q50raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinend,k=4)+
          s(southern_mnwinend,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnwinend, k = 4)`&&`s(northern_mnwinlen, k = 4)`)&&
              !(`s(southern_mnwinend, k = 4)`&&`s(southern_mnwinlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$MSW_Q50raw,
       Q50raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q50raw_preds-se.fit,ymax=Q50raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q50raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q50raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="MSW 50th Percentile Run Timing")+
  theme(text=element_text(size=14))->pmsw50th

pmsw50th

ggsave(filename="pmsw50th.png",plot=pmsw50th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(MSW_Q50raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="MSW 50th")+
      coord_cartesian(ylim=c(150,190))
    
  })
}

# 1
summary(gam(MSW_Q50raw~
              s(northern_mnwinlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(MSW_Q50raw~
              s(northern_mnwinlen,k=4)+
              s(southern_mnwinend,k=4),
            data=subdat,method="REML",select=TRUE))

# 3
summary(gam(MSW_Q50raw~
              s(northern_mnwinlen,k=4)+
              s(wscv,k=4),
            data=subdat,method="REML",select=TRUE))

# 4
summary(gam(MSW_Q50raw~
              s(northern_mnwinlen,k=4)+
              s(sprinfl,k=4),
            data=subdat,method="REML",select=TRUE))

# MSW 75th
# GAM
m1<-gam(MSW_Q75raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(northern_mnsumsrt,k=4)+
          s(southern_mnsumsrt,k=4)+
          s(northern_mnsumlen,k=4)+
          s(southern_mnsumlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnsumsrt, k = 4)`&&`s(northern_mnsumlen, k = 4)`)&&
              !(`s(southern_mnsumsrt, k = 4)`&&`s(southern_mnsumlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$MSW_Q75raw,
       Q75raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q75raw_preds-se.fit,ymax=Q75raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q75raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q75raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="MSW 75th Percentile Run Timing")+
  theme(text=element_text(size=14))->pmsw75th

pmsw75th

ggsave(filename="pmsw75th.png",plot=pmsw75th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(MSW_Q75raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="MSW 75th")+
      coord_cartesian(ylim=c(150,250))
    
  })
}

# 1
summary(gam(MSW_Q75raw~
              s(northern_mnwinlen,k=4)+
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(MSW_Q75raw~
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen),
            data=subdat,method="REML",select=TRUE))

# MSW 95th
# GAM
m1<-gam(MSW_Q95raw~
          s(northern_mnsumlen_lag,k=4)+
          s(northern_mnwinlen,k=4)+
          s(southern_mnwinlen,k=4)+
          s(northern_mnsumsrt,k=4)+
          s(southern_mnsumsrt,k=4)+
          s(northern_mnsumlen,k=4)+
          s(southern_mnsumlen,k=4)+
          s(wscv,k=4)+
          s(sprinfl,k=4)+
          s(Q50,k=4)+
          s(rivsumsrt,k=4),
        data=subdat,method="REML")

# model subset criteria
out<-dredge(m1,
            subset=!(`s(northern_mnwinlen, k = 4)`&&`s(southern_mnwinlen, k = 4)`)&&
              !(`s(northern_mnsumsrt, k = 4)`&&`s(northern_mnsumlen, k = 4)`)&&
              !(`s(southern_mnsumsrt, k = 4)`&&`s(southern_mnsumlen, k = 4)`),
            m.lim=c(1,4),extra=c("R^2"))

# keep all candidate models with delta AIC < d
d<-2
tab<-tibble(out[out$delta<d,])

# table formatting
newname<-c()
for(i in 1:length(colnames(tab))){
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])==0){
    newname<-c(newname,colnames(tab)[i])
  }
  if(length(codez$cleanname[which(colnames(tab)[i]==codez$modelname)])>0){
    newname<-c(newname,codez$cleanname[which(colnames(tab)[i]==codez$modelname)])
  }
}
colnames(tab)<-newname
tab$`(Intercept)`<-round(tab$`(Intercept)`,2)
tab$`R^2`<-round(tab$`R^2`,2)
tab$logLik<-round(tab$logLik,2)
tab$AICc<-round(tab$AICc,2)
tab$delta<-round(tab$delta,2)
tab$weight<-round(tab$weight,2)

# remove NAs
tab[sapply(tab,is.factor)]<-lapply(tab[sapply(tab,is.factor)],as.character)
tab<-replace(tab,is.na(tab),"")

tab

# model averaging, revised formula for averaging standard errors
MA.ests<-model.avg(out,subset=delta<d,revised.var=TRUE)

# plot model averaged coefficients
coef<-MA.ests$coefficients[2,-1]

# relative importance of variables
imp<-sw(MA.ests)

tibble(var=factor(rev(names(imp)),levels=rev(names(imp))),val=rev(as.numeric(imp[1:length(imp)])))%>%
  ggplot(aes(x=var,y=val))+
  geom_bar(stat="identity")+
  labs(x="",y="relative importance")+
  coord_flip()->varimp

varimp

# retrieve model subset
pred.parms<-get.models(out,subset=delta<d)

# predictions from each model
model.preds<-sapply(pred.parms,predict,newdata=subdat) 

# predictions from top models
fit<-c()
se.fit<-c()
for(j in 1:length(pred.parms)){
  predict(pred.parms[[j]],newdata=subdat,se.fit=TRUE)->predz
  fit<-cbind(fit,predz$fit)
  se.fit<-cbind(se.fit,predz$se.fit)
}

# average of model predictions weighted based on aic
# rescale subset of weights to equal 1
wts<-(out$weight[1:length(pred.parms)])/sum(out$weight[1:length(pred.parms)])
mod.ave.preds<-model.preds%*%wts

# average model predictions by aic weights
fit<-c(fit%*%wts)
se.fit<-c(se.fit%*%wts)

tibble(year=alldat$year,
       raw_obs=subdat$MSW_Q95raw,
       Q95raw_preds=fit,
       se.fit=2*se.fit)%>%
  ggplot()+
  geom_line(aes(x=year,y=raw_obs),linewidth=1)+
  geom_point(aes(x=year,y=raw_obs),pch=21,fill="black",alpha=0.7)+
  geom_ribbon(aes(x=year,ymin=Q95raw_preds-se.fit,ymax=Q95raw_preds+se.fit),alpha=0.5,fill="#E35612",color="transparent")+
  geom_line(aes(x=year,y=Q95raw_preds),color="#E35612",linewidth=1,alpha=0.8)+
  geom_point(aes(x=year,y=Q95raw_preds),pch=21,color="#E35612",fill="#E35612",alpha=0.7)+
  labs(x="Year",y="Numerical Day of Year",subtitle="MSW 95th Percentile Run Timing")+
  theme(text=element_text(size=14))->pmsw95th

pmsw95th

ggsave(filename="pmsw95th.png",plot=pmsw95th,dpi=300,width=7,height=3.5,units="in")

# variable names 
names(imp)%>%
  str_replace(",.*","")%>%
  str_sub(start=3,end=99)->v

# loop for variable effect plots
pltlst<-list()
for(i in 1:length(v)){
  
  # store in list
  pltlst[[i]]<-local({
    
    # variable
    var<-v[i]
    
    # sequence along variable range
    rng<-subdat%>%select(all_of(var))%>%range()
    envvar<-seq(rng[1],rng[2],length.out=100)
    
    # all other variables set to average values
    plotdata<-as.data.frame(lapply(lapply(subdat%>%select(-c(MSW_Q95raw,var)),mean),rep,length(envvar)))
    plotdata<-cbind(envvar,plotdata)
    colnames(plotdata)[1]<-var
    
    # predictions from top models
    fit<-c()
    se.fit<-c()
    for(j in 1:length(pred.parms)){
      predict(pred.parms[[j]],newdata=plotdata,se.fit=TRUE)->predz
      fit<-cbind(fit,predz$fit)
      se.fit<-cbind(se.fit,predz$se.fit)
    }
    
    # average model predictions by aic weights
    fit<-fit%*%wts
    se.fit<-se.fit%*%wts
    
    # save plot
    p1<-ggplot()+
      geom_line(aes(x=envvar,y=fit))+
      geom_ribbon(aes(x=envvar,ymin=fit-se.fit,ymax=fit+se.fit),alpha=0.5,fill="#F8766D")+
      labs(x=var,y="MSW 95th")+
      coord_cartesian(ylim=c(170,290))
    
  })
}

# 1
summary(gam(MSW_Q95raw~
              s(northern_mnwinlen,k=4)+
              s(Q50,k=4)+
              s(rivsumsrt,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# 2
summary(gam(MSW_Q95raw~
              s(northern_mnwinlen,k=4)+
              s(Q50,k=4)+
              s(southern_mnsumlen,k=4),
            data=subdat,method="REML",select=TRUE))

# Figure 4
gamtab<-read_csv("gamtableforfigure2xxxdev.csv")
astertab<-read_csv("gamtableforfigureasterisksEDFs2xxxdev.csv",col_types="ccccccccccccccc")

colnames(gamtab)<-c("percentile",
                    "Northern Summer Length Lag",
                    "Northern Winter End",
                    "Southern Winter End",
                    "Northern Winter Length",
                    "Southern Winter Length",
                    "Northern Summer Start",
                    "Southern Summer Start",
                    "Northern Summer Length",
                    "Southern Summer Length",
                    "Spring Inflection Point",
                    "Winter-Spring Center of Volume",
                    "Annual Center of Volume",
                    "River Summer Start",
                    "Dev. Expl. (# of Models)")

colnames(astertab)<-c("percentile",
                      "Northern Summer Length Lag",
                      "Northern Winter End",
                      "Southern Winter End",
                      "Northern Winter Length",
                      "Southern Winter Length",
                      "Northern Summer Start",
                      "Southern Summer Start",
                      "Northern Summer Length",
                      "Southern Summer Length",
                      "Spring Inflection Point",
                      "Winter-Spring Center of Volume",
                      "Annual Center of Volume",
                      "River Summer Start",
                      "Dev. Expl. (# of Models)")

# plot
gamtab%>%
  pivot_longer(-percentile)%>%
  left_join(.,astertab%>%pivot_longer(-percentile,values_to="aster"))%>%
  mutate(percentile=factor(percentile,levels=rev(gamtab$percentile)),
         name=factor(name,levels=colnames(gamtab)[-1]),)%>%
  filter(value>=0)%>%
  mutate(value=value*100)%>%
  filter(percentile!="75th_1SW"&percentile!="95th_1SW"&percentile!="75th_MSW"&percentile!="95th_MSW")%>%
  add_column(seaage=c(rep("1SW",30),rep("MSW",30),rep("1SW",24),rep("MSW",24)))%>%
  ggplot(aes(x=name,y=percentile))+
  geom_tile(aes(alpha=value),fill="#9ebcda")+          ##88419d
  geom_text(aes(label=aster),color="black",size=3)+ ##9ebcda
  scale_x_discrete(labels=wrap_format(10))+
  scale_y_discrete(labels=c(
    "5th_1SW"="5th",
    "25th_1SW"="25th",
    "50th_1SW"="50th",
    "75th_1SW_alt"="75th",
    "95th_1SW_alt"="95th",
    "5th_MSW"="5th",
    "25th_MSW"="25th",
    "50th_MSW"="50th",
    "75th_MSW_alt"="75th",
    "95th_MSW_alt"="95th"))+
  coord_cartesian(expand=FALSE)+
  facet_wrap(.~seaage,nrow=2,scales="free_y")+
  labs(x="",y="",alpha="Relative Importance %")+
  theme(legend.position="bottom",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_blank())->poi

poi

ggsave("Figure 4.pdf",poi,width=10,height=4,units="in",dpi=300)

# Supplemental Figure 1
runmetrics<-read_csv("runmetrics11142022.csv")

smoothlabs<-c("no smooth","5-day","10-day")
names(smoothlabs)<-c("raw","5day","10day")

perclabs<-c("5th","25th","50th","75th","95th")
names(perclabs)<-c("Q05","Q25","Q50","Q75","Q95")

# plots
runmetrics%>%
  pivot_longer(cols=Q05raw:Q9510day,names_to=c("percentile","smooth"),names_sep=3,values_to="yearday")%>%
  filter(year>=1983&year<2020)%>%
  ggplot(aes(x=year,y=yearday))+
  geom_point(aes(group=smooth,color=smooth),alpha=0.6,linewidth=2)+
  geom_smooth(aes(group=smooth,color=smooth),se=FALSE,method="lm",size=1)+
  labs(x="Year",y="Yearday",color="Observation smooths")+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"),labels=smoothlabs)+
  facet_wrap(.~seaage+percentile,scales="free_y",ncol=5,labeller=labeller(percentile=perclabs))+
  theme_bw()+
  theme(strip.background=element_blank(),legend.position="bottom",plot.margin=unit(c(5.5,10,5.5,5.5),"points"))->sf1

sf1

ggsave(filename="Supp Figure 1.pdf",plot=sf1,width=9,height=5,dpi=300)

# Supplemental Figure 2
# percentiles
runmetrics<-read_csv("runmetrics11142022.csv")
# raw
PNruntiming<-read_csv("PNruntiming11142022.csv")
# make Date object
PNruntiming$Date<-as.Date(PNruntiming$Date)

dispo<-ifelse(PNruntiming$Disposition=="CRAIG BROOK NAT. FISH HATCHERY-OOB","CRAIG BROOK","OTHER")

# daily counts
PNruntiming%>%
  mutate(dispo=dispo)%>%
  filter(seaage=="1SW")%>%
  group_by(Date,dispo)%>%
  count()->x

x%>%
  pivot_wider(names_from="dispo",values_from="n")->x

# full date time series
Date<-seq.Date(from=as.Date("1978-01-01"),to=as.Date("2019-12-31"),by="day")
fulldate<-tibble(Date)
fulldate<-left_join(fulldate,x)
fulldate<-fulldate[,-4]
fulldate$`CRAIG BROOK`<-ifelse(is.na(fulldate$`CRAIG BROOK`),0,fulldate$`CRAIG BROOK`)
fulldate$OTHER<-ifelse(is.na(fulldate$OTHER),0,fulldate$OTHER)

# ridge plot 
fulldate$yearday<-yday(fulldate$Date)
fulldate$year<-year(fulldate$Date)

#1SW
fulldate%>%
  group_by(year=year(Date))%>%
  mutate(runtotraw=cumsum(`CRAIG BROOK`))%>%
  mutate(run50raw=max(runtotraw)*0.50)%>%
  summarise(Q50raw=yday(Date[runtotraw>=run50raw][1]),tot=max(runtotraw))->oneSWcraigbrook50

fulldate%>%
  group_by(year=year(Date))%>%
  mutate(runtotraw=cumsum(`CRAIG BROOK`))%>%
  mutate(run05raw=max(runtotraw)*0.05,
         run25raw=max(runtotraw)*0.25,
         run50raw=max(runtotraw)*0.50,
         run75raw=max(runtotraw)*0.75,
         run95raw=max(runtotraw)*0.95)%>%
  summarise(Q05raw=yday(Date[runtotraw>=run05raw][1]),
            Q25raw=yday(Date[runtotraw>=run25raw][1]),
            Q50raw=yday(Date[runtotraw>=run50raw][1]),
            Q75raw=yday(Date[runtotraw>=run75raw][1]),
            Q95raw=yday(Date[runtotraw>=run95raw][1]))->oneSWcraigbrook

# raw
PNruntiming<-read_csv("PNruntiming11142022.csv")

# make Date object
PNruntiming$Date<-as.Date(PNruntiming$Date)
dispo<-ifelse(PNruntiming$Disposition=="CRAIG BROOK NAT. FISH HATCHERY-OOB","CRAIG BROOK","OTHER")

# daily counts
PNruntiming%>%
  mutate(dispo=dispo)%>%
  filter(seaage=="MSW")%>%
  group_by(Date,dispo)%>%
  count()->x

x%>%pivot_wider(names_from="dispo",values_from="n")->x

# full date time series
Date<-seq.Date(from=as.Date("1978-01-01"),to=as.Date("2019-12-31"),by="day")
fulldate<-tibble(Date)
fulldate<-left_join(fulldate,x)
fulldate<-fulldate[,-3]
fulldate$`CRAIG BROOK`<-ifelse(is.na(fulldate$`CRAIG BROOK`),0,fulldate$`CRAIG BROOK`)
fulldate$OTHER<-ifelse(is.na(fulldate$OTHER),0,fulldate$OTHER)

# ridge plot 
fulldate$yearday<-yday(fulldate$Date)
fulldate$year<-year(fulldate$Date)

#MSW
fulldate%>%
  group_by(year=year(Date))%>%
  mutate(runtotraw=cumsum(`CRAIG BROOK`))%>%
  mutate(run50raw=max(runtotraw)*0.50)%>%
  summarise(Q50raw=yday(Date[runtotraw>=run50raw][1]),tot=max(runtotraw))->MSWcraigbrook50

fulldate%>%
  group_by(year=year(Date))%>%
  mutate(runtotraw=cumsum(`CRAIG BROOK`))%>%
  mutate(run05raw=max(runtotraw)*0.05,
         run25raw=max(runtotraw)*0.25,
         run50raw=max(runtotraw)*0.50,
         run75raw=max(runtotraw)*0.75,
         run95raw=max(runtotraw)*0.95)%>%
  summarise(Q05raw=yday(Date[runtotraw>=run05raw][1]),
            Q25raw=yday(Date[runtotraw>=run25raw][1]),
            Q50raw=yday(Date[runtotraw>=run50raw][1]),
            Q75raw=yday(Date[runtotraw>=run75raw][1]),
            Q95raw=yday(Date[runtotraw>=run95raw][1]))->MSWcraigbrook

# color by group
# facet by seaage
perc.labs<-c("5th","25th","50th","75th","95th")
names(perc.labs)<-c("Q05raw","Q25raw","Q50raw","Q75raw","Q95raw")

runmetrics%>%
  select(c(year,Q05raw,Q25raw,Q50raw,Q75raw,Q95raw,seaage))%>%
  add_column(grp="All Returns")%>%
  bind_rows(.,oneSWcraigbrook%>%
              add_column(seaage="1SW",grp="Broodstock Returns"))%>%
  bind_rows(.,MSWcraigbrook%>%
              add_column(seaage="MSW",grp="Broodstock Returns"))%>%
  filter(year%in%c(1983:2019)&Q50raw!=1)%>%
  pivot_longer(cols=c(Q05raw,Q25raw,Q50raw,Q75raw,Q95raw),names_to="percentile",names_pattern="Q(.*)raw",values_to="yearday")%>%
  group_by(seaage,percentile)%>%
  do(tidy(lm(yearday~year*grp,.)))%>%
  mutate(p.value=round(p.value,4))%>%
  filter(term!="(Intercept)"&term!="grpBroodstock Returns")%>%
  select(-c(std.error,statistic))%>%
  pivot_wider(names_from=seaage,values_from=c(estimate,p.value))%>%
  relocate(estimate_MSW,.after=p.value_1SW)%>%
  filter(term!="year")%>%
  select(percentile,p.value_1SW,p.value_MSW)->pv

tibble(seaage=c(rep("1SW",5),rep("MSW",5)),
       percentile=rep(c("Q05raw","Q25raw","Q50raw","Q75raw","Q95raw"),2),
       pv=c(pv$p.value_1SW,pv$p.value_MSW),
       grp=rep("Broodstock Returns",10))->pvs

# plot percentiles together
runmetrics%>%
  select(c(year,Q05raw,Q25raw,Q50raw,Q75raw,Q95raw,seaage))%>%
  add_column(grp="All Returns")%>%
  bind_rows(.,oneSWcraigbrook%>%
              add_column(seaage="1SW",grp="Broodstock Returns"))%>%
  bind_rows(.,MSWcraigbrook%>%
              add_column(seaage="MSW",grp="Broodstock Returns"))%>%
  filter(year%in%c(1983:2019)&Q50raw!=1)%>%
  pivot_longer(cols=c(Q05raw,Q25raw,Q50raw,Q75raw,Q95raw),
               names_to="percentile",values_to="yearday")%>%
  ggplot(aes(x=year,y=yearday,group=interaction(percentile,grp),color=grp))+
  geom_line(alpha=0.2,linewidth=1)+
  geom_point(aes(fill=grp),color="black",shape=21,size=2,alpha=0.5)+
  geom_smooth(method="lm",se=FALSE)+
  labs(x="Year",y="Return Yearday",color="Group",fill="Group")+
  facet_wrap(.~seaage+percentile,nrow=2,labeller=labeller(percentile=perc.labs))+
  geom_text(data=pvs,mapping=aes(x=2016,y=288,label=round(pv,3)),color="black",fontface=c(rep("plain",3),"bold",rep("plain",6)))+
  theme(legend.position="bottom",strip.background=element_blank(),
        plot.margin=unit(c(5.5,10,5.5,5.5),"points"),text=element_text(size=13))->sf2

sf2

ggsave(filename="Supp Figure 2.pdf",plot=sf2,width=10,height=6,dpi=300)

# Supplemental Figure 3
# Run Timing data from Milford
dat<-read_csv("tblTendData2.csv")
dat<-dat%>%select(c(TendDate,SiteCode,WaterTemp))
colnames(dat)<-c("date","site","temperature")
dat$date<-as.Date(dat$date,format="%m/%d/%Y")
dat<-arrange(dat,date)

# filter out 1STILLW0.20
dat%>%filter(site!="1STILLW0.20")->dat

# daily mean
dat<-dat%>%
  group_by(date)%>%
  summarise(temperature=mean(temperature,na.rm=TRUE))

dat$temperature[which(dat$temperature=="NaN")]<-NA
dat<-dat%>%filter(!is.na(temperature))

# USGS temperature data
temp<-read_csv("PenobscotEddington.csv")
colnames(temp)<-c("date","temperature")

# date sequence
date<-tibble(date=seq.Date(as.Date("1978-01-01"),as.Date("2019-12-31"),"day"))

# merge run timing temperature to date sequence
riversites<-left_join(date,dat,by="date")
riversites$type<-NA
riversites$type[which(!is.na(riversites$temperature))]<-"Trap"

temp$date<-as.Date(temp$date)

# merge eddington temperature to riversites
riversites<-left_join(riversites,temp,by="date")
riversites%>%mutate(tmp=coalesce(temperature.x,temperature.y))->riversites
riversites$type[which(is.na(riversites$type)&!is.na(riversites$tmp))]<-"Gage"

# year, month, yearday
riversites<-riversites%>%mutate(year=year(date),month=month(date),yearday=yday(date))

# seasonal groupings
riversites$season<-ifelse(riversites$month%in%c(12,1,2),"djf",
                          ifelse(riversites$month%in%c(3:5),"mam",
                                 ifelse(riversites$month%in%c(6:8),"jja","son")))
riversites$season<-factor(riversites$season,levels=c("djf","mam","jja","son"))

# fill in remaining gaps with air temp model
airmod<-read_csv("Eddington_daily_obs.csv")

# kelvin to celsius, if negative make 0
airmod$daily_avg<-airmod$daily_avg-273.15
airmod$daily_avg<-ifelse(airmod$daily_avg<0,0,airmod$daily_avg)
colnames(airmod)<-c("daily_avg","date")

# merge airmod temperature to riversites
riversites<-left_join(riversites,airmod,by="date")
riversites%>%mutate(tmp=coalesce(tmp,daily_avg))->riversites
riversites$type[which(is.na(riversites$type)&!is.na(riversites$tmp))]<-"Air model"

# yearday temperature by year
riversites%>%
  filter(year%in%c(1983:2019))%>%
  filter(!is.na(tmp))%>%
  ggplot(aes(x=yearday,y=tmp,group=year,color=type))+
  geom_point(alpha=0.2,size=2)+
  labs(x="Yearday",y="Temperature C",color="Data Source")+
  scale_x_continuous(breaks=c(60,180,335))+
  facet_wrap(year~.,ncol=10)+
  theme_bw()+
  theme(legend.position="bottom",text=element_text(size=17),
        strip.background=element_blank())->sf3

sf3

ggsave(filename="Supp Figure 3.pdf",plot=sf3,width=11,height=7,dpi=300)

# Supplemental Figure 4
riversites%>%select(date,tmp,year,yearday)->riverTdaily

riverTdaily%>%
  filter(year%in%c(1982:2019))%>%
  group_by(month=month(date))%>%
  summarize(meantmp=mean(tmp,na.rm=T))%>%
  add_column(monthabb=month.abb)%>%
  ggplot(aes(x=month,y=meantmp))+
  geom_point(size=2)+
  geom_line()+
  geom_area(alpha=0.3)+
  labs(x="Month",y="Temperature C",subtitle="1982-2019")+
  scale_x_continuous(breaks=seq_along(month.name),labels=month.abb)->sf4

sf4

ggsave(filename="Supp Figure 4.pdf",plot=sf4,width=5,height=3,dpi=300)

# Supplemental Figure 5
yrz<-alldat$year
rundur1SW<-alldat$`1SW_Q95raw`-alldat$`1SW_Q05raw`
rundurMSW<-alldat$MSW_Q95raw-alldat$MSW_Q05raw

datz<-tibble(Year=rep(yrz,2),Duration=c(rundur1SW,rundurMSW),SeaAge=c(rep("1SW",37),rep("MSW",37)))

datz%>%
  ggplot(aes(x=Year,y=Duration,group=SeaAge,color=SeaAge))+
  geom_point()+
  geom_smooth(method="lm")->sf5

sf5

ggsave(filename="Supp Figure 5.pdf",plot=sf5,width=6,height=3,dpi=300)

# Supplemental Figure 6
alldat$rivsumlen<-alldat$rivsumend-alldat$rivsumsrt

met.labs<-c("Annual Center of Volume","River Summer End","Spring Inflection Point","River Summer Length","River Summer Start","Winter Spring Center of Volume")
names(met.labs)<-c("Q50","rivsumend","sprinfl","rivsumlen","rivsumsrt","wscv")

letz<-tibble(letz=c("A","B","C","D","E","F"),metric=c("wscv","Q50","sprinfl","rivsumsrt","rivsumend","rivsumlen"))

alldat%>%
  select(year,wscv,Q50,sprinfl,rivsumsrt,rivsumend,rivsumlen)%>%
  pivot_longer(-year,names_to="metric",values_to="yearday or number of days")%>%
  left_join(letz,by="metric")%>%
  mutate(metric=factor(metric,levels=c("wscv","Q50","sprinfl","rivsumsrt","rivsumend","rivsumlen")))->alldatriv

alldatriv%>%
  ggplot(aes(x=year,y=`yearday or number of days`))+
  geom_line()+
  stat_smooth(geom="line",method="gam",se=FALSE,linewidth=1.5,alpha=0.7)+
  geom_point(size=1.5,alpha=0.7)+
  scale_x_continuous(breaks=seq(1980,2020,5))+
  labs(x="Year",y="Yearday or Number of Days")+
  facet_wrap(.~metric,scales="free_y",labeller=labeller(metric=met.labs))+
  theme(strip.background=element_blank(),panel.grid.minor.x=element_blank(),text=element_text(size=15),
        plot.margin=margin(5.5,10,5.5,5.5),"points")->p

p+geom_text(data=alldatriv,aes(x=-Inf,y=Inf,label=letz),hjust=-1,vjust=1.5,inherit.aes=FALSE)->sf6

sf6

ggsave(filename="Supp Figure 6.pdf",plot=sf6,width=11.13,height=5.3,dpi=300)
