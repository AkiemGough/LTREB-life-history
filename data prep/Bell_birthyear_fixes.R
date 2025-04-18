
##CODE FROM QRQC

# Indiana LDW -------------------------------------------------------------

library(readxl)
library(tidyverse)
## Note on data read-in: read_excel was coercing some data to NA
## I saved the most complete and up to date spreadsheet as .csv and that is read in here
## read_excel still used for plot read-in because that is small and simple

# Poa autumnalis in TX ----------------------------------------------------
#hi this is a test
## read in raw POAU data from Miller Lab drive folder
tom_dir<-"G:/Shared drives/Miller Lab/LTREB/POAU/"
bell_dir<- "/Users/bell/Library/CloudStorage/GoogleDrive-ics2@rice.edu/Shared drives/Miller Lab/LTREB/POAU/"
#use_dir<-tom_dir 
use_dir<-bell_dir 


## These data are published on EDI here
## https://doi.org/10.6073/pasta/ea7db07a578fb030a173f37f76596b62 (Accessed 2024-03-26).
## read in data and apply some of the data transformations used above
#indiana<-read.csv("data prep/LDW_LTREB_20072022.csv")

indiana<-read.csv("/Users/bell/Documents/GitHub/LTREB-life-history/data prep/LDW_LTREB_20072022.csv")
str(indiana)
indiana %>% 
  rowwise %>% 
  mutate(mean_spike_t = mean(c_across(c(spike_a_t,spike_b_t,spike_c_t)),na.rm=T),
         mean_spike_t1 = mean(c_across(c(spike_a_t1,spike_b_t1,spike_c_t1)),na.rm=T))->indiana


## add age as the difference of year_t and year_recruit
indiana$age<-indiana$year_t-indiana$birth
## assign age as NA for original plants
indiana$age[indiana$origin_01==0]<-NA
## are there as many -1's as there are 0's?
#table(indiana$age) #no there are way more 0s
## this must be because the recruit data were managed differently in the early years
## all the -1s (first appearance in year t1) started in 2016 or later

## check that -1 ages are always cases with no data in year_t - FUCK
indiana[which(indiana$age==-1 & !is.na(indiana$size_t)),]->bad_plants
#indiana %>% filter(id %in% bad_plants$id) %>% View

## many of these have NA id, why?
indiana[which(is.na(indiana$id)),]->no_id

## apply rule that if size is non-NA and inf count is NA, then inf count should be zero
indiana$flw_count_t[!is.na(indiana$size_t) & is.na(indiana$flw_count_t)]<-0
indiana$flw_count_t1[!is.na(indiana$size_t1) & is.na(indiana$flw_count_t1)]<-0

## check that each time an individual appears it carries the same birth year
#indiana %>% 
  #group_by(id) %>% select(birth) %>% 
  #summarize(nbirths=length(unique(birth,na.rm=T))) %>% 
  #filter(nbirths>1)->birthyears


## some of the problems are likely numbers that we reused, which I could know
## if the plant died but the number reappeared later. In other cases it is likely 
## a data entry or copying error, and it is probably safe to assume the earlier year
## is the correct one, as long as everything else checks out
## I will go through these manually and hand-pick cases where the earliest
## year is not the correct birth year

## For now I want to keep moving, so I am just going to drop "problems".
## Losing 596 rows. Tom will return to this (-TM 4/9/2024)
#indiana_no_problems <- indiana %>% filter(!(id %in% birthyears$id[problems]))


###START OF MY CODE FIXING DATA
indiana %>% 
  group_by(id) %>% select(birth) %>% 
  summarize(nbirths=length(unique(birth,na.rm=T))) %>% 
  filter(nbirths>1)->birthyears # these are the plants that have more then one birth year
#can think on 2 ways that happens 1) repeted ID, 2) a plant can not be found 
#one year and is assumed dead and then is found the next year 

#organising data
sus_ids <- birthyears$id
indiana_sus <- indiana[which(indiana$id %in% sus_ids),]

#what percent of plants are E+ POSYFAG
length(unique(indiana_sus[indiana_sus$species== "POSY",]$id))/
  length(unique(indiana_sus$id))

 indiana_sus %>% 
  group_by(id) %>% 
  mutate( max_birth = max(birth),
          min_death = min(year_t1[surv_t1 == 0])) ->indiana_sus

#solution for #1
#split plants that have any recorded deaths greater then one of their birth years into 2 ids 

indiana_sus %>% 
  filter (max_birth > min_death) %>% #was it born after it died?
  ungroup() %>% 
  mutate(id_update = ifelse(max_birth > birth, paste0(id, "a"), paste0(id, "b"))) -> diff_plants 

#so these are all the plant whos tag got repeted. i added an a to the end of 
#the id if it was the first one and a b to the end of the id if it was the 
#second one 
length(unique(diff_plants$id))

indiana_sus %>% 
  filter (max_birth <= min_death) ->miss_reports #so that we get everything not in the diffrent plant data 

#for right now lets just go trough them 1 by one 

ids_list <- (unique(miss_reports$id))

  
  
#OK SO I guess we just go buy this one by one and see what to fix
  miss_reports[which(miss_reports$id == ids_list[1]),]
  #went back to the data. looks like we first saw this plant in 2017 but it 
  #was assumed to be new in 2016 but i am unsure so i am droping this 

miss_reports <- miss_reports[-which(miss_reports$id == ids_list[1]), ]


miss_reports[which(miss_reports$id == ids_list[2]),]
#had note "PEEL (could be original 2123) {or old recruit" so i do not
#trust this data b

miss_reports <- miss_reports[-which(miss_reports$id == ids_list[2]), ]


miss_reports[which(miss_reports$id == ids_list[3]),]
#no data from year_t = 2017 so year_t1=2017 must be birth year 
miss_reports[which(miss_reports$id == ids_list[3]),]$birth = 2017


miss_reports[which(miss_reports$id == ids_list[4]),]
#note: "peel could be old recruit" so drop for uncertanty 

miss_reports <- miss_reports[-which(miss_reports$id == ids_list[4]), ]


miss_reports[which(miss_reports$id == ids_list[5]),]
#size_t = 1 at year_t = 2012 so 2012 is brith year but there is some crazy stuff
#with the years that dont match so lets drop 

miss_reports <- miss_reports[-which(miss_reports$id == ids_list[5]), ]


#FESU
miss_reports[which(miss_reports$id == ids_list[6]),]
#looks like the birth is actually 2016?
#size_t = 1 at year_t = 2016
miss_reports[which(miss_reports$id == ids_list[6]),]$birth = 2016


miss_reports[which(miss_reports$id == ids_list[7]),]
#first data is size_t1 = 1 at year_t1 = 2016
miss_reports[which(miss_reports$id == ids_list[7]),]$birth = 2016


#to me it seems like the majority of mistakes are in 2016-2017 and 2012-2013


#POAL
miss_reports[which(miss_reports$id == ids_list[8]),]
#first observation is size_t = 2 and year t= 2012
#but there is also duplacites of the same year where the only thing that 
#changes is the recorded birth year 
#so 3 and 5 of miss_reports[which(miss_reports$id == ids_list[8]),] are bad and
#should be droped (where the year is repeted with wrong birth)
miss_reports <- miss_reports[-which(miss_reports$id == ids_list[8])[c(3,5)], ]
miss_reports[which(miss_reports$id == ids_list[8]),] #checking that it is right

miss_reports[which(miss_reports$id == ids_list[9]),]
#first report is size_t = 1 at year_t = 2011
#so birth year is 2011
miss_reports[which(miss_reports$id == ids_list[9]),]$birth <- 2011



#i know the code looks the same for all of these but i did go through the data
#manualy and all of these plants that were born in 2013 and in 2017 the birth changed to 2014

miss_reports[which(miss_reports$id == ids_list[10]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[10]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[11]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[11]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[12]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[12]),]$birth <- 2013


#POSY
miss_reports[which(miss_reports$id == ids_list[13]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[13]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[14]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[14]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[15]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[15]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[16]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[16]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[17]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[17]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[18]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[18]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[19]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[19]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[20]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[20]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[21]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[21]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[22]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[22]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[23]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[23]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[24]),]
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[24]),]$birth <- 2013

miss_reports[which(miss_reports$id == ids_list[25]),] #wow changed from 2013-2016!
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[25]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[26]),] #also 2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[26]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[27]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[27]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[28]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[28]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[29]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[29]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[30]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[30]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[31]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[31]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[32]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[32]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[33]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[33]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[34]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[33]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[35]),] #2016
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[35]),]$birth <- 2013

#i want to come back to to posy
#b/c i cant find it in the raw data 

miss_reports[which(miss_reports$id == ids_list[36]),] 
#size_t =1  in year_t = 2014
miss_reports[which(miss_reports$id == ids_list[36]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[37]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[37]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[38]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[38]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[39]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[39]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[40]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[40]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[41]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[41]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[42]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[42]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[43]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[43]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[44]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[44]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[45]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[45]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[46]),] 
#oh wow this one is diffrent!!! 
#went from 2009-2007
#size_t = 1 at year_t = 2009
#looks like 2 plants. one born in 2009 that dies in 2013 and then one
#born in ??? (but probably 2018 b/c only one year of dat) and died in 2019
# i cant find the second plant so for now i am going to drop it for now
#*maybe come back*
miss_reports[which(miss_reports$id == ids_list[46]),]
miss_reports <- miss_reports[-which(miss_reports$id == ids_list[46])[c(5)], ]


miss_reports[which(miss_reports$id == ids_list[47]),] 
#yeat_t = 2014 size_t = NA
#year_t1 = 2015 size_t1 = 1
miss_reports[which(miss_reports$id == ids_list[47]),]$birth <- 2015



miss_reports[which(miss_reports$id == ids_list[48]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[48]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[49]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[49]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[50]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[50]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[51]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[51]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[52]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[52]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[53]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[53]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[54]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[54]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[55]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[55]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[56]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[56]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[57]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[57]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[58]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[58]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[59]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[59]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[60]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[60]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[63]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[63]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[64]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[64]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[65]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[65]),]$birth <- 2014

miss_reports[which(miss_reports$id == ids_list[66]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[66]),]$birth <- 2014

miss_reports[which(miss_reports$id == ids_list[67]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[67]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[68]),] 
#size_t = 1 at year_t = 2014
miss_reports[which(miss_reports$id == ids_list[68]),]$birth <- 2014


miss_reports[which(miss_reports$id == ids_list[69]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[69]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[70]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[70]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[71]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[71]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[72]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[72]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[73]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[73]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[74]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[74]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[75]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[75]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[76]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[76]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[77]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[77]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[78]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[78]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[79]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[79]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[80]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[80]),]$birth <- 2013

##there is a kinda big gap between IDs

miss_reports[which(miss_reports$id == ids_list[81]),] 
#size_t = 2 at year_t = 2011
#2 tillers at birth year 
miss_reports[which(miss_reports$id == ids_list[81]),]$birth <- 2011


miss_reports[which(miss_reports$id == ids_list[82]),] 
#size_t = 2 at year_t = 2011
#2 tillers at birth year 
miss_reports[which(miss_reports$id == ids_list[82]),]$birth <- 2011


miss_reports[which(miss_reports$id == ids_list[83]),] 
#size_t = 1 at year_t = 2011
miss_reports[which(miss_reports$id == ids_list[83]),]$birth <- 2011


miss_reports[which(miss_reports$id == ids_list[84]),] 
#size_t = 1 at year_t = 2011
miss_reports[which(miss_reports$id == ids_list[82]),]$birth <- 2011



miss_reports[which(miss_reports$id == ids_list[85]),] 
#size_t = 1 at year_t = 2011
miss_reports[which(miss_reports$id == ids_list[85]),]$birth <- 2011


miss_reports[which(miss_reports$id == ids_list[86]),] 
#size_t = 1 at year_t = 2011
miss_reports[which(miss_reports$id == ids_list[86]),]$birth <- 2011


miss_reports[which(miss_reports$id == ids_list[87]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[87]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[88]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[88]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[89]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[89]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[90]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[90]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[91]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[91]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[92]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[92]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[93]),] 
#bad stuff is happing with this plant and i cant figure it out from 
#raw data so i am going to drop it 
miss_reports <- miss_reports[-which(miss_reports$id == ids_list[93]), ]



miss_reports[which(miss_reports$id == ids_list[94]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[94]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[95]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[95]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[96]),] 
#this looks like 2 plants and i cant tell what is what so i am just going to drop 
miss_reports <- miss_reports[-which(miss_reports$id == ids_list[96]), ]



miss_reports[which(miss_reports$id == ids_list[97]),] 
#same thing, 2 plant chaos so droping 
miss_reports <- miss_reports[-which(miss_reports$id == ids_list[97]), ]


miss_reports[which(miss_reports$id == ids_list[98]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[98]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[99]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[99]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[100]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[100]),]$birth <- 2013


#20_403- 20_409 are were all born in 2013 and were changed to 2014 in 2017
miss_reports[which(miss_reports$id == ids_list[101]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[101]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[102]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[102]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[103]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[103]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[104]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[104]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[105]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[105]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[106]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[106]),]$birth <- 2013



miss_reports[which(miss_reports$id == ids_list[107]),] 
#2 plants, need to drop 
miss_reports <- miss_reports[-which(miss_reports$id == ids_list[107]), ]


miss_reports[which(miss_reports$id == ids_list[108]),] 
#size_t = 1 at year_t = 2012
miss_reports[which(miss_reports$id == ids_list[108]),]$birth <- 2012


miss_reports[which(miss_reports$id == ids_list[109]),] 
#size_t = 1 at year_t = 2013
miss_reports[which(miss_reports$id == ids_list[109]),]$birth <- 2013


miss_reports[which(miss_reports$id == ids_list[110]),] 
#size_t1 = 1 at year_t1 = 2018
miss_reports[which(miss_reports$id == ids_list[110]),]$birth <- 2018



miss_reports[which(miss_reports$id == ids_list[111]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[111]),]$birth <- 2018


miss_reports[which(miss_reports$id == ids_list[112]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[112]),]$birth <- 2018


miss_reports[which(miss_reports$id == ids_list[113]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[113]),]$birth <- 2018



miss_reports[which(miss_reports$id == ids_list[114]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[114]),]$birth <- 2018


miss_reports[which(miss_reports$id == ids_list[115]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[115]),]$birth <- 2018



miss_reports[which(miss_reports$id == ids_list[116]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[116]),]$birth <- 2018



miss_reports[which(miss_reports$id == ids_list[117]),] 
#size_t = 1 at year_t = 2018
miss_reports[which(miss_reports$id == ids_list[117]),]$birth <- 2018


fixed_data<- full_join(diff_plants, miss_reports)
##if ID was updated, replace
fixed_data[!is.na(fixed_data$id_update),"id"]<-fixed_data$id_update[!is.na(fixed_data$id_update)]

##write this out to be merged back with the big data
write.csv(fixed_data,"data prep/birth_year_fixes.csv")
