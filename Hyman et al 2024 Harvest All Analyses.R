#------------------------------------------------------------------------------#
################################### Description ################################
#------------------------------------------------------------------------------#
# This file runs all analyses in Hyman et al 2024: Take it or leave it: 
# influence of environmental and management variables on harvest and discard 
# rates in a multispecies fishery; a case study with gag 
# (Mycteroperca microlepis)
#
# Required files can be downloaded directly from github.com via read.csv().
#
# The code below is annotated to explain to the user what each component 
# line executes.
#
# This file was written by A. Challen Hyman, PhD, on July 29th, 2024
#
#------------------------------------------------------------------------------#
##################################### Libraries ################################
#------------------------------------------------------------------------------#
## Syntax packages
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))

## Visualization packages
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(xtable))
suppressMessages(library(kableExtra))
options(width = 80)
suppressMessages(library(grid))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))

## Modeling
suppressMessages(library(brms))
suppressMessages(library(forecast))

## File reading
suppressMessages(library(readxl))
#------------------------------------------------------------------------------#
##################################### Housekeeping #############################
#------------------------------------------------------------------------------#
## Clear out old files in R
rm(list=ls(all=TRUE)) 

## Set working directory (set to wherever you store the files)
setwd("C:\\Users\\ichal\\Documents\\Hyman gag grouper models\\Hyman et al 2024b")

## Add helper functions from source file on Github
source("https://raw.githubusercontent.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/R-Files/Hyman%20et%20al%202024%20Harvest%20and%20Discard%20helper%20functions.R")

## set themes for plotting
theme_set(My_theme())

#------------------------------------------------------------------------------#
##################################### Load data ################################
#------------------------------------------------------------------------------#
## All MRIP data for analysis
### Note: for aggregation files see https://github.com/ChallenHymanPhD/Hyman-et-at-2024-Effort
All_Data <- read.csv("https://raw.githubusercontent.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/Data-Files/MRIP%20Effort%20and%20Catch%20data.csv")

## Mortality estimates from Sauls 2014 (https://doi.org/10.1016/j.fishres.2013.10.008)
Mortality <- read.csv("https://raw.githubusercontent.com/ChallenHymanPhD/Hyman-et-al-2024-Harvest/Data-Files/Sauls%202014%20Mortality.csv") 

## Effort model from Hyman et al 2024
download.file("https://raw.githubusercontent.com/ChallenHymanPhD/Hyman-et-at-2024-Effort/Data-Files/Effort_Model.rds","Effort_Model.rds", method = 'curl')
Effort_Model <- readRDS("Effort_Model.rds")

#------------------------------------------------------------------------------#
################################### Format data ################################
#------------------------------------------------------------------------------#
All_Data <- All_Data%>%group_by(year, Region)%>%                                ## Create lagged index of abundance for Harvest and discard models
  mutate(Index_Gag = CPUE_Gag_lag,
         CPUE_Gag_lag_annual = mean(CPUE_Gag_lag, na.rm = T))                   ## Mean is used to replace NA values


All_Data$Date <- as.Date(All_Data$Date)                                         ## Format dates
All_Data$Income <- All_Data$Income/1e4                                          ## Income variable is only for effort model (see https://doi.org/10.1016/j.fishres.2024.107136)

All_Data$Index_Gag <- All_Data$Index_Gag+0.001                                  ## Add small constant to lagged index to avoid non-finite values after log transformation

## Replace NA index values with region-and year-specific mean CPUE (only 5 instances total)
All_Data$Index_Gag[which(is.na(All_Data$Index_Gag))] <- 
  All_Data$CPUE_Gag_lag_annual[which(is.na(All_Data$Index_Gag))]

## Small correction to dataset (Gag season open in 4 county region in Big Bend in 2012 in addition to 2013-2022:
## https://www.woodsnwater.net/blogs/news/faq-2012-recreational-gulf-gag-grouper-season-changes#:~:text=Harvest%20and%20possession%20of%20gag,April%201%20through%20June%2030)
All_Data$M_Gag[which(All_Data$Region=="Panhandle" & 
                       All_Data$Date %in% as.Date(c("2012-04-01", 
                                                    "2012-05-01", 
                                                    "2012-06-01")))] <- 1

All_Data <- All_Data[which(All_Data$year<2024 & All_Data$year>2003),]           ## Subset data to desired years for analysis

## Removes four crazy outliers
All_Data <- All_Data[-which(All_Data$Region=="Panhandle" & 
                              All_Data$Date %in% c(as.Date("2005-04-01"), 
                                                   as.Date("2013-10-01"), 
                                                   as.Date("2004-10-01"), 
                                                   as.Date("2010-11-01"))),]


#------------------------------------------------------------------------------#
################################### Run models #################################
#------------------------------------------------------------------------------#
## Set training and testing data
set.seed(9876)
Test_data <- sample(1:nrow(All_Data), 50, replace = F)

Test <- All_Data[Test_data,]
Train <- All_Data[-Test_data,]

#------------------------------------------------------------------------------#
## Generally weakly informative priors for both models
Hyman_prior <- c(prior(normal(0,10), class = Intercept),
                 prior(normal(0,10), class = b),
                 prior(normal(0,10), class = Intercept, dpar = shape),
                 prior(normal(0,10), class = Intercept, dpar = hu),
                 prior(normal(0,10), class = b, dpar = hu))

#------------------------------------------------------------------------------#
## Harvest Model
Harvest_Model <- brm(                                                       
  bf(Harvest_Gag ~                                                              ## Models mean of gamma distribution
       Region*log(Index_Gag) +                                                  ## CPUE index of abundance in a given region and month from previous year
       Region*log(Juveniles) +                                                  ## Juvenile index of abundance (annual resolution for whole state)
       Region*sin1 +                                                            ## Harmonic terms (sin1 = annual sinusoidal term, sin2 = semi-annual sinusoidal term, etc) 
       Region*cos1 + 
       Region*sin2 + 
       Region*cos2 + 
       Region*SST +                                                             ## Sea-surface temperature
       Region*M_Gag +                                                           ## Fraction of month open to gag harest
       Region*M_RS +                                                            ## Fraction of month open to red snapper harest
       Region*log(Trips),                                                       ## log-effort                              
     hu ~ Region*ceiling(M_Gag),                                                ## Models hurdle component: Region-specific binary term denoting whether gag was open in a given month/region (at least one day open)
     shape ~ Region*sin1 +                                                      ## Models gamma shape parameter as functin of harmonic terms and gag index of abundance                                               
       Region*cos1 + 
       Region*sin2 + 
       Region*cos2 +
       Region*log(Index_Gag)), data = Train,                            
  family = hurdle_gamma(), 
  chains = 4, prior = Hyman_prior)

#------------------------------------------------------------------------------#
## Discard Model
Discard_Model <- brm(                                                       
  bf(Discard_Gag ~                                                              ## Models mean of gamma distribution
       Region*log(Index_Gag) +                                                  ## CPUE index of abundance in a given region and month from previous year
       Region*log(Juveniles) +                                                  ## Juvenile index of abundance (annual resolution for whole state)
       Region*sin1 +                                                            ## Harmonic terms (sin1 = annual sinusoidal term, sin2 = semi-annual sinusoidal term, etc) 
       Region*cos1 + 
       Region*sin2 + 
       Region*cos2 + 
       Region*SST +                                                             ## Sea-surface temperature
       Region*M_Gag +                                                           ## Fraction of month open to gag harvest
       Region*M_RS +                                                            ## Fraction of month open to red snapper harvest
       Region*log(Trips),                                                       ## log-effort                                                              
     hu ~ Region*log(Index_Gag),                                                ## Models hurdle component: Region-specific gag index of abundance    
     shape ~ Region*sin1 +                                                      ## Models gamma shape parameter as function of harmonic terms and gag index of abundance                                               
       Region*cos1 + 
       Region*sin2 + 
       Region*cos2 +
       Region*log(Index_Gag)), data = Train,                            
  family = hurdle_gamma(), 
  chains = 4, prior = Hyman_prior)

#------------------------------------------------------------------------------#
################################### Table 3 ####################################
#------------------------------------------------------------------------------#
## Calculate average gag weight for each month and region
Weight_Gag <- All_Data%>%
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(Weight_Gag, na.rm = T))%>%.[order(.$Region, .$month),3]%>%
  unlist()%>%as.vector()

## Combine discard mortality and aerage weight estimates into single data frame and format to create region-specific columns
Table_3 <- cbind(Mortality, Weight_Gag)%>%.[,c(2,1,4,3)]%>%
  pivot_wider(names_from = "Region", values_from = c(Weight_Gag, Mortality))
Table_3 <- Table_3[,c(1,2,4,3,5)] 

## Convert gag weight from Kg to lbs
Table_3$Weight_Gag_Panhandle <- Table_3$Weight_Gag_Panhandle*2.205
Table_3$Weight_Gag_Peninsula <- Table_3$Weight_Gag_Peninsula*2.205
Table_3[,2:5] <- round(Table_3[,2:5],3)
colnames (Table_3) <- c("Month", "Panhandle weight (lbs)", 
                        "Panhandle mortality (proportion)", 
                        "Peninsula weight (lbs)", 
                        "Peninsula mortality (proportion)")

print(xtable(Table_3),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)

#------------------------------------------------------------------------------#
################################### Figure 4 ###################################
#------------------------------------------------------------------------------#
## Out-of-sample Harvest predictions
Pred <- predict(Harvest_Model, Test, probs = c(0.1, 0.9), robust = T)%>%        ## Generate median out-of-sample estimates and upper/lower 80% CI using testing set
  as.data.frame()

## Append posterior predictions to testing dataset
Test$Median <- (Pred$Estimate)
Test$Min  <- (Pred$Q10)
Test$Max  <- (Pred$Q90)

## Determine if each withheld values are captured in OOS CI
Test$Status <- "Outside CI"
Test$Status[which(Test$Harvest_Gag >= Test$Min & 
                    Test$Harvest_Gag <= Test$Max)] <- "Within CI"

## Re-order by increasing value
Test_plot <- Test[order(Test$Harvest_Gag),]
Test_plot <-Test_plot%>%group_by(Region)%>%
  mutate(Rank = order(Harvest_Gag))
ggplot(Test_plot)+
  geom_errorbar(aes(x = Rank, ymin = Min, ymax = Max, col = Status), 
                show.legend = T, lwd = 1)+
  geom_point(aes(x = Rank, y = Median, col = Status), size = 3)+
  geom_point(aes(x = Rank, y = Harvest_Gag), alpha = 0.5, size = 3)+
  facet_wrap(~Region, scales = "free", ncol = 1)+
  scale_size(range=c(1,4))+
  scale_y_continuous(labels = function(x) format(x, scientific = F))+
  theme(strip.text.x = element_text(size = 14, 
                                    margin = margin(0.25,0,0.25,0, "cm")),
        strip.text.y = element_text(size = 14, 
                                    margin = margin(0,0.2,0,0.2, "cm")))+
  theme(legend.position = "top")+
  ylab("Harvests (number)")+labs(color = "Withheld observation:")

#------------------------------------------------------------------------------#
################################### Figure 5 ###################################
#------------------------------------------------------------------------------#
Pred <- predict(Harvest_Model, Train, probs = c(0.1, 0.9), robust = T)%>%       ## Generate median in-sample estimates and upper/lower 80% CI using training set
  as.data.frame()

## Append posterior predictions to Training dataset
Train$Median <- (Pred$Estimate)
Train$Min  <- (Pred$Q10)
Train$Max  <- (Pred$Q90)

## Plot results
ggplot(Train)+
  geom_point(aes(x = as.Date(Date), y = (Harvest_Gag)), alpha = 0.5, size = 3)+
  geom_ribbon(aes(x = Date, ymin = (Min), ymax = (Max)), fill ="#4575B4", 
              alpha = 0.3, show.legend = T)+
  geom_line(aes(x = as.Date(Date), y = (Harvest_Gag)))+
  geom_line(aes(x = Date, y = (Median)),col = "#4575B4",lwd = 1,  show.legend = T)+
  facet_wrap(~Region, scales = "free_y", ncol = 1)+ylab("Harvest (number)")+
  xlab("Date")+scale_size(range=c(1,4))+theme(legend.position = "top")+
  scale_y_continuous(labels = function(x) format(x, scientific = F))
#------------------------------------------------------------------------------#
################################### Table 4 ####################################
#------------------------------------------------------------------------------#
### Run function to get tables
All_Harvest_tables <- Hyman_tables(Harvest_Model, Type = "Harvest")
All_Discard_tables <- Hyman_tables(Discard_Model, Type = "Discard")
colnames(All_Harvest_tables[[1]]) <- c("Coefficient", "10%", "50%", 
                                       "90%","Coefficient", "10%", "50%", "90%")
rownames(All_Harvest_tables[[1]]) <- rep("", nrow(All_Harvest_tables[[1]]))
print(xtable(All_Harvest_tables[[1]]),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)
#------------------------------------------------------------------------------#
################################### Figure 6 ###################################
#------------------------------------------------------------------------------#
## Out-of-sample Discards predictions
Pred <- predict(Discard_Model, Test, probs = c(0.1, 0.9), robust = T)%>%        ## Generate median out-of-sample estimates and upper/lower 80% CI using testing set
  as.data.frame()

## Append posterior predictions to testing dataset
Test$Median <- (Pred$Estimate)
Test$Min  <- (Pred$Q10)
Test$Max  <- (Pred$Q90)

## Determine if each withheld values are captured in OOS CI
Test$Status <- "Outside CI"
Test$Status[which(Test$Discard_Gag >= Test$Min & 
                    Test$Discard_Gag <= Test$Max)] <- "Within CI"

## Re-order by increasing value
Test_plot <- Test[order(Test$Discard_Gag),]
Test_plot <-Test_plot%>%group_by(Region)%>%
  mutate(Rank = order(Discard_Gag))
ggplot(Test_plot)+
  geom_errorbar(aes(x = Rank, ymin = Min, ymax = Max, col = Status), 
                show.legend = T, lwd = 1)+
  geom_point(aes(x = Rank, y = Median, col = Status), size = 3)+
  geom_point(aes(x = Rank, y = Discard_Gag), alpha = 0.5, size = 3)+
  facet_wrap(~Region, scales = "free", ncol = 1)+
  scale_size(range=c(1,4))+
  scale_y_continuous(labels = function(x) format(x, scientific = F))+
  theme(strip.text.x = element_text(size = 14, 
                                    margin = margin(0.25,0,0.25,0, "cm")),
        strip.text.y = element_text(size = 14, 
                                    margin = margin(0,0.2,0,0.2, "cm")))+
  theme(legend.position = "top")+ylab("Discardss (number)")+
  labs(color = "Withheld observation:")

#------------------------------------------------------------------------------#
################################### Table 5 ####################################
#------------------------------------------------------------------------------#
colnames(All_Discard_tables[[1]]) <- c("Coefficient", "10%", "50%", "90%",
                                       "Coefficient", "10%", "50%", "90%")
rownames(All_Discard_tables[[1]]) <- rep("", nrow(All_Discard_tables[[1]]))
print(xtable(All_Discard_tables[[1]]),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)
#------------------------------------------------------------------------------#
################################### Figure 7 ###################################
#------------------------------------------------------------------------------#
Pred <- predict(Discard_Model, Train, probs = c(0.1, 0.9), robust = T)%>%       ## Generate median in-sample estimates and upper/lower 80% CI using training set
  as.data.frame()

## Append posterior predictions to Training dataset
Train$Median <- (Pred$Estimate)
Train$Min  <- (Pred$Q10)
Train$Max  <- (Pred$Q90)

## Plot results
ggplot(Train)+
  geom_point(aes(x = as.Date(Date), y = (Discard_Gag)), alpha = 0.5, size = 3)+
  geom_ribbon(aes(x = Date, ymin = (Min), ymax = (Max)), fill ="#4575B4", 
              alpha = 0.3, show.legend = T)+
  geom_line(aes(x = as.Date(Date), y = (Discard_Gag)))+
  geom_line(aes(x = Date, y = (Median)),col = "#4575B4",lwd = 1,  show.legend = T)+
  facet_wrap(~Region, scales = "free_y", ncol = 1)+ylab("Discards (number)")+
  xlab("Date")+scale_size(range=c(1,4))+theme(legend.position = "top")+
  scale_y_continuous(labels = function(x) format(x, scientific = F))

#------------------------------------------------------------------------------#
############################### Run simulation #################################
#------------------------------------------------------------------------------#
## Set starting values
Year <- 2024
Gag_Start <- c("No Season", "June 01", "September 01", "November 01")           ## Gag season start dates for counterfactual scenarios (including baseline)
ACT <-  333000*0.9*2.12                                                         ## ACL converted from SRFS units to MRIP-FES units with 10% buffer for ACT
RS_Season <- c("June 01", "July 31")                                            ## 2024 Florida Red Snapper season (https://myfwc.com/fishing/saltwater/recreational/snappers/#:~:text=Gulf%20Recreational%20Red%20Snapper,30.)
RS_Season2 <- c("September 01", "November 30")                                  ## Dates red snapper was reopened in Florida
RG_Season <- c("Jan 01", "June 30")                                             ## 2024 red grouper season (https://myfwc.com/news/all-news/red-grouper-524/)
RG_Season2 <- NULL                                                              ## Red grouper season did not re-open
N <- 1000                                                                       ## Number of iterations to run simulation

#------------------------------------------------------------------------------#
## Create fixed counterfactual data
Data <- All_Data[order(All_Data$Region, All_Data$Date),]                        ## Order Data
Start <- paste0(Year,'-01-01')                                                  ## Start
End <- paste0(Year,'-12-01')                                                    ## End

## Time and Region variables
Month <- rep(seq.Date(as.Date(Start), as.Date(End), by = 'month'),2)
Region <- sort(rep(unique(Data$Region),12))

## Environmental variables
## Sea-surface temperature
SST <- Data%>%
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(SST))%>%.[order(.$Region, .$month),3]%>%unlist()%>%as.vector()

## Wind (used effort model in Hyman et al. 2024)
Wind <-Data%>%
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(Wind))%>%.[order(.$Region, .$month),3]%>%unlist()%>%as.vector()

## Economic variables (used effort model in Hyman et al. 2024)
FIR <- Data$FIR[(nrow(Data)-11):nrow(Data)]
Sales <- Data$Sales[(nrow(Data)-11):nrow(Data)]

## Harmonic regression terms
## Set periods manually
per <- 365                                                                      ## Period 1 is annual (365 days in a year)
per2 <- 182.5                                                                   ## Period 2 is semi-annual (182.5 days every six months)

#------------------------------------------------------------------------------#
## Index variables
## Gag
CPUE_Gag <- Data%>%                                                             
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(CPUE_Gag))%>%.[order(.$Region, .$month),3]%>%unlist()%>%as.vector()

## Red Grouper
CPUE_RG  <- Data%>%                                                             
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(CPUE_RG))%>%.[order(.$Region, .$month),3]%>%unlist()%>%as.vector()

## Red Snapper
CPUE_RS  <- Data%>%                                                             
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(CPUE_RS))%>%.[order(.$Region, .$month),3]%>%unlist()%>%as.vector()  

## Juvenile gag index
Juveniles <- ((Data$Juveniles))[which(Data$year == Year-1)]                    

## Gag Weight
Weight_Gag <- Data%>%
  filter(year >2020)%>%
  group_by(month, Region)%>%
  summarize(mean(Weight_Gag, na.rm = T))%>%
  .[order(.$Region, .$month),3]%>%unlist()%>%as.vector()
Weight_Gag[which(is.na(Weight_Gag))] <- 0                                       ## Replace NA values with 0s (because gag has not been harvested in any closed month in last 3 years anyway)

## Gag seasonal live release mortality
LRM_Gag <- Mortality$Mortality

#------------------------------------------------------------------------------#
## Fixed management variables
## Red Snapper
M_RS <- Season_Month_Fraction(Year = Year,                                      ## Creates fraction of month open to red snapper harvest
                              Season = RS_Season,
                              Additional = RS_Season2,
                              Additional_KOD = "weekend")[[1]]$Open
S_RS <- Season_Month_Fraction(Year = Year,                                      ## Calculates number of days in red snapper season
                              Season = RS_Season,
                              Additional = RS_Season2,
                              Additional_KOD = "weekend")[[2]]
S_RSL <- M_RS*log(S_RS)                                                         ## Red snapper effort concentration variable

## Red Grouper (for CPUE variable in Hyman et al. 2024 effort model)
M_RG <- Season_Month_Fraction(Year = Year,                                      ## Creates fraction of month open to red grouper harvest
                              Season = RG_Season,
                              Additional = NULL,
                              Additional_KOD = "weekend")[[1]]$Open 
## Create data frames to hold all simulations
Complete_Harvest <- NULL
Complete_Discards <- NULL
Complete_Dead_Discards <- NULL
Complete_Landings <- NULL
for (Season_Start in Gag_Start){                                                ## For each of of the four season scenarios (including the baseline)
  for (d in 60:1){                                                              ## Begin at 60 day seasons and iteratively reduce until landings are at or below threshold
    if(Season_Start == "No Season"){                                            ## For baseline "no season" scenario, set fraction of month open to gag at 0
      M_Gag <- 0; S_Gag <- 1; 
      Gag_Season <- "No Season"} 
    else {
      Gag_Season <- c(Season_Start,paste(format(as.Date(Season_Start, "%B %d")  ## If not the baseline, add number of days in season (d) to season start value to get season length
                                                + d, "%B %d")))
      ## Obtain gag management variables based on simulated season
      M_Gag <- Season_Month_Fraction(Year = Year,                               ## Fraction of each month open to gag harvest
                                     Season = Gag_Season,Additional = NULL)[[1]]$Open        
      S_Gag <- Season_Month_Fraction(Year = Year,                               ## Season duration variable
                                     Season = Gag_Season,Additional = NULL)[[2]]
    }
    ## Construct counterfactual data frame
    Counterfactual <- data.frame(Wind = Wind,                                   ## For effort model
                                 SST = SST,                                     ## For harvest/discard models
                                 Date = Month,                                  ## For all models
                                 Region = Region,                               ## For all models
                                 Juveniles = Juveniles,                         ## For harvest/discard models
                                 FIR=FIR,                                       ## For effort model
                                 LRM = LRM_Gag,                                 ## To calculated dead discards
                                 Index_Gag = CPUE_Gag+0.001,                    ## For harvest/discard models; note small constant added to avoid 0s in log transformation
                                 CPUE_Gag_lag = CPUE_Gag*M_Gag,                 ## For effort model
                                 CPUE_RG_lag = CPUE_RG*M_RG,                    ## For effort model
                                 CPUE_RS_lag = CPUE_RS*M_RS,                    ## For effort model
                                 Weight_Gag = Weight_Gag,                       ## To convert monthly harvest in numbers to landings in pounds
                                 M_RG = M_RG,                                   ## For all models
                                 M_RS = M_RS,                                   ## For all models
                                 S_RSL = S_RSL,                                 ## For all models
                                 M_Gag = M_Gag,                                 ## For all models
                                 S_GagL = M_Gag*log(S_Gag),                     ## For all models
                                 Sales = Sales,                                 ## For effort model
                                 Scenario = Gag_Season)%>%
      mutate(Time = as.POSIXlt(Date)$yday,
             sin1 = sin(2*pi/per*Time),                                         ## Add harmonic terms (for all models)
             cos1 = cos(2*pi/per*Time),
             sin2 = sin(2*pi/per2*Time),
             cos2 = cos(2*pi/per2*Time)
      )
    
    ## Step 1) Extract posterior predictions for Effort
    Effort_preds <- posterior_epred(Effort_Model, Counterfactual, ndraws = N)   ## Compute reef effort expectation based on conditions in simulation
    
    ## Step 2) Extract model coefficients
    Harvest_terms <- as_draws_matrix(Harvest_Model)%>%as.data.frame()           ## Harvest model terms
    omega <- Harvest_terms[,-c(grep("b_shape_", colnames(Harvest_terms)),grep("b_hu_", colnames(Harvest_terms)))]%>%.[-ncol(.)]%>%.[,-ncol(.)]
    gamma <- Harvest_terms[,grep("b_hu_", colnames(Harvest_terms))]
    nu <- Harvest_terms[,grep("b_shape_", colnames(Harvest_terms))]
    
    Discard_terms <- as_draws_matrix(Discard_Model)%>%as.data.frame()           ## Discard model terms
    psi <- Discard_terms[,-c(grep("b_shape_", colnames(Discard_terms)),grep("b_hu_", colnames(Discard_terms)))]%>%.[-ncol(.)]%>%.[,-ncol(.)]
    zeta <- Discard_terms[,grep("b_hu_", colnames(Discard_terms))]
    xi <- Discard_terms[,grep("b_shape_", colnames(Discard_terms))]
    
    ## Data frames to hold values
    Discard_preds <- matrix(NA, nrow = N, ncol = nrow(Counterfactual))                  
    Harvest_preds <- matrix(NA, nrow = N, ncol = nrow(Counterfactual))
    Dead_Discard_preds <- matrix(NA, nrow = N, ncol = nrow(Counterfactual))    
    Removals_preds <- matrix(NA, nrow = N, ncol = nrow(Counterfactual))  
    Landings_preds <- matrix(NA, nrow = N, ncol = nrow(Counterfactual))  
    
    ## Loop to pull predictions from correlated model draws
    for (i in 1:N){                                                             ## For each of the 1000 simulations
      ## Model matrices for predictions
      Test_iter <- Counterfactual                                               ## Create temporary data frame for effort updating
      Test_iter$Trips <- (Effort_preds[i,])                                     ## Update design matrix with effort estimates from model predictions iteratively
      
      #--------------------------- Harvest ------------------------------------#
      omega_MM <- model.matrix(~Region*log(Index_Gag) +                         ## Design matrix for mean component of gag harvest based on fixed values, and set management conditions
                                 Region*log(Juveniles) + 
                                 Region*sin1 + 
                                 Region*cos1 + 
                                 Region*sin2 + 
                                 Region*cos2 + 
                                 Region*SST + 
                                 Region*M_Gag +
                                 Region*M_RS + 
                                 Region*log(Trips), data = Test_iter)
      gamma_MM <- model.matrix(~Region*M_Gag, data = Test_iter)                 ## Design matrix for hurdle component of gag HPUE based on fixed values, and set management conditions
      nu_MM <- model.matrix(~Region*sin1 + 
                              Region*cos1 + 
                              Region*sin2 + 
                              Region*cos2 +
                              Region*log(Index_Gag), data = Test_iter)          ## Design matrix for shape component of gag HPUE based on fixed values
      lambda <- exp(as.matrix(omega_MM)%*%t(as.matrix(omega[i,])))              ## Compute lambda (harvest nonzero mean) using posterior scan of coefficients and design matrix
      phi <- exp(as.matrix(nu_MM)%*%t(as.matrix(nu[i,])))                       ## Compute phi (harvest shape) using posterior scan of coefficients and design matrix
      delta <- plogis(as.matrix(gamma_MM)%*%t(as.matrix(gamma[i,])))            ## Compute delta (harvest hurdle) using posterior scan of coefficients and design matrix
      Harvest_preds[i,] <- (1-delta)*lambda;                                    ## Compute gag harvest expectation
      
      #--------------------------- Discards -----------------------------------#
      eta_MM <- model.matrix(~Region*log(Index_Gag) +                           ## Create discard model design matrix
                               Region*log(Juveniles) + 
                               Region*sin1 + 
                               Region*cos1 + 
                               Region*sin2 + 
                               Region*cos2 + 
                               Region*SST + 
                               Region*M_Gag +
                               Region*M_RS+
                               Region*log(Trips), data = Test_iter)                                       
      theta_MM <- model.matrix(~Region*log(Index_Gag), data = Test_iter)        ## Design matrix for hurdle component of gag discards based on fixed values, and set management conditions
      eta <- exp(as.matrix(eta_MM)%*%t(as.matrix(psi[i,])))                     ## Compute eta (discard nonzero mean) using posterior scan of coefficients and design matrix
      theta <- plogis(as.matrix(theta_MM)%*%t(as.matrix(zeta[i,])))             ## Compute theta (discard hurdle) using posterior scan of coefficients and design matrix
      Discard_preds[i,] <- (1-theta)*eta;                                       ## Compute gag discard expectation
      
      
      Dead_Discard_preds[i,] <- Discard_preds[i,]*Test_iter$LRM                 ## Compute expected dead discards as product of live releases and discard mortality
      Removals_preds[i,] <- Harvest_preds[i,] + Dead_Discard_preds[i,]          ## Compute expected total removals as sum of expected harvest and expected dead discards
      Landings_preds[i,] <- Harvest_preds[i,]*Test_iter$Weight_Gag*2.205/1.12   ## Compute numbers to landings in pounds and adjust for gutted weight 
    }
    
    ## Calculate pseudo-posterior distribution of cumulative landings
    Landings <- Landings_preds%>%t()%>%as.data.frame()
    Landings$Date <- as.Date(Counterfactual$Date)
    Landings$Region <-Counterfactual$Region
    Landings <- Landings%>%group_by(Date)%>%                                    ## Sum across regions
      summarize(across(1:N, list(sum)))
    Landings_cum <- Landings%>%                                                 ## compute cumulative sum for each of the 1000 scans
      summarize(across(2:(N+1), list(cumsum)))
    Landings_cum$Date <- as.Date(rep(Month[1:12],1))                            ## Add months (honestly not important)
    
    ## Obtain median and 80% CI for cumulative landings sums
    Landings_cum$Q10 <- apply(Landings_cum[,1:(N)], 1, function(x){quantile(x, 0.1, na.rm = T)})
    Landings_cum$Q50 <- apply(Landings_cum[,1:(N)], 1, function(x){quantile(x, 0.5, na.rm = T)})
    Landings_cum$Q90 <- apply(Landings_cum[,1:(N)], 1, function(x){quantile(x, 0.9, na.rm = T)})
    if(Landings_cum$Q50[12] <= ACT || Season_Start == "No Season"){             ## If median total landings fall below ACT, stop simulation for that season start date
      break                                                                     ## Otherwise reduce season by one day and repeat
    } 
  }
  
  #------------------------ Format estimates ----------------------------------#
  ## For each matrix of posterior predictive quantities, convert to data frame 
  ## and add counterfactual information on scenario, date, and region
  
  if(Season_Start == "No Season"){                                              ## For baseline scenario, set name without number of days in season (because it will be 0)
    Season_scenario <- "No season (baseline)"
  } else {
    Season_scenario <- paste0(Season_Start," (",d," days)")                     ## For all other scenarios, append number of days to season name
  }
  ## Harvest posterior predictive quantities
  Harvest <- Harvest_preds%>%t()%>%as.data.frame()
  Harvest$Date <- as.Date(Counterfactual$Date)                                  ## Add date
  Harvest$Region <-Counterfactual$Region                                        ## Add Region
  Harvest$Scenario <- Season_scenario                                           ## Add scenario name
  
  ## Discard posterior predictive quantities
  Discards <- Discard_preds%>%t()%>%as.data.frame()
  Discards$Date <- as.Date(Counterfactual$Date)                                 ## Add date
  Discards$Region <-Counterfactual$Region                                       ## Add Region
  Discards$Scenario <- Season_scenario                                          ## Add scenario name
  
  ## Dead discard posterior predictive quantities
  Dead_Discards <- Dead_Discard_preds%>%t()%>%as.data.frame()
  Dead_Discards$Date <- as.Date(Counterfactual$Date)                            ## Add date
  Dead_Discards$Region <-Counterfactual$Region                                  ## Add Region
  Dead_Discards$Scenario <- Season_scenario                                     ## Add scenario name
  
  ## Landings posterior predictive quantities
  Landings <- Landings_preds%>%t()%>%as.data.frame()
  Landings$Date <- as.Date(Counterfactual$Date)                                 ## Add date
  Landings$Region <-Counterfactual$Region                                       ## Add Region
  Landings$Scenario <- Season_scenario                                          ## Add scenario name
  
  #---------------- Bind to external objects to save --------------------------#
  Complete_Harvest <- rbind(Complete_Harvest, Harvest)
  Complete_Discards <- rbind(Complete_Discards, Discards)
  Complete_Dead_Discards <- rbind(Complete_Dead_Discards, Dead_Discards)
  Complete_Landings <- rbind(Complete_Landings, Landings)
}

#------------------------------------------------------------------------------#
Scenario_names <- unique(Complete_Dead_Discards$Scenario)
## Discards
All_Discard_differences <- NULL                                                 ## Stores all scan-wise differences in total discards
All_Discard_Scenarios <- NULL                                                   ## Stores monthly expected discard estimates for plotting
for (i in 2:length(Scenario_names)){                                            ## For each scenario (minus the baseline, which is used for comparisons)
  Baseline <- Complete_Discards[which(Complete_Dead_Discards$Scenario == Scenario_names[1]),] ## Extracts baseline data
  Scenario <- Complete_Discards[which(Complete_Dead_Discards$Scenario == Scenario_names[i]),] ## Extracts a given scenario's data
  Difference <- Scenario                                                        ## New data frame to store differences
  for (j in 1:N){                                                               ## For each of 1,000 posterior scans
    Difference[,j] <- Scenario[,j] - Baseline[,j]                               ## Calculate the difference between the scan from any given scenario and the baseline scenario (by month and region)
  }
  Panhandle <- Difference[which(Difference$Region=="Panhandle"),]               ## Subset for Panhandle
  Panhandle <- quantile(apply(Panhandle[,1:N],2,sum), c(0.1, 0.5, 0.9))         ## Sum all differences (for all regions and months for a given pairwise difference) and calculate summary statistics
  Panhandle <- t(as.data.frame(c(Panhandle, Scenario_names[i])))                ## Append Scenario name
  
  Peninsula <- Difference[which(Difference$Region=="Peninsula"),]               ## Subset for Peninsula
  Peninsula <- quantile(apply(Peninsula[,1:N],2,sum), c(0.1, 0.5, 0.9))         ## Sum all differences (for all regions and months for a given pairwise difference) and calculate summary statistics
  Peninsula <- t(as.data.frame(c(Peninsula, Scenario_names[i])))                ## Append Scenario name
  Difference <- rbind(Panhandle, Peninsula)%>%as.data.frame()                   ## Append Panhandle and Peninsula
  colnames(Difference) <- c("Q10", "Q50", "Q90", "Scenario")                    ## Column names
  Difference$Region <- c("Panhandle", "Peninsula")                              ## Region names
  All_Discard_differences <- rbind(All_Discard_differences, Difference)         ## Append to outside object to save information from loop                      
  ## New column for comparison between each scenario and baseline
  Scenario$Comparison <- "Season"                                               
  Baseline$Comparison <- "No season"
  ## Append each scenario's monthly estimates to the baseline for plotting
  Current_Scenario <- rbind(Scenario, Baseline)
  Current_Scenario$Scenario <- Scenario_names[i]                                ## Unique scenario identifier for plotting
  All_Discard_Scenarios <- rbind(All_Discard_Scenarios, Current_Scenario)       ## Append to outside object to save information from loop
}

All_Discard_differences <- as.data.frame(All_Discard_differences)               ## Convert to data frame

## Calculate summary statistics
All_Discard_differences$Q10 <- as.numeric(All_Discard_differences$Q10)
All_Discard_differences$Q50 <- as.numeric(All_Discard_differences$Q50)
All_Discard_differences$Q90 <- as.numeric(All_Discard_differences$Q90)

## Reorder factors
All_Discard_differences$Scenario <- as.factor(All_Discard_differences$Scenario)%>%
  factor(.,levels = unique(All_Discard_differences$Scenario))
All_Discard_differences$Category <- "Discards"                                  ## Unique discard identifier to distinguish from dead discards (below)

## Reorder factors
All_Discard_Scenarios$Scenario <- as.factor(All_Discard_Scenarios$Scenario)%>%
  factor(.,levels = unique(All_Discard_Scenarios$Scenario))

## Calculate summary statistics
All_Discard_Scenarios <- as.data.frame(All_Discard_Scenarios)
All_Discard_Scenarios$Q10 <- apply(All_Discard_Scenarios[,1:N], 1, function(x){quantile(x,0.1)})
All_Discard_Scenarios$Q50 <- apply(All_Discard_Scenarios[,1:N], 1, function(x){quantile(x,0.5)})
All_Discard_Scenarios$Q90 <- apply(All_Discard_Scenarios[,1:N], 1, function(x){quantile(x,0.9)})

#------------------------------------------------------------------------------#
## Dead Discards
All_Dead_Discard_differences <- NULL                                            ## Stores all scan-wise differences in total dead discards
All_Dead_Discard_Scenarios <- NULL                                              ## Stores monthly expected dead discard estimates for plotting
for (i in 2:length(Scenario_names)){                                            ## For each scenario (minus the baseline, which is used for comparisons)
  Baseline <- Complete_Dead_Discards[which(Complete_Dead_Discards$Scenario == Scenario_names[1]),] ## Extracts baseline data
  Scenario <- Complete_Dead_Discards[which(Complete_Dead_Discards$Scenario == Scenario_names[i]),] ## Extracts a given scenario's data
  Difference <- Scenario                                                        ## New data frame to store differences
  for (j in 1:N){                                                               ## For each of 1,000 posterior scans
    Difference[,j] <- Scenario[,j] - Baseline[,j]                               ## Calculate the difference between the scan from any given scenario and the baseline scenario (by month and region)
  }
  Panhandle <- Difference[which(Difference$Region=="Panhandle"),]               ## Subset for Panhandle
  Panhandle <- quantile(apply(Panhandle[,1:N],2,sum), c(0.1, 0.5, 0.9))         ## Sum all differences (for all regions and months for a given pairwise difference) and calculate summary statistics
  Panhandle <- t(as.data.frame(c(Panhandle, Scenario_names[i])))                ## Append Scenario name
  
  Peninsula <- Difference[which(Difference$Region=="Peninsula"),]               ## Subset for Peninsula
  Peninsula <- quantile(apply(Peninsula[,1:N],2,sum), c(0.1, 0.5, 0.9))         ## Sum all differences (for all regions and months for a given pairwise difference) and calculate summary statistics
  Peninsula <- t(as.data.frame(c(Peninsula, Scenario_names[i])))                ## Append Scenario name
  Difference <- rbind(Panhandle, Peninsula)%>%as.data.frame()                   ## Append Panhandle and Peninsula
  colnames(Difference) <- c("Q10", "Q50", "Q90", "Scenario")                    ## Column names
  Difference$Region <- c("Panhandle", "Peninsula")                              ## Region names
  All_Dead_Discard_differences <- rbind(All_Dead_Discard_differences, Difference)## Append to outside object to save information from loop                      
  ## New column for comparison between each scenario and baseline
  Scenario$Comparison <- "Season"                                               
  Baseline$Comparison <- "No season"
  ## Append each scenario's monthly estimates to the baseline for plotting
  Current_Scenario <- rbind(Scenario, Baseline)
  Current_Scenario$Scenario <- Scenario_names[i]                                ## Unique scenario identifier for plotting
  All_Dead_Discard_Scenarios <- rbind(All_Dead_Discard_Scenarios,Current_Scenario)## Append to outside object to save information from loop
}

All_Dead_Discard_differences <- as.data.frame(All_Dead_Discard_differences)     ## Convert to data frame

## Calculate summary statistics
All_Dead_Discard_differences$Q10 <- as.numeric(All_Dead_Discard_differences$Q10)
All_Dead_Discard_differences$Q50 <- as.numeric(All_Dead_Discard_differences$Q50)
All_Dead_Discard_differences$Q90 <- as.numeric(All_Dead_Discard_differences$Q90)

## Reorder factors
All_Dead_Discard_differences$Scenario <- as.factor(All_Dead_Discard_differences$Scenario)%>%
  factor(.,levels = unique(All_Dead_Discard_differences$Scenario))
All_Dead_Discard_differences$Category <- "Dead Discards"                        ## Unique discard identifier to distinguish from discards (above)

## Reorder factors
All_Dead_Discard_Scenarios$Scenario <- as.factor(All_Dead_Discard_Scenarios$Scenario)%>%
  factor(.,levels = unique(All_Dead_Discard_Scenarios$Scenario))

## Calculate summary statistics
All_Dead_Discard_Scenarios <- as.data.frame(All_Dead_Discard_Scenarios)
All_Dead_Discard_Scenarios$Q10 <- apply(All_Dead_Discard_Scenarios[,1:N], 1, function(x){quantile(x,0.1)})
All_Dead_Discard_Scenarios$Q50 <- apply(All_Dead_Discard_Scenarios[,1:N], 1, function(x){quantile(x,0.5)})
All_Dead_Discard_Scenarios$Q90 <- apply(All_Dead_Discard_Scenarios[,1:N], 1, function(x){quantile(x,0.9)})

## Combine discard differences and dead discard differences
All_differences <- rbind(All_Dead_Discard_differences, All_Discard_differences)

#------------------------------------------------------------------------------#
################################### Figure 8 ###################################
#------------------------------------------------------------------------------#
ggplot(All_Discard_Scenarios)+
  geom_ribbon(aes(Date, ymin = Q10, ymax = Q90, fill = Comparison), alpha = 0.3)+
  geom_line(aes(Date, y= Q50, col = Comparison), lwd = 1)+
  scale_x_date(date_breaks = "1 month", date_labels="%b")+
  facet_grid(rows = vars(Region), cols = vars(Scenario), scales = "free")+
  My_theme()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Expected gag discards (numbers)")

#------------------------------------------------------------------------------#
################################### Figure 9 ###################################
#------------------------------------------------------------------------------#
ggplot(All_Dead_Discard_Scenarios)+
  geom_ribbon(aes(Date, ymin = Q10, ymax = Q90, fill = Comparison), alpha = 0.3)+
  geom_line(aes(Date, y= Q50, col = Comparison), lwd = 1)+
  scale_x_date(date_breaks = "1 month", date_labels="%b")+
  facet_grid(rows = vars(Region), cols = vars(Scenario), scales = "free")+
  My_theme()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  ylab("Expected gag discards (numbers)")

#------------------------------------------------------------------------------#
################################### Figure 10 ##################################
#------------------------------------------------------------------------------#
All_differences$Category <- as.factor(All_differences$Category)%>%              ## Reorders levels to correct order
  factor(., levels = rev(unique(All_differences$Category)))
ggplot(All_differences)+
  geom_segment(aes(x=Scenario, xend = Scenario, y = Q10, yend = Q90), lwd = 3, 
               alpha = 0.5)+
  geom_point(aes(Scenario, y = Q50), size = 5)+
  facet_grid(rows = vars(Region), cols = vars(Category), scales = "free_x")+
  ylab("Expected change \n(numbers of fish)")+coord_flip()

#------------------------------------------------------------------------------#
################################### Figure S1 ##################################
#------------------------------------------------------------------------------#
## Supplemental figure diagnostic plots
N <- 100                                                                        ## Number of draws used in diagnostics
colors <- c("Predicted" = "#4575B4", "Observed" = "black")                      ## Color scheme for plots
Harvest_diag <- posterior_predict(Harvest_Model, Train, ndraws = N)%>%          ## Extract posterior predictive scans
  t()%>%as.data.frame()%>%cbind(.,Train$Region)%>%
  pivot_longer(., cols = c(1:N), names_to = "scan", values_to = "value")        ## Pivot to two columns (one for scan name, one for scan value for each region-month)
colnames(Harvest_diag)[1]<- "Region"

Harvest_diag$limit <- max(All_Data$Harvest_Gag)                                 ## Filter to avoid insanely high draws and obscure density plots
Harvest_diag <- Harvest_diag%>%filter(value < limit)

Dens <- Train                                                                   ## Plots kernel density plots of observed and predicted response values
Harvest_diagnostics_dens <- ggplot(Dens)+
  geom_line(data = Harvest_diag, aes(x = (value), group = scan, col = "Predicted"),
            stat = 'density', alpha = 0.1, show.legend = F)+
  geom_line(aes(x = Harvest_Gag, col = 'Observed'),stat = 'density', lwd = 1, 
            show.legend = F)+
  scale_color_manual(values = colors)+
  xlab("Gag harvest")+ylab("Density")+
  facet_wrap(~Region, scales = "free", ncol = 1)+#Supplemental_theme()+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "top")+labs(color = "Value:")


Dens$Predicted <- predict(Harvest_Model, Train, probs = c(0.1, 0.9), robust = T)%>%## Scatter plot for diagnostics (only median estimate for each region-monthis used)
  as.data.frame()%>%.[,1]

Harvest_diagnostics_scatter <- ggplot(Dens)+                                    ## Plot scatter plot diagnostics
  geom_point(aes(x = Predicted, y = Harvest_Gag), col = "#4575B4", alpha = 0.5)+
  geom_abline(lwd=1)+
  xlab("Predicted gag harvest")+ylab("Observed gag harvest")+
  facet_wrap(~Region, scales = "free", ncol = 1)+#Supplemental_theme()+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Combine scatter plot and density plot into single figure
ggarrange(Harvest_diagnostics_scatter,Harvest_diagnostics_dens, ncol = 2, labels = c("a)", "b)"))%>%
  annotate_figure(.,top = text_grob("Harvest model diagnostics", size = 18, family = 'serif'))

#------------------------------------------------------------------------------#
################################### Figure S2 ##################################
#------------------------------------------------------------------------------#
Discard_diag <- posterior_predict(Discard_Model, Train, ndraws = N)%>%          ## Extract posterior predictive scans
  t()%>%as.data.frame()%>%cbind(.,Train$Region)%>%
  pivot_longer(., cols = c(1:N), names_to = "scan", values_to = "value")        ## Pivot to two columns (one for scan name, one for scan value for each region-month)
colnames(Harvest_diag)[1]<- "Region"

Discard_diag$limit <- max(All_Data$Discard_Gag)                                 ## Filter to avoid insanely high draws and obscure density plots
Discard_diag <- Discard_diag%>%filter(value < limit)

Dens <- Train                                                                   ## Plots kernel density plots of observed and predicted response values
Discard_diagnostics_dens <- ggplot(Dens)+
  geom_line(data = Harvest_diag, aes(x = (value), group = scan, col = "Predicted"),
            stat = 'density', alpha = 0.1, show.legend = F)+
  geom_line(aes(x = Discard_Gag, col = 'Observed'),stat = 'density', lwd = 1, 
            show.legend = F)+
  scale_color_manual(values = colors)+
  xlab("Gag harvest")+ylab("Density")+
  facet_wrap(~Region, scales = "free", ncol = 1)+#Supplemental_theme()+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(legend.position = "top")+labs(color = "Value:")


Dens$Predicted <- predict(Discard_Model, Train, probs = c(0.1, 0.9), robust = T)%>%## Scatter plot for diagnostics (only median estimate for each region-monthis used)
  as.data.frame()%>%.[,1]

Discard_diagnostics_scatter <- ggplot(Dens)+                                    ## Plot scatter plot diagnostics
  geom_point(aes(x = Predicted, y = Discard_Gag), col = "#4575B4", alpha = 0.5)+
  geom_abline(lwd=1)+
  xlab("Predicted gag discards")+ylab("Observed gag discards")+
  facet_wrap(~Region, scales = "free", ncol = 1)+#Supplemental_theme()+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Combine scatter plot and density plot into single figure
ggarrange(Discard_diagnostics_scatter,Discard_diagnostics_dens, ncol = 2, labels = c("a)", "b)"))%>%
  annotate_figure(.,top = text_grob("Discard model diagnostics", size = 18, family = 'serif'))

#------------------------------------------------------------------------------#
################################### Table S3 ###################################
#------------------------------------------------------------------------------#
rownames(All_Harvest_tables[[2]]) <- rep("", nrow(All_Harvest_tables[[1]]))
print(xtable(All_Harvest_tables[[2]]),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)

#------------------------------------------------------------------------------#
################################### Table S4 ###################################
#------------------------------------------------------------------------------#
rownames(All_Discard_tables[[2]]) <- rep("", nrow(All_Discard_tables[[1]]))
print(xtable(All_Discard_tables[[2]]),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)
