#------------------------------------------------------------------------------#
################################### Description ################################
#------------------------------------------------------------------------------#
# This file produces user-defined functions employed in generating tables and 
# figures used in: Modeling effort in a multispecies recreational fishery; 
# influence of species-specific temporal closures, relative abundance, and 
# seasonality on angler-trips
#
# The code below is annotated to explain to the user what each component 
# line executes.
#
# This file was written by A. Challen Hyman, PhD, on July 29th, 2024
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
##################################### Functions ################################
#------------------------------------------------------------------------------#

# Simple helper functions
`%nin%` <- Negate(`%in%`)                                                       ## 'Not in' function

geometric.mean <- function(x){                                                  ## Geometric mean function
  exp(mean(log(na.omit(x))))
}
#------------------------------------------------------------------------------#
## Custom ggplot themes for plotting
Supplemental_theme <- function(){theme_bw(base_family = 'serif')%+replace%
    theme(axis.text = element_text(size = 12, color = 1, family = 'serif'),
          axis.title = element_text(size = 12, color = 1, family = 'serif'),
          strip.text.x = element_text(size = 14, margin = margin(0.25,0,0.25,0, "cm")),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)
    )}

My_theme <- function(){
  theme_bw(base_family = 'serif')%+replace%
    theme(axis.text = element_text(size = 16, color = 1, family = 'serif'),
          axis.title = element_text(size = 18, color = 1, family = 'serif'),
          strip.text = element_text(size = 14, margin = margin(0.2,0,0.2,0, "cm")),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18)
    )}

#------------------------------------------------------------------------------#
## Function to calculate fraction of months open and season length (for simulations)
Season_Month_Fraction <- function(
    Year = 2023,
    Season = c("June 16", "July 31"),
    Additional = c("September 01", "Nov 01"),
    Additional_KOD = 'weekend'
){
  if(is.null(Season)){stop("Please specify a season start and end date")}
  if(is.null(Year)){stop("Please specify a year")}
  
  ## This creates a data frame of all days in each month of the specified year
  All_Dates_Year <- as.Date(paste(c("January 01", "December 31"), Year), format = "%B %d %Y")
  All_Dates_Year <- data.frame(Date = seq.Date(All_Dates_Year[1], All_Dates_Year[2], by = 'day'))
  All_Dates_Year$Month <- month(All_Dates_Year$Date)
  Month_days <- as.data.frame(table(All_Dates_Year$Month))
  
  if(is.null(Additional)){                                                      ## If no additional days (e.g., red grouper)
    ## Create vector of season dates
    Dates <- as.Date(paste(Season, Year), format = "%B %d %Y")
    Date_vector <- data.frame(Date = seq.Date(Dates[1], Dates[2], by = 'day'))
    Season_length <- nrow(Date_vector)
  } else {                                                                      ## If additional days (e.g., red snapper)
    Dates_1 <- as.Date(paste(Season, Year), format = "%B %d %Y")
    Date_vector_1 <- data.frame(Date = seq.Date(Dates_1[1], Dates_1[2], by = 'day'))
    Dates_2 <- as.Date(paste(Additional, Year), format = "%B %d %Y")
    Date_vector_2 <- data.frame(Date = seq.Date(Dates_2[1], Dates_2[2], by = 'day'))
    if(tolower(Additional_KOD) == "weekday"){
      Date_vector_2 <- as.data.frame(Date_vector_2[which(paste(wday(Date_vector_2$Date, abbr = T, label = T)) %nin% c("Fri", "Sat", "Sun")),])
      colnames(Date_vector_2) <- "Date"
      ## Creates second season length variable
      Season_length <- sum(c(nrow(Date_vector_1), nrow(Date_vector_2)))
    }
    if(tolower(Additional_KOD) == "weekend"){
      Date_vector_2 <- as.data.frame(Date_vector_2[which(paste(wday(Date_vector_2$Date, abbr = T, label = T)) %in% c("Fri", "Sat", "Sun")),])
      colnames(Date_vector_2) <- "Date"
      ## Creates second season length variable
      Season_length <- sum(c(nrow(Date_vector_1), nrow(Date_vector_2)))
    }
    Date_vector <- rbind(Date_vector_1, Date_vector_2)
  }
  All_Dates_Year$Open <- 0
  All_Dates_Year$Open[which(All_Dates_Year$Date %in% Date_vector$Date)] <- 1
  All_Dates_Year$day <- ifelse(day(All_Dates_Year$Date) < 15, "01", "01")
  All_Dates_Year <- All_Dates_Year%>%
    mutate(Date2 = as.Date(paste(year(Date), month(Date), day, sep = "-")))
  Fraction_Open <- All_Dates_Year%>%group_by(Date2)%>%summarize(Open = mean(Open))
  return(list(Fraction_Open, Season_length))                                    ## First item returns vector of fractions of each month open
}        
#------------------------------------------------------------------------------#
### Function to generate tables for harvest and discard parameters
Hyman_tables <- function(Model, Type = "Discard"){
  Hyman_preds <- as_draws_matrix(Model)%>%                                      ## Extract terms from each model
    as.data.frame()%>%.[,-ncol(.)]%>%.[,-ncol(.)]
  Mean_terms <- c(-grep("shape", colnames(Hyman_preds)),                        ## Mean component
                  -grep("hu", colnames(Hyman_preds)))
  Hurdle_terms <- grep("hu", colnames(Hyman_preds))                             ## Hurdle component
  Shape_terms <- grep("shape", colnames(Hyman_preds))                           ## Shape component
  
  ## Names of each term in correct order
  ### Mean
  Mean_names <- c("PH", "PN",
                  paste0("PH:Index"), paste0("PN:Index"),
                  paste0("PH:Juvenile"), paste0("PN:Juvenile"),
                  paste0("PH:Open_{Gag}"), paste0("PN:Open_{Gag}"),
                  paste0("PH:Open_{RS}"), paste0("PN:Open_{RS}"),
                  "PH:F", "PN:F",
                  "PH:Temp", "PN:Temp",
                  "PH:sin_{12}",
                  "PN:sin_{12}", 
                  "PH:cos_{12}",
                  "PN:cos_{12}",
                  "PH:sin_{6}", 
                  "PN:sin_{6}", 
                  "PH:cos_{6}",
                  "PN:cos_{6}")
  
  ## Hurdle structure changes depending on harvest or discard model
  if(Type == "Discard"){
    Hurdle_names <- c("PH", "PN",
                      paste0("PH:Index"), paste0("PN:Index"))
  } else {
    Hurdle_names <- c("PH", "PN",
                      paste0("PH:Season_{Gag}"), paste0("PN:Season_{Gag}"))
  }
  
  ### Shape
  Shape_names <-c("PH", "PN",
                  "PH:sin_{12}",
                  "PN:sin_{12}", 
                  "PH:cos_{12}",
                  "PN:cos_{12}",
                  "PH:sin_{6}", 
                  "PN:sin_{6}", 
                  "PH:cos_{6}",
                  "PN:cos_{6}",
                  paste0("PH:Index"), paste0("PN:Index"))
  
  #-------------------- Supplemental table --------------------#
  Hyman_supplemental_tables <-apply(Hyman_preds, 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})
  
  ## Mean
  Mean_supp <- Hyman_supplemental_tables[,Mean_terms]
  Mean_supp <-Mean_supp[,c(
    1,2,                                                                       
    grep("Index", colnames(Mean_supp)),
    grep("Juvenile", colnames(Mean_supp)),
    grep("M_Gag", colnames(Mean_supp)),
    grep("M_RS", colnames(Mean_supp)),
    grep("Trips", colnames(Mean_supp)),
    grep("SST", colnames(Mean_supp)),
    grep("sin1", colnames(Mean_supp)),
    grep("cos1", colnames(Mean_supp)),
    grep("sin2", colnames(Mean_supp)),
    grep("cos2", colnames(Mean_supp))
  )]%>%t()%>%round(.,3)
  
  ##Hurdle
  Hurdle_supp <- Hyman_supplemental_tables[,Hurdle_terms]
  Hurdle_supp <-Hurdle_supp[,c(
    1,2,                                                                       
    grep("Index", colnames(Hurdle_supp)),
    grep("M", colnames(Hurdle_supp))                                                              
  )]%>%t()%>%round(.,3)
  
  ##Shape
  Shape_supp <- Hyman_supplemental_tables[,Shape_terms]
  Shape_supp <-Shape_supp[,c(
    1,2,
    grep("sin1", colnames(Shape_supp)),
    grep("cos1", colnames(Shape_supp)),
    grep("sin2", colnames(Shape_supp)),
    grep("cos2", colnames(Shape_supp)),
    grep("Index", colnames(Shape_supp))                                                           
  )]%>%t()%>%round(.,3)
  
  ## parameter values based on harvest or discard model
  if(Type == "Discard"){
    Mean_params <- paste0("$\\psi_{",0:(nrow(Mean_supp)-1),"}$")
    Hurdle_params <- paste0("$\\zeta_{",0:(nrow(Hurdle_supp)-1),"}$")
    Shape_params <- paste0("$\\xi_{",0:(nrow(Shape_supp)-1),"}$")
    
    Mean_supp <- cbind(c("$\\eta$", rep("", (nrow(Mean_supp)-1))),Mean_names,Mean_params,Mean_supp)
    Hurdle_supp <- cbind(c("$\\theta$", rep("", (nrow(Hurdle_supp)-1))),Hurdle_names,Hurdle_params,Hurdle_supp)
    Shape_supp <- cbind(c("$\\tau$", rep("", (nrow(Shape_supp)-1))),Shape_names,Shape_params,Shape_supp)
  } else {
    Mean_params <- paste0("$\\omega_{",0:(nrow(Mean_supp)-1),"}$")
    Hurdle_params <- paste0("$\\nu_{",0:(nrow(Hurdle_supp)-1),"}$")
    Shape_params <- paste0("$\\gamma_{",0:(nrow(Shape_supp)-1),"}$")
    
    Mean_supp <- cbind(c("$\\lambda$", rep("", (nrow(Mean_supp)-1))),Mean_names,Mean_params,Mean_supp)
    Hurdle_supp <- cbind(c("$\\delta$", rep("", (nrow(Hurdle_supp)-1))),Hurdle_names,Hurdle_params,Hurdle_supp)
    Shape_supp <- cbind(c("$\\phi$", rep("", (nrow(Shape_supp)-1))),Shape_names,Shape_params,Shape_supp)
  }
  ## Column names consistent for each component
  colnames(Mean_supp)[1:3] <- 
    colnames(Hurdle_supp)[1:3] <- 
    colnames(Shape_supp)[1:3] <- c("Component", "Predictor", "Regression Coefficient")
  
  ## Bind to single data frame
  Supplemental_table <- rbind(Mean_supp, Hurdle_supp, Shape_supp)
  rownames(Supplemental_table) <- NULL
  
  ## Add 'significance' notation
  Signif <- ifelse(sign(as.numeric(Supplemental_table[,4])) == sign(as.numeric(Supplemental_table[,6])), "*", "")
  Supplemental_table[,2] <- paste0(Supplemental_table[,2], Signif)
  
  #------------------------ Main table ------------------------#
  ## Mean
  Mean_main <- Hyman_preds[,Mean_terms]
  Mean_main <- Mean_main[,c(
    1,2,                                                                       
    grep("Index", colnames(Mean_main)),
    grep("Juvenile", colnames(Mean_main)),
    grep("M_Gag", colnames(Mean_main)),
    grep("M_RS", colnames(Mean_main)),
    grep("Trips", colnames(Mean_main)),
    grep("SST", colnames(Mean_main)),
    grep("sin1", colnames(Mean_main)),
    grep("cos1", colnames(Mean_main)),
    grep("sin2", colnames(Mean_main)),
    grep("cos2", colnames(Mean_main))
  )]
  
  ## Marginal effects for Peninsula
  for (i in seq(1, (ncol(Mean_main)-1), 2)){
    Mean_main[,i+1] <- Mean_main[,i+1] + Mean_main[,i]
  }
  
  ## Separate by region
  PH <- c(-grep("Peninsula", colnames(Mean_main)))
  PN <- grep("Peninsula", colnames(Mean_main))
  
  ## Generate summary statistics
  Mean_main_PH <- Mean_main[,PH]%>%apply(., 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})%>%t()%>%round(.,3)
  Mean_main_PN <- Mean_main[,PN]%>%apply(., 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})%>%t()%>%round(.,3)
  
  ## Append significance
  Signif_PH <- ifelse(sign(as.numeric(Mean_main_PH[,1])) == sign(as.numeric(Mean_main_PH[,3])), "*", "")
  Signif_PN <- ifelse(sign(as.numeric(Mean_main_PN[,1])) == sign(as.numeric(Mean_main_PN[,3])), "*", "")
  
  ## Generate parameter names based on model
  if(Type == "Discard"){
    PH_names <- paste0("$\\psi_{", gsub("\\$", "", Mean_names[PH]), "}$", Signif_PH)
    PN_names <- paste0("$\\psi_{", gsub("\\$", "", Mean_names[PN]), "}$", Signif_PN)
  } else {
    PH_names <- paste0("$\\omega_{", gsub("\\$", "", Mean_names[PH]), "}$", Signif_PH)
    PN_names <- paste0("$\\omega_{", gsub("\\$", "", Mean_names[PN]), "}$", Signif_PN)
  }
  
  ## Append columns to generate mean table
  Mean_main_table <- cbind(PH_names, Mean_main_PH,PN_names, Mean_main_PN)
  
  ## Hurdle
  Hurdle_main <- Hyman_preds[,Hurdle_terms]
  Hurdle_main <-Hurdle_main[,c(
    1,2,                                                                       
    grep("Index", colnames(Hurdle_main)),
    grep("M", colnames(Hurdle_main))                                                               
  )]
  
  ## Marginal effects for Peninsula
  for (i in seq(1, (ncol(Hurdle_main)-1), 2)){
    Hurdle_main[,i+1] <- Hurdle_main[,i+1] + Hurdle_main[,i]
  }
  
  ## Separate by region
  PH <- c(-grep("Peninsula", colnames(Hurdle_main)))
  PN <- grep("Peninsula", colnames(Hurdle_main))
  
  ## Generate summary statistics
  Hurdle_main_PH <- Hurdle_main[,PH]%>%apply(., 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})%>%t()%>%round(.,3)
  Hurdle_main_PN <- Hurdle_main[,PN]%>%apply(., 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})%>%t()%>%round(.,3)
  
  ## Append significance
  Signif_PH <- ifelse(sign(as.numeric(Hurdle_main_PH[,1])) == sign(as.numeric(Hurdle_main_PH[,3])), "*", "")
  Signif_PN <- ifelse(sign(as.numeric(Hurdle_main_PN[,1])) == sign(as.numeric(Hurdle_main_PN[,3])), "*", "")
  
  ## Generate parameter names based on model
  if(Type == "Discard"){
    PH_names <- paste0("$\\theta_{", gsub("\\$", "", Hurdle_names[PH]), "}$", Signif_PH)
    PN_names <- paste0("$\\theta_{", gsub("\\$", "", Hurdle_names[PN]), "}$", Signif_PN)
  } else {
    PH_names <- paste0("$\\delta_{", gsub("\\$", "", Hurdle_names[PH]), "}$", Signif_PH)
    PN_names <- paste0("$\\delta_{", gsub("\\$", "", Hurdle_names[PN]), "}$", Signif_PN)
  }
  
  ## Append columns to generate mean table
  Hurdle_main_table <- cbind(PH_names, Hurdle_main_PH,PN_names, Hurdle_main_PN)
  
  ## Shape
  Shape_main <- Hyman_preds[,Shape_terms]
  Shape_main <-Shape_main[,c(
    1,2,
    grep("sin1", colnames(Shape_main)),
    grep("cos1", colnames(Shape_main)),
    grep("sin2", colnames(Shape_main)),
    grep("cos2", colnames(Shape_main)),
    grep("Index", colnames(Shape_main))                                                               
  )]
  
  ## Marginal effects for Peninsula
  for (i in seq(1, (ncol(Shape_main)-1), 2)){
    Shape_main[,i+1] <- Shape_main[,i+1] + Shape_main[,i]
  }
  
  ## Separate by region
  PH <- c(-grep("Peninsula", colnames(Shape_main)))
  PN <- grep("Peninsula", colnames(Shape_main))
  
  ## Generate summary statistics
  Shape_main_PH <- Shape_main[,PH]%>%apply(., 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})%>%t()%>%round(.,3)
  Shape_main_PN <- Shape_main[,PN]%>%apply(., 2, function(x){quantile(x,c(0.1, 0.5, 0.9))})%>%t()%>%round(.,3)
  
  ## Append significance
  Signif_PH <- ifelse(sign(as.numeric(Shape_main_PH[,1])) == sign(as.numeric(Shape_main_PH[,3])), "*", "")
  Signif_PN <- ifelse(sign(as.numeric(Shape_main_PN[,1])) == sign(as.numeric(Shape_main_PN[,3])), "*", "")
  
  ## Generate parameter names based on model
  if(Type == "Discard"){
    PH_names <- paste0("$\\tau_{", gsub("\\$", "", Shape_names[PH]), "}$", Signif_PH)
    PN_names <- paste0("$\\tau_{", gsub("\\$", "", Shape_names[PN]), "}$", Signif_PN)
  } else {
    PH_names <- paste0("$\\phi_{", gsub("\\$", "", Shape_names[PH]), "}$", Signif_PH)
    PN_names <- paste0("$\\phi_{", gsub("\\$", "", Shape_names[PN]), "}$", Signif_PN)
  }
  
  ## Append columns to generate mean table
  Shape_main_table <- cbind(PH_names, Shape_main_PH,PN_names, Shape_main_PN)
  
  ## Append mean, hurdle, and shape tables to create main table
  Main_Table <- rbind(Mean_main_table,
                      Hurdle_main_table,
                      Shape_main_table)
  colnames(Main_Table)[c(1, 4)] <- c("Panhandle", "Peninsula")
  return(list(Main_Table, Supplemental_table))
}
