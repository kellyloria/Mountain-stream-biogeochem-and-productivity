#### Possible nutrient data infill from EGRET package

# https://github.com/DOI-USGS/EGRET/blob/main/R/readUserDaily.r

# https://cran.r-project.org/web/packages/EGRET/vignettes/EGRET.html


library(EGRET)
library(dataRetrieval)

############################
# Check flow history data:

############################
# Gather discharge data:
siteNumber <- "01491000" #Choptank River near Greensboro, MD
startDate <- "2007-09-30" #Gets earliest date
endDate <- "2011-09-30"
# Gather sample data:
parameter_cd<-"00631" #5 digit USGS code
Sample <- readNWISSample(siteNumber,parameter_cd,
                         startDate=as.Date("2007-09-30"),endDate=as.Date("2011-09-30"))
#Gets earliest date from Sample record:
#This is just one of many ways to assure the Daily record
#spans the Sample record
startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:

# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "Choptank River at Greensboro, MD"

# Merge discharge with sample data:
eList <- mergeReport(INFO, Daily, Sample)
############################

savePath<-"/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/egret_chem_infil/" 
saveResults(savePath, eList)


############################
# Check sample data:
boxConcMonth(eList)
boxQTwice(eList)
plotConcTime(eList)
plotConcQ(eList)
multiPlotDataOverview(eList)

####
bg_nuts <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/NS_chem_dat_nh4_24.rds") %>%
  filter(location=="stream") %>%
  mutate(NO3_mgL_dl = if_else(NO3_mgL_dl < 0.0015, 0.0015, NO3_mgL_dl))%>%
  mutate(NH4_mgL_dl=if_else(NH4_mgL_dl < 0.001, 0.001, NH4_mgL_dl))%>%
  mutate(PO4_ugL_dl=if_else(PO4_ugL_dl < 0.201, 0.201, PO4_ugL_dl)) %>%
  mutate(PO4_ugL_dl=if_else(DOC_mgL_dl < 0.25, 0.125, DOC_mgL_dl)) %>%
  group_by(site, shore, date, depth, location, substrate) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


##### pivot wider 
bg_nuts_wide <- bg_nuts %>%
  dplyr::select(site, date, substrate, NO3_mgL_dl, NH4_mgL_dl, PO4_ugL_dl, DOC_mgL_dl, pH_infill) %>%
  group_by(site, date,substrate) %>% 
  summarise(across(c(NO3_mgL_dl, NH4_mgL_dl, PO4_ugL_dl, DOC_mgL_dl, pH_infill), mean, na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = substrate, values_from = c(NO3_mgL_dl, NH4_mgL_dl, PO4_ugL_dl, DOC_mgL_dl, pH_infill), names_sep = "_")

######
BWL_NH4 <- bg_nuts_wide %>%
  filter(site=="BWL") %>%
  dplyr::select(date, NH4_mgL_dl_sw)


###########################
# Run WRTDS model:
nrow(eList$Sample)

eList$INFO$minNumObs <- nrow(eList$Sample)

eList$INFO$minNumObs <-nrow(eList$Sample)



modelEstimation(eList, windowY = 7, windowQ = 2, windowS = 0.5,
                minNumObs = 50, minNumUncen = 25, edgeAdjust = TRUE, verbose = TRUE,
                run.parallel = FALSE)

############################
# 1. Model Setup and Execution

# windowY = 7: The model considers a 7-year moving window for trends in time.
# windowQ = 2: A 2-year window for flow (Q), meaning that discharge trends are smoothed

# windowS = 0.5: A seasonal smoothing window of 0.5 years (6 months).
# minNumObs = 50: At least 50 observations are required for estimation.
# minNumUncen = 25: At least 25 uncensored (non-below-detection) observations are required.
# edgeAdjust = TRUE: Adjusts for potential biases at the start and end of the time series.

# 2. Cross-validation Process
# The model first runs estCrossVal, which validates the regression model's predictions by withholding subsets of the data and checking how well the model can predict them.

# 3. Survival Regression for Trend Estimation:
# The model then runs estSurfaces using a survival regression, a statistical approach to handle censored data (e.g., values below detection limits).
############################
#Check model results:

# After running modelEstimation(), the estimated concentration and flux values are stored in the eList object:
# Extract Estimated Concentrations:

eList_mod <- modelEstimation(eList, minNumObs = 20, minNumUncen = 10)

eList_mod$Daily
summary(eList$Daily)

# ConcDay: Estimated daily concentration (mg/L).
# FluxDay: Estimated daily flux (kg/day).

# Extract Estimated Trends
tableChange(eList)
summary(eList$Sample$DecYear)
?tableChange()

tableChange(eList, fluxUnit = 2, yearPoints = c(2009, 2010, 2012))


# 2. Plot Model Results
plotConcTime(eList)



#Require Sample + INFO:
plotConcTimeDaily(eList) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList) # Plot Modeled Flux over Time

plotConcHist(eList) # Flow-Normalized Concentration Trends

plotConcPred(eList)
plotFluxPred(eList)
plotResidPred(eList) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.

plotResidTime(eList) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList)
boxConcThree(eList) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList)
plotFluxHist(eList)


nrow(eList$Sample)

# Check how many observations exist for different flow conditions
table(eList$Sample$Q)

# Check how many observations exist per year
table(format(eList$Sample$Date, "%Y"))

qBottom <- 0.1
qTop <- 20

plotConcQSmooth(eList, date1, date2, date3, qBottom, qTop, concMax=2, legendTop=0.8, windowY = 3)

plotConcQSmooth(eList, date1, date2, date3, qBottom, qTop, concMax=2, legendTop=0.8, windowY = 5)


eList$Daily <- subset(eList$Daily, Date >= as.Date("2007-09-30") & Date <= as.Date("2011-09-30"))
eList$Sample <- subset(eList$Sample, Date >= as.Date("2007-09-30") & Date <= as.Date("2011-09-30"))

date1 <- as.Date("2008-09-01")
date2 <- as.Date("2009-09-01")
date3 <- as.Date("2010-09-01")
qBottom <- 0.2
qTop <- 10
plotConcQSmooth(eList, date1, date2, date3, qBottom, qTop, concMax=2, legendTop=0.8)


# Multi-line plots:
date1 <- as.Date("2008-09-01")
date2 <- as.Date("2009-09-0")
date3 <- as.Date("2010-09-01")
qBottom<-0.2
qTop<-10
plotConcQSmooth(eList, date1, date2, date3, qBottom, qTop, 
                concMax=2,legendTop = 0.8)
q1 <- 2
q2 <- 10
q3 <- 20
centerDate <- "07-01"
yearEnd <- 2007
yearStart <- 2012
plotConcTimeSmooth(eList, q1, q2, q3, centerDate, yearStart, yearEnd, legendTop = 0.7)

# Multi-plots:
fluxBiasMulti(eList)


plotContours(eList, yearStart,yearEnd,qBottom,qTop, 
             contourLevels = clevel,qUnit=2)

plotDiffContours(eList, yearStart,yearEnd,
                 qBottom,qTop,maxDiff,qUnit=2)

# modify this for your own computer file structure:
savePath<-"/Users/rhirsch/Desktop/" 
saveResults(savePath, eList)



### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod[["Daily"]]$Date <- as.Date(eList_mod[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_agg <- eList_mod[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_agg, by = "Date", all.x = TRUE)

# Create the plot
logQ_plt3 <- ggplot(combined_data, aes(x = Date)) +
  # Sample data as black points
  geom_point(aes(y = ConcAve), color = "black", shape = 21, alpha = 0.9) + 
  
  # Modeled data as red triangles
  geom_point(aes(y = ConcDay_mod), color = "#D62828", shape = 17, alpha = 0.9) + 
  
  # Adding titles and labels
  labs(
    title = "Comparison of Modeled vs. Sampled Nitrogen Concentrations",
    x = "Date",
    y = "Nitrogen Concentration (mg/L)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate date labels for better visibility

# Print the plot
print(logQ_plt3)



########

siteNo_GB <- "10336730"
siteNo_BW <- "10336660"

flow_BWL <- flow_df%>%
  filter(site=="BWL") %>%
  dplyr::select(date, Q_m)

INFO_BWL <- readNWISInfo(siteNo_BW,"00060")
INFO_BWL$shortName <- "BW Lower"


############################
# MY sites
############################

####### TDN 
############################
# Gather discharge data:
siteNumber <- "10336660" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-09-30"
# Gather sample data:
parameter_cd<-"00631" #5 digit USGS code
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
#Gets earliest date from Sample record:
#This is just one of many ways to assure the Daily record
#spans the Sample record
startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:

# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "BWL"

# Merge discharge with sample data:
eList_BWL <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_BWL)
boxQTwice(eList_BWL)
plotConcTime(eList_BWL)
plotConcQ(eList_BWL)
multiplot <- multiPlotDataOverview(eList_BWL)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWL$Sample)

eList_BWL$INFO$minNumObs <- nrow(eList_BWL$Sample)


eList_mod_BWL <- modelEstimation(eList_BWL, windowY = 3, windowQ = 1, windowS = 0.5,
                minNumObs = 70, minNumUncen = 25, edgeAdjust = TRUE, verbose = TRUE,
                run.parallel = FALSE)


eList_mod_BWL$Daily
summary(eList_mod_BWL$Daily)

plotConcTime(eList_mod_BWL)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWL) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_BWL) # Plot Modeled Flux over Time

plotConcHist(eList_mod_BWL) # Flow-Normalized Concentration Trends

plotConcPred(eList_mod_BWL)
plotFluxPred(eList_mod_BWL)
plotResidPred(eList_mod_BWL) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWL) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.

plotResidTime(eList_mod_BWL) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWL)
boxConcThree(eList_mod_BWL) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWL)
plotFluxHist(eList_mod_BWL)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWL[["Daily"]]$Date <- as.Date(eList_mod_BWL[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWL_agg <- eList_mod_BWL[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_BWL_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt3 <- ggplot(combined_data, aes(x = Date)) +
  # Sample data as black points
  geom_point(aes(y = ConcAve, color = "Sampled TIN"), shape = 21, alpha = 0.9) + 
  geom_line(aes(y = ConcAve, color = "Sampled TIN")) +
  
  # Modeled data as red triangles
  geom_point(aes(y = ConcDay_mod, color = "Modeled TIN"), shape = 2, alpha = 0.9) + 
  geom_line(aes(y = ConcDay_mod, color = "Modeled TIN"), alpha = 0.9) + 
  
  # Adding titles and labels
  labs(
    title = "Comparison of Modeled vs. Sampled Nitrogen Concentrations",
    x = "Date",
    y = "TIN Nitrogen Concentration (mg/L)"
  ) +
  scale_color_manual(values = c("Sampled TIN" = "black", "Modeled TIN" = "#D62828")) +  # Custom colors for the legend
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +  # Set breaks every 4 months
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate date labels for better visibility

# Print the plot
print(logQ_plt3)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample_TS.png", plot = logQ_plt3, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_BWL_TIN.csv")

###

############################
# MY sites + my chem 
############################
##### NO3?  
############################
# Gather discharge data:
siteNumber <- "10336660" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-09-30"
# Gather sample data:
parameter_cd<-"00631" #5 digit USGS code
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
Sample$source <- "nwis"


# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
BWL_NO3$Date <- as.Date(BWL_NO3$date)
BWL_NO3$ConcAve <- (BWL_NO3$NO3_mgL_dl_sw)


# Add the "source" column to each dataset
Sample <- Sample %>%
  mutate(source = "USGS")

BWL_NO3 <- BWL_NO3 %>%
  mutate(source = "KAL")

# Perform the full join by matching the date columns
Sample <- full_join(Sample, BWL_NO3, by = c("Date", "source", "ConcAve"))


# Create the plot
logQ_plt3 <- ggplot(Sample, aes(x = Date, colour = source, shape=source)) +
  # Sample data as black points
  geom_point(aes(y = ConcAve), alpha = 0.9) + 
  # Adding titles and labels
  labs(
    title = "Comparison of Modeled vs. Sampled Nitrogen Concentrations",
    x = "Date",
    y = "TIN Nitrogen Concentration (mg/L)"
  ) +
 # scale_color_manual(values = c("Sampled TIN" = "black", "Modeled TIN" = "#D62828")) +  # Custom colors for the legend
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +  # Set breaks every 4 months
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate date labels for better visibility

# Print the plot
print(logQ_plt3)



startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:

# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "BWL"

# Merge discharge with sample data:
eList_BWL <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_BWL)
boxQTwice(eList_BWL)
plotConcTime(eList_BWL)
plotConcQ(eList_BWL)
multiplot <- multiPlotDataOverview(eList_BWL)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWL$Sample)

eList_BWL$INFO$minNumObs <- nrow(eList_BWL$Sample)

eList_mod_BWL$Daily
summary(eList_mod_BWL$Daily)

plotConcTime(eList_mod_BWL)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWL) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_BWL) # Plot Modeled Flux over Time


plotConcPred(eList_mod_BWL)
plotFluxPred(eList_mod_BWL)
plotResidPred(eList_mod_BWL) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWL) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_BWL) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_BWL) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWL)
boxConcThree(eList_mod_BWL) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWL)
plotFluxHist(eList_mod_BWL)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWL[["Daily"]]$Date <- as.Date(eList_mod_BWL[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWL_agg <- eList_mod_BWL[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_BWL_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt3 <- ggplot(combined_data, aes(x = Date)) +
  # Sample data as black points
  geom_point(aes(y = ConcAve, color = "Sampled TIN",  shape = source), alpha = 0.9) + 
  geom_line(aes(y = ConcAve, color = "Sampled TIN")) +
  
  # Modeled data as red triangles
  geom_point(aes(y = ConcDay_mod, color = "Modeled TIN"), shape = 2, alpha = 0.9) + 
  geom_line(aes(y = ConcDay_mod, color = "Modeled TIN"), alpha = 0.9) + 
  
  # Adding titles and labels
  labs(
    title = "Comparison of Modeled vs. Sampled Nitrogen Concentrations",
    x = "Date",
    y = "N Nitrogen Concentration (mg/L)"
  ) +
  scale_color_manual(values = c("Sampled TIN" = "black", "Modeled TIN" = "#D62828")) +  # Custom colors for the legend
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +  # Set breaks every 4 months
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate date labels for better visibility

# Print the plot
print(logQ_plt3)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_NO3_plot_sample_TS.png", plot = logQ_plt3, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_BWL_NO3.csv")

###
###


############################
##### NH4?  
############################
# Gather discharge data:
siteNumber <- "10336660"
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-09-30"
# Gather sample data:
parameter_cd<-"00608" #5 digit USGS code
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)

# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
BWL_NH4$Date <- as.Date(BWL_NH4$date)
BWL_NH4$ConcAve <- (BWL_NH4$NH4_mgL_dl_sw)

BWL_NH4$ConcLow <- (BWL_NH4$NH4_mgL_dl_sw)
BWL_NH4$ConcHigh <- (BWL_NH4$NH4_mgL_dl_sw)

str(Sample)

# Add the "source" column to each dataset
Sample <- Sample %>%
  mutate(source = "USGS")

BWL_NH4 <- BWL_NH4 %>%
  mutate(source = "KAL")

# Perform the full join by matching the date columns
Sample <- full_join(Sample, BWL_NH4, by = c("Date", "source", 
                                            "ConcAve", "ConcLow", "ConcHigh"))

Sample$Julian <- as.numeric(Sample$Date - as.Date("1960-01-01"))
Sample$Day <- lubridate::yday(Sample$Date)
Sample$DecYear <- lubridate::year(Sample$Date) + (lubridate::yday(Sample$Date) - 1) / 365.25
Sample$SinDY <- sin(2 * pi * Sample$Julian / 365.25)
Sample$CosDY <- cos(2 * pi * Sample$Julian / 365.25)
Sample$Month <- month(Sample$Date)
Sample$dateTime <- month(Sample$Date)
Sample$CharacteristicName <- "Ammonia and ammonium"
Sample$USGSPCode <- 00608
Sample$ActivityMediaName <- "Water"
Sample$ActivityStartDateTime <- as.POSIXct(paste(Sample$Date, "13:00:00"), tz = "America/Los_Angeles")
Sample$ResultSampleFractionText <- "Filtered field and/or lab"
Sample$ResultStatusIdentifier <- "Accepted"
Sample$ResultValueTypeName <- "Actual"
Sample$date <- as.Date(Sample$Date)
Sample <- water_year(Sample)
Sample$waterYear <- c(Sample$water_year)
Sample$MonthSeq <- (year(Sample$Date) - 1800) * 12 + month(Sample$Date)
Sample$Uncen <- 1


# Create the plot
logQ_plt3 <- ggplot(Sample, aes(x = Date, colour = source, shape=source)) +
  # Sample data as black points
  geom_point(aes(y = ConcAve), alpha = 0.9) + 
  # Adding titles and labels
  labs(
    title = "Comparison of Modeled vs. Sampled Nitrogen Concentrations",
    x = "Date",
    y = "TIN Nitrogen Concentration (mg/L)"
  ) +
  # scale_color_manual(values = c("Sampled TIN" = "black", "Modeled TIN" = "#D62828")) +  # Custom colors for the legend
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +  # Set breaks every 4 months
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate date labels for better visibility

# Print the plot
print(logQ_plt3)

Sample <- Sample %>%
  dplyr::select("Date", "ConcLow", "ConcHigh", "Uncen", "ConcAve",
                "Julian", "Month", "Day", "DecYear", "waterYear",
                "MonthSeq", "SinDY", "CosDY", "dateTime", "CharacteristicName",
                "USGSPCode", "ActivityStartDateTime", "ActivityMediaSubdivisionName",
                "ActivityMediaName", "ResultSampleFractionText", "ResultStatusIdentifier",
                "ResultValueTypeName")

startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:

Daily$Julian <- as.numeric(Daily$Date - as.Date("1960-01-01"))
Daily$Day <- lubridate::yday(Daily$Date)
Daily$DecYear <- lubridate::year(Daily$Date) + (lubridate::yday(Daily$Date) - 1) / 365.25
# Sample$SinDY <- sin(2 * pi * Sample$Julian / 365.25)
# Sample$CosDY <- cos(2 * pi * Sample$Julian / 365.25)
Daily$Month <- month(Daily$Date)
# Sample$dateTime <- month(Sample$Date)
# Sample$CharacteristicName <- "Ammonia and ammonium"
# Sample$USGSPCode <- 00608
# Sample$ActivityMediaName <- "Water"
# Sample$ActivityStartDateTime <- as.POSIXct(paste(Sample$Date, "13:00:00"), tz = "America/Los_Angeles")
# Sample$ResultSampleFractionText <- "Filtered field and/or lab"
# Sample$ResultStatusIdentifier <- "Accepted"
# Sample$ResultValueTypeName <- "Actual"
# Sample$date <- as.Date(Sample$Date)
# Sample <- water_year(Sample)
# Sample$waterYear <- c(Sample$water_year)
Daily$MonthSeq <- (year(Daily$Date) - 1800) * 12 + month(Daily$Date)


# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd, interactive=FALSE)
INFO$shortName <- "BWL"
# 
# INFO <- data.frame(
#   shortName = "BWL",
#   siteName = "BWL Site Name",
#   siteType = "Non-USGS",
#   drainageArea = 29.00787,
#   parameterUnits = "mg/l as N"
# )

Daily$Date <- as.Date(Daily$Date)
Sample$Date <- as.Date(Sample$Date)


Daily$ConcAve[is.nan(Daily$ConcAve)] <- NA

if (any(Sample$ConcLow == 0, na.rm = TRUE)) {
  stop("modelEstimation cannot be run with 0 values in ConcLow.")
}

INFO <- data.frame(
  siteName = "BWL",
  CharacteristicName = "Ammonia and ammonium",
  USGSPCode = "00608",
  stringsAsFactors = FALSE
)

eList_BWL <- mergeReport(INFO = INFO, Daily = Daily, Sample = Sample)


# Merge discharge with sample data

# Check the structure of the eList to confirm it's correctly formatted
str(eList_BWL)

# Merge Sample and Daily by Date (keeping all Sample records)
edat_BWL <- merge(Sample, Daily, by = "Date", all.x = TRUE)

# Create a list that includes INFO, the merged data, and optionally the original Daily and Sample
eList_BWL <- list(INFO = INFO, Daily = Daily, Sample = Sample, eList = eList_BWL)


############################


#### YOU broke things that you are fixing here!!!
eList

eList_mod_BWL <- EGRET::setUpEstimation(eList=eList_BWL, windowY = 3, windowQ = 1, windowS = 0.5,
                                        minNumObs = 70, minNumUncen = 25, verbose = TRUE)


eList1 <- setUpEstimation(eList = eList_BWL, windowY = 3, 
                         windowQ = 1, windowS = 0.5, minNumObs = 30, 
                         minNumUncen = 10, edgeAdjust = TRUE, verbose = TRUE)
if (verbose) 
  cat("\n first step running estCrossVal may take about 1 minute")
Sample1 <- estCrossVal(eList_BWL$Daily$DecYear[1], eList_BWL$Daily$DecYear[length(eList_BWL$Daily$DecYear)], 
                       eList_BWL$Sample, windowY = windowY, windowQ = windowQ, windowS = windowS, 
                       minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust = edgeAdjust, 
                       verbose = verbose)
eList$Sample <- Sample1
if (verbose) 
  cat("\nNext step running  estSurfaces with survival regression:\n")
surfaces1 <- estSurfaces(eList, windowY = windowY, windowQ = windowQ, 
                         windowS = windowS, minNumObs = minNumObs, minNumUncen = minNumUncen, 
                         edgeAdjust = edgeAdjust, verbose = verbose, run.parallel = run.parallel)
eList$surfaces <- surfaces1
Daily1 <- estDailyFromSurfaces(eList)
eList$Daily <- Daily1
checkSurfaceSpan(eList)
return(eList)
}


############################
# Check sample data:
boxConcMonth(edat_BWL)
boxQTwice(eList_BWL)
plotConcTime(eList_BWL)
plotConcQ(eList_BWL)
multiPlotDataOverview(eList_BWL)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWL$Sample)

eList_BWL$INFO$minNumObs <- nrow(eList_BWL$Sample)


eList_mod_BWL <- modelEstimation(eList_BWL, windowY = 3, windowQ = 1, windowS = 0.5,
                                 minNumObs = 30, minNumUncen = 15, edgeAdjust = TRUE, verbose = TRUE,
                                 run.parallel = FALSE)


eList_mod_BWL$Daily
summary(eList_mod_BWL$Daily)

plotConcTime(eList_mod_BWL)


# Check for missing values in Daily and Sample data frames
summary(eList_BWL$Daily)
summary(eList_BWL$Sample)

# Remove rows with NA values (if necessary)
eList_BWL$Daily <- na.omit(eList_BWL$Daily)
eList_BWL$Sample <- na.omit(eList_BWL$Sample)



######

# Step 1: Create the eList_BWL list object
eList_BWL <- list(
  Daily = data.frame(
    Date = edat_BWL$Date,
    Q = edat_BWL$Q,
    LogQ = log(edat_BWL$Q),
    Qualifier = edat_BWL$Qualifier
  ),
  Sample = data.frame(
    Date = edat_BWL$Date,
    ConcLow = edat_BWL$ConcLow,
    ConcHigh = edat_BWL$ConcHigh,
    ConcAve = edat_BWL$ConcAve
  ),
  INFO = data.frame(
    siteName = "BWL", # Replace with your site name
    CharacteristicName = unique(INFO$paramShortName)[1], # Assuming it’s constant
    USGSPCode = unique(INFO$paramNumber)[1], # Assuming it’s constant
    stringsAsFactors = FALSE
  )
)

# Check the column names of the Daily and Sample data frames
colnames(eList_BWL$Daily)
colnames(eList_BWL$Sample)
colnames(eList_BWL$INFO)

INFO$stringsAsFactors <- "FALSE"

# Check if eList_BWL exists and is a valid list
str(eList_BWL)

# Run modelEstimation with adjusted arguments
eList_mod_BWL <- modelEstimation(
  eList_BWL,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 35,  # Minimum number of observations
  minNumUncen = 25,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

# Check the result to see if it worked
str(eList_mod_BWL)


#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWL) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_BWL) # Plot Modeled Flux over Time


plotConcPred(eList_mod_BWL)
plotFluxPred(eList_mod_BWL)
plotResidPred(eList_mod_BWL) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWL) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_BWL) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_BWL) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWL)
boxConcThree(eList_mod_BWL) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWL)
plotFluxHist(eList_mod_BWL)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWL[["Daily"]]$Date <- as.Date(eList_mod_BWL[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWL_agg <- eList_mod_BWL[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_BWL_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt3 <- ggplot(combined_data, aes(x = Date)) +
  # Sample data as black points
  geom_point(aes(y = ConcAve, color = "Sampled N"), alpha = 0.9) + 
  geom_line(aes(y = ConcAve, color = "Sampled N")) +
  
  # Modeled data as red triangles
  geom_point(aes(y = ConcDay_mod, color = "Modeled N"), shape = 2, alpha = 0.9) + 
  geom_line(aes(y = ConcDay_mod, color = "Modeled N"), alpha = 0.9) + 
  
  # Adding titles and labels
  labs(
    title = "Comparison of Modeled vs. Sampled Nitrogen Concentrations",
    x = "Date",
    y = "N Nitrogen Concentration (mg/L)"
  ) +
  scale_color_manual(values = c("Sampled N" = "black", "Modeled N" = "#D62828")) +  # Custom colors for the legend
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +  # Set breaks every 4 months
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate date labels for better visibility

# Print the plot
print(logQ_plt3)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_NH4_plot_sample_TS.png", plot = logQ_plt3, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_NH4.csv")

######
###### 

## Alright BWU?



