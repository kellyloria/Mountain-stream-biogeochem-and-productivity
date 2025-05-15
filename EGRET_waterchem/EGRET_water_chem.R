#### Possible nutrient data infill from EGRET package

# https://github.com/DOI-USGS/EGRET/blob/main/R/readUserDaily.r

# https://cran.r-project.org/web/packages/EGRET/vignettes/EGRET.html


library(EGRET)


############################
savePath<-"/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/egret_chem_infil/" 
saveResults(savePath, eList)

############################
## KAL sample data:
############################

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

BWL_NO3 <- bg_nuts_wide %>%
  filter(site=="BWL") %>%
  dplyr::select(date, NO3_mgL_dl_sw)

BWL_NH4 <- bg_nuts_wide %>%
  filter(site=="BWL") %>%
  dplyr::select(date, NH4_mgL_dl_sw)

BWU_NO3 <- bg_nuts_wide %>%
  filter(site=="BWU") %>%
  dplyr::select(date, NO3_mgL_dl_sw)

BWU_NH4 <- bg_nuts_wide %>%
  filter(site=="BWU") %>%
  dplyr::select(date, NH4_mgL_dl_sw)

GBL_NO3 <- bg_nuts_wide %>%
  filter(site=="GBL") %>%
  dplyr::select(date, NO3_mgL_dl_sw)

GBL_NH4 <- bg_nuts_wide %>%
  filter(site=="GBL") %>%
  dplyr::select(date, NH4_mgL_dl_sw)


GBU_NO3 <- bg_nuts_wide %>%
  filter(site=="GBU") %>%
  dplyr::select(date, NO3_mgL_dl_sw)

GBU_NH4 <- bg_nuts_wide %>%
  filter(site=="GBU") %>%
  dplyr::select(date, NH4_mgL_dl_sw)

###########################
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
## MY sites + my chem 
## NO3?  
############################
# Gather discharge data:
siteNumber <- "10336660" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-09-30"
# Gather sample data:
parameter_cd<-"00631" # nitrate + nitrite as nitrogen is 00631, and for nitrate as nitrogen alone, it is 00618.
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

Sample_BWL <- Sample_C1 %>%
  dplyr::select(!c(date.x, date.y, NO3_mgL_dl_sw))


# Ensure date columns are in Date format
Sample_BWL$Date <- as.Date(Sample_BWL$dateTime)
BWL_NO3$Date <- (BWL_NO3$date)
BWL_NO3$ConcAve <- (BWL_NO3$NO3_mgL_dl_sw)
BWL_NO3$ConcLow <- (BWL_NO3$NO3_mgL_dl_sw)
BWL_NO3$ConcHigh <- (BWL_NO3$NO3_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample_BW <- full_join(Sample, BWL_NO3, by = c("Date", "source", 
                                                "ConcAve", "ConcLow", "ConcHigh"))


Sample <- Sample_BW %>%
  filter(!is.na(Date))

Sample$Julian <- as.numeric(Sample$Date - as.Date("1960-01-01"))
Sample$Day <- lubridate::yday(Sample$Date)
Sample$DecYear <- lubridate::year(Sample$Date) + (lubridate::yday(Sample$Date) - 1) / 365.25
Sample$SinDY <- sin(2 * pi * Sample$Julian / 365.25)
Sample$CosDY <- cos(2 * pi * Sample$Julian / 365.25)
Sample$Month <- month(Sample$Date)
Sample$dateTime <- month(Sample$Date)
Sample$CharacteristicName <- "Inorganic nitrogen (nitrate and nitrite)"
Sample$USGSPCode <- 00631
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


Sample <- Sample %>%
  dplyr::select("Date", "ConcLow", "ConcHigh", "Uncen", "ConcAve",
                "Julian", "Month", "Day", "DecYear", "waterYear",
                "MonthSeq", "SinDY", "CosDY", "dateTime", "CharacteristicName",
                "USGSPCode", "ActivityStartDateTime", "ActivityMediaSubdivisionName",
                "ActivityMediaName", "ResultSampleFractionText", "ResultStatusIdentifier",
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


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

Daily$Q <- (Daily$Q*0.69)
Daily$LogQ <- log(Daily$Q)


# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "BWL"

# Merge discharge with sample data:
eList_BWL1 <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_BWL1)
boxQTwice(eList_BWL1)
plotConcTime(eList_BWL1)
plotConcQ(eList_BWL1)
multiplot <- multiPlotDataOverview(eList_BWL1)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWL1$Sample)
eList_BWL1$INFO$minNumObs <- nrow(eList_BWL1$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_BWL1 <- modelEstimation(
  eList_BWL1,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 100,  # Minimum number of observations
  minNumUncen = 50,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

# Check the result to see if it worked
str(eList_mod_BWL1)



eList_mod_BWL1$Daily
summary(eList_mod_BWL1$Daily)

plotConcTime(eList_mod_BWL1)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWL1) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_BWL1) # Plot Modeled Flux over Time


plotConcPred(eList_mod_BWL1)
plotFluxPred(eList_mod_BWL1)
plotResidPred(eList_mod_BWL1) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWL1) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_BWL1) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_BWL1) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWL1)
boxConcThree(eList_mod_BWL1) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWL1)
plotFluxHist(eList_mod_BWL1)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWL1[["Daily"]]$Date <- as.Date(eList_mod_BWL1[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWL_agg <- eList_mod_BWL1[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_BWL_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt <- ggplot(combined_data, aes(x = Date)) +
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
print(logQ_plt)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_NO3_plot_sample_TS1.png", plot = logQ_plt, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_BWL_NO31.csv")

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
Daily$Q <- (Daily$Q*0.69)
Daily$LogQ <- log(Daily$Q)

# Gather site and parameter information:
Daily$Julian <- as.numeric(Daily$Date - as.Date("1960-01-01"))
Daily$Day <- lubridate::yday(Daily$Date)
Daily$DecYear <- lubridate::year(Daily$Date) + (lubridate::yday(Daily$Date) - 1) / 365.25
Daily$Month <- month(Daily$Date)
Daily$MonthSeq <- (year(Daily$Date) - 1800) * 12 + month(Daily$Date)


# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "BWL"

# Merge discharge with sample data:
eList_BWLa <- mergeReport(INFO, Daily, Sample)

# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd, interactive=FALSE)
INFO$shortName <- "BWL"

############################



############################
# Check sample data:
boxConcMonth(eList_BWLa)
boxQTwice(eList_BWLa)
plotConcTime(eList_BWLa)
plotConcQ(eList_BWLa)
multiPlotDataOverview(eList_BWLa)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWLa$Sample)

eList_BWLa$INFO$minNumObs <- nrow(eList_BWLa$Sample)


eList_mod_BWLa <- modelEstimation(eList_BWLa, windowY = 3, windowQ = 1, windowS = 0.5,
                                 minNumObs = 50, minNumUncen = 25, edgeAdjust = TRUE, verbose = TRUE,
                                 run.parallel = FALSE)


eList_mod_BWLa$Daily
summary(eList_mod_BWLa$Daily)

plotConcTime(eList_mod_BWLa)


# Check for missing values in Daily and Sample data frames
summary(eList_BWLa$Daily)
summary(eList_BWLa$Sample)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWLa) # Plot Modeled Concentration over Time
plotFluxTimeDaily(eList_mod_BWLa) # Plot Modeled Flux over Time
plotConcPred(eList_mod_BWLa)
plotFluxPred(eList_mod_BWLa)
plotResidPred(eList_mod_BWLa) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWLa) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_BWLa) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_BWLa) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWLa)
boxConcThree(eList_mod_BWLa) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWLa)
plotFluxHist(eList_mod_BWLa)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWLa[["Daily"]]$Date <- as.Date(eList_mod_BWLa[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWL_agga <- eList_mod_BWLa[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_BWL_agga, by = "Date", all.x = TRUE)


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
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_BWL_NH4.csv")



###########################
### BWU
###########################
# Gather discharge data:
siteNumber <- "10336660" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-09-30"
# Gather sample data:
parameter_cd<-"00631" #5 digit USGS code
# Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)

Sample_BWU <- Sample_C1 %>%
  dplyr::select(!c(date.x, date.y, NO3_mgL_dl_sw))

BWU_NO3 <- BWU_NO3 %>%
  mutate(source = "KAL")

# Ensure date columns are in Date format
Sample_BWU$Date <- as.Date(Sample_BWU$dateTime)
BWU_NO3$Date <- (BWU_NO3$date)
BWU_NO3$ConcAve <- (BWU_NO3$NO3_mgL_dl_sw)
BWU_NO3$ConcLow <- (BWU_NO3$NO3_mgL_dl_sw)
BWU_NO3$ConcHigh <- (BWU_NO3$NO3_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample <- full_join(Sample_BWU, BWU_NO3, by = c("Date", "source", 
                                            "ConcAve", "ConcLow", "ConcHigh"))


Sample <- Sample %>%
  filter(!is.na(Date))

Sample$Julian <- as.numeric(Sample$Date - as.Date("1960-01-01"))
Sample$Day <- lubridate::yday(Sample$Date)
Sample$DecYear <- lubridate::year(Sample$Date) + (lubridate::yday(Sample$Date) - 1) / 365.25
Sample$SinDY <- sin(2 * pi * Sample$Julian / 365.25)
Sample$CosDY <- cos(2 * pi * Sample$Julian / 365.25)
Sample$Month <- month(Sample$Date)
Sample$dateTime <- month(Sample$Date)
Sample$CharacteristicName <- "Inorganic nitrogen (nitrate and nitrite)"
Sample$USGSPCode <- 00631
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
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:
Daily$Q <- c(Daily$Q*0.444)
Daily$LogQ <-c(log(Daily$Q))
Daily$Julian <- as.numeric(Daily$Date - as.Date("1960-01-01"))
Daily$Day <- lubridate::yday(Daily$Date)
Daily$DecYear <- lubridate::year(Daily$Date) + (lubridate::yday(Daily$Date) - 1) / 365.25
Daily$Month <- month(Daily$Date)
Daily$MonthSeq <- (year(Daily$Date) - 1800) * 12 + month(Daily$Date)

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
INFO$shortName <- "BWU"

# Merge discharge with sample data:
eList_BWU <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_BWU)
boxQTwice(eList_BWU)
plotConcTime(eList_BWU)
plotConcQ(eList_BWU)
multiplot <- multiPlotDataOverview(eList_BWU)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWU$Sample)

eList_BWU$INFO$minNumObs <- nrow(eList_BWU$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_BWU <- modelEstimation(
  eList_BWU,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 90,  # Minimum number of observations
  minNumUncen = 25,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

eList_BWU$Daily
summary(eList_BWU$Daily)

plotConcTime(eList_BWU)




#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWU) # Plot Modeled Concentration over Time
plotFluxTimeDaily(eList_mod_BWU) # Plot Modeled Flux over Time


plotConcPred(eList_mod_BWU)
plotFluxPred(eList_mod_BWU)
plotResidPred(eList_mod_BWU) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWU) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_BWU) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_BWU) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWU)
boxConcThree(eList_mod_BWU) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWU)
plotFluxHist(eList_mod_BWU)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWU[["Daily"]]$Date <- as.Date(eList_mod_BWU[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWU_agg <- eList_mod_BWU[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE)
            )

# Now merge the aggregated model data with the sample data
combined_data_BWU <- merge(Sample, eList_mod_BWU_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt3 <- ggplot(combined_data_BWU, aes(x = Date)) +
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


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWU_NO3_plot_sample_TS.png", plot = logQ_plt3, width = 8, height = 4, units = "in")
# write.csv(combined_data_BWU, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_BWU_NO3.csv")



### BWU NH4

###########################
### BWU
###########################
# Gather discharge data:
siteNumber <- "10336660" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-09-30"
# Gather sample data:
parameter_cd<-"00608" 
SampleB <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
Sample$source <- "nwis"

BWU_NH4 <- BWU_NH4 %>%
  mutate(source = "KAL")

# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
BWU_NH4$Date <- (BWU_NH4$date)
BWU_NH4$ConcAve <- (BWU_NH4$NH4_mgL_dl_sw)
BWU_NH4$ConcLow <- (BWU_NH4$NH4_mgL_dl_sw)
BWU_NH4$ConcHigh <- (BWU_NH4$NH4_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample <- full_join(Sample, BWU_NH4, by = c("Date", "source", 
                                                "ConcAve", "ConcLow", "ConcHigh"))


Sample <-BWU_NH4  %>%
  filter(!is.na(Date))

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
Sample$ActivityMediaSubdivisionName <- NA

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
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:
Daily$Q <- c(Daily$Q*0.444)
Daily$LogQ <-c(log(Daily$Q))
Daily$Julian <- as.numeric(Daily$Date - as.Date("1960-01-01"))
Daily$Day <- lubridate::yday(Daily$Date)
Daily$DecYear <- lubridate::year(Daily$Date) + (lubridate::yday(Daily$Date) - 1) / 365.25
Daily$Month <- month(Daily$Date)
Daily$MonthSeq <- (year(Daily$Date) - 1800) * 12 + month(Daily$Date)

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
INFO$shortName <- "BWU"

# Merge discharge with sample data:
eList_BWUa <- mergeReport(INFO, Daily, Sample)
############################

############################
# Check sample data:
boxConcMonth(eList_BWUa)
boxQTwice(eList_BWUa)
plotConcTime(eList_BWUa)
plotConcQ(eList_BWUa)
multiplot <- multiPlotDataOverview(eList_BWUa)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_BWUa$Sample)

eList_BWUa$INFO$minNumObs <- nrow(eList_BWUa$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_BWUa <- modelEstimation(
  eList_BWUa,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 22,  # Minimum number of observations
  minNumUncen = 15,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

eList_BWUa$Daily
summary(eList_BWUa$Daily)

plotConcTime(eList_BWUa)



#Require Sample + INFO:
plotConcTimeDaily(eList_mod_BWUa) # Plot Modeled Concentration over Time
plotFluxTimeDaily(eList_mod_BWUa) # Plot Modeled Flux over Time

plotConcPred(eList_mod_BWUa)
plotFluxPred(eList_mod_BWUa)
plotResidPred(eList_mod_BWUa) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_BWUa) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_BWUa) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_BWUa) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_BWUa)
boxConcThree(eList_mod_BWUa) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_BWUa)
plotFluxHist(eList_mod_BWUa)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_BWUa[["Daily"]]$Date <- as.Date(eList_mod_BWUa[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_BWU_agga <- eList_mod_BWUa[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE)
  )

# Now merge the aggregated model data with the sample data
combined_data_BWU <- merge(Sample, eList_mod_BWU_agga, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt <- ggplot(combined_data_BWU, aes(x = Date)) +
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
print(logQ_plt)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWU_NH4_plot_sample_TS.png", plot = logQ_plt, width = 8, height = 4, units = "in")
# write.csv(combined_data_BWU, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_BWU_NH4.csv")

#############################################
#############################################

############################
## GBU
## NO3  
############################
# Gather discharge data:
siteNumber <- "10336730" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-10-30"
# Gather sample data:
parameter_cd<-"00631" # nitrate + nitrite as nitrogen is 00631, and for nitrate as nitrogen alone, it is 00618.
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
Sample$source <- "nwis"


# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
GBU_NO3$Date <- as.Date(GBU_NO3$date)
GBU_NO3$ConcAve <- (GBU_NO3$NO3_mgL_dl_sw)


# Add the "source" column to each dataset
Sample <- Sample %>%
  mutate(source = "USGS")

GBU_NO3 <- GBU_NO3 %>%
  mutate(source = "KAL")

Sample_BWL <- Sample_C1 %>%
  dplyr::select(!c(date.x, date.y, NO3_mgL_dl_sw))


# Ensure date columns are in Date format
Sample_BWL$Date <- as.Date(Sample_BWL$dateTime)
GBU_NO3$Date <- (GBU_NO3$date)
GBU_NO3$ConcAve <- (GBU_NO3$NO3_mgL_dl_sw)
GBU_NO3$ConcLow <- (GBU_NO3$NO3_mgL_dl_sw)
GBU_NO3$ConcHigh <- (GBU_NO3$NO3_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample_BW <- full_join(Sample, GBU_NO3, by = c("Date", "source", 
                                               "ConcAve", "ConcLow", "ConcHigh"))


Sample <- GBU_NO3 %>%
  filter(!is.na(Date))

Sample$Julian <- as.numeric(Sample$Date - as.Date("1960-01-01"))
Sample$Day <- lubridate::yday(Sample$Date)
Sample$DecYear <- lubridate::year(Sample$Date) + (lubridate::yday(Sample$Date) - 1) / 365.25
Sample$SinDY <- sin(2 * pi * Sample$Julian / 365.25)
Sample$CosDY <- cos(2 * pi * Sample$Julian / 365.25)
Sample$Month <- month(Sample$Date)
Sample$dateTime <- month(Sample$Date)
Sample$CharacteristicName <- "Inorganic nitrogen (nitrate and nitrite)"
Sample$USGSPCode <- 00631
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
Sample$ActivityMediaSubdivisionName <- NA

Sample <- Sample %>%
  dplyr::select("Date", "ConcLow", "ConcHigh", "Uncen", "ConcAve",
                "Julian", "Month", "Day", "DecYear", "waterYear",
                "MonthSeq", "SinDY", "CosDY", "dateTime", "CharacteristicName",
                "USGSPCode", "ActivityStartDateTime", "ActivityMediaSubdivisionName",
                "ActivityMediaName", "ResultSampleFractionText", "ResultStatusIdentifier",
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


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
INFO$shortName <- "GBU"

# Merge discharge with sample data:
eList_GBU <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_GBU)
boxQTwice(eList_GBU)
plotConcTime(eList_GBU)
plotConcQ(eList_GBU)
multiplot <- multiPlotDataOverview(eList_GBU)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_GBU$Sample)
eList_GBU$INFO$minNumObs <- nrow(eList_GBU$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_GBU <- modelEstimation(
  eList_GBU,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 35,  # Minimum number of observations
  minNumUncen = 20,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

# Check the result to see if it worked
str(eList_mod_GBU)



eList_mod_GBU$Daily
summary(eList_mod_GBU$Daily)

plotConcTime(eList_mod_GBU)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_GBU) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_GBU) # Plot Modeled Flux over Time


plotConcPred(eList_mod_GBU)
plotFluxPred(eList_mod_GBU)
plotResidPred(eList_mod_GBU) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_GBU) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_GBU) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_GBU) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_GBU)
boxConcThree(eList_mod_GBU) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_GBU)
plotFluxHist(eList_mod_GBU)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_GBU[["Daily"]]$Date <- as.Date(eList_mod_GBU[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_GBU_agg <- eList_mod_GBU[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_GBU_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt <- ggplot(combined_data, aes(x = Date)) +
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
print(logQ_plt)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_GBU_NO3_plot_sample_TS.png", plot = logQ_plt, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_GBU_NO3.csv")

###
###

############################
## GBL
## NH4  
############################
# Gather discharge data:
siteNumber <- "10336730" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-10-30"
# Gather sample data:
parameter_cd<-"00608" # nitrate + nitrite as nitrogen is 00631, and for nitrate as nitrogen alone, it is 00618.
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
Sample$source <- "nwis"


# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
GBU_NH4$Date <- as.Date(GBU_NH4$date)
GBU_NH4$ConcAve <- (GBU_NH4$NO3_mgL_dl_sw)


# Add the "source" column to each dataset
Sample <- Sample %>%
  mutate(source = "USGS")

GBU_NH4 <- GBU_NH4 %>%
  mutate(source = "KAL")

Sample_BWL <- Sample_C1 %>%
  dplyr::select(!c(date.x, date.y, NO3_mgL_dl_sw))


# Ensure date columns are in Date format
Sample_BWL$Date <- as.Date(Sample_BWL$dateTime)
GBU_NH4$Date <- (GBU_NH4$date)
GBU_NH4$ConcAve <- (GBU_NH4$NH4_mgL_dl_sw)
GBU_NH4$ConcLow <- (GBU_NH4$NH4_mgL_dl_sw)
GBU_NH4$ConcHigh <- (GBU_NH4$NH4_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample_BW <- full_join(Sample, GBU_NH4, by = c("Date", "source", 
                                               "ConcAve", "ConcLow", "ConcHigh"))


Sample <- GBU_NH4 %>%
  filter(!is.na(Date))

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
Sample$ActivityMediaSubdivisionName <- NA

Sample <- Sample %>%
  dplyr::select("Date", "ConcLow", "ConcHigh", "Uncen", "ConcAve",
                "Julian", "Month", "Day", "DecYear", "waterYear",
                "MonthSeq", "SinDY", "CosDY", "dateTime", "CharacteristicName",
                "USGSPCode", "ActivityStartDateTime", "ActivityMediaSubdivisionName",
                "ActivityMediaName", "ResultSampleFractionText", "ResultStatusIdentifier",
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

startDate <- min(as.character(Sample$Date)) 
# Gather discharge data:
Daily <- readNWISDaily(siteNumber,"00060",startDate,endDate)
# Gather site and parameter information:



# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "GBU"

# Merge discharge with sample data:
eList_GBUb <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_GBUb)
boxQTwice(eList_GBUb)
plotConcTime(eList_GBUb)
plotConcQ(eList_GBUb)
multiplot <- multiPlotDataOverview(eList_GBUb)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_GBUb$Sample)
eList_GBUb$INFO$minNumObs <- nrow(eList_GBUb$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_GBUb <- modelEstimation(
  eList_GBUb,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 35,  # Minimum number of observations
  minNumUncen = 20,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

# Check the result to see if it worked
str(eList_mod_GBUb)



eList_mod_GBUb$Daily
summary(eList_mod_GBUb$Daily)

plotConcTime(eList_mod_GBUb)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_GBUb) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_GBUb) # Plot Modeled Flux over Time


plotConcPred(eList_mod_GBUb)
plotFluxPred(eList_mod_GBUb)
plotResidPred(eList_mod_GBUb) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_GBUb) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_GBUb) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_GBUb) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_GBUb)
boxConcThree(eList_mod_GBUb) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_GBUb)
plotFluxHist(eList_mod_GBUb)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_GBUb[["Daily"]]$Date <- as.Date(eList_mod_GBUb[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_GBU_aggb<- eList_mod_GBUb[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_GBU_aggb, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt <- ggplot(combined_data, aes(x = Date)) +
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
print(logQ_plt)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_GBU_NH4_plot_sample_TS.png", plot = logQ_plt, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_GBU_NH4.csv")




############################
## GBU
## NO3  
############################

# Gather discharge data:
siteNumber <- "10336730" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-10-30"
# Gather sample data:
parameter_cd<-"00631" # nitrate + nitrite as nitrogen is 00631, and for nitrate as nitrogen alone, it is 00618.
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
Sample$source <- "nwis"


# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
GBU_NH4$Date <- as.Date(GBU_NH4$date)
GBU_NH4$ConcAve <- (GBU_NH4$NO3_mgL_dl_sw)


# Add the "source" column to each dataset
Sample <- Sample %>%
  mutate(source = "USGS")

GBU_NH4 <- GBU_NH4 %>%
  mutate(source = "KAL")

Sample_BWL <- Sample_C1 %>%
  dplyr::select(!c(date.x, date.y, NO3_mgL_dl_sw))


# Ensure date columns are in Date format
Sample_BWL$Date <- as.Date(Sample_BWL$dateTime)
GBU_NH4$Date <- (GBU_NH4$date)
GBU_NH4$ConcAve <- (GBU_NH4$NO3_mgL_dl_sw)
GBU_NH4$ConcLow <- (GBU_NH4$NO3_mgL_dl_sw)
GBU_NH4$ConcHigh <- (GBU_NH4$NO3_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample_BW <- full_join(Sample, GBU_NH4, by = c("Date", "source", 
                                               "ConcAve", "ConcLow", "ConcHigh"))


Sample <- GBU_NH4 %>%
  filter(!is.na(Date))

Sample$Julian <- as.numeric(Sample$Date - as.Date("1960-01-01"))
Sample$Day <- lubridate::yday(Sample$Date)
Sample$DecYear <- lubridate::year(Sample$Date) + (lubridate::yday(Sample$Date) - 1) / 365.25
Sample$SinDY <- sin(2 * pi * Sample$Julian / 365.25)
Sample$CosDY <- cos(2 * pi * Sample$Julian / 365.25)
Sample$Month <- month(Sample$Date)
Sample$dateTime <- month(Sample$Date)
Sample$CharacteristicName <- "Inorganic nitrogen (nitrate and nitrite)"
Sample$USGSPCode <- 00631
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
Sample$ActivityMediaSubdivisionName <- NA

Sample <- Sample %>%
  dplyr::select("Date", "ConcLow", "ConcHigh", "Uncen", "ConcAve",
                "Julian", "Month", "Day", "DecYear", "waterYear",
                "MonthSeq", "SinDY", "CosDY", "dateTime", "CharacteristicName",
                "USGSPCode", "ActivityStartDateTime", "ActivityMediaSubdivisionName",
                "ActivityMediaName", "ResultSampleFractionText", "ResultStatusIdentifier",
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


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

Daily$Q <- (Daily$Q*0.812)
Daily$LogQ <- log(Daily$Q)


# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "GBL"

# Merge discharge with sample data:
eList_GBU <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_GBU)
boxQTwice(eList_GBU)
plotConcTime(eList_GBU)
plotConcQ(eList_GBU)
multiplot <- multiPlotDataOverview(eList_GBU)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_GBU$Sample)
eList_GBU$INFO$minNumObs <- nrow(eList_GBU$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_GBU <- modelEstimation(
  eList_GBU,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 35,  # Minimum number of observations
  minNumUncen = 20,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

# Check the result to see if it worked
str(eList_mod_GBU)



eList_mod_GBU$Daily
summary(eList_mod_GBU$Daily)

plotConcTime(eList_mod_GBU)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_GBU) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_GBU) # Plot Modeled Flux over Time


plotConcPred(eList_mod_GBU)
plotFluxPred(eList_mod_GBU)
plotResidPred(eList_mod_GBU) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_GBU) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_GBU) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_GBU) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_GBU)
boxConcThree(eList_mod_GBU) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_GBU)
plotFluxHist(eList_mod_GBU)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_GBU[["Daily"]]$Date <- as.Date(eList_mod_GBU[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_GBU_agg <- eList_mod_GBU[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_GBU_agg, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt <- ggplot(combined_data, aes(x = Date)) +
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
print(logQ_plt)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_GBU_NH4_plot_sample_TS.png", plot = logQ_plt, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_GBU_NH4.csv")

###
###

############################
## GBL
## NH4  
############################
# Gather discharge data:
siteNumber <- "10336730" #Choptank River near Greensboro, MD
startDate <- "2020-09-30" #Gets earliest date
endDate <- "2024-10-30"
# Gather sample data:
parameter_cd<-"00608" # nitrate + nitrite as nitrogen is 00631, and for nitrate as nitrogen alone, it is 00618.
Sample <- readNWISSample(siteNumber,parameter_cd,startDate,endDate)
Sample$source <- "nwis"


# Ensure date columns are in Date format
Sample$Date <- as.Date(Sample$dateTime)
GBU_NH4$Date <- as.Date(GBU_NH4$date)
GBU_NH4$ConcAve <- (GBU_NH4$NO3_mgL_dl_sw)


# Add the "source" column to each dataset
Sample <- Sample %>%
  mutate(source = "USGS")

GBU_NH4 <- GBU_NH4 %>%
  mutate(source = "KAL")

Sample_BWL <- Sample_C1 %>%
  dplyr::select(!c(date.x, date.y, NO3_mgL_dl_sw))


# Ensure date columns are in Date format
Sample_BWL$Date <- as.Date(Sample_BWL$dateTime)
GBU_NH4$Date <- (GBU_NH4$date)
GBU_NH4$ConcAve <- (GBU_NH4$NH4_mgL_dl_sw)
GBU_NH4$ConcLow <- (GBU_NH4$NH4_mgL_dl_sw)
GBU_NH4$ConcHigh <- (GBU_NH4$NH4_mgL_dl_sw)


# Perform the full join by matching the date columns
Sample_BW <- full_join(Sample, GBU_NH4, by = c("Date", "source", 
                                               "ConcAve", "ConcLow", "ConcHigh"))


Sample <- GBU_NH4 %>%
  filter(!is.na(Date))

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
Sample$ActivityMediaSubdivisionName <- NA

Sample <- Sample %>%
  dplyr::select("Date", "ConcLow", "ConcHigh", "Uncen", "ConcAve",
                "Julian", "Month", "Day", "DecYear", "waterYear",
                "MonthSeq", "SinDY", "CosDY", "dateTime", "CharacteristicName",
                "USGSPCode", "ActivityStartDateTime", "ActivityMediaSubdivisionName",
                "ActivityMediaName", "ResultSampleFractionText", "ResultStatusIdentifier",
                "ResultValueTypeName") %>%
  group_by(Date, CharacteristicName, ActivityStartDateTime, 
           ActivityMediaSubdivisionName, ActivityMediaName, ResultSampleFractionText,
           ResultStatusIdentifier,ResultValueTypeName) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


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

Daily$Q <- (Daily$Q*0.812)
Daily$LogQ <- log(Daily$Q)


# Here user must input some values:
INFO<- readNWISInfo(siteNumber,parameter_cd)
INFO$shortName <- "GBL"

# Merge discharge with sample data:
eList_GBUa <- mergeReport(INFO, Daily, Sample)
############################


############################
# Check sample data:
boxConcMonth(eList_GBUa)
boxQTwice(eList_GBUa)
plotConcTime(eList_GBUa)
plotConcQ(eList_GBUa)
multiplot <- multiPlotDataOverview(eList_GBUa)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_BWL_TIN_plot_sample.png", plot = multiplot, width = 5.25, height = 6, units = "in")



###########################
# Run WRTDS model:
nrow(eList_GBUa$Sample)
eList_GBUa$INFO$minNumObs <- nrow(eList_GBUa$Sample)


# Run modelEstimation with adjusted arguments
eList_mod_GBUa <- modelEstimation(
  eList_GBUa,
  windowY = 3,  # Adjust year window as needed
  windowQ = 1,  # Adjust flow window as needed
  windowS = 0.5,  # Adjust stage window as needed
  minNumObs = 35,  # Minimum number of observations
  minNumUncen = 20,  # Minimum uncensored observations
  edgeAdjust = TRUE,  # Edge adjustment flag
  verbose = TRUE,  # Verbosity flag
  run.parallel = FALSE  # Parallel processing flag
)

# Check the result to see if it worked
str(eList_mod_GBUa)



eList_mod_GBUa$Daily
summary(eList_mod_GBUa$Daily)

plotConcTime(eList_mod_GBUa)

#Require Sample + INFO:
plotConcTimeDaily(eList_mod_GBUa) # Plot Modeled Concentration over Time

plotFluxTimeDaily(eList_mod_GBUa) # Plot Modeled Flux over Time


plotConcPred(eList_mod_GBUa)
plotFluxPred(eList_mod_GBUa)
plotResidPred(eList_mod_GBUa) # Plots model residuals (differences between observed and predicted values) against predicted concentrations.
plotResidQ(eList_mod_GBUa) # Plots residuals against discharge (Q) to assess whether model errors vary with flow.
plotConcHist(eList_mod_GBUa) # Flow-Normalized Concentration Trends

plotResidTime(eList_mod_GBUa) # Plots residuals over time to check for temporal trends in model errors.
boxResidMonth(eList_mod_GBUa)
boxConcThree(eList_mod_GBUa) # Plots boxplots of concentration grouped by flow conditions (low, mid, high).


#Require Daily + INFO:
plotConcHist(eList_mod_GBUa)
plotFluxHist(eList_mod_GBUa)

### for me 
# Merge the sample and modeled data
# Convert Sample Date to Date format (if not already)
Sample$Date <- as.Date(Sample$Date)

# Aggregate modeled data to match the sample data's Date range
eList_mod_GBUa[["Daily"]]$Date <- as.Date(eList_mod_GBUa[["Daily"]][["Date"]])  # Ensure the modeled data has a Date column
eList_mod_GBU_agga<- eList_mod_GBUa[["Daily"]] %>%
  group_by(Date) %>%
  summarize(ConcDay_mod = mean(ConcDay, na.rm = TRUE),
            FluxDay_mod = mean(FluxDay, na.rm = TRUE))

# Now merge the aggregated model data with the sample data
combined_data <- merge(Sample, eList_mod_GBU_agga, by = "Date", all.x = TRUE)


# Create the plot
logQ_plt <- ggplot(combined_data, aes(x = Date)) +
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
print(logQ_plt)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/EGRET_GBU_NH4_plot_sample_TS.png", plot = logQ_plt, width = 8, height = 4, units = "in")
# write.csv(combined_data, "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/EGRET_GBU_NH4.csv")


