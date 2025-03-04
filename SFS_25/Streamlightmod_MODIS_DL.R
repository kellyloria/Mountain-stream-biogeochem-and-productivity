
## Folling vinegette form https://psavoy.github.io/StreamLight/articles/2%20Download%20and%20process%20MODIS%20LAI.html

setwd("~/Documents/Publications/CH1 biogeochem linkages")


#Set the download location (add your own directory)
working_dir <- "~/Documents/Publications/CH1 biogeochem linkages/data/Temp_NLDAS"

save_dir <- "~/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/Stream_light_driver"
## ======================================== 
##  Temporary example of to get NLDAS headers: 
# #Download NLDAS data at NC_NHC
# NLDAS_DL(
#   save_dir = working_dir,
#   Site_ID = "NC_NHC",
#   Lat = 35.9925, 
#   Lon = -79.0460, 
#   startDate = "2017-01-01"
# )
# 
# #Process the downloaded data
# NLDAS_processed <- NLDAS_proc(
#   read_dir = working_dir, 
#   Site_IDs = "NC_NHC"
# )

## ======================================== 
#Use the devtools packge to install StreamLightUtils
#devtools::install_github("psavoy/StreamLightUtils")
#devtools::install_github("psavoy/StreamLight")
library("StreamLightUtils")
library("StreamLight")
library("streamMetabolizer")

## ======================================== 
## BW dat:
## where GBU cord - is actually GBL here:
MOD_unpack <- AppEEARS_unpack_QC(
  zip_file = "MODIS_streamDat.zip", 
  zip_dir = "./data", 
  c("BWL", "BWU", "GBU")
)

MOD_processed <- AppEEARS_proc(
  unpacked_LAI = MOD_unpack,  
  fit_method = "Gu", 
  plot = TRUE
)


## ======================================== 
## GBU dat:
MOD_unpack_GBU <- AppEEARS_unpack_QC(
  zip_file = "MODIS_GB_streamDat.zip", 
  zip_dir = "./data", 
  c("GBL","GBU")
)

MOD_processed <- AppEEARS_proc(
  unpacked_LAI = MOD_unpack_GBU,  
  fit_method = "Gu", 
  plot = TRUE
)

### Read in light data from 
##  Make timesireies empty dataframe.

# Create a sequence of dates
date_seq <- seq(from = as.POSIXct("2021-03-20"), to = as.POSIXct("2024-10-10"), by = "15 min")

# Convert the sequence to a dataframe
date_df <- data.frame(datetime = date_seq)
date_df$site <- "BWL"
date_df$catch <- "BW"
date_df1 <- data.frame(datetime = date_seq)
date_df1$site <- "GBL"
date_df1$catch <- "GB"
date_df2 <- data.frame(datetime = date_seq)
date_df2$site <- "GBU"
date_df2$catch <- "GB"
date_df3 <- data.frame(datetime = date_seq)
date_df3$site <- "BWU"
date_df3$catch <- "BW"

date_datf <- rbind(date_df,date_df1,date_df2, date_df3)


## A. the server par file. 
## B. inflill missing data with stream metabolizer. 

light_dat <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/NLDAS/stream/24stream_NLDAS_light.rds") %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") 

light_dat_BWL <- light_dat %>% 
  dplyr::filter(site =="BWL") %>%
  mutate(par= c(light* 2.114))

light_dat_BWU <- light_dat %>% 
  dplyr::filter(site =="BWU") %>%
  mutate(par= c(light* 2.114))

light_dat_GBU <- light_dat %>% 
  dplyr::filter(site =="GBU") %>%
  mutate(par= c(light* 2.114))

light_dat_GBL <- light_dat %>% 
  dplyr::filter(site =="GBL") %>%
  mutate(par= c(light* 2.114))


NLDAS_dat <- rbind(light_dat_BWL, light_dat_BWU, light_dat_GBU)


## fill in light from stream metabolizer
BWL_dat <- date_df %>%
  mutate(solar.time = calc_solar_time(datetime, -120.164335),
         calclight = calc_light(solar.time, 39.10754, -120.164335)) #,light_in = ifelse(is.na(par), calclight, par))


BWU_dat <- date_df3 %>%
  mutate(solar.time = calc_solar_time(datetime, -120.1957),
         calclight = calc_light(solar.time, 39.1052, -120.1957))

## Real upper 
GBU_dat <- date_df2 %>%
  mutate(solar.time = calc_solar_time(datetime, -119.937181),
         calclight = calc_light(solar.time, 39.08679, -119.937181))
# 
# 
# GBU_dat <- date_df2 %>%
#   mutate(solar.time = calc_solar_time(datetime, -119.930379),
#          calclight = calc_light(solar.time, 39.088427, -119.930379))
# 



SM_light <- rbind(BWL_dat,BWU_dat,GBU_dat)


light_dat<- SM_light%>%
  left_join(NLDAS_dat, by = c("site", "datetime")) %>%
  mutate(light_in = ifelse(is.na(par), calclight, par))


TS_plot <- ggplot(light_dat, aes(x = datetime, y = light_in, color = site, shape = site)) +
  geom_point(size = 3, alpha = 0.95, position = position_dodge(width = 0.3)) +
  scale_color_viridis_d(option = "viridis") +
  geom_point(size = 2, alpha = 0.7) +
  scale_shape_manual(values = c(15, 0, 17, 2)) +
  theme_classic() + facet_grid(catch~.)


### To create NLDAS_processed object 
## We need col for Year, DOY, Hour, and SW
### to par.to.sw.base(par, coeff=0.473)
NLDAS_proc <- light_dat %>%
  mutate(Year = year(datetime),
         DOY = yday(datetime), 
         Hour = hour(datetime),
         SW = light_in *0.473) %>%
  group_by(site, Year, DOY, Hour) %>%
  summarize(SW = mean(SW, na.rm=T)) %>%
  dplyr::rename(Site_ID ="site") %>%
  mutate(
    local_time = as.POSIXct(paste(Year, DOY, Hour), format = "%Y %j %H", tz = "America/Los_Angeles"),
    offset = -8,  # Pacific Standard Time (adjust for daylight saving if needed)
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    SW_inc = SW  # Ensure this column exists (NLDAS incoming shortwave radiation)
  ) %>%
  dplyr::select(Site_ID,Year, DOY, Hour, SW)

  ##

str(NLDAS_proc)

NLDAS_proc <- as.data.frame(NLDAS_proc)

NLDAS_processed <- split(NLDAS_proc, NLDAS_proc$Site_ID)

NLDAS_processed <- lapply(NLDAS_processed, function(x) x[, !colnames(x) %in% "Site_ID"])


str(MOD_processed)

## step 2 make driver file 
# make_driver(site_locs, NLDAS_processed, MOD_processed, write_output, save_dir)

site_locs <- read.csv("/Users/kellyloria/Documents/Old mac transfer files 20241229/Downloads/TahoeCords_MODIS.csv") 

site_locs  <- site_locs%>%
  dplyr::filter(!site =="GBL") %>%
  dplyr::rename(Site_ID ="site", Lat="lat", Lon = "long")

# Convert to an sf object with EPSG:4326 (WGS84)


site_locs <- st_as_sf(site_locs, coords = c("Lat", "Lon"), crs = 4326)
site_locs$epsg_crs <- 4326


names(MOD_processed) <- c("BWL","BWU", "GBU")

#site_locs <- as.data.frame(site_locs)


# Function to trim MOD_processed to match NLDAS_processed
trim_MOD_to_NLDAS <- function(NLDAS_processed, MOD_processed) {
  trimmed_MOD <- list()
  
  for (site in names(NLDAS_processed)) {
    if (site %in% names(MOD_processed)) {
      # Extract unique Year and DOY from NLDAS for the site
      nldas_dates <- unique(NLDAS_processed[[site]][, c("Year", "DOY")])
      
      # Filter MOD_processed to match NLDAS time range
      trimmed_MOD[[site]] <- MOD_processed[[site]] %>%
        dplyr::semi_join(nldas_dates, by = c("Year", "DOY"))
    }
  }
  
  return(trimmed_MOD)
}

# Apply the function
MOD_trimmed <- trim_MOD_to_NLDAS(NLDAS_processed, MOD_processed)

# Check the structure after trimming
str(MOD_trimmed)



MOD_trimmed2 <- lapply(MOD_trimmed, function(x) {
  x %>%
    mutate(
      jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
      LAI = zoo::na.approx(Lai, rule = 2))  # Interpolate missing values if needed
})

single_site_driver <- make_driver(site_locs,  NLDAS_processed,
                                  MOD_trimmed)




tahoe_stream_driver <-make_driver(site_locs, 
                                   NLDAS_processed, 
                                   MOD_trimmed2,  
                                   write_output = FALSE,
                                   save_dir = NULL)


str(site_locs)

str(NLDAS_processed)

str(MOD_trimmed2)

# Extract Year and DOY from jday in MOD_trimmed2
MOD_trimmed2 <- lapply(MOD_trimmed2, function(df) {
  df$Year <- as.integer(substr(df$jday, 1, 4)) # Extract first 4 digits as Year
  df$DOY <- as.integer(df$DOY)                 # Ensure DOY is integer
  df$LAI <- as.numeric(df$LAI)                 # Ensure LAI is numeric
  df <- df[, c("Year", "DOY", "LAI")]          # Keep relevant columns
  return(df)
})


source("./make_driver_mod.R")
tahoe_stream_driver <- make_driver(site_locs, 
                                   NLDAS_processed, 
                                   MOD_trimmed2,  
                                   write_output = FALSE,
                                   save_dir = NULL)


lapply(NLDAS_processed, function(x) summary(x))

lapply(MOD_trimmed2, function(x) summary(x))


str(site_locs)
str(NLDAS_processed)
str(MOD_trimmed2)

test_merge_BWL <- merge(NLDAS_processed$BWL, MOD_trimmed2$BWL, by = c("Year", "DOY"), all = TRUE)

test_merge_BWU <- merge(NLDAS_processed$BWU, MOD_trimmed2$BWU, by = c("Year", "DOY"), all = TRUE)

test_merge_GBU <- merge(NLDAS_processed$GBU, MOD_trimmed2$GBU, by = c("Year", "DOY"), all = TRUE)


test_merge_BWL2<- test_merge_BWL%>%
  mutate(
    local_time = as.POSIXct(paste(Year, DOY, Hour), format = "%Y %j %H", tz = "America/Los_Angeles"),
    offset = -8,  # Adjust if necessary based on daylight saving
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    SW_inc = SW  # Ensure correct shortwave radiation column exists
  ) %>%
 dplyr:: select(local_time, offset, jday, DOY, Hour, SW_inc, LAI)



test_merge_BWU2<- test_merge_BWU%>%
  mutate(
    local_time = as.POSIXct(paste(Year, DOY, Hour), format = "%Y %j %H", tz = "America/Los_Angeles"),
    offset = -8,  # Adjust if necessary based on daylight saving
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    SW_inc = SW  # Ensure correct shortwave radiation column exists
  ) %>%
  dplyr:: select(local_time, offset, jday, DOY, Hour, SW_inc, LAI)





test_merge_GBU2<- test_merge_GBU%>%
  mutate(
    local_time = as.POSIXct(paste(Year, DOY, Hour), format = "%Y %j %H", tz = "America/Los_Angeles"),
    offset = -8,  # Adjust if necessary based on daylight saving
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    SW_inc = SW  # Ensure correct shortwave radiation column exists
  ) %>%
  dplyr:: select(local_time, offset, jday, DOY, Hour, SW_inc, LAI)


str(test_merge_BWL)

extract_height(
  Site_ID = site_locs[Site_ID, "BWL"], 
  Lat = site_locs[Lat, "39.10754"],
  Lon = site_locs[Lon, "120.1647"],
  site_crs = site_locs[site_crs, "epsg_crs"]
)

extract_height(
  Site_ID ="BWL", 
  Lat = 39.10754,
  Lon = -120.1647,
  site_crs = "epsg_crs"
)

#Load the example driver file for NC_NHC


### YOU ARE HERE :: 
# Run the model
BWL_SL_modeled <- stream_light(
  test_merge_BWL2, 
  Lat = 39.1052, 
  Lon = -120.1957, 
  channel_azimuth = 290, 
  bottom_width = 8.1, 
  BH = 1, 
  BS = 100, 
  WL = 0.62, 
  TH = 23, 
  overhang = 2.3, 
  overhang_height = NA, 
  x_LAD = 1
)

BWL_SL_modeled$site <- "BWL"

# saveRDS(BWL_SL_modeled, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/CH1_StreamLightMod_BWL.rds")




# Run the model
BWU_SL_modeled <- stream_light(
  test_merge_BWU2, 
  Lat = 39.1052, 
  Lon = -120.1957, 
  channel_azimuth = 315, 
  bottom_width = 5.8, 
  BH = 1, 
  BS = 100, 
  WL = 1.2, 
  TH = 23, 
  overhang = 2.3, 
  overhang_height = NA, 
  x_LAD = 1
)

BWU_SL_modeled$site <- "BWU"


# saveRDS(BWU_SL_modeled, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/CH1_StreamLightMod_BWU.rds")


### Get Glenbrook: 
# Run the model
GBL_SL_modeled <- stream_light(
  test_merge_GBU2, 
  Lat = 39.1052, 
  Lon = -120.1957, 
  channel_azimuth = 90, 
  bottom_width = 1.67, 
  BH = 1, 
  BS = 100, 
  WL = 0.53, 
  TH = 23, 
  overhang = 2.3, 
  overhang_height = NA, 
  x_LAD = 1
)


GBL_SL_modeled$site <- "GBL"

# saveRDS(GBL_SL_modeled, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/CH1_StreamLightMod_GBL.rds")


light_mod_df <- rbind(BWL_SL_modeled, BWU_SL_modeled, GBU_SL_modeled, GBL_SL_modeled)

GBU_SL_modeled <- stream_light(
  test_merge_GBU2, 
  Lat = 39.1052, 
  Lon = -120.1957, 
  channel_azimuth = 85, 
  bottom_width = 1.5, 
  BH = 1, 
  BS = 100, 
  WL = 0.19, 
  TH = 23, 
  overhang = 2.3, 
  overhang_height = NA, 
  x_LAD = 1
)
GBU_SL_modeled$site <- "GBU"


# saveRDS(GBU_SL_modeled, file = "/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/CH1_StreamLightMod_GBU.rds")


light_mod_df<- light_mod_df%>%
  mutate(date=as.Date(local_time))
light_mod_df<-water_year(light_mod_df)

light_mod_df_day <- light_mod_df%>%
  group_by(site, date, water_year)%>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(DOY = yday(date),
         year = year(date))
  

TS_plot <- ggplot(light_mod_df_day%>%dplyr::filter(water_year<2025), aes(x = DOY, y = PAR_surface, color = site, shape = site)) +
  geom_point(size = 2, alpha = 0.5, position = position_dodge(width = 0.3)) +
 # scale_color_viridis_d(option = "viridis") +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = siteC_colors) +
  scale_shape_manual(values = c(15, 0, 17, 2)) +
  ylab(expression(PAR~at~stream~surface~(mu~mol~m^-2~s^-1))) +
  theme_classic() + facet_grid(year~.)

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/Sup_Figures/Stream_light_TS.png", plot = TS_plot, width = 7.75, height = 8.5, units = "in")


#########
########


TS_plot <- ggplot(GBU_SL_modeled, aes(x = local_time, y = PAR_surface)) +
  geom_point(size = 3, alpha = 0.95, position = position_dodge(width = 0.3)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_shape_manual(values = c(15, 0, 17, 2)) +
  theme_classic() 


TS_plot <- ggplot(BWU_SL_modeled, aes(x = local_time, y = PAR_surface)) +
  geom_point(size = 3, alpha = 0.95, position = position_dodge(width = 0.3)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_shape_manual(values = c(15, 0, 17, 2)) +
  theme_classic() 




######################
library(dplyr)
library(lubridate)

# Format NLDAS data
NLDAS_BWL <- NLDAS_processed$BWL %>%
  mutate(
    local_time = as.POSIXct(paste(Year, DOY, Hour), format = "%Y %j %H", tz = "America/Los_Angeles"),
    offset = -8,  # Pacific Standard Time (adjust for daylight saving if needed)
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    SW_inc = SW  # Ensure this column exists (NLDAS incoming shortwave radiation)
  ) %>%
  dplyr::select(local_time, offset, jday, DOY, Hour, SW_inc)

# Format MODIS data (interpolating to daily LAI values if needed)
MODIS_BWL <- MOD_processed$BWL %>%
  mutate(
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    LAI = zoo::na.approx(Lai, rule = 2)  # Interpolate missing values if needed
  ) %>%
  dplyr::select(jday, DOY, LAI)


driver_BWL <- NLDAS_BWL %>%
  left_join(MODIS_BWL %>% dplyr::select(jday, LAI), by = "jday") %>%
  dplyr::select(local_time, offset, jday, DOY, Hour, SW_inc, LAI)

str(MODIS_BWL)

str(NLDAS_BWL)


site_BWL  <- site_locs%>%
  dplyr::filter(site =="BWL") %>%
  dplyr::rename(Site_ID ="site", Lat="lat", Lon = "long")

# Convert to an sf object with EPSG:4326 (WGS84)
site_BWL <- st_as_sf(site_BWL, coords = c("Lat", "Lon"), crs = 4326)
site_BWL$epsg_crs <- 4326

names(MOD_processed) <- c("BWL","BWU", "GBU")


site_BWL <-  site_locs%>%filter(Site_ID=="BWL")

single_site_driver <- make_driver(
  site_locs = site_BWL,  # Filter for BWL site
  NLDAS_processed = NLDAS_BWL,  
  MOD_processed = MODIS_BWL,  
  write_output = TRUE, 
  save_dir = save_dir
)



MODIS_BWL$Site_ID <- site_locs$Site_ID[match(
  paste(MODIS_BWL$Lat, MODIS_BWL$Lon),
  paste(site_locs$Lat, site_locs$Lon)
)]


NLDAS_BWL$Site_ID <- site_locs$Site_ID[match(
  paste(NLDAS_BWL$Lat, NLDAS_BWL$Lon),
  paste(site_locs$Lat, site_locs$Lon)
)]



NLDAS_BWL <- NLDAS_processed$BWL %>%
  mutate(
    local_time = as.POSIXct(paste(Year, DOY, Hour), format = "%Y %j %H", tz = "America/Los_Angeles"),
    offset = -8,  # Adjust if necessary based on daylight saving
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    SW_inc = SW  # Ensure correct shortwave radiation column exists
  ) %>%
  select(local_time, offset, jday, DOY, Hour, SW_inc)

MODIS_BWL <- MOD_processed$BWL %>%
  mutate(
    jday = as.integer(paste0(Year, sprintf("%03d", DOY))),  # YYYYDDD format
    LAI = zoo::na.approx(Lai, rule = 2)  # Interpolate missing values if needed
  ) %>%
  select(jday, DOY, LAI)


single_site_driver <- make_driver(
  site_locs = site_BWL,  # Filter for BWL site
  NLDAS_processed = (NLDAS_BWL),  # Pass the formatted NLDAS data
  MOD_processed = (MODIS_BWL),  # Pass the formatted MODIS data
  write_output = TRUE,  # Set to TRUE if you want to save the driver file
  save_dir = save_dir  # Define the save directory
)



## Gu methods prevents the interpolation of the final ~30 days of 2016 in nwis11044000
## fill in the missing LAI_proc data with a final simple linear interpolation
MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date == "2017-01-01"),]$LAI_proc <- 0.2
nrow(MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date >= "2016-12-02" & MOD_processed$nwis11044000$Date <= "2017-01-01"),])


interp_LAI <- rev(approx(x = c(MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date == "2016-12-02"),]$LAI_proc,
                               MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date == "2017-01-01"),]$LAI_proc),
                         y = c(0,0),
                         n = 31)$x)
## replace
MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date >= "2016-12-02" & MOD_processed$nwis11044000$Date <= "2017-01-01"),]$LAI_proc <- interp_LAI
View(MOD_processed$nwis11044000) ## scroll to bottom to confirm
