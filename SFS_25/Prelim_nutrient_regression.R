lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2", "ggpubr", "ggpattern",
         "lme4", "lmerTest", "MuMIn", "PerformanceAnalytics", "car"), require, character.only=T)
site_colors <- c(
  "BWL" = "#054fb9",
  "BWU" = "#97b5c2",
  "GBL" = "#DD6E42",
  "GBU" = "#d9cb7c"
)


catch_colors <- c(
  "BW" = "#054fb9",
  "GB" = "#DD6E42"
)

# Create a sequence of dates
date_seq <- seq(from = as.Date("2021-03-20"), to = as.Date("2024-10-10"), by = "day")

# Convert the sequence to a dataframe
date_df <- data.frame(date = date_seq)
date_df$site <- "BWL"
date_df$catch <- "BW"
date_df1 <- data.frame(date = date_seq)
date_df1$site <- "GBL"
date_df1$catch <- "GB"
date_df2 <- data.frame(date = date_seq)
date_df2$site <- "GBU"
date_df2$catch <- "GB"
date_df3 <- data.frame(date = date_seq)
date_df3$site <- "BWU"
date_df3$catch <- "BW"

date_datf <- rbind(date_df,date_df1,date_df2, date_df3)

### fxns:
## Fxn for water year:
water_year <- function(data) {
  # Ensure the `date` column exists
  if (!"date" %in% names(data)) {
    stop("The dataframe must contain a column named 'date'.")
  }
  
  # Convert the `date` column to Date format (if not already)
  data$date <- as.Date(data$date)
  
  # Add the water_year column
  data$water_year <- with(data, {
    year <- as.numeric(format(date, "%Y"))
    month <- as.numeric(format(date, "%m"))
    ifelse(month >= 10, year + 1, year)
  })
  
  return(data)
}
## Function to perform ANOVA and post-hoc test within each site
comparisons <- list(
  c("2021", "2022"),
  c("2021", "2023"),
  c("2022", "2023"),
  c("2023", "2024")
)

####
date_datf <- water_year(date_datf)

covariat_dat <- readRDS("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/data/CH1_covariate_dat.rds") %>%
  dplyr::select("date", "Site","bulk.density","AFDM_mgg","AFDM_mgcm2", "Chla_ugL_Q", "Pheo_ugL_Q")
summary(covariat_dat)

covariat_dat <- date_datf %>%
  left_join(covariat_dat, by=c("date", "site"="Site"))

bg_nuts <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/NS_chem_dat_nh4_24.rds") %>%
  #mutate(date = as.Date(date, format="%m/%d/%y")) %>%
  filter(location=="stream")

WQ_dat <- read.csv("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/Stream_Lake_YSI_WaterQuality.csv")%>%
  mutate(date = as.Date(date, format="%m/%d/%y")) %>%
  dplyr::select("site", "date","pH") %>%
  dplyr:: group_by(site, date) %>%
  dplyr::summarise(pH=mean(pH, na.rm=T))

covariat_datq <- covariat_dat%>%
  left_join(bg_nuts[ , c("site", "substrate", "date", "NO3_mgL_dl", "NH4_mgL_dl", "PO4_ugL_dl", "DOC_mgL_dl")], by=c("date", "site"))

covariat_datq <- covariat_datq%>%
  left_join(WQ_dat,  by=c("date", "site"))


## bring in SPC data 
### missing 2024 BW SPC
SPC_BWL <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWL_SPCv2.rds") %>%
  filter(datetime>as.POSIXct("2021-04-29 24:00:00"))
SPC_BWL$site <- "BWL"
SPC_BWU <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWU_SPCv2.rds")
SPC_BWU$site <- "BWU"
SPC_GBL <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_GBL_SPCv2.rds")
SPC_GBL$site <- "GBL"
SPC_GBU <- readRDS("/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_GBU_SPCv2.rds")
SPC_GBU$site <- "GBU"

SPC_dat <- rbind(SPC_BWL, SPC_BWU, SPC_GBL, SPC_GBU)
SPC_dat <- SPC_dat %>%
  mutate(date =as.Date(datetime)) 

SPC_dat_day <- SPC_dat %>%
  dplyr::group_by(site, date) %>%
  dplyr::summarise(
    wt= mean(wtr, na.rm=T),
    wt_sd= sd(wtr, na.rm=T),
    SPC_m=mean(SPC, na.rm=T),
    SPC_sd=sd(SPC, na.rm=T))

hist(SPC_dat_day$wt_sd)
hist(SPC_dat_day$SPC_sd)

SPC_datD <- SPC_dat_day %>%
  filter(wt_sd<4.1 & SPC_sd <35)
  
hist(SPC_datD$wt_sd)
hist(SPC_datD$SPC_sd)
  
covariat_datq <- covariat_datq%>%
  left_join(SPC_datD, by=c("date", "site"))

### 
## Bring in WQ data : 


## BWL  
BWL_dat <- readRDS("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/input_data/25_BWL_DO_flag_sat_light_v2.rds") %>%
  filter(maintenance_flag<1 & fouling_flag<1 & solar.time > as.POSIXct("2021-04-22 00:00:00") & solar.time < as.POSIXct("2024-08-01 00:00:00") & 
           do.obs>3.85 & wtr>0.01) %>%
  mutate(site="BWL")
summary(BWL_dat) 

## BWU 
BWU_dat <- readRDS("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/input_data/25_BWU_DO_flag_sat_light_v2.rds") %>%
  filter(maintenance_flag<1 & fouling_flag<1 & solar.time > as.POSIXct("2021-06-29 00:00:00") & solar.time < as.POSIXct("2024-08-01 00:00:00") & 
           do.obs>3.85 & wtr>0.01)  %>%
  mutate(site="BWU")
summary(BWU_dat)

## GBL 
GBL_dat <- readRDS("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/input_data/25_GBL_DO_flag_sat_light_v2.rds") %>%
  filter(maintenance_flag<1 & fouling_flag<1 & solar.time > as.POSIXct("2021-03-12 00:00:00") & solar.time < as.POSIXct("2024-08-01 00:00:00") & 
           do.obs>3.85 & wtr>0.01)  %>%
  mutate(site="GBL")
summary(GBL_dat)

## GBU 
GBU_dat <- readRDS("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/input_data/25_GBU_DO_flag_sat_light_v2.rds") %>%
  filter(fouling_flag<1 & solar.time > as.POSIXct("2021-03-12 00:00:00") & solar.time < as.POSIXct("2024-08-01 00:00:00") & 
           do.obs>3.85 & wtr>0.01)  %>%
  mutate(site="GBU")
summary(GBU_dat)


miniDOT_dat <- rbind(BWL_dat, BWU_dat, GBL_dat, GBU_dat)

miniDOT_dat_day <- miniDOT_dat%>%
  mutate(date =as.Date(solar.time))%>%
  dplyr::group_by(site, date) %>%
  dplyr::summarise(
    do.obs_m=mean(do.obs, na.rm=T),
    do.obs_sd=sd(do.obs, na.rm=T),
    wtr_m=mean(wtr, na.rm=T),
    wtr_sd=sd(wtr, na.rm=T),
    Q_m= mean(dischargeCMS, na.rm=T),
    Q_sd= sd(dischargeCMS, na.rm=T),
    v_m=mean(v, na.rm=T),
    v_sd=sd(v, na.rm=T),
    w_m=mean(w, na.rm=T),
    w_sd=sd(w, na.rm=T),
    depth_m=mean(depth, na.rm=T),
    depth_sd=sd(depth, na.rm=T))

### N-uptake metrics

n_up <- read.csv("/Users/kellyloria/Documents/UNR/Ncycle/BTC_N_Table_2025_figure.csv")%>%
  mutate(date = as.Date(date, format="%m/%d/%y")) 

covariat_datq <- covariat_datq%>%
  left_join(n_up, by=c("date", "site"))

covariat_datq <- covariat_datq%>%
  left_join(miniDOT_dat_day, by=c("date", "site"))

##================================
## Load metabolism model data
##================================
BWL_met_binned <- read_csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/final_output_2025/BWL_TVL_flagged_flow_excluded_2025-01-20_daily.csv") %>%
  filter(GPP_daily_Rhat<1.05)%>%
  filter(GPP_97.5pct>0)%>%
  filter(ER_daily_Rhat<1.05) %>%
  filter(ER_2.5pct<0)%>%
  filter(K600_daily_Rhat<1.05) %>%
  mutate(
    GPP_mean = ifelse(GPP_mean < 0, 0.001, GPP_mean)) %>%
  filter(ER_mean < 0.01) %>%
  mutate(NEP_daily_mean = GPP_daily_mean + ER_daily_mean) %>%
  dplyr::select(date, GPP_mean, ER_mean, NEP_daily_mean, K600_daily_mean,
                GPP_2.5pct, GPP_97.5pct, 
                ER_2.5pct, ER_97.5pct,
                K600_daily_2.5pct, K600_daily_predlog_97.5pct,
                GPP_daily_Rhat, ER_daily_Rhat, K600_daily_Rhat) %>% 
  mutate(site="BWL",
         catch="BW")

summary(BWL_met_binned)

GBL_met_binned <- read_csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/final_output_2025/GBL_Full_2025-02-02_daily.csv") %>%
  filter(GPP_daily_Rhat<1.05)%>%
  filter(GPP_97.5pct>0)%>%
  filter(ER_daily_Rhat<1.05) %>%
  filter(ER_2.5pct<0)%>%
  filter(K600_daily_Rhat<1.05) %>%
  mutate(
    GPP_mean = ifelse(GPP_mean < 0, 0.001, GPP_mean)) %>%
  filter(ER_mean < 0.01) %>%
  mutate(NEP_daily_mean = GPP_daily_mean + ER_daily_mean) %>%
  dplyr::select(date, GPP_mean, ER_mean, NEP_daily_mean, K600_daily_mean,
                GPP_2.5pct, GPP_97.5pct, 
                ER_2.5pct, ER_97.5pct,
                K600_daily_2.5pct, K600_daily_predlog_97.5pct,
                GPP_daily_Rhat, ER_daily_Rhat, K600_daily_Rhat) %>%
  mutate(site="GBL",
         catch="GB")


BWU_met_binned <- read_csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/final_output_2025/BWU_Full_2025-02-04_daily.csv") %>%
  filter(GPP_daily_Rhat<1.05)%>%
  filter(GPP_97.5pct>0)%>%
  filter(ER_daily_Rhat<1.05) %>%
  filter(ER_2.5pct<0)%>%
  filter(K600_daily_Rhat<1.05) %>%
  mutate(
    GPP_mean = ifelse(GPP_mean < 0, 0.001, GPP_mean)) %>%
  filter(ER_mean < 0.01) %>%
  mutate(NEP_daily_mean = GPP_daily_mean + ER_daily_mean) %>%
  dplyr::select(date, GPP_mean, ER_mean, NEP_daily_mean, K600_daily_mean,
                GPP_2.5pct, GPP_97.5pct, 
                ER_2.5pct, ER_97.5pct,
                K600_daily_2.5pct, K600_daily_predlog_97.5pct,
                GPP_daily_Rhat, ER_daily_Rhat, K600_daily_Rhat) %>%
  mutate(site="BWU",
         catch="BW")

GBU_met_binned <- read_csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/final_output_2025/GBU_Full_2025-02-03_daily.csv") %>%
  filter(GPP_daily_Rhat<1.05)%>%
  filter(GPP_97.5pct>0)%>%
  filter(ER_daily_Rhat<1.05) %>%
  filter(ER_2.5pct<0)%>%
  filter(K600_daily_Rhat<1.05) %>%
  mutate(
    GPP_mean = ifelse(GPP_mean < 0, 0.001, GPP_mean)) %>%
  filter(ER_mean < 0.01) %>%
  mutate(NEP_daily_mean = GPP_daily_mean + ER_daily_mean) %>%
  dplyr::select(date,GPP_mean, ER_mean, NEP_daily_mean, K600_daily_mean,
                GPP_2.5pct, GPP_97.5pct, 
                ER_2.5pct, ER_97.5pct,
                K600_daily_2.5pct, K600_daily_predlog_97.5pct,
                GPP_daily_Rhat, ER_daily_Rhat, K600_daily_Rhat) %>%
  mutate(site="GBU",
         catch="GB")

metab_dat <- rbind(BWL_met_binned, BWU_met_binned, GBL_met_binned, GBU_met_binned)



covariat_datq <- covariat_datq%>%
  left_join(metab_dat, by=c("date", "site", "catch"))


hist(covariat_datq$AFDM_mgcm2)
hist(covariat_datq$Chla_ugL_Q)
hist(covariat_datq$AFDM_mgg)

covariat_datq$year <- year(covariat_datq$date)
covariat_datq$yday <- yday(covariat_datq$date)
covariat_datq$month <- month(covariat_datq$date)

## 
names(covariat_datq)
str(covariat_datq)
cor_data <- covariat_datq[, c(5, 10:14, 17, 33, 37)]
str(cor_data)

chart.Correlation(cor_data, histogram=TRUE, pch=19)

covariat_datq$TIN_mgL <- c(covariat_datq$NO3_mgL_dl + covariat_datq$NH4_mgL_dl)

####
#### GLMS
####

covariat_datq_BW <- covariat_datq%>%
  filter(site=="BWL" | site=="BWU")

biomass_mod_bw <- lmer(log(AFDM_mgg+1)~ 
                         scale(NO3_mgL_dl* 1000) + 
                         scale(NH4_mgL_dl* 1000) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(wtr_m) + 
                         scale(Q_m) + 
                         (1|water_year), data=covariat_datq_BW)

summary(biomass_mod_bw)
vif(biomass_mod_bw)
hist(residuals(biomass_mod_bw))
r.squaredGLMM(biomass_mod_bw)


library(broom.mixed)
# Extract fixed effects and confidence intervals
fixed_effects <- broom.mixed::tidy(biomass_mod_bw, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects1 <- fixed_effects[fixed_effects$term != "(Intercept)", ]
fixed_effects1 <- fixed_effects1[order(fixed_effects1$estimate), ]  # Order by estimate
fixed_effects1$Catchment <- "BW"


### gB
covariat_datq_GB <- covariat_datq%>%
  filter(site=="GBL" | site=="GBU")

biomass_mod_gb <- lmer(log(AFDM_mgg+1)~ 
                         scale(NO3_mgL_dl* 1000) + 
                         scale(NH4_mgL_dl* 1000) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(wtr_m) + 
                         scale(Q_m) + 
                         (1|water_year), data=covariat_datq_GB)

summary(biomass_mod_gb)
vif(biomass_mod_gb)
hist(residuals(biomass_mod_gb))
r.squaredGLMM(biomass_mod_gb)


fixed_effects_gb <- broom.mixed::tidy(biomass_mod_gb, effects = "fixed", conf.int = TRUE)

# Filter to exclude intercept and arrange predictors
fixed_effects_gb1 <- fixed_effects_gb[fixed_effects_gb$term != "(Intercept)", ]
fixed_effects_gb1 <- fixed_effects_gb1[order(fixed_effects_gb1$estimate), ]  # Order by estimate
fixed_effects_gb1$Catchment <- "GB"

fxd_eff <- rbind(fixed_effects_gb1, fixed_effects1)

fxd_eff1 <- fxd_eff %>%
  mutate(Significance = case_when(
    p.value <= 0.055 ~ "sig",
    p.value > 0.055 & p.value <= 0.09  ~ "marg.",
    p.value > 0.09  ~ "not sig."))


term_labels <- c(
  "scale(NO3_mgL_dl * 1000)" = expression(NO[3]~(µg~L^-1)),
  "scale(PO4_ugL_dl)" = expression(PO[4]~(µg~L^-1)),
  "scale(DOC_mgL_dl)" = expression(DOC~(mg~L^-1)),
  "scale(NH4_mgL_dl * 1000)" = expression(NH[4]~(µ~L^-1)),
  "scale(wtr_m)" = expression(Temp~(degree~C)),
  "scale(Q_m)" = expression(Q~(m^3~s^-1))
)

# Add an offset column based on Catchment
fxd_eff1 <- fxd_eff1 %>%
  mutate(y_position = as.numeric(reorder(term, estimate)) + ifelse(Catchment == "BW", 0.2, 0))

Biomass_coeff_plot <-ggplot(fxd_eff1, aes(x = estimate, y = reorder(term, estimate), color = Catchment, shape = Significance)) +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 0.85, position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_shape_manual(values = c(1, 4, 19)) +
  scale_color_manual(values = catch_colors) +
  geom_vline(xintercept = 0) +
  labs(
    x = "Scaled effect size",
    y = "Predictors",
    title = "Epilithic biomass"
  ) +
  scale_y_discrete(labels = term_labels) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )


Biomass_coeff_plot


##########
## CHL-a
#########
chla_mod_bw <- lmer(log(Chla_ugL_Q+1)~ 
                         scale(NO3_mgL_dl* 1000) + 
                         scale(NH4_mgL_dl* 1000) + 
                         scale(PO4_ugL_dl)+ 
                         scale(DOC_mgL_dl)+
                         scale(wtr_m) + 
                         scale(Q_m) + 
                         (1|water_year), data=covariat_datq_BW)

summary(chla_mod_bw)
vif(chla_mod_bw)
hist(residuals(chla_mod_bw))
r.squaredGLMM(chla_mod_bw)



# Extract fixed effects and confidence intervals
fixed_effects_chl <- broom.mixed::tidy(chla_mod_bw, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects_chl1 <- fixed_effects_chl[fixed_effects_chl$term != "(Intercept)", ]
fixed_effects_chl1 <- fixed_effects_chl1[order(fixed_effects_chl1$estimate), ]  # Order by estimate
fixed_effects_chl1$Catchment <- "BW"

## GB
chla_mod_gb <- lmer(log(Chla_ugL_Q+1)~ 
                      scale(NO3_mgL_dl* 1000) + 
                      scale(NH4_mgL_dl* 1000) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(wtr_m) + 
                      scale(Q_m) + 
                      (1|water_year), data=covariat_datq_GB)

summary(chla_mod_gb)
vif(chla_mod_gb)
hist(residuals(chla_mod_gb))
r.squaredGLMM(chla_mod_gb)

fixed_effects_chlg <- broom.mixed::tidy(chla_mod_gb, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects_chlg1 <- fixed_effects_chlg[fixed_effects_chlg$term != "(Intercept)", ]
fixed_effects_chlg1 <- fixed_effects_chlg1[order(fixed_effects_chlg1$estimate), ]  # Order by estimate
fixed_effects_chlg1$Catchment <- "GB"


fxd_eff_chla <- rbind(fixed_effects_chlg1, fixed_effects_chl1)

fxd_eff_chla1 <- fxd_eff_chla %>%
  mutate(Significance = case_when(
    p.value <= 0.055 ~ "sig.",
    p.value > 0.055 & p.value <= 0.09  ~ "marg.",
    p.value > 0.09  ~ "not sig."))


term_labels_chla <- c(
  "scale(NO3_mgL_dl * 1000)" = expression(NO[3]~(µg~L^-1)),
  "scale(PO4_ugL_dl)" = expression(PO[4]~(µg~L^-1)),
  "scale(DOC_mgL_dl)" = expression(DOC~(mg~L^-1)),
  "scale(NH4_mgL_dl * 1000)" = expression(NH[4]~(µ~L^-1)),
  "scale(wtr_m)" = expression(Temp~(degree~C)),
  "scale(Q_m)" = expression(Q~(m^3~s^-1))
)

# Add an offset column based on Catchment
fxd_eff_chla1 <- fxd_eff_chla1 %>%
  mutate(y_position = as.numeric(reorder(term, estimate)) + ifelse(Catchment == "BW", 0.2, 0))


chla_coeff_plot <-ggplot(fxd_eff_chla1, aes(x = estimate, y = reorder(term, estimate), color = Catchment, shape = Significance)) +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 0.85, position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_shape_manual(values = c(1,4, 19)) +
  scale_color_manual(values = catch_colors) +
  geom_vline(xintercept = 0) +
  labs(
    x = "Scaled effect size",
    y = "Predictors",
    title = "Epilithic chl-a"
  ) +
  scale_y_discrete(labels = term_labels_chla) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )
chla_coeff_plot

#########

##########
## GPP 
#########
hist(covariat_datq_BW$GPP_mean)
hist(log(covariat_datq_BW$GPP_mean+1))

gpp_mod_bw <- lmer(log(GPP_mean+1)~ 
                      scale(NO3_mgL_dl* 1000) + 
                      scale(NH4_mgL_dl* 1000) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(wtr_m) + 
                      scale(Q_m) + 
                      (1|water_year), data=covariat_datq_BW)

summary(gpp_mod_bw)
vif(gpp_mod_bw)
hist(residuals(gpp_mod_bw))
r.squaredGLMM(gpp_mod_bw)


# library(lme4)
# gpp_mod_bw <- glmer(GPP_mean ~ 
#                       scale(NO3_mgL_dl * 1000) + 
#                       scale(NH4_mgL_dl * 1000) + 
#                       scale(PO4_ugL_dl) + 
#                       scale(DOC_mgL_dl) + 
#                       scale(wtr_m) + 
#                       scale(Q_m) + 
#                       (1 | water_year), 
#                     family = Gamma(link = "log"), 
#                     data = covariat_datq_BW)
# summary(gpp_mod_bw)



library(glmmTMB)
gpp_mod_bw <- glmmTMB(GPP_mean ~ 
                        scale(NO3_mgL_dl * 1000) + 
                        scale(NH4_mgL_dl * 1000) + 
                        scale(PO4_ugL_dl) + 
                        scale(DOC_mgL_dl) + 
                        scale(wtr_m) + 
                        scale(Q_m) + 
                        (1 | water_year), 
                      family = gaussian(link = "log"), 
                      data = covariat_datq_BW)


# Extract fixed effects and confidence intervals
fixed_effects_nep <- broom.mixed::tidy(gpp_mod_bw, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects_nep <- fixed_effects_nep[fixed_effects_nep$term != "(Intercept)", ]
fixed_effects_nep <- fixed_effects_nep[order(fixed_effects_nep$estimate), ]  # Order by estimate
fixed_effects_nep$Catchment <- "BW"

## GB
hist(covariat_datq_GB$NEP_daily_mean)
hist(log(covariat_datq_GB$NEP_daily_mean))


gpp_mod_gb <- lmer(log(GPP_mean+1) ~ 
                      scale(NO3_mgL_dl* 1000) + 
                      scale(NH4_mgL_dl* 1000) + 
                      scale(PO4_ugL_dl)+ 
                      scale(DOC_mgL_dl)+
                      scale(wtr_m) + 
                      scale(Q_m) + 
                      (1|water_year), data=covariat_datq_GB)

summary(gpp_mod_gb)
vif(gpp_mod_gb)
hist(residuals(gpp_mod_gb))
r.squaredGLMM(gpp_mod_gb)





library(glmmTMB)
gpp_mod_gb <- glmmTMB(GPP_mean ~ 
                        scale(NO3_mgL_dl * 1000) + 
                        scale(NH4_mgL_dl * 1000) + 
                        scale(PO4_ugL_dl) + 
                        scale(DOC_mgL_dl) + 
                        scale(wtr_m) + 
                        scale(Q_m) + 
                        (1 | water_year), 
                      family = gaussian(link = "log"), 
                      data = covariat_datq_GB)


fixed_effects_nep1 <- broom.mixed::tidy(gpp_mod_gb, effects = "fixed", conf.int = TRUE)
# Filter to exclude intercept and arrange predictors
fixed_effects_nep1 <- fixed_effects_nep1[fixed_effects_nep1$term != "(Intercept)", ]
fixed_effects_nep1 <- fixed_effects_nep1[order(fixed_effects_nep1$estimate), ]  # Order by estimate
fixed_effects_nep1$Catchment <- "GB"


fxd_eff_nep <- rbind(fixed_effects_nep, fixed_effects_nep1)

fxd_eff_nep1 <- fxd_eff_nep %>%
  mutate(Significance = case_when(
    p.value <= 0.055 ~ "sig.",
    p.value > 0.055 & p.value <= 0.09  ~ "marg.",
    p.value > 0.09  ~ "not sig."))


term_labels_nep <- c(
  "scale(NO3_mgL_dl * 1000)" = expression(NO[3]~(µg~L^-1)),
  "scale(PO4_ugL_dl)" = expression(PO[4]~(µg~L^-1)),
  "scale(DOC_mgL_dl)" = expression(DOC~(mg~L^-1)),
  "scale(NH4_mgL_dl * 1000)" = expression(NH[4]~(µ~L^-1)),
  "scale(wtr_m)" = expression(Temp~(degree~C)),
  "scale(Q_m)" = expression(Q~(m^3~s^-1))
)

# Add an offset column based on Catchment
fxd_eff_nep1 <- fxd_eff_nep1 %>%
  mutate(y_position = as.numeric(reorder(term, estimate)) + ifelse(Catchment == "BW", 0.2, 0))


nep_coeff_plot <-ggplot(fxd_eff_nep1, aes(x = estimate, y = reorder(term, estimate), color = Catchment, shape = Significance)) +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, alpha = 0.85, position = position_dodge(width = 0.3)) +
  theme_bw() +
  scale_shape_manual(values = c(1,4, 19)) +
  scale_color_manual(values = catch_colors) +
  geom_vline(xintercept = 0) +
  labs(
    x = "Scaled effect size",
    y = "Predictors",
    title = "GPP"
  ) +
  scale_y_discrete(labels = term_labels_nep) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )








### big grid 
coeff_grid <- ggarrange(
  Biomass_coeff_plot,
  chla_coeff_plot,
  nep_coeff_plot,
  ncol = 1, nrow = 3,
  common.legend = TRUE, 
  legend = "bottom")


coeff_grid <- ggarrange(
  Biomass_coeff_plot + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)),
  chla_coeff_plot + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)),
  nep_coeff_plot + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 3)),
  labels=c("a", "b", "c"),
  
  ncol = 1, nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/Fig4_coeff_ecol_CH1_grid.png", plot = coeff_grid, width = 3.5, height = 9, units = "in")



## =======================
## Mean WY values 

metab_dat_wy <- water_year(metab_dat)

metab_dat_wy1 <- metab_dat_wy %>%
  filter(water_year < 2025) %>%
  mutate(mon = month(date))%>%
  filter(mon %in% c(6:10))%>%
 # filter( == "BWL" | site =="GBL") %>%
  group_by(site, water_year) %>%
  dplyr:: summarise(
    GPP = mean(GPP_mean, na.rm=T), 
    GPP_sd = sd(GPP_mean, na.rm=T), 
    ER = mean(ER_mean, na.rm=T), 
    ER_sd = sd(ER_mean, na.rm=T),
    NEP = mean(NEP_daily_mean, na.rm=T), 
    NEP_sd = sd(NEP_daily_mean, na.rm=T))



gpp_plot <- metab_dat_wy1 %>%
  ggplot(aes(y = GPP, x= water_year, col = site))+
  geom_point(position = position_dodge(width = 0.1), alpha = .7, size = 1.3) + 
  geom_line(alpha = .5) +
  geom_pointrange(aes(ymin = (GPP-GPP_sd), ymax = (GPP+GPP_sd)),  position = position_dodge(width = 0.1)) +
  theme_bw() +
  scale_color_manual(values = site_colors)

er_plot <-metab_dat_wy1 %>%
  ggplot(aes(y = ER, x= water_year, col = site))+
  geom_point(position = position_dodge(width = 0.1), alpha = .7, size = 1.3) +
  geom_line(alpha = .5) +
  geom_pointrange(aes(ymin = (ER-ER_sd), ymax = (ER+ER_sd)),position = position_dodge(width = 0.1)) +
  theme_bw() +
  scale_color_manual(values = site_colors)

NEP_plot <- metab_dat_wy1 %>%
  ggplot(aes(y = NEP, x= water_year, col = site))+
  geom_point(position = position_dodge(width = 0.1), alpha = .7, size = 1.3) +
  geom_line(alpha = .5) +
  geom_pointrange(aes(ymin = (NEP-NEP_sd), ymax = (NEP+NEP_sd)), position = position_dodge(width = 0.1)) +
  theme_bw() +
  scale_color_manual(values = site_colors)

metab_grid <- ggarrange(
  gpp_plot,
  er_plot,
  NEP_plot,
  ncol = 1, nrow = 3,
  common.legend = TRUE, 
  legend = "bottom")


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/WY_metab_CH1_grid.png", plot = metab_grid, width = 4.5, height = 6, units = "in")


## % changes 
BWL_dry_gpp <- (2.69 + 1.90) /2
BWL_wet_gpp <- 0.442
(BWL_dry_gpp- BWL_wet_gpp)/BWL_dry_gpp

BWL_dry_er <- (12.0 + 11.1) /2
BWL_wet_er <- 8.89

(BWL_dry_er- BWL_wet_er)/BWL_dry_er

BWL_dry_nep <- (9.36 + 9.30) /2
BWL_wet_nep <- 8.62

(BWL_dry_nep- BWL_wet_nep)/BWL_dry_nep


GBL_dry_gpp <- (0.0117 + 0.0121) /2
GBL_wet_gpp <- 0.0311
(GBL_wet_gpp- GBL_dry_gpp)/GBL_wet_gpp


GBL_dry_gpp <- (2.50 + 5.59) /2
GBL_wet_gpp <- 10.2
(GBL_wet_gpp- GBL_dry_gpp)/GBL_wet_gpp


GBL_dry_gpp <- (3.06 + 6.00) /2
GBL_wet_gpp <- 10.5
(GBL_wet_gpp- GBL_dry_gpp)/GBL_wet_gpp




nut_dat_wy1 <- covariat_datq %>%
  filter(water_year < 2025) %>%
  mutate(mon = month(date))%>%
  filter(mon %in% c(6:10))%>%
  group_by(site, water_year) %>%
  dplyr:: summarise(
    NO3 = mean(NO3_mgL_dl, na.rm=T), 
    NO3_sd= sd(NO3_mgL_dl, na.rm=T),
    NH4= mean(NH4_mgL_dl, na.rm=T),
    NH4_sd= sd(NH4_mgL_dl, na.rm=T))




nut_plot <- nut_dat_wy1 %>%
  ggplot(aes(y = NO3, x= water_year, col = site))+
  geom_point(position = position_dodge(width = 0.1), alpha = .7, size = 1.3) + 
  geom_line(alpha = .5) +
  #ylim(0,0.05) +
 geom_pointrange(aes(ymin = (NO3-NO3_sd), ymax = (NO3+NO3_sd)),  position = position_dodge(width = 0.1), alpha=.5) +
  theme_bw() +
  scale_color_manual(values = site_colors)

nut1_plot <- nut_dat_wy1 %>%
  ggplot(aes(y = NH4, x= water_year, col = site))+
  geom_point(position = position_dodge(width = 0.1), alpha = .7, size = 1.3) + 
  geom_line(alpha = .5) +
  geom_pointrange(aes(ymin = (NH4-NH4_sd), ymax = (NH4+NH4_sd)),  position = position_dodge(width = 0.1)) +
  theme_bw() +
  scale_color_manual(values = site_colors)



nut_grid <- ggarrange(
  nut_plot,
  nut1_plot,
  ncol = 1, nrow = 2,
  common.legend = TRUE, 
  legend = "bottom")


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/WY_nh4_CH1_grid.png", plot = nut_grid, width = 4.5, height = 4, units = "in")




GBL_dry_gpp <- (0.0593 + 0.0368 + 0.0183 + 0.0446) /4
GBL_wet_gpp <- (0.0486 + 0.0468) /2
(GBL_wet_gpp- GBL_dry_gpp)/GBL_wet_gpp

# dry to wet
# NO3 16% higher 

GBL_dry_gpp <- (0.0161 + 0.0371 + 0.00934 + 0.0366) /4
GBL_wet_gpp <- (0.0215 + 0.00703) /2
(GBL_dry_gpp-GBL_wet_gpp)/GBL_dry_gpp

# dry to wet 
# NH4 42% lower 


GBL_dry_gpp <- (0.0346 + 0.0302 + 0.0094 + 0.0229) /4
GBL_wet_gpp <- (0.0381 + 0.0313) /2
(GBL_wet_gpp- GBL_dry_gpp)/GBL_wet_gpp
# 30% higher NO3 at BW

(GBL_dry_gpp-GBL_wet_gpp)/GBL_dry_gpp

GBL_dry_gpp <- (0.0111 + 0.0191 + 0.0106 +0.0198) /4
GBL_wet_gpp <-  (0.0133 + 0.00679) /2
(GBL_dry_gpp-GBL_wet_gpp)/GBL_dry_gpp

# 33% lower NH3 at BW


#### ===================================
## Not quick run through of lms for each   
##    reach for uptake experiments 
### ===================================
library(dplyr)
library(broom) # For tidy summaries of model outputs

dat <- covariat_datq%>%
  filter(site=="GBU" & method == "NH3")

# Filter data by site and method
site_methods <- list(
  #list(site = "GBL", method = "NO3")
  # list(site = "BWL", method = "NH3")
  # list(site = "GBL", method = "NO3"),
  # list(site = "GBL", method = "NH3")
  # list(site = "BWU", method = "NO3"),
  # list(site = "BWU", method = "NH3"),
  # list(site = "GBU", method = "NO3")
   list(site = "GBU", method = "NH3")
)

# Variables to model
covariates <- c("scale(AFDM_mgcm2)", 
                "scale(Chla_ugL_Q)", 
                "scale(Pheo_ugL_Q)", 
                "scale(AFDM_mgg)", 
                "scale(NO3_mgL_dl)", 
                "scale(NH4_mgL_dl)",
                "scale(GPP_mean)",
                "scale(ER_mean)")

# Placeholder for results
results <- data.frame()

# Loop through sites, methods, and covariates
for (site_method in site_methods) {
  site <- site_method$site
  method <- site_method$method
  
  # Filter data
  data_filtered <- dat 
  # Check if there are data for this site and method
  if (nrow(data_filtered) == 0) next
  
  for (covariate in covariates) {
    formula <- as.formula(paste("Uadd_ug_L_min ~", covariate))
    model <- lm(formula, data = data_filtered)
    
    # Extract summary stats
    model_summary <- tidy(model) %>%
      filter(term != "(Intercept)") %>% # Exclude intercept
      mutate(
        site = site,
        method = method,
        df = glance(model)$df.residual # Degrees of freedom
      )
    
    # Append results
    results <- bind_rows(results, model_summary)
  }
}

# Ensure dplyr functions are explicitly called to avoid conflicts
final_table <- results %>%
  mutate(estimate= round(estimate,2),
         std.error =round(std.error,2),
         p.value = round(p.value, 3)
  ) %>%
  dplyr::select(site, method, term, estimate, std.error, df, p.value) %>%
  dplyr::rename(
    Covariate = term,
    Estimate = estimate,
    `Std. Error` = std.error,
    `Degrees of Freedom` = df,
    `p-value` = p.value
  )


# View or export table
print(final_table)
# Write table to CSV if needed
# write.csv(final_table, "/Users/kellyloria/Documents/UNR/Ncycle/Ncycle_linear_model_drivers/BWL_NO3_Uadd_model.csv", row.names = FALSE)
getwd()


NO3_covar <- covariat_datq%>%
  filter(method=="NH3" & site=="BWL")

BW_no3_mod <- lm(Vf_m_min~ scale(_mean), data = NO3_covar)
summary(BW_no3_mod)


NO3_covar <- covariat_datq%>%
  filter(method=="NH3" & site=="GBU")

BW_no3_mod <- lm(Vf_m_min~ scale(ER_mean), data = NO3_covar)
summary(BW_no3_mod)

####
## Quick summary of results by site:

n_up_sum <-n_up%>%
  filter(success=="yes") %>%
  dplyr::group_by(method, site) %>%
  dplyr::summarise(
    Uadd_ug_L = mean(Uadd_ug_L_min, na.rm=T),
    Uadd_ug_L_MIN = min(Uadd_ug_L_min, na.rm=T),
    Uadd_ug_L_MAX = max(Uadd_ug_L_min, na.rm=T),
    sw = mean(sw_m, na.rm=T),
    sw_MIN = min(sw_m, na.rm=T),
    sw_MAX = max(sw_m, na.rm=T),
    Vf_m = mean(Vf_m_min, na.rm=T),
    Vf_m_MIN = min(Vf_m_min, na.rm=T),
    Vf_m_MAX = max(Vf_m_min, na.rm=T))

ggplot(n_up_sum, aes(x=site, y = Uadd_ug_L, color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")
  

nup_sum_plot <-ggplot(n_up_sum, aes(x = site, y = Uadd_ug_L, color = site, shape = method)) +
  geom_point(size = 3, alpha = 0.95, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Uadd_ug_L_MIN, ymax = Uadd_ug_L_MAX),  width = 0.05, height = 0, alpha = 0.75, position = position_dodge(width = 0.3)) +
  theme_bw() + ylab(expression(U[t]~(µg~L^-1~min^-1)))+
  xlab(NULL) +
  #scale_shape_manual(values = c(4, 19)) +
  scale_color_manual(values = site_colors) 

nvf_sum_plot <-ggplot(n_up_sum, aes(x = site, y = Vf_m, color = site, shape = method)) +
  geom_point(size = 3, alpha = 0.95, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Vf_m_MIN, ymax = Vf_m_MAX),  width = 0.05, height = 0, alpha = 0.75, position = position_dodge(width = 0.3)) +
  theme_bw() + ylab(expression(V[f]~(µg~m^-1~min^-1)))+
  xlab(NULL) +
  #scale_shape_manual(values = c(4, 19)) +
  scale_color_manual(values = site_colors) 

nsw_sum_plot <-ggplot(n_up_sum, aes(x = site, y = sw, color = site, shape = method)) +
  geom_point(size = 3, alpha = 0.95, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = sw_MIN, ymax = sw_MAX),  width = 0.05, height = 0, alpha = 0.75, position = position_dodge(width = 0.3)) +
  theme_bw() + ylab(expression(S[w]~(m)))+
  xlab(NULL) +
  #scale_shape_manual(values = c(4, 19)) +
  scale_color_manual(values = site_colors) 



### big grid 
Nup_grid <- ggarrange(
  nup_sum_plot,
  nvf_sum_plot,
  nsw_sum_plot,
  ncol = 1, nrow = 3,
  common.legend = TRUE, 
  legend = "right")

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages/supp\ figures/Nuptake_CH1_grid.png", plot = Nup_grid, width = 3.25, height = 7, units = "in")



# Combine plots with a shared legend
Nup_grid <- ggarrange(
  nup_sum_plot,
  nvf_sum_plot,
  nsw_sum_plot,
  ncol = 1, nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)


Nup_grid <- Nup_grid + 
  theme(legend.box = "vertical",        # Arrange legends vertically
        legend.text = element_text(size = 10),  # Adjust legend text size if needed
        legend.title = element_text(size = 12)) +
  guides(
    color = guide_legend(nrow = 1),    # Put colors on one row
    shape = guide_legend(nrow = 1)    # Put shapes on another row
  )

##################
##################
##################
hist(covariat_datq$AFDM_mgg)
hist(log(covariat_datq$AFDM_mgg)+1)
ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtr_m, y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(AFDM_mgg+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) + 
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm", se=F)

## biomass
hist(covariat_datq$AFDM_mgcm2)
hist(log(covariat_datq$AFDM_mgcm2+1))

ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtr_m, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)+
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm", se=F)

ggplot(covariat_datq, aes(x=NH4_mgL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm", se=F)

ggplot(covariat_datq, aes(x=PO4_ugL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=DOC_mgL_dl, y = log(AFDM_mgcm2+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

####
# chla
## biomass
hist(covariat_datq$Chla_ugL_Q)
hist(log(covariat_datq$Chla_ugL_Q+1))

ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtemp, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)+
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=NH4_mgL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=PO4_ugL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=DOC_mgL_dl, y = log(Chla_ugL_Q+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

###

####
# chla
## biomass

covariat_datq$biom_chla <- (covariat_datq$Chla_ugL_Q)/(covariat_datq$AFDM_mgcm2)
hist(covariat_datq$biom_chla)
hist(log(covariat_datq$biom_chla+1))

ggplot(covariat_datq, aes(x=log(Q_m+1), y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=wtemp, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=SPC_m, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)

ggplot(covariat_datq, aes(x=NO3_mgL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors)+
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=NH4_mgL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=PO4_ugL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

ggplot(covariat_datq, aes(x=DOC_mgL_dl, y = log(biom_chla+1), color =site)) + 
  geom_point(size=2, alpha=0.5) + theme_bw() + scale_color_manual(values = site_colors) +
  facet_grid(.~substrate, scales = "free") + geom_smooth(method="lm")

### 
# fill in NA wtemp with temps from SPC 
covariat_datq <- covariat_datq  %>%
  mutate(wtemp = ifelse(is.na(wtemp), wt, wtemp))
  
str(covariat_datq)



### big grid 
up_grid3 <- ggarrange(
  Q_plot,
  temp_plot,
  SPC_plot,
  pH_plot, # not done
  NO3_plot,
  NH4_plot,
  PO4_plot,
  DOC_plot,
  OM_plot,
  biom_plot,
  chla_plot,
  chla_plot,
  ncol = 4, nrow = 3,
  common.legend = TRUE, 
  legend = "bottom")

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/Draft_figure3_CH1_grid.png", plot = up_grid3, width = 14, height = 10, units = "in")



### ================== 
## estimate n-demand 
# Compute median values for each site

covariat_nitrogen <- covariat_datq %>%
  filter(substrate=="sw") %>%
  group_by(site) %>%
  arrange(date) %>%  # Ensure chronological order
  mutate(v_m = zoo::na.locf(v_m, na.rm = FALSE),  # Fill downwards
         w_m = zoo::na.locf(w_m, na.rm = FALSE)) %>%
  ungroup()

# 
# covariat_nitrogen1 <- covariat_nitrogen %>%
#   group_by(site) %>%
#   arrange(date) %>%  # Ensure chronological order
#   mutate(v_m = zoo::na.locf(v_m, na.rm = FALSE, fromLast=T),  # Fill downwards
#          w_m = zoo::na.locf(w_m, na.rm = FALSE, fromLast=T)) %>%
#   ungroup()

# Constants
ra <- 0.5  # Autotrophic respiration coefficient (Hall & Tank 2003)
C_Nauto <- 16  # Autotrophic C:N ratio (Stelzer & Lamberti 2001)
C_Nhetero <- 20  # Heterotrophic C:N ratio (Hall & Tank 2003)
HGE <- 0.05  # Heterotrophic Growth Efficiency (Hall & Tank 2003)

# Calculate components of nitrogen demand
covariat_nitrogen2 <- covariat_nitrogen %>%
  filter(substrate=="sw")%>%
  mutate(
    # Autotrophic respiration (raGPP)
    raGPP = GPP_mean * ra,
    # Autotrophic assimilation of N
    Auto_N_assim = GPP_mean / C_Nauto,
    # Heterotrophic respiration (Rh)
    Rh = ER_mean - raGPP,
    # Heterotrophic assimilation of N
    Hetero_N_assim = (Rh * HGE) / C_Nhetero,
    # Total nitrogen demand
    Ndemand = Auto_N_assim + Hetero_N_assim,
    # calculate reach length in m
    reachL = c((v_m*w_m*86.4)/K600_daily_mean),
   # Calculate no3 supply 
   NO3_supply = c((86400*Q_m*NO3_mgL_dl)/(w_m*reachL)),
   # Calculate nh3 supply 
   NH4_supply = c((86400*Q_m*NH4_mgL_dl)/(w_m*reachL)))

hist(covariat_nitrogen2$reachL)
hist(covariat_nitrogen2$Auto_N_assim)
hist(covariat_nitrogen2$Rh)
hist(covariat_nitrogen2$Hetero_N_assim)
hist(covariat_nitrogen2$Ndemand)
hist(covariat_nitrogen2$NO3_supply)
hist(covariat_nitrogen2$NH4_supply)

## need to calculate reach length
# NO3_mgL_w = N
# NO3_supply = c(86400*Q_m*NO3_mgL_dl/w_m*L)

# Check results
head(covariat_datq1 %>% select(date, site, GPP_mean, ER_mean, raGPP, Auto_N_assim, Rh, Hetero_N_assim, Ndemand))
names(covariat_datq1)