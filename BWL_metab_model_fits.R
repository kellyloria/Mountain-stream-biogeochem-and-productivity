##===============================================================================
## Created  01/09/2025 by KAL
##===============================================================================
## Evaluate different stream metabolism model fits based on different data filtering, 
##    aggregation levels, or priors on GPP or ER. 

packages <- c("dplyr", "scales", "ggplot2", "lubridate", "tidyr")
lapply(packages, library, character.only = TRUE)


## create an empty vector of times 
# Define the start and end dates
start_date <- as.Date("2021-06-30")
end_date <- as.Date("2024-08-02")
# Create a sequence of datetime values in 15-minute intervals
date_seq <- seq(from = start_date, to = end_date, by = "1 day")
# Convert to a dataframe
date_df <- data.frame(date = date_seq)

BWL_prior_1208g <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_measured_prior_NLDAS_daily241208.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d")) %>%
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))

BWL_prior_1208g$info <- "Filtered out biofouled days + strong GPP and ER prior"
str(BWL_prior_1208g)

BWL_filtered_1211g <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_measured_prior_noGPP_prior_NLDAS_daily241211.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d")) %>%
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
BWL_filtered_1211g$info <- "Filtered + no prior on ER/GPP"


BWL_ufiltered_1213g <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_measured_prior_noGPP_prior_no_filter_NLDAS_daily241213.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d")) %>%
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
BWL_ufiltered_1213g$info <- "Ufiltered + no prior on ER/GPP"


BWL_hr_0107g <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_m_prior_filter_NLDAS_daily250107.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d")) %>%
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
BWL_hr_0107g$info <- "hourly + no prior on ER/GPP"



all_BWL <- rbind(BWL_prior_1208g, BWL_filtered_1211g, BWL_ufiltered_1213g, BWL_hr_0107g)

str(all_BWL)


library(ggplot2)

# Create the density plot with overlapping transparency
gpp_outcomes<-ggplot(all_BWL, aes(x = GPP_mean, fill = info, color = info)) +
  geom_density(alpha = 0.4, size = 0.7) +  # Use alpha for transparency
  labs(
    title = "",
    x= expression(GPP~(g~O[2]~m^-2~d^-1)), 
    y = "Density"
  ) +
  theme_bw() +
  geom_vline(xintercept = 0)+
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

# Create the density plot with overlapping transparency
ER_outcomes<-ggplot(all_BWL, aes(x = ER_mean, fill = info, color = info)) +
  geom_density(alpha = 0.4, size = 0.7) +  # Use alpha for transparency
  labs(
    title = "BWL",
    x= expression(ER~(g~O[2]~m^-2~d^-1)), 
    y = "Density"
  ) +
  theme_bw() +
  geom_vline(xintercept = 0)+
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )


library(ggpubr)

GBL_modfit_grid <- ggarrange(ER_outcomes,
                             gpp_outcomes,
                          ncol = 2, nrow = 1,
                          common.legend = TRUE, 
                          legend = "bottom")


## ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/BWL_metab_fits_grid_v2.png", plot = GBL_modfit_grid, width = 10.5, height = 5, units = "in")

############
## DO plots
############

BWL_prior_1208 <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_measured_prior_NLDAS_mod_and_obs.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d")) %>%
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
info <- "Filtered out biofouled days + strong GPP and ER prior"
str(BWL_prior_1208)
BWL_prior_1208 <- BWL_prior_1208 %>%
  filter(!is.na(DO.mod), !is.na(DO.obs))  
rmse_values <- BWL_prior_1208 %>%
  # group_by(y_var) %>%
  summarize(RMSE = sqrt(mean((DO.mod - DO.obs)^2, na.rm = TRUE)))

BWL_prior_1208 <- date_df%>%
  left_join(BWL_prior_1208, by=c("date"))

BWL_prior_1208 <- BWL_prior_1208 %>%
  left_join(BWL_prior_1208g, by=c("date"))

# Create a dataframe with shading regions based on GPP_mean < 0
shading_data <- BWL_prior_1208 %>%
  mutate(shade = GPP_mean < 0) %>%
  group_by(group = cumsum(c(0, diff(shade)) != 0)) %>%  # Group consecutive TRUE/FALSE regions
  filter(shade) %>%  # Keep only the TRUE (GPP_mean < 0) regions
  summarize(
    xmin = min(date),
    xmax = max(date)
  )

# Create the plot with background shading
plot1 <- ggplot(BWL_prior_1208, aes(x = date, y = DO.obs)) +
  # Add shaded regions
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_line(color = "blue") +
  geom_line(aes(y = DO.mod), color = "red") +
  scale_x_date(labels = date_format("%m-%y"),
                   breaks = date_breaks("120 days"))+
  theme_bw() +
  labs(title = info, x = "Date (month-year)", y = "DO Value") +
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_prior_1208$date)),
                                    y = min(na.omit(BWL_prior_1208$DO.obs)) + 0.4,
                                    label = paste("RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)

plot1


footnote<- "Shaded areas indicate GPP is negative or ER is postive"

# Create the plot with background shading
gpp_plot1 <- ggplot(BWL_prior_1208, aes(x = date, y = GPP_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymin = GPP_2.5pct, ymax = GPP_97.5pct),
              fill = "#377a46",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#377a46") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  ylim(-3.25, 3.25) +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(GPP~(g~O[2]~m^-2~d^-1))) + 
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_prior_1208$date)),
                                    y = max(na.omit(BWL_prior_1208$GPP_97.5pct)) + 0.5,
                                    label = paste("DO RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)
gpp_plot1

# Create the plot with background shading
er_plot1 <- ggplot(BWL_prior_1208, aes(x = date, y = ER_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymax = ER_daily_2.5pct, ymin = ER_daily_97.5pct),
              fill = "#b36a22",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#b36a22") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(ER~(g~O[2]~m^-2~d^-1))) 
  
er_plot1




#####
BWL_filtered_1211 <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_measured_prior_noGPP_prior_NLDAS_mod_and_obs.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d"))  %>% 
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
info <- "Filtered + no prior on ER/GPP"
BWL_filtered_1211<- BWL_filtered_1211 %>%
  filter(!is.na(DO.mod), !is.na(DO.obs))
rmse_values <- BWL_filtered_1211 %>%
  # group_by(y_var) %>%
  summarize(RMSE = sqrt(mean((DO.mod - DO.obs)^2, na.rm = TRUE)))

BWL_filtered_1211 <- BWL_filtered_1211 %>%
  left_join(BWL_filtered_1211g, by=c("date"))

str(BWL_filtered_1211)
# Create a dataframe with shading regions based on GPP_mean < 0
shading_data <- BWL_filtered_1211 %>%
  mutate(shade = GPP_mean < -0.1) %>%
  group_by(group = cumsum(c(0, diff(shade)) != 0)) %>%  # Group consecutive TRUE/FALSE regions
  filter(shade) %>%  # Keep only the TRUE (GPP_mean < 0) regions
  summarize(
    xmin = min(date),
    xmax = max(date)
  )

BWL_filtered_1211 <- date_df%>%
  left_join(BWL_filtered_1211, by=c("date"))

plot2<-ggplot(BWL_filtered_1211, aes(x = date, y = DO.obs)) +
  # Add shaded regions
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_line(color = "blue") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  geom_line(aes(y = DO.mod), color = "red") + theme_bw() +
  labs(title =  info, x = "Date (month-year)", y = "DO Value") +
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_filtered_1211$date)), 
                                    y = min(na.omit(BWL_filtered_1211$DO.obs)+0.4), 
                                    label = paste("RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1) +
  annotate("text", 
           x = mean(BWL_filtered_1211$date), y = min(na.omit(BWL_filtered_1211$DO.obs)), 
           label = "Shaded areas indicate GPP is negative or ER is positive",
           #hjust = 1, vjust = -1, 
           size = 3, color = "black")






# Create the plot with background shading
gpp_plot2 <- ggplot(BWL_filtered_1211, aes(x = date, y = GPP_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymin = GPP_2.5pct, ymax = GPP_97.5pct),
              fill = "#377a46",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#377a46") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  ylim(-3.25, 3.25) +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(GPP~(g~O[2]~m^-2~d^-1))) + 
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_prior_1208$date)),
                                    y = max(na.omit(BWL_prior_1208$GPP_97.5pct)) + 0.5,
                                    label = paste("DO RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)+
  annotate("text", 
           x = mean(BWL_filtered_1211$date), y = min(na.omit(BWL_filtered_1211$GPP_daily_mean-5)), 
           label = "Shaded areas indicate GPP is negative or ER is positive",
           #hjust = 1, vjust = -1, 
           size = 3, color = "black")


gpp_plot2

# Create the plot with background shading
er_plot2 <- ggplot(BWL_filtered_1211, aes(x = date, y = ER_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymax = ER_daily_2.5pct, ymin = ER_daily_97.5pct),
              fill = "#b36a22",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#b36a22") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(ER~(g~O[2]~m^-2~d^-1))) 

er_plot2




#####
BWL_ufiltered_1213 <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_measured_prior_noGPP_prior_no_filter_NLDAS_mod_and_obs.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d"))%>% 
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
info <- "Unfiltered + no prior on ER/GPP"
BWL_ufiltered_1213<- BWL_ufiltered_1213 %>%
  filter(!is.na(DO.mod), !is.na(DO.obs))
rmse_values <- BWL_ufiltered_1213 %>%
 # group_by(y_var) %>%
  summarize(RMSE = sqrt(mean((DO.mod - DO.obs)^2, na.rm = TRUE)))


BWL_ufiltered_1213 <- BWL_ufiltered_1213 %>%
  left_join(BWL_ufiltered_1213g, by=c("date"))

str(BWL_ufiltered_1213)
# Create a dataframe with shading regions based on GPP_mean < 0
shading_data <- BWL_ufiltered_1213 %>%
  mutate(shade = GPP_mean < -0.1) %>%
  group_by(group = cumsum(c(0, diff(shade)) != 0)) %>%  # Group consecutive TRUE/FALSE regions
  filter(shade) %>%  # Keep only the TRUE (GPP_mean < 0) regions
  summarize(
    xmin = min(date),
    xmax = max(date)
  )


BWL_ufiltered_1213 <- date_df%>%
  left_join(BWL_ufiltered_1213, by=c("date"))


plot3 <- ggplot(BWL_ufiltered_1213, aes(x = date, y = DO.obs)) +
  geom_line(color = "blue") +
  # Add shaded regions
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_line(aes(y = DO.mod), color = "red") + theme_bw() +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  labs(title = info, x = "Date (month-year)", y = "DO Value") +
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_ufiltered_1213$date)), 
                                    y = min(na.omit(BWL_ufiltered_1213$DO.obs) +0.4), 
                                    label = paste("RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)





# Create the plot with background shading
gpp_plot3 <- ggplot(BWL_ufiltered_1213, aes(x = date, y = GPP_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymin = GPP_2.5pct, ymax = GPP_97.5pct),
              fill = "#377a46",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#377a46") + 
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  ylim(-3.25, 3.25) +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(GPP~(g~O[2]~m^-2~d^-1))) + 
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_prior_1208$date)),
                                    y = max(na.omit(BWL_prior_1208$GPP_97.5pct)) + 0.5,
                                    label = paste("DO RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)
gpp_plot3

# Create the plot with background shading
er_plot3 <- ggplot(BWL_ufiltered_1213, aes(x = date, y = ER_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymax = ER_daily_2.5pct, ymin = ER_daily_97.5pct),
              fill = "#b36a22",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#b36a22") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(ER~(g~O[2]~m^-2~d^-1))) 

er_plot3






#####
BWL_ufiltered_0107 <- read.csv("/Users/kellyloria/Documents/Publications/2024_stream_metab_output/BWL/BWL_K600_m_prior_filter_NLDAS_0107_mod_and_obs.csv")%>%
  mutate(date= as.Date(date, format="%Y-%m-%d")) %>%
  filter(date> as.Date("2021-06-30") & date< as.Date("2024-08-01"))
info <- "hourly + no prior on ER/GPP"
BWL_ufiltered_0107<- BWL_ufiltered_0107 %>%
  filter(!is.na(DO.mod), !is.na(DO.obs))
rmse_values <- BWL_ufiltered_0107 %>%
  # group_by(y_var) %>%
  summarize(RMSE = sqrt(mean((DO.mod - DO.obs)^2, na.rm = TRUE)))


BWL_ufiltered_0107 <- BWL_ufiltered_0107 %>%
  left_join(BWL_hr_0107g, by=c("date"))

str(BWL_ufiltered_0107)
# Create a dataframe with shading regions based on GPP_mean < 0
shading_data <- BWL_ufiltered_0107 %>%
  mutate(shade = GPP_mean < -0.1) %>%
  group_by(group = cumsum(c(0, diff(shade)) != 0)) %>%  # Group consecutive TRUE/FALSE regions
  filter(shade) %>%  # Keep only the TRUE (GPP_mean < 0) regions
  summarize(
    xmin = min(date),
    xmax = max(date)
  )


BWL_ufiltered_0107 <- date_df%>%
  left_join(BWL_ufiltered_0107, by=c("date"))


plot4 <- ggplot(BWL_ufiltered_0107, aes(x = date, y = DO.obs)) +
  geom_line(color = "blue") +
  # Add shaded regions
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_line(aes(y = DO.mod), color = "red") + theme_bw() +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  labs(title = info, x = "Date (month-year)", y = "DO Value") +
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_ufiltered_1213$date)), 
                                    y = min(na.omit(BWL_ufiltered_1213$DO.obs) +0.4), 
                                    label = paste("RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)





# Create the plot with background shading
gpp_plot4 <- ggplot(BWL_ufiltered_0107, aes(x = date, y = GPP_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymin = GPP_2.5pct, ymax = GPP_97.5pct),
              fill = "#377a46",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#377a46") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  ylim(-3.25, 3.25) +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(GPP~(g~O[2]~m^-2~d^-1))) + 
  geom_text(data = rmse_values, aes(x = max(na.omit(BWL_prior_1208$date)),
                                    y = max(na.omit(BWL_prior_1208$GPP_97.5pct)) + 0.5,
                                    label = paste("DO RMSE:", round(RMSE, 3))),
            inherit.aes = FALSE, hjust = 1, vjust = 1)
gpp_plot4

# Create the plot with background shading
er_plot4 <- ggplot(BWL_ufiltered_0107, aes(x = date, y = ER_daily_mean)) +
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4, inherit.aes = FALSE) +
  # Plot lines for DO.obs and DO.mod
  geom_ribbon(aes(ymax = ER_daily_2.5pct, ymin = ER_daily_97.5pct),
              fill = "#b36a22",
              linetype = 0, alpha = 0.5) +
  geom_line(color = "#b36a22") +
  scale_x_date(labels = date_format("%m-%y"),
               breaks = date_breaks("120 days"))+
  theme_bw() +
  geom_hline(yintercept = 0)+
  labs(title = info, x = "Date (month-year)", y= expression(ER~(g~O[2]~m^-2~d^-1))) 

er_plot4


library(ggpubr)

GBL_dofit_grid <- ggarrange(plot1,
                             plot2,
                             plot3,
                            plot4,
                             ncol = 1, nrow = 4,
                             common.legend = TRUE, 
                             legend = "bottom")


# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/BWL_DO_fits_gridv3.png", plot = GBL_dofit_grid, width = 8.5, height = 9.5, units = "in")







GBL_metabfit_grid <- ggarrange(gpp_plot1,
                               gpp_plot2,
                               gpp_plot3,
                               gpp_plot4,
                            er_plot1,
                            er_plot2,
                            er_plot3,
                            er_plot4,
                            ncol = 4, nrow = 2,
                            common.legend = TRUE, 
                            legend = "bottom")

# ggsave("/Users/kellyloria/Documents/Publications/CH1\ biogeochem\ linkages\ /supp\ figures/BWL_metab_fits_grid_rmse_v2.png", plot = GBL_metabfit_grid, width = 20, height = 5, units = "in")



############
## Metab plots
############


