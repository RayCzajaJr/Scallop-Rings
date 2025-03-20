
library(dplyr)
library(lubridate)
library(MuMIn)
library(betareg)
library(ggridges)
library(lme4)
library(stats)
library(car)
library(viridis) 
library(DHARMa)
library(ggeffects)
library(effects)
library(mgcv)
library(DescTools)



####### OBJECTIVE: Compare mean/median sizes at sites where we have multiple yrs of data 



# Rename the DataFrame
sizedistscallops <- X_Master_ring_total_shell_ht_22Feb2025

unique_locations <- unique(sizedistscallops$Location)
print(unique_locations)

# Create the new 'location' column with grouped labels
sizedistscallops <- sizedistscallops %>%
  mutate(location = case_when(
    Location %in% c("East Marion", "East Marion (NE of Bay Ave)") ~ "East Marion",
    Location %in% c("Southold (Off Cedar Beach)", "Southold Bay (off Cedar Beach)", "Southold Bay (N central)") ~ "Southold - Cedar Beach",
    Location %in% c(
                    "NW Harbor (East side)", "NW Harbor (Off Split Rock)", 
                    "NW Harbor (S of Alewife Creek)", "NW Harbor (South of Alewife Creek)") ~ "NW Harbor - E Side",
    Location %in% c("NW Harbor (Barcelona Point)",  "NW Harbor (Barcelona Pt)", "NWHarbor (N of Mile Hill Road)") ~ "NW Harbor - Barcelona ",
    Location == "Flanders" ~ "Flanders",
    Location == "Noyack Bay (EW Side)" ~ "Noyack Bay",
    Location %in% c("Hallock Bay", "Hallock Bay (Outside Narrow River)", "Hallock Bay (central flats)") ~ "Hallock Bay",
    Location %in% c("Hog Neck", "Hog Neck Bay (SE corner)") ~ "Hog Neck",
    Location %in% c("OH Harbor", "OH Harbor (North)", "OH Harbor North", "OH Harbor - N") ~ "OH Harbor",
    Location %in% c("Robins Island (West central side)", "Robins Island (W side)") ~ "Robins Island",
    Location %in% c("Shelter Island (Hay Beach)") ~ "Shelter Island - H Beach",
    Location %in% c("Shelter Island (East Side)", "Shelter Island (NE Side)") ~ "Shelter Island - E Side",
    TRUE ~ Location # Keep original value if no match
  ))

unique_locations <- unique(sizedistscallops$location)
print(unique_locations)

sizedistscallops$sample <- paste(sizedistscallops$Location, sizedistscallops$date, sep = "_")

unique_sample <- unique(sizedistscallops$sample)
print(unique_sample)

# Convert date column to Date format
sizedistscallops <- sizedistscallops %>%
  mutate(date = mdy(date),  
         year = year(date)) 


# Make and apply Function to get the first full week of oct for a given year
get_first_full_week_oct <- function(year) {
  # Generate all October dates for the given year
  october_days <- seq(ymd(paste0(year, "-10-01")), ymd(paste0(year, "-10-31")), by = "day")
  
  # Find the first Monday in October
  first_monday <- october_days[wday(october_days) == 2][1]  # Monday is 2 (Sunday = 1)
  
  # If there's no Monday found (shouldn't happen), return empty
  if (is.na(first_monday)) return(as.Date(character()))
  
  # Define the full week (Monday to Sunday)
  week_range <- seq(first_monday, first_monday + 6, by = "day")
  
  return(week_range)
}

first_full_weeks_oct <- do.call(rbind, lapply(unique(sizedistscallops$year), function(y) {
  data.frame(year = y, date = get_first_full_week_oct(y))
}))

first_full_weeks_oct$date <- as.Date(first_full_weeks_oct$date) 

sizedistscallops_oct <- sizedistscallops %>%
  inner_join(first_full_weeks_oct, by = c("year", "date"))

# Make and apply Function to get the first full week of nov for a given year
get_first_full_week_nov <- function(year) {
  # Generate all October dates for the given year
  november_days <- seq(ymd(paste0(year, "-11-01")), ymd(paste0(year, "-11-30")), by = "day")
  
  # Find the first Monday in October
  first_monday <- november_days[wday(november_days) == 2][1]  # Monday is 2 (Sunday = 1)
  
  # If there's no Monday found (shouldn't happen), return empty
  if (is.na(first_monday)) return(as.Date(character()))
  
  # Define the full week (Monday to Sunday)
  week_range <- seq(first_monday, first_monday + 6, by = "day")
  
  return(week_range)
}

first_full_weeks_nov <- do.call(rbind, lapply(unique(sizedistscallops$year), function(y) {
  data.frame(year = y, date = get_first_full_week_nov(y))
}))

first_full_weeks_nov$date <- as.Date(first_full_weeks_nov$date) 

sizedistscallops_nov <- sizedistscallops %>%
  inner_join(first_full_weeks_nov, by = c("year", "date"))


# Convert year to character for analyses 
sizedistscallops_nov <- sizedistscallops_nov %>%
  mutate(Year = as.character(Year))

sizedistscallops_oct <- sizedistscallops_oct %>%
  mutate(Year = as.character(Year))

sizedistscallops<- sizedistscallops %>%
  mutate(Year = as.character(Year))


# Ridge plot for shell and ring heights by site and year (two seperte months for shell heights)
novtotalhtyear<-ggplot(sizedistscallops_nov, aes(x = `Total Ht`, y = Year, fill = Year)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BurgYl") + 
  theme_bw() +  
  labs(
       x = "Total Height (Nov)",
       y = "Year") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(53, NA))

novtotalhtloc<-ggplot(sizedistscallops_nov, aes(x = `Total Ht`, y = location, fill = location)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BluGrn") + 
  theme_bw() +  
  labs(
       x = "Total Height (Nov)",
       y = "Site") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(53, NA))

octtotalhtyear<-ggplot(sizedistscallops_oct, aes(x = `Total Ht`, y = Year, fill = Year)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BurgYl") + 
  theme_bw() +  
  labs(
       x = "Total Height (Oct)",
       y = "Year") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(46, 88))

octtotalhtloc<-ggplot(sizedistscallops_oct, aes(x = `Total Ht`, y = location, fill = location)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BluGrn") + 
  theme_bw() +  
  labs(
       x = "Total Height (Oct)",
       y = "Site") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(46, 88))

ridgeplots<-grid.arrange(octtotalhtloc, octtotalhtyear, novtotalhtloc, novtotalhtyear,
                                 ncol = 2, nrow = 2)

ggsave("ridgeplots.tiff",ridgeplots, dpi = 300, bg = "white",
       width = 24,
       height = 24,
       units = "cm")

ringhtloc<-ggplot(sizedistscallops, aes(x = `Ring Ht`, y = location, fill = location)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BluGrn") + 
  theme_bw() +  
  labs(
       x = "Ring Height",
       y = "Site") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 80))

# Boxplot since there are too many years for a nice ridge plot
ringhtyr<-ggplot(sizedistscallops, aes(x = as.factor(Year), y = `Ring Ht`, fill = as.factor(Year))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +  
  scale_fill_carto_d(palette = "BurgYl") +  
  theme_bw() +  
  labs(
       x = "Year",  
       y = "Ring Height") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 78))

ringhts<-grid.arrange(ringhtloc, ringhtyr,
                         ncol = 1, nrow = 2)

ggsave("ringhts.tiff",ringhts, dpi = 300, bg = "white",
       width = 20,
       height = 24,
       units = "cm")


# Prelim analyses for differences in shell and ring ht by year and location.... may wanna switch to a KS test?
moct <- aov(sizedistscallops_oct$`Total Ht`~ location + Year + location * Year, data=sizedistscallops_oct)
summary(moct)
res1 <- simulateResiduals(moct)
plot(res1)
mnov<- aov(sizedistscallops_nov$`Total Ht`~ location + Year + location * Year, data=sizedistscallops_nov)
summary(mnov)
mring <- aov(sizedistscallops$`Ring Ht`~ location + Year + location * Year, data=sizedistscallops)
summary(mring)


####### OBJECTIVE: Examine relationship between growth after the winter vs ring size (ie size reached before the winter)


# Add a month column to each dataframe before merging
sizedistscallops_oct <- sizedistscallops_oct %>%
  mutate(Month = "October")

sizedistscallops_nov <- sizedistscallops_nov %>%
  mutate(Month = "November")

# Merge the two monthly dataframes into one
sizedistscallops_octandnov <- bind_rows(sizedistscallops_oct, sizedistscallops_nov)

colnames(sizedistscallops_octandnov) <- gsub(" ", "_", colnames(sizedistscallops_octandnov))

m <- lmer(Growthafterfirstwinter~ Ring_Ht + (1 | location) + (1 | Year), 
          data = sizedistscallops_octandnov)
Anova(m)
r.squaredGLMM(m)
res1 <- simulateResiduals(m)
plot(res1)

efct <- effect("Ring_Ht", mod = m, xlevels = list(Ring_Ht = seq(0, 72, length.out = 100)))

efct <- as.data.frame(efct)

growthpostwinter<-ggplot() +
  geom_point(data=sizedistscallops_octandnov, aes(x = Ring_Ht, y =Growthafterfirstwinter), alpha = 0.4, color = "black")+  
  geom_line(data = efct, aes(x = Ring_Ht, y = fit), color = "darkseagreen4", size=1) +
  geom_ribbon(data = efct, aes(x = Ring_Ht, ymin = lower, ymax = upper),  fill = "darkseagreen4", alpha = 0.3) +
  theme_bw() +
  labs(
       x = "Ring Height",
       y = "Growth after 1st Winter") +
  theme(plot.title = element_text(hjust = 0.5))  

ggsave("growthpostwinter.tiff",growthpostwinter, dpi = 300, bg = "white",
       width = 18,
       height = 16,
       units = "cm")




####### OBJECTIVE: Examine relationship between % ripe in Oct/Nov vs small ring adults the next yea


scallopringsummarydf <- X_2024_master_scallop_ring_paper_summary_data_highlighting_revised15Mar2025_data_UNCHANGED_4

# Create column for days post 9-30
scallopringsummarydf$Date <- as.Date(scallopringsummarydf$Date)
scallopringsummarydf$dayspostsep30 <- as.numeric(scallopringsummarydf$Date - as.Date(paste0(format(scallopringsummarydf$Date, "%Y"), "-09-30")))

# Convert percent to decimal and remove NAS for percent ripe
scallopringsummarydf <- scallopringsummarydf[!is.na(scallopringsummarydf$orangegonad_percentripe), ]
scallopringsummarydf$orangegonad_percentripe <- scallopringsummarydf$orangegonad_percentripe / 100

# Convert percent to decimal and remove NAS for percent small rings
scallopringsummarydf <- scallopringsummarydf[!is.na(scallopringsummarydf$percentsmallrings_lessthan20mm), ]
scallopringsummarydf$percentsmallrings_lessthan20mm <- scallopringsummarydf$percentsmallrings_lessthan20mm / 100

# Remove months after the fall
scallopringsummarydf <- scallopringsummarydf[scallopringsummarydf$dayspostsep30 <= 60, ]

# Get annnaul means with lagged predictor
scallopsummarydf_annualmeans <- scallopringsummarydf %>%
  mutate(year = format(Date, "%Y")) %>%
  group_by(year) %>%
  summarise(
    mean_orangegonad = mean(orangegonad_percentripe, na.rm = TRUE),
    mean_percentsmallrings = mean(percentsmallrings_lessthan20mm, na.rm = TRUE)
  ) %>%
  ungroup()

scallopsummarydf_annualmeans <- scallopsummarydf_annualmeans %>%
  arrange(year) %>%
  mutate(lagged_orangegonad = lag(mean_orangegonad)) 

scallopsummarydf_annualmeans <- na.omit(scallopsummarydf_annualmeans)

plot(scallopsummarydf_annualmeans$lagged_orangegonad, scallopsummarydf_annualmeans$mean_percentsmallrings)

# Make beta regression model and plot
mygam <- gam(mean_percentsmallrings~ lagged_orangegonad, family=betar(link="logit"), data = scallopsummarydf_annualmeans)
summary(mygam)
res1 <- simulateResiduals(mygam)
plot(res1)
min <- min(scallopsummarydf_annualmeans$lagged_orangegonad)
max <- max(scallopsummarydf_annualmeans$lagged_orangegonad)
new.x <- expand.grid(lagged_orangegonad = seq(min, max, length.out = 1000))
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, y = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = y - 1.96 * SE, upr = y + 1.96 * SE) # calculating the 95% confidence interval
addThese <- rename(addThese, mean_percentsmallrings = y)
RipePlot<-ggplot(scallopsummarydf_annualmeans, aes(x = lagged_orangegonad, y = mean_percentsmallrings)) +
  geom_point(size =2.5, alpha = .75)+
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity',color="darkseagreen4")+
  theme_bw() +
  ylab("Percent Small Rings the Next Year")+
  xlab("Percent Ripe Scallops")+
  theme(text = element_text(size=10)) +
  theme(panel.background = element_blank())
RipePlot

ggsave("RipePlot.tiff",RipePlot, dpi = 300, bg = "white",
       width = 14,
       height = 16,
       units = "cm")

### Create scallopsummarydf_annualmeans with annual means SEPERATE SITES

scallopsummarydf_annualmeans <- scallopringsummarydf %>%
  mutate(year = format(Date, "%Y")) %>%
  group_by(year, Site) %>%  # Group by both year and Site
  summarise(
    mean_orangegonad = mean(orangegonad_percentripe, na.rm = TRUE),
    mean_percentsmallrings = mean(percentsmallrings_lessthan20mm, na.rm = TRUE)
  ) %>%
  ungroup()

scallopsummarydf_annualmeans <- scallopsummarydf_annualmeans %>%
  arrange(year) %>%
  mutate(lagged_orangegonad = lag(mean_orangegonad)) # Create lagged predictor

scallopsummarydf_annualmeans <- na.omit(scallopsummarydf_annualmeans)

plot(scallopsummarydf_annualmeans$lagged_orangegonad, scallopsummarydf_annualmeans$mean_percentsmallrings)

# Create the new 'location' column with grouped labels
scallopsummarydf_annualmeans<- scallopsummarydf_annualmeans %>%
  mutate(Site = case_when(
    Site %in% c("NW Harbor (Off Split Rock)", "NW Harbor - S of Alewife Creek", "NW Harbor - Barcelona Neck", "NW Harbor - N of Mile Hill Rd") ~ "NW Harbor",
    Site %in% c("Southold - off Cedar Beach", "Southold Bay - N central hole") ~ "Southold",
    Site %in% c("Robin's Island - W side", "Robin's Island - W side (dredged)") ~ "Robins",
    Site %in% c("Shelter Island - N Tip  (Off Hay Beach)", "Shelter Island - E side") ~ "Shelter Island",
    TRUE ~ Site # Keep original value if no match
  ))

mygam <- gam(mean_percentsmallrings~ lagged_orangegonad, family=betar(link="logit"), data = scallopsummarydf_annualmeans)
summary(mygam) 

#model has trouble converging when trying to keep sites/bays separate, so prob just go with previous model where we used annual means



####### OBJECTIVE: Examine relationship between the number of days past the start of fall and fecundity



scallopringsummarydf <- X_2024_master_scallop_ring_paper_summary_data_highlighting_revised15Mar2025_data_UNCHANGED_4

scallopringsummarydf$Date <- as.Date(scallopringsummarydf$Date)
scallopringsummarydf$dayspostsep30 <- as.numeric(scallopringsummarydf$Date - as.Date(paste0(format(scallopringsummarydf$Date, "%Y"), "-09-30")))

scallopringsummarydf <- scallopringsummarydf[!is.na(scallopringsummarydf$orangegonad_percentripe), ]
scallopringsummarydf$orangegonad_percentripe <- scallopringsummarydf$orangegonad_percentripe / 100


mygam <- gam(orangegonad_percentripe~ dayspostsep30, family=betar(link="logit"), data = scallopringsummarydf)
summary(mygam)
min <- min(scallopringsummarydf$dayspostsep30)
max <- max(scallopringsummarydf$dayspostsep30)
new.x <- expand.grid(dayspostsep30 = seq(min, max, length.out = 1000))
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, y = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = y - 1.96 * SE, upr = y + 1.96 * SE) # calculating the 95% confidence interval
addThese <- rename(addThese, orangegonad_percentripe = y)
RipePlotDaysPost<-ggplot(scallopringsummarydf, aes(x = dayspostsep30, y = orangegonad_percentripe )) +
  geom_point(size =2.5, alpha = .75)+
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity',color="darkseagreen4")+
  theme_bw() +
  ylab("Percent Ripe Scallops")+
  xlab("Days Post Sept 30")+
  theme(text = element_text(size=10)) +
  theme(panel.background = element_blank())
RipePlotDaysPost

ggsave("RipePlotDaysPost.tiff",RipePlotDaysPost, dpi = 300, bg = "white",
       width = 14,
       height = 14,
       units = "cm")



