### Script to analyze and plot DO concentrations and loads in the Upper Llano
## First section will plot concentrations at project sites and use linear regression
## to assess trends per TCEQ standards

### Second part will estimate streamflows at each station using the DAR approach

### Part three will use GAMs to estimate trends in concentrations and loads at each station

library(readr)
library(dataRetrieval)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(smwrBase)
library(mgcv)
library(nlme)
library(rloadest)
library(hrbrthemes)
update_geom_font_defaults(font_rc)

##Old code. See below for updated
# # Download Data -----------------------------------------------------------
# 
# ## Retrieve water quality info from storet
# LR1 <- readWQPqw(siteNumbers = c("TCEQMAIN-17471"), parameterCd = "Oxygen")
# LR2 <- readWQPqw(siteNumbers = c("TCEQMAIN-21489"), parameterCd = "Oxygen")
# NL1 <- readWQPqw(siteNumbers = c("TCEQMAIN-17425"), parameterCd = "Oxygen")
# NL2 <- readWQPqw(siteNumbers = c("TCEQMAIN-21548"), parameterCd = "Oxygen")
# SL1 <- readWQPqw(siteNumbers = c("TCEQMAIN-18197"), parameterCd = "Oxygen")
# 
# ## Retrieve discharge from NWIS
# USGSLR <- readNWISdv("08150000", "00060", "2001-01-01", "2016-12-31")
# USGSLR <- renameNWISColumns(USGSLR)
# 
# USGSNL <- readNWISdv("08148500", "00060", "2001-01-01", "2016-12-31")
# USGSNL <- renameNWISColumns(USGSNL)
# 
# USGSSL <- readNWISdv("08149900", "00060", "2001-01-01", "2016-12-31")
# USGSSL <- renameNWISColumns(USGSSL)
# 
# ## Clean storet data
# 
# LR1 <- select(LR1, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue, Unit = ResultMeasure.MeasureUnitCode)
# LR2 <- select(LR2, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue, Unit = ResultMeasure.MeasureUnitCode)
# NL1 <- select(NL1, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue, Unit = ResultMeasure.MeasureUnitCode)
# NL2 <- select(NL2, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue, Unit = ResultMeasure.MeasureUnitCode)
# SL1 <- select(SL1, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue, Unit = ResultMeasure.MeasureUnitCode)
# 
# #need to rename station to the names in the report
# cframe <- bind_rows(LR1, LR2, NL1, NL2, SL1)
# cframe <- cframe %>%
#   mutate(Location = case_when(
#     Location == "TCEQMAIN-17471" ~ "LR1",
#     Location == "TCEQMAIN-21489" ~ "LR2",
#     Location == "TCEQMAIN-17425" ~ "NL1",
#     Location == "TCEQMAIN-21548" ~ "NL2",
#     Location == "TCEQMAIN-18197" ~ "SL1"
#   ))
# cframe$ndate <- dectime(cframe$Date)
# cframe <- as.data.frame(cframe)
# cframe <- cframe[cframe$Unit == "mg/l",]
# 
# 
# #fit linear model for each station
# lm.lr1 <- lm(Value ~ ndate, data = cframe[cframe$Location=="LR1",])
# plot(lm.lr1)
# a <- summary(lm.lr1)
# lm.lr2 <- lm(Value ~ ndate, data = cframe[cframe$Location=="LR2",])
# plot(lm.lr2)
# summary(lm.lr2)
# lm.nl1 <- lm(Value ~ ndate, data = cframe[cframe$Location=="NL1",])
# plot(lm.nl1)
# summary(lm.nl1)
# lm.nl2 <- lm(Value ~ ndate, data = cframe[cframe$Location=="NL2",])
# plot(lm.nl2)
# summary(lm.nl2)
# lm.sl1 <- lm(Value ~ ndate, data = cframe[cframe$Location=="SL1",])
# plot(lm.sl1)
# summary(lm.sl1)
# 
# a <- data.frame(Location = as.factor(c("LR1", "LR2", "NL1", "NL2", "SL1")), # verbose and not R-like but works
#                 pvalue = c(round(coef(summary(lm.lr1))[2,4],3), 
#                            round(coef(summary(lm.lr2))[2,4],3),
#                            round(coef(summary(lm.nl1))[2,4],3),
#                            round(coef(summary(lm.nl2))[2,4],3),
#                            round(coef(summary(lm.sl1))[2,4],3)),
#                 tstat = c(coef(summary(lm.lr1))[2,3], 
#                           coef(summary(lm.lr2))[2,3],
#                           coef(summary(lm.nl1))[2,3],
#                           coef(summary(lm.nl2))[2,3],
#                           coef(summary(lm.sl1))[2,3]),
#                 Date = as.Date(c("2001-1-1","2010-1-1","2010-1-1","2001-1-1","2001-1-1")),
#                 Value = c(10,10,10,10,10))
# 
# p <- ggplot() + 
#   geom_point(data = cframe, aes(Date, Value), color = "dodgerblue", alpha = 0.75) + 
#   stat_smooth(data = cframe, aes(Date, Value), method = "lm") + 
#   geom_label(data = a, aes(Date, Value, label = paste("atop(italic(T)==",round(tstat,3),",italic(p)==",pvalue,")")), 
#              family = "serif", size = 3, parse = T, hjust = "left",vjust = "bottom", alpha = 0.5) +
#   geom_hline(yintercept = 5, linetype = 2) +
#   facet_wrap(~Location) +
#   background_grid() +
#   ylab("Dissolved Oxygen (mg/L)")
# p
# ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DOConcentration.png", width = 6.5, height = 4, units = "in")



##try with CRP data
CRPField <- read_delim("RawData/CRPField.txt", 
                       "|", escape_double = FALSE, col_types = cols(`End Date` = col_date(format = "%m/%d/%Y")), 
                       trim_ws = TRUE)

## Retrieve discharge from NWIS
USGSLR <- readNWISdv("08150000", "00060", "2001-01-01", "2016-12-31")
USGSLR <- renameNWISColumns(USGSLR)

USGSNL <- readNWISdv("08148500", "00060", "2001-01-01", "2016-12-31")
USGSNL <- renameNWISColumns(USGSNL)

USGSSL <- readNWISdv("08149900", "00060", "2001-01-01", "2016-12-31")
USGSSL <- renameNWISColumns(USGSSL)

#clean
CRPField <- select(CRPField, Date = `End Date`, Location = `Station ID`, Value = `Value`, Parameter = `Parameter Code`, ParameterDesc = `Parameter Description`)
CRPField <- CRPField[CRPField$Parameter == "00300",]
cframe <- CRPField %>%
  mutate(Location = case_when(
    Location == "17471" ~ "LR1",
    Location == "21489" ~ "LR2",
    Location == "17425" ~ "NL1",
    Location == "21548" ~ "NL2",
    Location == "18197" ~ "SL1"
  ))
cframe$ndate <- dectime(cframe$Date)
cframe <- as.data.frame(cframe)


#fit linear model for each station
lm.lr1 <- lm(Value ~ ndate, data = cframe[cframe$Location=="LR1",])
plot(lm.lr1)
summary(lm.lr1)
lm.lr2 <- lm(Value ~ ndate, data = cframe[cframe$Location=="LR2",])
plot(lm.lr2)
summary(lm.lr2)
lm.nl1 <- lm(Value ~ ndate, data = cframe[cframe$Location=="NL1",])
plot(lm.nl1)
summary(lm.nl1)
lm.nl2 <- lm(Value ~ ndate, data = cframe[cframe$Location=="NL2",])
plot(lm.nl2)
summary(lm.nl2)
lm.sl1 <- lm(Value ~ ndate, data = cframe[cframe$Location=="SL1",])
plot(lm.sl1)
summary(lm.sl1)



a <- data.frame(Location = as.factor(c("LR1", "LR2", "NL1", "NL2", "SL1")), # verbose and not R-like but works
                pvalue = c(round(coef(summary(lm.lr1))[2,4],3), 
                           round(coef(summary(lm.lr2))[2,4],3),
                           round(coef(summary(lm.nl1))[2,4],3),
                           round(coef(summary(lm.nl2))[2,4],3),
                           round(coef(summary(lm.sl1))[2,4],3)),
                tstat = c(coef(summary(lm.lr1))[2,3], 
                          coef(summary(lm.lr2))[2,3],
                          coef(summary(lm.nl1))[2,3],
                          coef(summary(lm.nl2))[2,3],
                          coef(summary(lm.sl1))[2,3]),
                Date = as.Date(c("2001-1-1","2010-1-1","2010-1-1","2001-1-1","2001-1-1")),
                Value = c(10,10,10,10,10))

p <- ggplot() + 
  geom_point(data = cframe, aes(Date, Value), color = "dodgerblue", alpha = 0.75) + 
  stat_smooth(data = cframe, aes(Date, Value), method = "lm") + 
  #geom_label(data = a, aes(Date, Value, label = paste("atop(italic(T)==",round(tstat,3),",italic(p)==",pvalue,")")), 
  #           family = "serif", size = 3, parse = T, hjust = "left",vjust = "bottom", alpha = 0.5) +
  geom_hline(yintercept = 5, linetype = 2) +
  facet_wrap(~Location) +
  background_grid() +
  ylab("Dissolved Oxygen (mg/L)") +
  theme_ipsum_rc()
p
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DOConcentration.png", width = 6.5, height = 4, units = "in")

# Estimate Streamflow -----------------------------------------------------

# assign values of phi based on streamflow perentile

phi <- function(data, ...){
  d <- data %>%
    mutate(phi = case_when(
      cumdist >= 0 & cumdist < 0.06 ~ 0.885,
      cumdist >= 0.06 & cumdist < 0.10 ~ 0.886,
      cumdist >= 0.10 & cumdist < 0.18 ~ 0.887,
      cumdist >= 0.18 & cumdist < 0.24 ~ 0.888,
      cumdist >= 0.24 & cumdist < 0.30 ~ 0.889,
      cumdist >= 0.30 & cumdist < 0.34 ~ 0.890,
      cumdist >= 0.34 & cumdist < 0.36 ~ 0.891,
      cumdist >= 0.36 & cumdist < 0.38 ~ 0.892,
      cumdist >= 0.38 & cumdist < 0.40 ~ 0.893,
      cumdist >= 0.40 & cumdist < 0.42 ~ 0.894,
      cumdist >= 0.42 & cumdist < 0.44 ~ 0.895,
      cumdist >= 0.44 & cumdist < 0.46 ~ 0.897,
      cumdist >= 0.46 & cumdist < 0.48 ~ 0.899,
      cumdist >= 0.48 & cumdist < 0.50 ~ 0.902,
      cumdist >= 0.50 & cumdist < 0.52 ~ 0.905,
      cumdist >= 0.52 & cumdist < 0.54 ~ 0.908,
      cumdist >= 0.54 & cumdist < 0.56 ~ 0.912,
      cumdist >= 0.56 & cumdist < 0.58 ~ 0.916,
      cumdist >= 0.58 & cumdist < 0.60 ~ 0.920,
      cumdist >= 0.60 & cumdist < 0.62 ~ 0.924,
      cumdist >= 0.62 & cumdist < 0.64 ~ 0.927,
      cumdist >= 0.64 & cumdist < 0.66 ~ 0.930,
      cumdist >= 0.66 & cumdist < 0.68 ~ 0.932,
      cumdist >= 0.68 & cumdist < 0.70 ~ 0.934,
      cumdist >= 0.70 & cumdist < 0.76 ~ 0.935,
      cumdist >= 0.76 & cumdist < 0.80 ~ 0.934,
      cumdist >= 0.80 & cumdist < 0.84 ~ 0.933,
      cumdist >= 0.84 & cumdist < 0.86 ~ 0.931,
      cumdist >= 0.86 & cumdist < 0.88 ~ 0.927,
      cumdist >= 0.88 & cumdist < 0.90 ~ 0.920,
      cumdist >= 0.90 & cumdist < 0.92 ~ 0.906,
      cumdist >= 0.92 & cumdist < 0.94 ~ 0.890,
      cumdist >= 0.94 & cumdist < 0.95 ~ 0.865,
      cumdist >= 0.95 & cumdist < 0.96 ~ 0.850,
      cumdist >= 0.96 & cumdist < 0.97 ~ 0.830,
      cumdist >= 0.97 & cumdist < 0.98 ~ 0.806,
      cumdist >= 0.98 & cumdist < 0.99 ~ 0.773,
      cumdist >= 0.99 & cumdist < 1.00 ~ 0.737,
      cumdist == 1.00 ~ 0.700
    ))
  d
}

# DARs used:
# lr = Q[usgslr]
# lr2 = Q[usgslr] * (4,789/4,809)^phi
# nr1 = Q[usgsnl] * (2,372/2,332)^phi
# nr2 = Q[usgsnl] * (2,342/2,332)^phi
# sr1 = Q[usgssl] * (1,941/2,277)^phi

# Calculate flow exceedance and phi at each gage
USGSLR <- USGSLR %>%
  arrange(Flow) %>%
  mutate(cumdist = cume_dist(-Flow))
USGSLR <- phi(USGSLR)

USGSNL <- USGSNL %>%
  arrange(Flow) %>%
  mutate(cumdist = cume_dist(-Flow))
USGSNL <- phi(USGSNL)

USGSSL <- USGSSL %>%
  arrange(Flow) %>%
  mutate(cumdist = cume_dist(-Flow))
USGSSL <- phi(USGSSL)

# create plot and report of USGS streamflow data
fframe <- bind_rows(USGSLR, USGSNL, USGSSL)

usgsplot <- ggplot(fframe) +
  geom_line(aes(Date, Flow)) +
  facet_wrap(~site_no, ncol = 1) +
  background_grid() +
  ylab("Daily Discharge (cfs)")
usgsplot
# add columns to calculate flows at each station using DAR^phi then gather into tidy format

fframe <- fframe %>%
  mutate(
    LR1 = case_when(
      site_no == "08150000" ~ Flow),
    LR2 = case_when(
      site_no == "08150000" ~ Flow * (4789/4809)^phi),
    NL1 = case_when(
      site_no == "08148500" ~ Flow * (2372/2332)^phi),
    NL2 = case_when(
      site_no == "08148500" ~ Flow * (2342/2332)^phi),
    SL1 = case_when(
      site_no == "08149900" ~ Flow * (1941/2277)^phi )) %>%
  gather(key = station, value = estflow, LR1, LR2, NL1, NL2, SL1 ) %>%
  filter(!is.na(estflow))

fframe <- fframe %>%
  left_join(cframe, by = c("Date" = "Date", "station" = "Location"))



##import temperature

CRPField <- read_delim("RawData/CRPField.txt", 
                       "|", escape_double = FALSE, col_types = cols(`End Date` = col_date(format = "%m/%d/%Y")), 
                       trim_ws = TRUE)
#clean
CRPField <- select(CRPField, Date = `End Date`, Location = `Station ID`, Value = `Value`, Parameter = `Parameter Code`, ParameterDesc = `Parameter Description`)
CRPField <- CRPField[CRPField$Parameter == "00010",]
cframe <- CRPField %>%
  mutate(Location = case_when(
    Location == "17471" ~ "LR1",
    Location == "21489" ~ "LR2",
    Location == "17425" ~ "NL1",
    Location == "21548" ~ "NL2",
    Location == "18197" ~ "SL1"
  ))
cframe$ndate <- dectime(cframe$Date)
cframe <- as.data.frame(cframe)

fframe <- fframe %>%
  left_join(cframe, by = c("Date" = "Date", "station" = "Location"))
fframe <- fframe %>%
  rename(DissolvedOxygen = Value.x)
fframe <- fframe %>%
  rename(Temperature = Value.y)
fframe$lflow <- log(fframe$estflow + 0.001)


##import total nitrate nitrogen
## need to figure out
CRPField <- read_delim("RawData/CRPNutrient.txt", 
                       "|", escape_double = FALSE, col_types = cols(`End Date` = col_date(format = "%m/%d/%Y")), 
                       trim_ws = TRUE)
#clean
CRPField <- select(CRPField, Date = `End Date`, Location = `Station ID`, Value = `Value`, Parameter = `Parameter Code`, ParameterDesc = `Parameter Description`)
CRPField <- as.data.frame(CRPField)
CRPField <- CRPField[CRPField$Parameter == "00630",]
cframe <- CRPField %>%
  mutate(Location = case_when(
    Location == "17471" ~ "LR1",
    Location == "21489" ~ "LR2",
    Location == "17425" ~ "NL1",
    Location == "21548" ~ "NL2",
    Location == "18197" ~ "SL1"
  ))
cframe$ndate <- dectime(cframe$Date)
cframe <- as.data.frame(cframe)

fframe <- fframe %>%
  left_join(cframe, by = c("Date" = "Date", "station" = "Location"))

fframe <- fframe %>%
  rename(TN = Value)

##import total TP
## need to figure out
CRPField <- read_delim("RawData/CRPNutrient.txt", 
                       "|", escape_double = FALSE, col_types = cols(`End Date` = col_date(format = "%m/%d/%Y")), 
                       trim_ws = TRUE)
#clean
CRPField <- select(CRPField, Date = `End Date`, Location = `Station ID`, Value = `Value`, Parameter = `Parameter Code`, ParameterDesc = `Parameter Description`)
CRPField <- as.data.frame(CRPField)
CRPField <- CRPField[CRPField$Parameter == "00665",]
cframe <- CRPField %>%
  mutate(Location = case_when(
    Location == "17471" ~ "LR1",
    Location == "21489" ~ "LR2",
    Location == "17425" ~ "NL1",
    Location == "21548" ~ "NL2",
    Location == "18197" ~ "SL1"
  ))
cframe$ndate <- dectime(cframe$Date)
cframe <- as.data.frame(cframe)

fframe <- fframe %>%
  left_join(cframe, by = c("Date" = "Date", "station" = "Location"))

fframe <- fframe %>%
  rename(TP = Value)
### TP values are nearly all at the detectible limit. Not much data we can get from this, remove from analysis

##import total cholorphyll-a
## need to figure out
CRPField <- read_delim("RawData/CRPNutrient.txt", 
                       "|", escape_double = FALSE, col_types = cols(`End Date` = col_date(format = "%m/%d/%Y")), 
                       trim_ws = TRUE)
#clean
CRPField <- select(CRPField, Date = `End Date`, Location = `Station ID`, Value = `Value`, Parameter = `Parameter Code`, ParameterDesc = `Parameter Description`)
CRPField <- as.data.frame(CRPField)
CRPField <- CRPField[CRPField$Parameter == "70953" | CRPField$Parameter == "32211",]
cframe <- CRPField %>%
  mutate(Location = case_when(
    Location == "17471" ~ "LR1",
    Location == "21489" ~ "LR2",
    Location == "17425" ~ "NL1",
    Location == "21548" ~ "NL2",
    Location == "18197" ~ "SL1"
  ))
cframe$ndate <- dectime(cframe$Date)
cframe <- as.data.frame(cframe)

fframe <- fframe %>%
  left_join(cframe, by = c("Date" = "Date", "station" = "Location"))

fframe <- fframe %>%
  rename(CA = Value)

### chlorophyll-a values are nearly all at the detectible limit. Not much data we can get from this, remove from analysis

library(GGally)


### some custom functions #####
my_custom_cor <- function(data, mapping, color = I("grey50"), sizeRange = c(2, 4), ...) {
  
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ct <- cor.test(x,y, method = "kendall")
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]
  
  # since we can't print it to get the strsize, just use the max size range
  cex <- max(sizeRange)
  
  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }
  
  # plot the cor value
  ggally_text(
    label = as.character(rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = I(percent_of_range(cex * abs(r), sizeRange)),
    color = color,
    ...
  ) + 
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = sig, 
      size = I(cex),
      color = color,
      ...
    ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        color = color, 
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}

sub <- ggpairs(fframe[fframe$station=="LR1",], columns = c(10, 14, 18, 19, 27),
               lower = list(continuous = wrap("smooth", color = "dodgerblue", alpha = 0.5)),
               diag = list(continuous = "barDiag"),
               upper = list(continuous = my_custom_cor),
               columnLabels = c("DO", "Temp", "log(Q)", "TN", "Chl-A"),
               axisLabels = "none",
               switch = "both")

  
sub <- sub +
  #theme_ipsum_rc() +
  theme(strip.background = element_blank())
sub
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DOLR1_corr.png", width = 3.25, height = 3.25, units = "in")

sub2 <- ggpairs(fframe[fframe$station=="LR2",], columns = c(10, 14, 18, 19, 27),
               lower = list(continuous = wrap("smooth", color = "dodgerblue", alpha = 0.5)),
               diag = list(continuous = "barDiag"),
               upper = list(continuous = my_custom_cor),
               columnLabels = c("DO", "Temp", "log(Q)", "TN", "Chl-A"),
               axisLabels = "none",
               switch = "both")


sub2 <- sub2 + theme(strip.background = element_blank())
sub2
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DOLR2_corr.png", width = 3.25, height = 3.25, units = "in")


sub3 <- ggpairs(fframe[fframe$station=="NL1",], columns = c(10, 14, 18),
                lower = list(continuous = wrap("smooth", color = "dodgerblue", alpha = 0.5)),
                diag = list(continuous = "barDiag"),
                upper = list(continuous = my_custom_cor),
                columnLabels = c("DO", "Temp", "log(Q)"),
                axisLabels = "none",
                switch = "both")


sub3 <- sub3 + theme(strip.background = element_blank())
sub3
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DONL1_corr.png", width = 3.25, height = 3.25, units = "in")


sub4 <- ggpairs(fframe[fframe$station=="NL2",], columns = c(10, 14, 18),
                lower = list(continuous = wrap("smooth", color = "dodgerblue", alpha = 0.5)),
                diag = list(continuous = "barDiag"),
                upper = list(continuous = my_custom_cor),
                columnLabels = c("DO", "Temp", "log(Q)"),
                axisLabels = "none",
                switch = "both")


sub4 <- sub4 + theme(strip.background = element_blank())
sub4
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DONL2_corr.png", width = 3.25, height = 3.25, units = "in")


sub5 <- ggpairs(fframe[fframe$station=="SL1",], columns = c(10, 14, 18),
                lower = list(continuous = wrap("smooth", color = "dodgerblue", alpha = 0.5)),
                diag = list(continuous = "barDiag"),
                upper = list(continuous = my_custom_cor),
                columnLabels = c("DO", "Temp", "log(Q)"),
                axisLabels = "none",
                switch = "both")


sub5 <- sub5 + theme(strip.background = element_blank())

sub5
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DOSL1_corr.png", width = 3.25, height = 3.25, units = "in")



### Fit GAM to DO concentrations at LR2 and NL1
library(mgcv)
library(nlme)

subset.fframe <- fframe[fframe$station == "LR2" | fframe$station == "NL1",]
subset.fframe$yday <- lubridate::yday(subset.fframe$Date)
subset.fframe$year <- lubridate::year(subset.fframe$Date)
subset.fframe$month <- as.numeric(lubridate::month(subset.fframe$Date))

lr2.subset.fframe <- subset.fframe[subset.fframe$station=="LR2",]
nl1.subset.fframe <- subset.fframe[subset.fframe$station=="NL1",]


do.model1.lr2 <- gamm(DissolvedOxygen ~ s(lflow) + s(Temperature) + s(year) + s(month, bs = "cc", k = 4),
                      data = lr2.subset.fframe)
plot(do.model1.lr2$gam)
plot(do.model1.lr2$lme)
summary(do.model1.lr2$gam)



do.model1.nl1 <- gamm(DissolvedOxygen ~ s(lflow) + s(Temperature) + s(year) + s(month, bs = "cc", k = 4),
                      data = nl1.subset.fframe)
plot(do.model1.nl1$gam)
plot(do.model1.nl1$lme)
summary(do.model1.nl1$gam)


### predict long term contributions

pdat <- data_frame(year = 2001:2015,
                   lflow = median(nl1.subset.fframe$lflow, na.rm = TRUE),
                   Temperature = mean(nl1.subset.fframe$Temperature, na.rm = TRUE),
                   month = 6
                   )
p <- predict.gam(do.model1.nl1$gam, type = "terms", se.fit = TRUE, newdata = pdat)
pdat$fit <- p$fit[,3] + mean(nl1.subset.fframe$DissolvedOxygen, na.rm = TRUE)
pdat$lower <- pdat$fit-p$se[,3]
pdat$upper <- pdat$fit+p$se[,3]
pdat$station <- "NL1"

pdat2 <- data_frame(year = 2001:2015,
                   lflow = median(lr2.subset.fframe$lflow, na.rm = TRUE),
                   Temperature = mean(lr2.subset.fframe$Temperature, na.rm = TRUE),
                   month = 6
)
p2 <- predict.gam(do.model1.lr2$gam, type = "terms", se.fit = TRUE, newdata = pdat)
pdat2$fit <- p2$fit[,3] + mean(lr2.subset.fframe$DissolvedOxygen, na.rm = TRUE)
pdat2$lower <- pdat2$fit-p2$se[,3]
pdat2$upper <- pdat2$fit+p2$se[,3]
pdat2$station <- "LR2"

pdat <- bind_rows(pdat, pdat2)

ggplot(pdat) +
  geom_line(aes(x=year, y=fit, color = station)) +
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper, fill = station), alpha = 0.5) +
  theme_ipsum_rc() + scale_color_ipsum() + scale_fill_ipsum() +
  labs(x = "Year", y = "DO (mg/L)") +
  theme(legend.title = element_blank(), legend.position = "bottom")
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/DOGAMM.png", width = 6, height = 4, units = "in")

## now plot changes in DO concentration as function of temp, month and year

min.temp <- min(do.model1.nl1$gam$model$Temperature)
max.temp <- max(do.model1.nl1$gam$model$Temperature)
pred.temps <- seq(from = min.temp, to = max.temp, by = 1)

min.date <- as.Date("2001-01-01")
max.date <- as.Date("2016-12-30")

pdat <- expand.grid(Date = seq.Date(min.date, max.date, by = "month"),
                    Temperature = pred.temps,
                    lflow = median(nl1.subset.fframe$lflow, na.rm = TRUE),
                    KEEP.OUT.ATTRS = TRUE)
pdat <- merTools:::stripAttributes(pdat)

pdat$month <- as.numeric(lubridate::month(pdat$Date))
pdat$monthlab <- lubridate::month(pdat$Date, label = TRUE, abbr = FALSE)
#pdat$yday <- yday(pdat$Date)
pdat$year <- year(pdat$Date)

## need to constrain temps to those observed per month
lim.vals <- group_by(data.frame(do.model1.nl1$gam$model), month) %>%
  summarize(Low = quantile(Temperature, 0.05, na.rm = TRUE),
            High = quantile(Temperature, 0.95, na.rm = TRUE))
lim.vals$month <- as.integer(lim.vals$month)
# merge limits with months
pdat <- left_join(pdat, lim.vals, by = 'month')

# reduce data
sel_vec <- with(pdat, 
                Temperature >= Low &
                  Temperature <= High
)
pdat <- pdat[sel_vec, !names(pdat) %in% c('Low', 'High')]
pdat <- arrange(pdat, year, month)

pdat <- pdat %>%
  filter(!is.na(lflow))

pred.do <- predict.gam(do.model1.nl1$gam, newdata = pdat)
pdat$fit <- pred.do

ggplot(pdat) +
  geom_line(aes(x = Temperature, y = fit, color = year, group = year),size=0.4) +
  facet_wrap(~monthlab) +
  geom_hline(yintercept = 5, linetype = 2) +
  viridis::scale_colour_viridis(direction = -1, option = "B") +
  labs(x = expression("Temperature "( degree*C)), y = "Concentration (mg/L)",title=NULL, subtitle = "NL1 Dissolved Oxygen") +
  theme_ipsum_rc() +
  theme(legend.position = "right",
        panel.spacing.x = unit(.2, "in"),
        panel.spacing.y = unit(.1, "in")) +
  guides(colour = guide_colourbar("Year", barwidth = 0.5))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/NL1DO_Temp.png", width = 6, height = 8, units = "in")


## now plot changes in DO concentration as function of temp, month and year

min.temp <- min(do.model1.lr2$gam$model$Temperature)
max.temp <- max(do.model1.lr2$gam$model$Temperature)
pred.temps <- seq(from = min.temp, to = max.temp, by = 1)

min.date <- as.Date("2001-01-01")
max.date <- as.Date("2016-12-30")

pdat <- expand.grid(Date = seq.Date(min.date, max.date, by = "month"),
                    Temperature = pred.temps,
                    lflow = median(lr2.subset.fframe$lflow, na.rm = TRUE),
                    KEEP.OUT.ATTRS = TRUE)
pdat <- merTools:::stripAttributes(pdat)

pdat$month <- as.numeric(lubridate::month(pdat$Date))
pdat$monthlab <- lubridate::month(pdat$Date, label = TRUE, abbr = FALSE)
#pdat$yday <- yday(pdat$Date)
pdat$year <- year(pdat$Date)

## need to constrain temps to those observed per month
lim.vals <- group_by(data.frame(do.model1.lr2$gam$model), month) %>%
  summarize(Low = quantile(Temperature, 0.05, na.rm = TRUE),
            High = quantile(Temperature, 0.95, na.rm = TRUE))
lim.vals$month <- as.integer(lim.vals$month)
# merge limits with months
pdat <- left_join(pdat, lim.vals, by = 'month')

# reduce data
sel_vec <- with(pdat, 
                Temperature >= Low &
                  Temperature <= High
)
pdat <- pdat[sel_vec, !names(pdat) %in% c('Low', 'High')]
pdat <- arrange(pdat, year, month)

pdat <- pdat %>%
  filter(!is.na(lflow))

pred.do <- predict.gam(do.model1.lr2$gam, newdata = pdat)
pdat$fit <- pred.do

ggplot(pdat) +
  geom_line(aes(x = Temperature, y = fit, color = year, group = year),size=0.4) +
  facet_wrap(~monthlab, ncol = 2) +
  geom_hline(yintercept = 5, linetype = 2) +
  viridis::scale_colour_viridis(direction = -1, option = "B") +
  labs(x = expression("Temperature "( degree*C)), y = "Concentration (mg/L)",title=NULL, subtitle = "LR2 Dissolved Oxygen") +
  theme_ipsum_rc() +
  theme(legend.position = "right",
        panel.spacing.x = unit(.2, "in"),
        panel.spacing.y = unit(.1, "in")) +
  guides(colour = guide_colourbar("Year", barwidth = 0.5))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/LR2DO_Temp.png", width = 6, height = 8, units = "in")
