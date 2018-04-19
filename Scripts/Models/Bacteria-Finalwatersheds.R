### Script to analyze and plot bacteria concentrations and loads in the Upper Llano
## First section will plot concentrations at project sites and use linear regression
## to assess trends per TCEQ standards

### Second part will estimate streamflows at each station using the DAR approach

### Part three will use GAMs to estimate trends in concentrations and loads at each station

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

# Download Data -----------------------------------------------------------

## Retrieve water quality info from storet
LR1 <- readWQPqw(siteNumbers = c("TCEQMAIN-17471"), parameterCd = "Escherichia coli")
LR2 <- readWQPqw(siteNumbers = c("TCEQMAIN-21489"), parameterCd = "Escherichia coli")
NL1 <- readWQPqw(siteNumbers = c("TCEQMAIN-17425"), parameterCd = "Escherichia coli")
NL2 <- readWQPqw(siteNumbers = c("TCEQMAIN-21548"), parameterCd = "Escherichia coli")
SL1 <- readWQPqw(siteNumbers = c("TCEQMAIN-18197"), parameterCd = "Escherichia coli")

## Retrieve discharge from NWIS
USGSLR <- readNWISdv("08150000", "00060", "2001-01-01", "2016-12-31")
USGSLR <- renameNWISColumns(USGSLR)

USGSNL <- readNWISdv("08148500", "00060", "2001-01-01", "2016-12-31")
USGSNL <- renameNWISColumns(USGSNL)

USGSSL <- readNWISdv("08149900", "00060", "2001-01-01", "2016-12-31")
USGSSL <- renameNWISColumns(USGSSL)



## Clean storet data
LR1 <- LR1[LR1$ResultMeasure.MeasureUnitCode!="hours",]
LR2 <- LR2[LR2$ResultMeasure.MeasureUnitCode!="hours",]
NL1 <- NL1[NL1$ResultMeasure.MeasureUnitCode!="hours",]
NL2 <- NL2[NL2$ResultMeasure.MeasureUnitCode!="hours",]
SL1 <- SL1[SL1$ResultMeasure.MeasureUnitCode!="hours",]

LR1 <- select(LR1, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue)
LR2 <- select(LR2, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue)
NL1 <- select(NL1, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue)
NL2 <- select(NL2, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue)
SL1 <- select(SL1, Date = ActivityStartDate, Location = MonitoringLocationIdentifier, Value = ResultMeasureValue)



# Concentration Trends ----------------------------------------------------

#need to rename station to the names in the report
cframe <- bind_rows(LR1, LR2, NL1, NL2, SL1)
cframe <- cframe %>%
  mutate(Location = case_when(
    Location == "TCEQMAIN-17471" ~ "LR1",
    Location == "TCEQMAIN-21489" ~ "LR2",
    Location == "TCEQMAIN-17425" ~ "NL1",
    Location == "TCEQMAIN-21548" ~ "NL2",
    Location == "TCEQMAIN-18197" ~ "SL1"
  ))
cframe$ndate <- dectime(cframe$Date)

#fit linear model for each station
lm.lr1 <- lm(log(Value) ~ ndate, data = cframe[cframe$Location=="LR1",])
plot(lm.lr1)
a <- summary(lm.lr1)
lm.lr2 <- lm(log(Value) ~ ndate, data = cframe[cframe$Location=="LR2",])
plot(lm.lr2)
summary(lm.lr2)
lm.nl1 <- lm(log(Value) ~ ndate, data = cframe[cframe$Location=="NL1",])
plot(lm.nl1)
summary(lm.nl1)
lm.nl2 <- lm(log(Value) ~ ndate, data = cframe[cframe$Location=="NL2",])
plot(lm.nl2)
summary(lm.nl2)
lm.sl1 <- lm(log(Value) ~ ndate, data = cframe[cframe$Location=="SL1",])
plot(lm.sl1)
summary(lm.sl1)

#plot linear regressions for each station
# tstat for date = coef(summary(model))[2,3]
# pvalue for date = coef(summary(model))[2,4]
# create a dataframe of t stats and p for each station
# a <- data.frame(Location = as.factor(c("LR1", "LR2", "NL1", "NL2", "SL1")), # verbose and not R-like but works
#                 pvalue = c(paste0("p=",round(coef(summary(lm.lr1))[2,4],3)), 
#                            paste0("p=",round(coef(summary(lm.lr2))[2,4],3),"*"),
#                            paste0("p=",round(coef(summary(lm.nl1))[2,4],3)),
#                            paste0("p=",round(coef(summary(lm.nl2))[2,4],3)),
#                            paste0("p=",round(coef(summary(lm.sl1))[2,4],3))),
#                 tstat = c(coef(summary(lm.lr1))[2,3], 
#                           coef(summary(lm.lr2))[2,3],
#                           coef(summary(lm.nl1))[2,3],
#                           coef(summary(lm.nl2))[2,3],
#                           coef(summary(lm.sl1))[2,3]),
#                 Date = as.Date("2014-1-1"),
#                 Value = 4500)

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
                Value = c(500,500,500,500,500))

p <- ggplot() + 
  geom_point(data = cframe, aes(Date, Value), alpha = 0.75) + 
  stat_smooth(data = cframe, aes(Date, Value), method = "lm") + 
  #geom_label(data = a, aes(Date, Value, label = paste("atop(italic(T)==",round(tstat,3),",italic(p)==",pvalue,")")), 
  #          family = "serif", size = 3, parse = T, hjust = "left",vjust = "bottom", alpha = 0.5) +
  geom_hline(yintercept = 126, linetype = 2) +
  facet_wrap(~Location) +
  scale_y_log10() +
  background_grid() +
  ylab(expression(paste(italic("E. coli"), " (cfu/100 mL)"))) +
  theme_ipsum_rc() + scale_color_ipsum()
p
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/ecoliConcentrationRC.png", width = 6.5, height = 4, units = "in")


#paste("italic(T)==",round(tstat,3),"~italic(p) ==",pvalue)

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
fframe.sum <- fframe %>%
  group_by(site_no) %>%
  summarise(min = min(Flow), max = max(Flow), mean = mean(Flow), median = median(Flow))
fframe.sum


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

#plot FDCs for each station

ggplot(filter(fframe, estflow != 0)) + 
  geom_line(aes(cumdist,estflow)) +
  scale_y_log10(labels = scales::comma) + scale_x_continuous(labels = scales::percent) + facet_wrap(~station) +
  ylab("Q (cfs)") + xlab("Duration Interval") +
  background_grid() +
  theme_ipsum_rc() +
  theme(axis.text.x = element_text(size = 9))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/fdc.png", width = 6.5, height = 4, units = "in")


# calculate and plot LDCs


fframe <- fframe %>%
  left_join(cframe, by = c("Date" = "Date", "station" = "Location"))


# allowable load = cfs * 126/100 mpn/ml * mL/cubic feet * sec/day = mpn/day  
fframe$AllowableLoad <- fframe$estflow * 126/100 * 28316.8 * 86400 / 1000000
fframe$MeasuredLoad <- fframe$estflow * fframe$Value/100 * 28316.8 * 86400 / 1000000

ggplot(filter(fframe, estflow != 0)) + 
  geom_line(aes(cumdist,AllowableLoad)) +
  geom_point(aes(cumdist,MeasuredLoad), color = "dodgerblue", alpha = 0.75) +
  scale_y_log10(labels = scales::scientific) + scale_x_continuous(labels = scales::percent) + facet_wrap(~station) +
  ylab(expression(paste(italic("E. coli"), " (million cfu/Day)"))) + 
  xlab("Duration Interval") +
  background_grid() +
  theme_ipsum_rc() +
  theme(axis.text.x = element_text(size = 9))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/ldc.png", width = 6.5, height = 4, units = "in")

## export the data
mdfLR2 <- fframe[fframe$station=="LR2",] %>%
  select(Date = Date,
         Flow = estflow,
         Value = Value) %>%
  arrange(Date)

write.csv(mdfLR2,"~/Data-Analysis-Projects/UpperLlano-WaterQuality/OutputFiles/LR2Data.csv")

mdfNL1 <- fframe[fframe$station=="NL1",] %>%
  select(Date = Date,
         Flow = estflow,
         Value = Value) %>%
  arrange(Date)
write.csv(mdfNL1,"~/Data-Analysis-Projects/UpperLlano-WaterQuality/OutputFiles/NL1Data.csv")

##### WRTDS Models
library(EGRET)
library(EGRETci)
library(fields)
library(loadflex)

plotFlowConc <- function(eList, month = c(1:12), years = NULL, col_vec = c('red', 'green', 'blue'), ylabel = NULL, xlabel = NULL, alpha = 1, size = 1,  allflo = FALSE, ncol = NULL, grids = TRUE, scales = NULL, interp = 4, pretty = TRUE, use_bw = TRUE, fac_nms = NULL, ymin = 0){
  
  localDaily <- getDaily(eList)
  localINFO <- getInfo(eList)
  localsurfaces <- getSurfaces(eList)
  
  # plot title
  toplab <- with(eList$INFO, paste(shortName, paramShortName, sep = ', '))
  
  # flow, date info for interpolation surface
  LogQ <- seq(localINFO$bottomLogQ, by=localINFO$stepLogQ, length.out=localINFO$nVectorLogQ)
  year <- seq(localINFO$bottomYear, by=localINFO$stepYear, length.out=localINFO$nVectorYear)
  jday <- 1 + round(365 * (year - floor(year)))
  surfyear <- floor(year)
  surfdts <- as.Date(paste(surfyear, jday, sep = '-'), format = '%Y-%j')
  surfmos <- as.numeric(format(surfdts, '%m'))
  surfday <- as.numeric(format(surfdts, '%d'))
  
  # interpolation surface
  ConcDay <- localsurfaces[,,3]
  
  # convert month vector to those present in data
  month <- month[month %in% surfmos]
  if(length(month) == 0) stop('No observable data for the chosen month')
  
  # salinity/flow grid values
  flo_grd <- LogQ
  
  # get the grid
  to_plo <- data.frame(date = surfdts, year = surfyear, month = surfmos, day = surfday, t(ConcDay))
  
  # reshape data frame, average by year, month for symmetry
  to_plo <- to_plo[to_plo$month %in% month, , drop = FALSE]
  names(to_plo)[grep('^X', names(to_plo))] <- paste('flo', flo_grd)
  to_plo <- tidyr::gather(to_plo, 'flo', 'res', 5:ncol(to_plo)) %>% 
    mutate(flo = as.numeric(gsub('^flo ', '', flo))) %>% 
    select(-day)
  
  # smooth the grid
  if(!is.null(interp)){
    
    to_interp <- to_plo
    to_interp <- ungroup(to_interp) %>% 
      select(date, flo, res) %>% 
      tidyr::spread(flo, res)
    
    # values to pass to interp
    dts <- to_interp$date
    fit_grd <- select(to_interp, -date)
    flo_fac <- length(flo_grd) * interp
    flo_fac <- seq(min(flo_grd), max(flo_grd), length.out = flo_fac)
    yr_fac <- seq(min(dts), max(dts), length.out = length(dts) *  interp)
    to_interp <- expand.grid(yr_fac, flo_fac)
    
    # bilinear interpolation of fit grid
    interps <- interp.surface(
      obj = list(
        y = flo_grd,
        x = dts,
        z = data.frame(fit_grd)
      ), 
      loc = to_interp
    )
    
    # format interped output
    to_plo <- data.frame(to_interp, interps) %>% 
      rename(date = Var1, 
             flo = Var2, 
             res = interps
      ) %>% 
      mutate(
        month = as.numeric(format(date, '%m')), 
        year = as.numeric(format(date, '%Y'))
      )
    
  }
  
  # subset years to plot
  if(!is.null(years)){
    
    to_plo <- to_plo[to_plo$year %in% years, ]
    to_plo <- to_plo[to_plo$month %in% month, ]
    
    if(nrow(to_plo) == 0) stop('No data to plot for the date range')
    
  }
  
  # summarize so no duplicate flos for month/yr combos
  to_plo <- group_by(to_plo, year, month, flo) %>% 
    summarize(res = mean(res, na.rm = TRUE)) %>% 
    ungroup
  
  # axis labels
  if(is.null(ylabel))
    ylabel <- localINFO$paramShortName
  if(is.null(xlabel))
    xlabel <- expression(paste('Discharge in ', m^3, '/s'))
  
  # constrain plots to salinity/flow limits for the selected month
  if(!allflo){
    
    #min, max flow values to plot
    lim_vals<- group_by(data.frame(localDaily), Month) %>% 
      summarize(
        Low = quantile(LogQ, 0.05, na.rm = TRUE),
        High = quantile(LogQ, 0.95, na.rm = TRUE)
      )
    
    # month flo ranges for plot
    lim_vals <- lim_vals[lim_vals$Month %in% month, ]
    lim_vals <- rename(lim_vals, month = Month)
    
    # merge limits with months
    to_plo <- left_join(to_plo, lim_vals, by = 'month')
    to_plo <- to_plo[to_plo$month %in% month, ]
    
    # reduce data
    sel_vec <- with(to_plo, 
                    flo >= Low &
                      flo <= High
    )
    to_plo <- to_plo[sel_vec, !names(to_plo) %in% c('Low', 'High')]
    to_plo <- arrange(to_plo, year, month)
    
  }
  
  # months labels as text
  mo_lab <- data.frame(
    num = seq(1:12), 
    txt = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
  )
  mo_lab <- mo_lab[mo_lab$num %in% month, ]
  to_plo$month <- factor(to_plo$month, levels =  mo_lab$num, labels = mo_lab$txt)
  
  # reassign facet names if fac_nms is provided
  if(!is.null(fac_nms)){
    
    if(length(fac_nms) != length(unique(to_plo$month))) 
      stop('fac_nms must have same lengths as months')
    
    to_plo$month <- factor(to_plo$month, labels = fac_nms)
    
  }
  
  # convert discharge to arithmetic scale
  to_plo$flo <- exp(to_plo$flo)*35.31467
  
  # make plot
  p <- ggplot(to_plo, aes(x = flo, y = res, group = year)) + 
    facet_wrap(~month, ncol = ncol, scales = scales)
  
  # set lower limit for y-axis if applicable
  lims <- coord_cartesian(ylim = c(ymin, max(to_plo$res, na.rm = TRUE)))
  if(!is.null(scales)){
    if(scales == 'free_x') p <- p + lims
  } else {
    p <- p + lims
  }
  
  
  # return bare bones if FALSE
  if(!pretty) return(p + geom_line())
  
  # get colors
  cols <- col_vec
  
  # use bw theme
  if(use_bw) p <- p + theme_bw()
  
  # log scale breaks
  brks <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
  
  p <- p + 
    geom_line(size = size, aes(colour = year), alpha = alpha) +
    scale_y_continuous(ylabel, expand = c(0, 0)) +
    scale_x_log10(xlabel, expand = c(0, 0), breaks = brks) +
    theme(
      legend.position = 'top', 
      axis.text.x = element_text(size = 8), 
      axis.text.y = element_text(size = 8)
    ) +
    scale_colour_gradientn('Year', colours = cols) +
    guides(colour = guide_colourbar()) + 
    ggtitle(toplab)
  
  # remove grid lines
  if(!grids) 
    p <- p + 
    theme(      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
  
}


path <- here("/OutputFiles/LR2Data.csv")
LR2 <- read_csv(path)

### Subset data to calibration data and transform concentration
LR2.calib <- LR2 %>%
  filter(!is.na(Value)) %>%
  mutate(Ecoli = as.numeric(Value),
         logEcoli = log(Ecoli))

LR2.fit <- LR2 %>%
  mutate(Ecoli = as.numeric(Value),
         logEcoli = log(Ecoli))


meta <- metadata(constituent = "Ecoli", flow = "Flow", dates = "Date", conc.units = "col/100mL", 
                 flow.units = "cfs", load.units = "million_colonies", load.rate.units = "colonies/day", 
                 site.name = "Llano River, Junction", site.id = "2")

### user flow needs two columns, date and value
### EGRET reads from user file, so export as temp csv

temp <- data_frame(Date = LR2$Date,
                   QDaily = LR2$Flow)
path <- tempdir()

write_csv(x = temp, path = paste0(path,"\\tempflow.csv"))

Daily <- readUserDaily(filePath = path, fileName = "tempflow.csv", separator = ",", interactive = F)

LR2.calib <- filter(LR2.calib, Value < 7000)
temp <- data_frame(Date = LR2.calib$Date,
                   remark = "",
                   value = LR2.calib$Value)
write_csv(x = temp, path = paste0(path,"\\tempecoli.csv"))

Sample <- readUserSample(filePath = path, fileName = "tempecoli.csv", separator = ",")

meta <- updateMetadata(metadata = meta,
                       consti.name = "E.coli")
INFO <- loadflex:::convertToEGRETInfo(meta)

eList <- mergeReport(INFO, Daily, Sample)

newFlux <- new("fluxUnit", shortName = "Mcol/day",
               unitFactor = 10,
               unitName = "Mcol/day",
               unitExpress = expression("Flux in millions of colonies / day"),
               unitExpressTiny = expression("Flux (Mcol/day)"),
               unitUSGS = "Millions of colonies per day",
               shortCode = 14)


eList <- setPA(eList, paStart = 1, paLong = 12, window = 4)
eList <- modelEstimation(eList, minNumObs = 30, minNumUncen = 0, windowY = 4, windowQ = 1, windowS = 10)

fluxBiasMulti(eList,qUnit = 1, fluxUnit = newFlux)


library(cowplot)
p <- plotFlowConc(eList, years = seq(2001,2016), scales = "fixed", alpha = 0.75, size = .4, use_bw = FALSE, xlabel = "Discharge (cfs)", ncol = 3)
p + scale_colour_viridis_c(direction = -1, option = "B") + scale_y_log10() +
  geom_hline(yintercept = 126, linetype = 2) +
  labs(y = "Concentration (cfu/100mL)",title=NULL, subtitle = "LR2 E. coli") +
  theme_ipsum_rc() +
  theme(legend.position = "right",
        panel.spacing.x = unit(.2, "in"),
        panel.spacing.y = unit(.1, "in")) +
  guides(colour = guide_colourbar("Year", barwidth = 0.5))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/LR2FlowConc.png", width = 6, height = 8, units = "in")

d <- tableResults(eList, fluxUnit = newFlux)
d

conc <- ggplot(d) +
  geom_line(aes(Year, `FN Conc [mg/L]`)) +
  geom_point(aes(Year, `Conc [mg/L]`)) +
  #geom_ribbon(data = CIAnnualResults, aes(as.numeric(Year),ymin = FNConcLow, ymax = FNConcHigh), alpha = 0.25) +
  scale_y_log10(breaks = c(1, 10, 100), limits = c(1,500)) +
  theme_ipsum_rc() +
  labs(y = "Concentration (cfu/100mL)", subtitle = "a")

fl <- ggplot(d) +
  geom_line(aes(Year, `FN Flux [Mcol/day]`)) +
  geom_point(aes(Year, `Flux [Mcol/day]`)) +
  #geom_ribbon(data = CIAnnualResults, aes(as.numeric(Year), ymin = FNFluxLow*10, ymax = FNFluxHigh*10), alpha = 0.25) +
  scale_y_log10(limits = c(.1, 1E8)) +
  theme_ipsum_rc() +
  labs(y = "Load (million colonies/day)", subtitle = "b")

path <- here("/OutputFiles/NL1Data.csv")
NL1 <- read_csv(path)

NL1.calib <- NL1 %>%
  filter(!is.na(Value)) %>%
  mutate(Ecoli = as.numeric(Value),
         logEcoli = log(Ecoli))

NL1.fit <- NL1 %>%
  mutate(Ecoli = as.numeric(Value),
         logEcoli = log(Ecoli))


meta <- metadata(constituent = "Ecoli", flow = "Flow", dates = "Date", conc.units = "col/100mL", 
                 flow.units = "cfs", load.units = "million_colonies", load.rate.units = "colonies/day", 
                 site.name = "NL, Junction", site.id = "1")

temp <- data_frame(Date = NL1$Date,
                   QDaily = NL1$Flow)
path <- tempdir()

write_csv(x = temp, path = paste0(path,"\\tempflow.csv"))

DailyNL <- readUserDaily(filePath = path, fileName = "tempflow.csv", separator = ",", interactive = F)

NL1.calib <- filter(NL1.calib, Flow >= 1 )
temp <- data_frame(Date = NL1.calib$Date,
                   remark = "",
                   value = NL1.calib$Value)
write_csv(x = temp, path = paste0(path,"\\tempecoli.csv"))

SampleNL <- readUserSample(filePath = path, fileName = "tempecoli.csv", separator = ",")

metaNL <- updateMetadata(metadata = meta,
                         consti.name = "E.coli")
INFONL <- loadflex:::convertToEGRETInfo(metaNL)
eListNL <- mergeReport(INFONL, DailyNL, SampleNL)

eListNL <- setPA(eListNL, paStart = 1, paLong = 12, window = 4)
eListNL <- modelEstimation(eListNL, minNumObs = 20, minNumUncen = 0, windowY = 4, windowQ = 1, windowS = 4)

fluxBiasMulti(eListNL,qUnit = 1, fluxUnit = newFlux)

p1 <- plotFlowConc(eListNL, years = seq(2001,2016), scales = "fixed", alpha = 0.75, size = .4, use_bw = FALSE, xlabel = "Discharge (cfs)", ncol = 3)
p1 + scale_colour_viridis_c(direction = -1, option = "B") + scale_y_log10(breaks = c(10, 100, 500)) +
  geom_hline(yintercept = 126, linetype = 2) +
  labs(x = "Discharge (cfs)", y = "Concentration (cfu/100mL)",title=NULL, subtitle = "NL1 E. coli") +
  scale_x_log10(breaks = c(.1, 1,10,100), labels = c("0.1", "1", "10", "100")) +
  theme_ipsum_rc() +
  theme(legend.position = "right",
        panel.spacing.x = unit(.2, "in"),
        panel.spacing.y = unit(.1, "in")) +
  guides(colour = guide_colourbar("Year", barwidth = 0.5))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/NL1FlowConc.png", width = 6, height = 8, units = "in")

dNL <- tableResults(eListNL, qUnit = 1, fluxUnit = newFlux)
dNL


concNL <- ggplot(dNL) +
  geom_line(aes(Year, `FN Conc [mg/L]`)) +
  geom_point(aes(Year, `Conc [mg/L]`)) +
  scale_y_log10(breaks = c(1, 10, 100), limits = c(1,500)) +
  theme_ipsum_rc() +
  labs(y = "Concentration (cfu/100mL)", subtitle = "c")

flNL <- ggplot(dNL) +
  geom_line(aes(Year, `FN Flux [Mcol/day]`)) +
  geom_point(aes(Year, `Flux [Mcol/day]`)) +
  scale_y_log10(limits = c(.1, 1E8)) +
  theme_ipsum_rc() +
  labs(y = "Load (million colonies/day)", subtitle = "d")

cowplot::plot_grid(conc, fl, concNL, flNL, rel_widths = c(1,1))
ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/WRTDSOutput.png", width = 6, height = 6, units = "in")
