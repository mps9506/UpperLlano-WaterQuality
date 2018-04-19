
# Load Packages -----------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(dataRetrieval)



# Import Data -------------------------------------------------------------

# Stations 12212, 16701, 17009, 17425, 18197, 21263, 21264, 21266, 21267, 21268, 21269, 21270, 21271, 21272,
# 21283, 21289, 21548
# 
# # Import Instantaneous Flow
# LlanoFlow <- read_delim("~/Data-Analysis-Projects/UpperLlano-WaterQuality/RawData/LlanoFlow.txt", 
#                         "|", escape_double = FALSE, trim_ws = TRUE, 
#                         col_types = cols(`Segment`= col_character(),
#                                          `Station ID`= col_character(),
#                                          `End Date` = col_date(format = "%m/%d/%Y")))
# LlanoFlow <- LlanoFlow[LlanoFlow$`Parameter Code`=="00061",]
# 
# # Import Bacteria
# LlanoBact <- read_delim("~/Data-Analysis-Projects/UpperLlano-WaterQuality/RawData/LlanoBacteria.txt",
#                         "|", escape_double = FALSE, trim_ws = TRUE,
#                         col_types = cols(`Segment`= col_character(),
#                                          `Station ID`= col_character(),
#                                          `End Date` = col_date(format = "%m/%d/%Y")))
# LlanoBact <- LlanoBact[LlanoBact$`Parameter Code` == "31699" | LlanoBact$`Parameter Code` == "31648",]

Bacteria <- readWQPqw(siteNumbers = c("TCEQMAIN-12212", "TCEQMAIN-16701", "TCEQMAIN-17009", "TCEQMAIN-17425", "TCEQMAIN-18197", 
                    "TCEQMAIN-21263", "TCEQMAIN-21264", "TCEQMAIN-21266", "TCEQMAIN-21267", "TCEQMAIN-21268", "TCEQMAIN-21269", 
                    "TCEQMAIN-21270", "TCEQMAIN-21271", "TCEQMAIN-21272", "TCEQMAIN-21283", "TCEQMAIN-21289", "TCEQMAIN-21489", 
                    "TCEQMAIN-21548"),
                    parameterCd = "Escherichia coli")
Bacteria <- Bacteria[Bacteria$ResultMeasure.MeasureUnitCode != "hours",]

Bacteria %>%
  group_by(MonitoringLocationIdentifier) %>%
  summarize(n = n())


# #Field
# LlanoField <- read_delim("~/Data-Analysis-Projects/UpperLlano-WaterQuality/RawData/LlanoField.txt", 
#                          "|", escape_double = FALSE, trim_ws = TRUE,
#                          col_types = cols(`Segment`= col_character(),
#                                           `Station ID`= col_character(),
#                                           `End Date` = col_date(format = "%m/%d/%Y")))
# LlanoGrabDO <- LlanoField[LlanoField$`Parameter Code`=="00300",]
# 
# #Nutrients
# LlanoNutrient <- read_delim("~/Data-Analysis-Projects/UpperLlano-WaterQuality/RawData/LlanoNutrient.txt", 
#                             "|", escape_double = FALSE, trim_ws = TRUE,
#                             col_types = cols(`Segment`= col_character(),
#                                              `Station ID`= col_character(),
#                                              `End Date` = col_date(format = "%m/%d/%Y")))
# LlanoChloro <- LlanoNutrient[LlanoNutrient$`Parameter Code`=="70953",]
# LlanoKjel <- LlanoNutrient[LlanoNutrient$`Parameter Code`=="00625",]
# LlanoPhos <- LlanoNutrient[LlanoNutrient$`Parameter Code`=="00665",]

# import USGS discharge data from 08149900, 08150000, and 08148500

discharge08149900 <- readNWISdv("08149900", "00060", "2000-01-01", "2017-07-30")
discharge08149900 <- renameNWISColumns(discharge08149900)
discharge08150000 <- readNWISdv("08150000", "00060", "2000-01-01", "2017-07-30")
discharge08150000 <- renameNWISColumns(discharge08150000)
discharge08148500 <- readNWISdv("08148500", "00060", "2000-01-01", "2017-07-30")
discharge08148500 <- renameNWISColumns(discharge08148500)

# Graphically explore -----------------------------------------------------


# Instantaneous Flow
ggplot(LlanoFlow) +
  geom_point(aes(x = `End Date`, y = Value, color = `Station ID`))

# Bacteria
ggplot(Bacteria) +
  geom_point(aes(x = ActivityStartDate, y = ResultMeasureValue, color = MonitoringLocationIdentifier)) + scale_y_log10()

# Grab DO
ggplot(LlanoGrabDO) +
  geom_point(aes(x = `End Date`, y = Value, color = `Station ID`))

# Chloro
ggplot(LlanoChloro) +
  geom_point(aes(x = `End Date`, y = Value, color = `Station ID`))

# Kjel
ggplot(LlanoKjel) +
  geom_point(aes(x = `End Date`, y = Value, color = `Station ID`))

# Phos
ggplot(LlanoPhos) +
  geom_point(aes(x = `End Date`, y = Value, color = `Station ID`))

# Need to decide how to handle censored data. Appears there are multiple lower detection limits.


ggplot(discharge08149900) + geom_line(aes(Date, Flow)) + scale_y_log10()
ggplot(discharge08150000) + geom_line(aes(Date, Flow)) + scale_y_log10()
ggplot(discharge08148500) + geom_line(aes(Date, Flow)) + scale_y_log10()

# Data Munging ------------------------------------------------------------

##TEMP NOTE
## CHANGE THE STATIONS USED IN THE ANALYSIS. REMOVE THE SPRINGS, GROUP THE STATIONS BY AUID
## IT APPEARS THAT INSTANANEOUS FLOW DATA WAS TAKEN ON DIFFERENT DAYS THAT SAMPLES COLLECTED BY TTU
## THEREFORE, CAN'T DO FLOW NORMALIZED REGRESSION ON DATA COLLECTED FOR THE WPP
## PLAN ON USING SIMPLE MANN-kENDALL TREND TEST ON DATA GROUPED BY AU
## FLOW NORMALIZED REGRESSION WILL BE DONE ON DATA COLLECTED BY LCRA ONLY. THIS CONSISTS OF 3-4 SITES IN TOTAL.
## DUE TO PROXIMITY TO GAUGE STATIONS, CONSIDER USING MEAN DAILY FLOWS FROM USGS GAUGES FOR THAT DATA


# join instant flow data to water quality data via dplyr

#bacteria + flow
LlanoBact <- LlanoBact %>%
  left_join(LlanoFlow, by = "RFA(Sample Set ID)/Tag_id")

ggplot(LlanoBact) +
  geom_point(aes(x = Value.y, y = Value.x, color = `Station ID.x`)) + scale_y_log10()
