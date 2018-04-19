# notes ==================================================================
# This script downloads water quality data and
# evaluates changes in concentrations and loads
# in the North Llano, South Llano, and the Confluence


# North Llano SWQM Stations: 21268, 21267, 21266, 21264, 21263,
# 21283 (USGS08150000), 21548, 17425
# South llano SWQM stations: 16701, 18197, 21271, 21269, 21270,
# USGS08149900
# Confluence: 21489, USGS08148500, 17471

# import libraries ========================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(dataRetrieval)
library(cowplot)
library(lubridate)
library(smwrBase)

# Custom Functions


# rank =====================================================================

# take a dataframe of daily mean streamflow values and generate the flow
# duration curve

# input should be a dataframe with
# column 0 = date in yyyy-mm-dd format
# column 2 = daily mean flow in cfs
#' Flow duration curve
#'
#' Calculates the flow exceedance probabilities from an input dataframe with daily mean streamflow values.

rank <- function(data, filter, ...){
  
  if(filter == TRUE) { #if filter = TRUE, calculate the cum dist, then remove the zeros. Handy for plotting
    d <- data %>%
      arrange(.[[4]]) %>%
      mutate(cumdist = cume_dist(.[[4]])) %>%
      filter(.[[4]] > 0)
  } else { #filter == FALSE, don't remove zero's, needed for calculating the flow exceedance probabilities
    d <- data %>%
      arrange(.[[4]]) %>%
      mutate(cumdist = cume_dist(.[[4]]))
  }
  
  d
  
}


# phi =======================================================================

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

# wtr_yr =========================================================================

# function to return water year

wtr_yr <- function(dates, start_month=9) {
  # Convert dates into POSIXlt
  dates.posix = as.POSIXlt(dates)
  # Year offset
  offset = ifelse(dates.posix$mon >= start_month - 1, 1, 0)
  # Water year
  adj.year = dates.posix$year + 1900 + offset
  # Return the water year
  adj.year
}


#function to return day of water year

wtr_day <- function(dates, start_month=9) {
  dates.posix = as.POSIXlt(dates)
  offset = ifelse(dates.posix$mon >= start_month -1, 1, 0)
  
}


#############################################
## Functions for derivatives of GAM models ##
#############################################
Deriv <- function(mod, n = 200, eps = 1e-7, newdata) {
  if(isTRUE(all.equal(class(mod), "list")))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  # number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    ## Xi <- Xp * 0 ##matrix(0, nrow = Xp.r, ncol = Xp.c)
    ## J <- bs.dims[i]
    ## Xi[,(i-1) * J + 1:J + 1] <- Xp[,(i-1) * J + 1:J +1]
    ## df <- Xi %*% coef(mod)
    ## df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    ## lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  return(lD)
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  ##term <- term[match(term, term.labs)]
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  ## if(is.na(term))
  ##     stop("'term' not a valid model term.")
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- length(object$gamModel$y) - sum(object$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  ## tVal <- qt(1 - (alpha/2), object$gamModel$df.residual)
  for(i in seq_along(term)) {
    upr <- object[[term[i]]]$deriv + tVal * object[[term[i]]]$se.deriv
    lwr <- object[[term[i]]]$deriv - tVal * object[[term[i]]]$se.deriv
    res[[term[i]]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(is.na(Term)))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  ## tVal <- qt(1 - (alpha/2), x$gamModel$df.residual)
  residual.df <- length(x$gamModel$y) - sum(x$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")[Term]
    names(xlab) <- xlab
  }
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  CI <- confint(x, term = term, alpha = alpha)
  for(i in seq_along(term)) {
    ## for(i in seq_len(l)) {
    upr <- CI[[term[i]]]$upper
    lwr <- CI[[term[i]]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,term[i]], x[[term[i]]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[term[i]], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,term[i]], rev(x$eval[,term[i]])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,term[i]], upr, lty = "dashed")
      lines(x$eval[,term[i]], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 1)
      S <- signifD(x[[term[i]]]$deriv, x[[term[i]]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,term[i]], S$incr, lwd = lwd, col = "#ca0020")
      lines(x$eval[,term[i]], S$decr, lwd = lwd, col = "#0571b0")
    } else {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}

# Import Data ====================================================================

# Step 1: Use the DAR approach to estimate flows at respective stations based on USGS stations within respective assessment units
# Approach for Texas Streams documented in Asquith, Roussel, and Vrabel (2006)

# import USGS discharge data from 08149900, 08150000, and 08148500

discharge08149900 <- readNWISdv("08149900", "00060", "2000-01-01", "2017-07-30")
discharge08149900 <- renameNWISColumns(discharge08149900)
discharge08150000 <- readNWISdv("08150000", "00060", "2000-01-01", "2017-07-30")
discharge08150000 <- renameNWISColumns(discharge08150000)
discharge08148500 <- readNWISdv("08148500", "00060", "2000-01-01", "2017-07-30")
discharge08148500 <- renameNWISColumns(discharge08148500)


#rank streamflows
discharge08148500 <- rank(discharge08148500, filter = FALSE)
discharge08149900 <- rank(discharge08149900, filter = FALSE)
discharge08150000 <- rank(discharge08150000, filter = FALSE)



# exponent (phi) based on Asquith, Roussel, and Vrabel (2006)

discharge08148500 <- phi(discharge08148500)
discharge08149900 <- phi(discharge08149900)
discharge08150000 <- phi(discharge08150000)


## make a new column for flows at each station, then convert to long (narrow) format
discharge08148500$usgs_08148500 <- discharge08148500$Flow
discharge08148500$`TCEQMAIN-21489` <- discharge08148500$Flow * (1849/1854)^discharge08148500$phi
discharge08148500$`TCEQMAIN-17471` <- discharge08148500$Flow * (1856/1854)^discharge08148500$phi

discharge08149900$usgs_08149900 <- discharge08149900$Flow
discharge08149900$`TCEQMAIN-21270` <- discharge08149900$Flow * (868/879)^discharge08149900$phi
discharge08149900$`TCEQMAIN-21269` <- discharge08149900$Flow * (812/879)^discharge08149900$phi
discharge08149900$`TCEQMAIN-21271` <- discharge08149900$Flow * (783/879)^discharge08149900$phi
discharge08149900$`TCEQMAIN-18197` <- discharge08149900$Flow * (750/879)^discharge08149900$phi
discharge08149900$`TCEQMAIN-16701` <- discharge08149900$Flow * (502/879)^discharge08149900$phi

discharge08150000$usgs_08150000 <- discharge08150000$Flow
discharge08150000$`TCEQMAIN-17425` <- discharge08150000$Flow * (916/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21548` <- discharge08150000$Flow * (904/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21283` <- discharge08150000$Flow * (900/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21263` <- discharge08150000$Flow * (701/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21264` <- discharge08150000$Flow * (641/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21266` <- discharge08150000$Flow * (437/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21267` <- discharge08150000$Flow * (422/914)^discharge08150000$phi
discharge08150000$`TCEQMAIN-21268` <- discharge08150000$Flow * (416/914)^discharge08150000$phi

tidy08148500 <- discharge08148500 %>%
  select(-agency_cd, -site_no, -Flow, -Flow_cd) %>%
  gather(key = Station, value = Flow, usgs_08148500:`TCEQMAIN-17471`)

tidy08149900 <- discharge08149900 %>%
  select(-agency_cd, -site_no, -Flow, -Flow_cd) %>%
  gather(key = Station, value = Flow, usgs_08149900:`TCEQMAIN-16701`)

tidy08150000 <- discharge08150000 %>%
  select(-agency_cd, -site_no, -Flow, -Flow_cd) %>%
  gather(key = Station, value = Flow, usgs_08150000:`TCEQMAIN-21268`)

# import bacteria water quality monitoring data and join to above dataframe
Bacteria <- readWQPqw(siteNumbers = c("TCEQMAIN-17471", "TCEQMAIN-12212", "TCEQMAIN-16701", "TCEQMAIN-17009", "TCEQMAIN-17425", "TCEQMAIN-18197", 
                                      "TCEQMAIN-21263", "TCEQMAIN-21264", "TCEQMAIN-21266", "TCEQMAIN-21267", "TCEQMAIN-21268", "TCEQMAIN-21269", 
                                      "TCEQMAIN-21270", "TCEQMAIN-21271", "TCEQMAIN-21272", "TCEQMAIN-21283", "TCEQMAIN-21289", "TCEQMAIN-21489", 
                                      "TCEQMAIN-21548"),
                      parameterCd = "Escherichia coli")
Bacteria <- Bacteria[Bacteria$ResultMeasure.MeasureUnitCode != "hours",]

#check
Bacteria %>%
  group_by(MonitoringLocationIdentifier) %>%
  summarize(n = n())
Bacteria <- Bacteria %>%
  select(ActivityStartDate, MonitoringLocationIdentifier, ResultMeasureValue)

tidy08148500 <- tidy08148500 %>%
  left_join(Bacteria, by = c("Date" = "ActivityStartDate", "Station" = "MonitoringLocationIdentifier"))

#since sampling data goes back further than USGS data, need to use a full join, but keep out stations not in the South llano
tidy08149900 <- tidy08149900 %>%
  full_join(Bacteria[Bacteria$MonitoringLocationIdentifier %in% unique(tidy08149900$Station),], 
            by = c("Date" = "ActivityStartDate", "Station" = "MonitoringLocationIdentifier"))

tidy08150000 <- tidy08150000 %>%
  left_join(Bacteria, by = c("Date" = "ActivityStartDate", "Station" = "MonitoringLocationIdentifier"))

# confluence
ecoliA <- ggplot(tidy08148500) + geom_point(aes(Date, ResultMeasureValue), color = "dodgerblue", alpha = 0.8) + 
  geom_smooth(aes(Date, ResultMeasureValue), method ="lm") + 
  scale_y_log10() +
  theme_bw() +
  ylab("E. coli (MPN/100mL)") + ggtitle("Llano River Confluence")
ecoliA
# South
ecoliB <- ggplot(tidy08149900) + geom_point(aes(Date, ResultMeasureValue), color = "dodgerblue", alpha = 0.8) + 
  geom_smooth(aes(Date, ResultMeasureValue), method ="lm") + 
  scale_y_log10() +
  theme_bw() +
  ylab("E. coli (MPN/100mL)") + ggtitle("South Llano River")
ecoliB
# north
ecoliC <- ggplot(tidy08150000) + geom_point(aes(Date, ResultMeasureValue), color = "dodgerblue", alpha = 0.8) + 
  geom_smooth(aes(Date, ResultMeasureValue), method ="lm") + 
  scale_y_log10() +
  theme_bw() +
  ylab("E. coli (MPN/100mL)") + ggtitle("North Llano River")
ecoliC

plot_grid(ecoliA, ecoliB, ecoliC, labels = "AUTO")

#ggsave("~/Data-Analysis-Projects/UpperLlano-WaterQuality/Figures/ecoliConcentration.png", width = 6.5, height = 4, units = "in")

summary(lm(log(tidy08150000$ResultMeasureValue) ~ tidy08150000$Date))
summary(lm(log(tidy08149900$ResultMeasureValue) ~ tidy08149900$Date))
summary(lm(log(tidy08148500$ResultMeasureValue) ~ tidy08148500$Date))

