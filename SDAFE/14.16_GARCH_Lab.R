
library(tidyverse)
library(ggfortify)

# 14.16 R Lab

# 14.16.1 Fitting GARCH Models

# Run the following code to load the data set TbGdpPi.csv, which has three
# variables: the 91-day T-bill rate, the log of real GDP, and the inflation rate.
# In this lab you will use only the T-bill rate.

   # set working directory
setwd("~/Dropbox/w2019/cis696/SDAFE/datasets")
TbGdpPi = read.csv("TbGdpPi.csv", header=TRUE)
   # r = the 91-day treasury bill rate
   # y = the log of real GDP
   # pi = the inflation rate
TbGdpPi = ts(TbGdpPi, start = 1955, freq = 4)
Tbill = TbGdpPi[,1]
Tbill.diff = diff(Tbill)

# Problem 1 
# Plot both Tbill and Tbill.diff. Use both time series and ACF
# plots. Also, perform ADF and KPSS tests on both series. Which series do you
# think are stationary? Why? What types of heteroskedasticity can you see in
# the Tbill.diff series?

plot(Tbill)
acf(Tbill)

plot(Tbill.diff)
acf(Tbill.diff)

library(tseries)
   # Ha: stationarity
adf.test(Tbill)         # (fail to reject H0)
adf.test(Tbill.diff)    # (reject H0)

   # H0: stationarity
kpss.test(Tbill)        # (reject H0)
kpss.test(Tbill.diff)   # (fail to reject H0)

# ANSWER
# Tbill.diff is stationary, but Tbill is not. Tbill.diff reverts to zero, wheras Tbill does not. 

# In the following code, the variable Tbill can be used if you believe that series
# is stationary. Otherwise, replace Tbill by Tbill.diff. This code will fit an
# ARMA+
#    GARCH model to the series.

library(rugarch)
arma.garch.norm = ugarchspec(mean.model=list(armaOrder=c(1,0)), 
                             variance.model=list(garchOrder=c(1,1)))
Tbill.arma.garch.norm = ugarchfit(data=Tbill, spec=arma.garch.norm)
Tbill.diff.arma.garch.norm = ugarchfit(data=Tbill.diff, spec=arma.garch.norm)
show(Tbill.arma.garch.norm)
show(Tbill.diff.arma.garch.norm)

# Problem 2 (a) Which ARMA+GARCH model is being fit? Write down the
# model using the same parameter names as in the R output.

# ANSWER
# This is an sGARCH(1,1), ARFIMA(1,0,0) model. AR(1) MA(0) / GARCH(1,1) model.

# (b) What are the estimates of each of the parameters in the model?

# ANSWER
# Parameter estimates:  mu      0.018556
#                       ar1     0.139619
#                       omega   0.016962
#                       alpha1  0.383654
#                       beta1   0.615346

# Next, plot the residuals (ordinary or raw) and standardized residuals in various
# ways using the code below. The standardized residuals are best for checking
# the model, but the residuals are useful to see if there are GARCH effects in
# the series.

res = ts(residuals(Tbill.arma.garch.norm, standardize=FALSE),
         start = 1955, freq = 4)

res.std = ts(residuals(Tbill.arma.garch.norm, standardize=TRUE),
             start = 1955, freq = 4)
par(mfrow=c(2,3))
plot(res)
acf(res)
acf(res^2)
plot(res.std)
acf(res.std)
acf(res.std^2)

# Problem 3 (a) Describe what is plotted by acf(res). What, if anything,
# does the plot tell you about the fit of the model?
#    (b) Describe what is plotted by acf(res^2). What, if anything, does the plot
# tell you about the fit of the model?
#    (c) Describe what is plotted by acf(res_std^2). What, if anything, does the
# plot tell you about the fit of the model?
#    (d) Is there anything noteworthy in the figure produced by the command
# plot(res.std)?

# ANSWER
# Not sure.

# Problem 4 Now find an ARMA+GARCH model for the series diff.log.Tbill, which we will define as 
# diff(log(Tbill)). Do you see any advantages
# of working with the differences of the logarithms of the T-bill rate, rather than
# with the difference of Tbill as was done earlier?

diff.log.Tbill <- diff(log(Tbill))

plot(diff.log.Tbill)
acf(diff.log.Tbill)           # auto-correlation looks good.
   # HA: stationary
adf.test(diff.log.Tbill)      # (reject H0)
   # H0: stationary
kpss.test(diff.log.Tbill)     # (fail to reject H0)

diff.log.Tbill.spec <- ugarchspec(mean.model=list(armaOrder=c(1,0)), 
                                  variance.model=list(garchOrder=c(1,1)))

diff.log.Tbill.fit  <- ugarchfit(data=diff.log.Tbill, spec=diff.log.Tbill.spec)

show(diff.log.Tbill.fit)

res = ts(residuals(diff.log.Tbill.fit, standardize=FALSE),
         start = 1955, freq = 4)

res.std = ts(residuals(diff.log.Tbill.fit, standardize=TRUE),
             start = 1955, freq = 4)
par(mfrow=c(2,3))
plot(res)
acf(res)
acf(res^2)
plot(res.std)
acf(res.std)
acf(res.std^2)

# 14.16.2 The GARCH-in-Mean (GARCH-M) Model

GPRO = read.table("GPRO.csv")
garchm = ugarchspec(mean.model=list(armaOrder=c(0,0),
                                    archm=T,archpow=1),
                                    variance.model=list(garchOrder=c(1,1)))
GPRO.garchm = ugarchfit(garchm, data=GPRO)
show(GPRO.garchm)

# Problem 5

fitted(GPRO.garchm)
























































































































