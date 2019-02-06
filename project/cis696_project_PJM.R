
library(xts)
library(readr)
library(keras)
library(tseries)
library(rugarch)
library(tidyverse)
library(ggfortify)
library(PerformanceAnalytics)


## PRELIMINARIES
oldw <- getOption("warn")
start_time <- proc.time()
set.seed(1337)

   # use GMT TZ to avoid DST dupes
Sys.setenv(TZ = "GMT")

   # '2-yr time series, before and after blackout'
BEG_DT <- as.character("2002-08-01 01:00")
END_DT <- as.character("2004-08-31 24:00")
all_time_idx <- seq(from = as.POSIXct(BEG_DT),
                    to = as.POSIXct(END_DT),
                    by = "hour")

   # '1-yr before blackout' time index
BEG_DT <- as.character("2002-08-01 01:00")
END_DT <- as.character("2003-08-14 24:00")
before_time_idx <- seq(from = as.POSIXct(BEG_DT),
                       to = as.POSIXct(END_DT),
                       by = "hour")

   # '1-yr after blackout' time index
BEG_DT <- as.character("2003-08-15 01:00")
END_DT <- as.character("2004-08-31 24:00")
after_time_idx <- seq(from = as.POSIXct(BEG_DT),
                      to = as.POSIXct(END_DT),
                      by = "hour")

   # 'blackout month' (bomo) time index (August 2003)
BEG_DT <- as.character("2003-08-01 01:00")
END_DT <- as.character("2003-08-31 24:00")
bomo_time_idx <- seq(from = as.POSIXct(BEG_DT),
                     to = as.POSIXct(END_DT),
                     by = "hour")

   # '5-days during blackout' time index (08/13 - 08/15, 2003)
BEG_DT <- as.character("2003-08-12 01:00")
END_DT <- as.character("2003-08-16 24:00")
during_time_idx <- seq(from = as.POSIXct(BEG_DT),
                       to = as.POSIXct(END_DT),
                       by = "hour")

## END PRELIMINARIES


## FUNCTIONS

    # df2xts - convert ts dataframe to xts time-series
    # params: ts_df - dataframe
    #
df2xts <- function(ts_df) {
    # remove 1st three columns
    tmp.df <- ts_df[-c(1,2,3)]
    tmp.vec <- as.numeric(as.character(unlist(t(tmp.df))))
   
    if(length(tmp.vec) == 18288) { tidx = all_time_idx }        # 2-yr time series 
    else if(length(tmp.vec) == 9096) { tidx = before_time_idx } # 1-yr before blackout
    else if(length(tmp.vec) == 9192) { tidx = after_time_idx }  # 1-yr after blackout
    else if(length(tmp.vec) == 744) { tidx = bomo_time_idx }    # blackout month
    else if(length(tmp.vec) == 120) { tidx = during_time_idx }  # 5-days during bo
    
    # create time-series object
    tmp.xts <- xts(tmp.vec, order.by = tidx)
    colnames(tmp.xts)[1] <- "Price"
    
    return(tmp.xts)
}

   # ens - eight-number summary
   #  min, Q1, med, Q3, max, mean, var, sd
   #
   # param: ts_obj - time series-like object, where 'ts_obj' could be:
   #  xts   time-series
   #  tbl   tibble 
   #  df    data.frame
   #
ens <- function(ts_obj) {
   df1  <- apply(ts_obj, 2, fivenum)
   row1 <- apply(ts_obj, 2, mean, na.rm=TRUE)
   row2 <- apply(ts_obj, 2, var, na.rm=TRUE)
   row3 <- apply(ts_obj, 2, sd, na.rm=TRUE)
   df2  <- as.data.frame(rbind(row1, rbind(row2, row3)))
   
   out_df <- rbind(df1, df2)
   rownames(out_df) <- c("min", "Q1", "med", "Q3", "max", "mean", "var", "sd")
   
   return(out_df)
}

   # subset_zone - extract data for a specific zone from master datasets
   #
subset_zone <- function(rt.data, da.data, zone) {
        # subset for zone
    rt.all <- subset(rt.data, Name == zone)
    da.all <- subset(da.data, Name == zone)
        # subset before and after RT & DA data for zone
    rt.before <- subset(rt.data, Name == zone & Date <= 20030814)
    rt.after  <- subset(rt.data, Name == zone & Date  > 20030814) 
    da.before <- subset(da.data, Name == zone & Date <= 20030814)
    da.after  <- subset(da.data, Name == zone & Date  > 20030814)
        # subset 31 days centered on the blackout (August 2003)
    rt.bomo <- subset(rt.data, Name == zone & Date >= 20030801 & Date <= 20030831)
    da.bomo <- subset(da.data, Name == zone & Date >= 20030801 & Date <= 20030831)
        # subset 5 days centered on the blackout (August 12-16, 2003)
    rt.during <- subset(rt.data, Name == zone & Date >= 20030812 & Date <= 20030816)
    da.during <- subset(da.data, Name == zone & Date >= 20030812 & Date <= 20030816)
        # collect in a list
    df.list <- list(rt.all, da.all,
                    rt.before, rt.after,
                    da.before, da.after,
                    rt.bomo, da.bomo,
                    rt.during, da.during)
    
    return(df.list)
}

   # hyper-parameter search for GARCH models
   # setup vectors for hyper-parameter search
   #
hyper_param_search <- function(ds_list) { 
   results_df   <- data.frame()
   
   vm_vec <- c("sGARCH", "gjrGARCH", "csGARCH")
   dm_vec <- c("std", "sstd")
   
      # initialize dataset index
   idx = 1
      # mark rugarch loop start time
   loop_start_time <- proc.time()
   count = 1
   for (ds in ds_list) {
         # get dataset name
      dataset <- names(ds_list)[idx]
      for (vm in vm_vec) {
            # skip eGARCH for non-standardized diff datasets
         if ( (dataset == "rt_diff" || dataset == "da_diff") && vm == "eGARCH") { next }
         for (dm in dm_vec) {
            for (alpha in c(1, 2)) {
               for (beta in c(1, 2)) {
                  for (AR in c(0, 1)) {
                     for (MA in c(0, 1)) {
                        spec <- ugarchspec(
                           variance.model = list(model = vm, garchOrder = c(alpha, beta)),
                           mean.model = list(armaOrder = c(AR, MA)),
                           distribution.model = dm
                        )
                        
                        fit_start <- proc.time()
                        fit <- ugarchfit(
                           data = ds,
                           spec = spec
                        )
                        fit_elapsed <- proc.time() - fit_start
                        
                        v.model        <- fit@model[["modeldesc"]]$vmodel
                        alpha          <- fit@model[["modelinc"]][["alpha"]]
                        beta           <- fit@model[["modelinc"]][["beta"]]
                        AR             <- fit@model[["modelinc"]][["ar"]]
                        MA             <- fit@model[["modelinc"]][["ma"]]
                        distribution   <- fit@model[["modeldesc"]]$distribution
                        mIC            <- mean(infocriteria(fit))
                        log_likelihood <- likelihood(fit)
                        elapsed        <- fit_elapsed[[3]]
                        
                        row <- data.frame(dataset,
                                          v.model, alpha, beta,
                                          AR, MA, distribution, 
                                          mIC, log_likelihood,
                                          elapsed)
                        
                        results_df <- rbind(results_df, row)
                        
                        print(row)
                        print(count)
                        count <- count + 1
                     }
                  }
               }
            }
         }
      }
         # increment dataset index
      idx <- idx + 1
   }
      # print elapsed time
   print(proc.time() - loop_start_time)
   
   return(results_df)
}

## END FUNCTIONS


## DATA PREP

    # set the working directory
setwd("~/Dropbox/w2019/cis696/project/pjmdata")

    # import PJM wholesale electricy price data
    # (files generated by 'subset_m.zsh')
    # real-time and day-ahead
    # Aug 2002 - Aug 2004

   # skip first row "Start of LMP Data,,,'
blackout_rt = read_csv('blackout_rt_data.csv', skip=1, col_types="dcclldddddddddddddddddddddddd")
blackout_da = read_csv('blackout_da_data.csv', skip=1, col_types="dcclldddddddddddddddddddddddd")

    # COLUMN GUIDE
    # Date, Name, Voltage(Type): columns 01, 02, 03
    # TotalLMP - Locational Marginal Price (hourly): columns 06 - 29

    # subset data, LMP only
blackout_rt = blackout_rt[ , c(01:03, 06:29) ]
blackout_da = blackout_da[ , c(01:03, 06:29) ]

    # set meaningful column names
colnames(blackout_rt) = c("Date", "Name", "Type", 
                          "0100", "0200", "0300", "0400", "0500", 
                          "0600", "0700", "0800", "0900", "1000",
                          "1100", "1200", "1300", "1400", "1500",
                          "1600", "1700", "1800", "1900", "2000",
                          "2100", "2200", "2300", "2400")
colnames(blackout_da) = c("Date", "Name", "Type",
                          "0100", "0200", "0300", "0400", "0500",
                          "0600", "0700", "0800", "0900", "1000",
                          "1100", "1200", "1300", "1400", "1500",
                          "1600", "1700", "1800", "1900", "2000",
                          "2100", "2200", "2300", "2400")

   # rt.all, da.all, rt.before, rt.after, da.before, da.after, rt.bomo, da.bomo, rt.during, da.during
df.list <- subset_zone(blackout_rt, blackout_da, "PJM")

   # convert rt.all and da.all to xts
rt_xts <- df2xts(df.list[[1]])
da_xts <- df2xts(df.list[[2]])

   # generate log-returns
   # real-time
rt_ret <- CalculateReturns(rt_xts, method = "log")
which(is.na(rt_ret))
rt_ret[1] <- 0
rt_ret[!is.finite(rt_ret)] <- NA
rt_ret <- na.approx(rt_ret)
   # day-ahead
da_ret <- CalculateReturns(da_xts, method = "log")
which(is.na(da_ret))
da_ret[1] <- 0
da_ret[!is.finite(da_ret)] <- NA
da_ret <- na.approx(da_ret)

   # diff-log (log-returns)
   # equivalent to CalculateReturns (log returns) above
   # not needed but left in as comments for reference
# rt_diff.log <- diff(log(rt_xts))
# which(is.na(rt_diff.log))
# rt_diff.log[1] <- 0
# rt_diff.log[!is.finite(rt_diff.log)] <- NA
# rt_diff.log <- na.approx(rt_diff.log)


   # generate volatility time-series for 2 candidate fits found by hyper_param_search()
   #
   # dataset   v.model  alpha beta AR MA  dist  mIC            log_likelihood elapsed
   # RUN #2
   # da_ret	   gjrGARCH	2	   2	   1	1	sstd	-1.155528836	10593.14828	   28.991
   # rt_ret	   gjrGARCH	1	   2	   1	1	sstd	0.623543476	   -5678.3942	   18.875

   # da_ret
spec.da_ret <- ugarchspec(
   variance.model = list(model = "gjrGARCH", garchOrder = c(2,2)),
   mean.model = list(armaOrder = c(1,1)),
   distribution.model = "sstd"
)
fit.da_ret <- ugarchfit(
   data = da_ret,
   spec = spec.da_ret
)
   # rt_ret
spec.rt_ret <- ugarchspec(
   variance.model = list(model = "gjrGARCH", garchOrder = c(1,2)),
   mean.model = list(armaOrder = c(1,1)),
   distribution.model = "sstd"
)
fit.rt_ret <- ugarchfit(
   data = rt_ret,
   spec = spec.rt_ret
)

   # tests for time series stationarity
   #
   # library(tseries)
   # ADF (Augmented Dickey-Fuller) test for stationarity
   #     HA: stationary
   # KPSS (Kwiatkowski-Phillips-Schmidt-Shin) test for stationarity
   #     H0: stationary

options(warn = -1)
   # ADF/KPSS on raw RT price [obj. must be ts() ]
adf.test(as.ts(rt_xts))    # H0: non-stationary; p-val 0.01, reject H0 
adf.test(as.ts(rt_ret))    # H0: non-stationary; p-val 0.01, reject H0
kpss.test(as.ts(rt_xts))   # H0: stationary; p-val 0.01, reject H0
kpss.test(as.ts(rt_ret))   # H0: stationary; p-val 0.10, fail to reject H0
   # ADF/KPSS on raw DA price [obj. must be ts() ]
adf.test(as.ts(da_xts))    # H0: non-stationary; p-val 0.01, reject H0 
adf.test(as.ts(da_ret))    # H0: non-stationary; p-val 0.01, reject H0
kpss.test(as.ts(da_xts))   # H0: stationary; p-val 0.01, reject H0
kpss.test(as.ts(da_ret))   # H0: stationary; p-val 0.10, fail to reject H0
options(warn = oldw)

   # plot stationary price time series
   # and volatility
# autoplot(da_ret)
# autoplot(sigma(fit.da_ret))
# autoplot(rt_ret)
# autoplot(sigma(fit.rt_ret))
   # store the volatility time series
da_ret.vol_xts  <- sigma(fit.da_ret)
rt_ret.vol_xts  <- sigma(fit.rt_ret)

   # build datasets
   # real-time
realtime_xts <- rt_xts
realtime_xts$Log.Return <- rt_ret
# realtime_xts$Volatility <- c(0, sigma(fit.rt_ret))
realtime_xts$Volatility <- sigma(fit.rt_ret)
   # day-ahead
dayahead_xts <- da_xts
dayahead_xts$Log.Return <- da_ret
# dayahead_xts$Volatility <- c(0, sigma(fit.da_ret))
dayahead_xts$Volatility <- sigma(fit.da_ret)
   # datasets as tibbles
realtime.tbl <- as_tibble(realtime_xts)
dayahead.tbl <- as_tibble(dayahead_xts)

## END DATA PREP

## MAIN

   # set up zone vector. NOTE: 
   # GPU   stopped reporting 2003-04-30
   # COMED started reporting 2004-05-01
# zone_vec <- c("PJM", "PSEG", "PECO", "PPL",
#               "BGE", "JCPL", "PENELEC", "METED",
#               "PEPCO", "AECO", "DPL", "APS", "RECO")
 
   # model search results
# dataset_list <- list(rt_ret = rt_ret, da_ret = da_ret) 
# search_results.df <- hyper_param_search(dataset_list)


   # RECURRENT NEURAL NET

   # normalize training data
   # 60 %
train_prop <- 0.6    # training proportion
val_prop   <- 0.2    # validation proportion
test_prop  <- 0.2    # testing proportion

   # here switch between real-time and day-ahead
data <- data.matrix(realtime.tbl)
# data <- data.matrix(dayahead.tbl)

   # scale the data
   # find mean and sd of training data
   # see Chollet/Allaire pg. 195 &
   # https://www.listendata.com/2017/04/how-to-standardize-variable-in-regression.html
train_data <- data[1:(nrow(data) * train_prop), ]
train_mean <- apply(train_data, 2, mean)
train_sd   <- apply(train_data, 2, sd)
data <- scale(data, center=train_mean, scale=train_sd)
   # this scales both price and volatility!!
   # do we want that?

   # generator yielding time-series samples and targets
generator <- function(data,
                      lookback,
                      delay,
                      min_index,
                      max_index,
                      #shuffle = FALSE,
                      shuffle = TRUE,
                      batch_size = 24,
                      step = 3) {
   if (is.null(max_index))
      max_index <- nrow(data) - delay - 1
   i <- min_index + lookback
   function() {
      if (shuffle) {
         rows <- sample(c((min_index + lookback):max_index), size = batch_size)
      } else {
         if (i + batch_size >= max_index)
            i <<- min_index + lookback
         rows <- c(i:min(i + batch_size, max_index))
         i <<- i + length(rows)
      }
      samples <- array(0, dim = c(length(rows),
                                  lookback / step,
                                  dim(data)[[-1]]))
      targets <- array(0, dim = c(length(rows)))
      for (j in 1:length(rows)) {
         indices <- seq(rows[[j]] - lookback, rows[[j]],
                        length.out = dim(samples)[[2]])
         samples[j, , ] <- data[indices, ]
         targets[[j]] <- data[rows[[j]] + delay, 3]   # response variable
      }
      return(list(samples, targets))
   }
}

   # prepare training, validation, and test generators
lookback <- 240
step <- 3
delay <- 24
batch_size <- 24

train_min_index = 1
train_max_index = as.integer(nrow(data) * train_prop)
train_gen <- generator(
   data,
   lookback = lookback,
   delay = delay,
   min_index = train_min_index,
   max_index = train_max_index,
   shuffle = TRUE,
   step = step,
   batch_size = batch_size
)

val_min_index = train_max_index + 1
val_max_index = train_max_index + as.integer((nrow(data) * val_prop))
val_gen = generator(
   data,
   lookback = lookback,
   delay = delay,
   min_index = val_min_index,
   max_index = val_max_index,
   step = step,
   batch_size = batch_size
)

test_min_index = val_max_index + 1
test_max_index = NULL
test_gen <- generator(
   data,
   lookback = lookback,
   delay = delay,
   min_index = test_min_index,
   max_index = test_max_index,
   step = step,
   batch_size = batch_size
)

val_steps <- as.integer((val_max_index - val_min_index - lookback) / batch_size)
test_steps <- as.integer((nrow(data) - test_min_index - lookback) / batch_size)

   # common-sense baseline MSE
naive_mse <- function() {
   batch_mses <- c()
   for (step in 1:val_steps) {
      c(samples, targets) %<-% val_gen()
      preds <-
         samples[, dim(samples)[[2]], 3]     # 3rd dimension is response
      mse <- mean((preds - targets) ^ 2)
      batch_mses <- c(batch_mses, mse)
   }
   print(paste("MSE: ", mean(batch_mses)))
}
naive_mse()

elapsed_time <- proc.time() - start_time
print(elapsed_time)

   # dense model
model <- keras_model_sequential() %>%
   layer_flatten(input_shape = c(lookback / step, dim(data)[-1])) %>%
   layer_dense(units = 32, activation = "relu") %>%
   layer_dense(units = 1)
model %>% compile(optimizer = optimizer_rmsprop(),
                  loss = "mse")
history <- model %>% fit_generator(
   train_gen,
   steps_per_epoch = 500,
   epochs = 20,
   validation_data = val_gen,
   validation_steps = val_steps
)
## END MAIN