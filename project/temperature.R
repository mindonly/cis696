# temperature forecasting example
# Deep Learning with R

# download and uncompress the data
dir.create("~/Dropbox/w2019/cis696/project/jena_climate", recursive = TRUE)
download.file(
   "https://s3.amazonaws.com/keras-datasets/jena_climate_2009_2016.csv.zip",
   "~/Dropbox/w2019/cis696/project/jena_climate/jena_climate_2009_2016.csv.zip"
)
unzip(
   "~/Dropbox/w2019/cis696/project/jena_climate/jena_climate_2009_2016.csv.zip",
   exdir = "~/Dropbox/w2019/cis696/project/jena_climate"
)

# 6.28 inspect the data
library(tibble)
library(readr)

data_dir <- "~/Dropbox/w2019/cis696/project/jena_climate"
fname <- file.path(data_dir, "jena_climate_2009_2016.csv")
data <- read_csv(fname)

glimpse(data)

# 6.29 plot the timeseries
library(ggplot2)

ggplot(data, aes(x = 1:nrow(data), y = `T (degC)`)) + geom_line()

# 6.30 plot the first 10 days
ggplot(data[1:1440,], aes(x = 1:1440, y = `T (degC)`)) + geom_line()

# Listing 6.31 Converting the data into a floating-point matrix
data <- data.matrix(data[,-1])

# Listing 6.32 Normalizing the data
train_data <- data[1:200000,]
mean <- apply(train_data, 2, mean)
std <- apply(train_data, 2, sd)
data <- scale(data, center = mean, scale = std)

# Listing 6.33 Generator yielding timeseries samples and their targets
generator <- function(data, lookback, delay, min_index, max_index,
                      shuffle = FALSE, batch_size = 128, step = 6) {
   if (is.null(max_index)) max_index <- nrow(data) - delay - 1
   i <- min_index + lookback
   function() {
      if (shuffle) {
         rows <- sample(c((min_index+lookback):max_index), size = batch_size)
      } else {
         if (i + batch_size >= max_index)
            i <<- min_index + lookback
         rows <- c(i:min(i+batch_size, max_index))
         i <<- i + length(rows)
      }
      samples <- array(0, dim = c(length(rows),
                                  lookback / step,
                                  dim(data)[[-1]]))
      targets <- array(0, dim = c(length(rows)))
      for (j in 1:length(rows)) {
         indices <- seq(rows[[j]] - lookback, rows[[j]],
                        length.out = dim(samples)[[2]])
         samples[j,,] <- data[indices,]
         targets[[j]] <- data[rows[[j]] + delay,2]
      }
      list(samples, targets)
   }
}

# Listing 6.34 Preparing the training, validation, and test generators
library(keras)

lookback <- 1440
step <- 6
delay <- 144
batch_size <- 128
train_gen <- generator(
   data,
   lookback = lookback,
   delay = delay,
   min_index = 1,
   max_index = 200000,
   shuffle = TRUE,
   step = step,
   batch_size = batch_size
)
val_gen = generator(
   data,
   lookback = lookback,
   delay = delay,
   min_index = 200001,
   max_index = 300000,
   step = step,
   batch_size = batch_size
)
test_gen <- generator(
   data,
   lookback = lookback,
   delay = delay,
   min_index = 300001,
   max_index = NULL,
   step = step,
   batch_size = batch_size
)
val_steps <- (300000 - 200001 - lookback) / batch_size
test_steps <- (nrow(data) - 300001 - lookback) / batch_size

# Listing 6.35 Computing the common-sense baseline MAE
evaluate_naive_method <- function() {
   batch_maes <- c()
   for (step in 1:val_steps) {
      c(samples, targets) %<-% val_gen()
      preds <- samples[,dim(samples)[[2]],2]
      mae <- mean(abs(preds - targets))
      batch_maes <- c(batch_maes, mae)
   }
   print(mean(batch_maes))
}
evaluate_naive_method()

# Listing 6.36 Converting the MAE back to a Celsius error
celsius_mae <- 0.29 * std[[2]]

# Listing 6.37 Training and evaluating a densely connected model
library(keras)

model <- keras_model_sequential() %>%
   layer_flatten(input_shape = c(lookback / step, dim(data)[-1])) %>%
   layer_dense(units = 32, activation = "relu") %>%
   layer_dense(units = 1)
model %>% compile(
   optimizer = optimizer_rmsprop(),
   loss = "mae"
)
history <- model %>% fit_generator(
   train_gen,
   steps_per_epoch = 500,
   epochs = 20,
   validation_data = val_gen,
   validation_steps = val_steps
)