source('analysis/train_models.R')
source('analysis/load_data.R')
#library(ggplot2)

train_stacked_models(indels, indel_calls)
