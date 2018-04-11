
library(ggplot2)

# =========================
#       GGPLOT2 Tests
# =========================
# df = read.csv("C:/Users/Telma/Desktop/test.csv")
# head(df)
# df_zwf1 = df[c('ZWF1_real_flux', 'ZWF1_sim_flux')]
# df_adh3 = df[c('ADH3_real_flux', 'ADH3_sim_flux')]
# 
# 
# ggplot(data = df_zwf1, aes(x, y)) + 
#   geom_point(color = ) +
#   geom_point(data = df_adh3, colour='red') +
# 
# ggplot() + 
# geom_point(data = df_adh3, aes(x = ADH3_sim_flux, y = ADH3_real_flux, color = 'ADH3')) +
# geom_point(data = df_zwf1, aes(x = ZWF1_sim_flux, y = ZWF1_real_flux, color = 'ZWF1')) +
# xlab('Simulated flux') +
# ylab('Real flux') +
# geom_abline(intercept = 0, color = 'slategrey', size = 0.5, linetype = 'dashed')



# p <- ggplot() + 
#   xlab('Simulated flux') +
#   ylab('Real flux') +
#   geom_abline(intercept = 0, color = 'slategrey', size = 0.5, linetype = 'dashed')
# 
# p = p + geom_point(data = fba_datasets$d7_adh3_fba, aes(x = eval(parse(text = 'fba_datasets$d7_adh3_fba[,3]')), y = fba_datasets$d7_adh3_fba[,2], color = 'ADH3'))
# # p = p + geom_text(label = as.character(fba_datasets$d7_adh3_fba$X))
# # p = p + geom_text(label = rownames(fba_datasets$d7_adh3_fba), nudge_x = 0.25, nudge_y = 0.25, check_overlap = T)
# p = p + labs(x = '', colour = 'Experiment')
# p
# 
# p = p + geom_point(data = fba_datasets$d7_zwf1_fba, aes(x = fba_datasets$d7_zwf1_fba[,3], y = fba_datasets$d7_zwf1_fba[,2], color = 'ZWF1'))
# p = p + geom_point(data = fba_datasets$d13_ferm_fba, aes(x = fba_datasets$d13_ferm_fba[,3], y = fba_datasets$d13_ferm_fba[,2], color = 'FERM'))
# p = p + geom_point(data = fba_datasets$d13_oxi_fba, aes(x = fba_datasets$d13_oxi_fba[,3], y = fba_datasets$d13_oxi_fba[,2], color = 'OXI'))
# p = p + geom_point(data = fba_datasets$d13_rferm_fba, aes(x = fba_datasets$d13_rferm_fba[,3], y = fba_datasets$d13_rferm_fba[,2], color = 'RFERM'))
# p = p + geom_point(data = fba_datasets$d14_hxk2_fba, aes(x = fba_datasets$d14_hxk2_fba[,3], y = fba_datasets$d14_hxk2_fba[,2], color = 'HXK2'))
# p
# 

# ======================
#    Prepare the data
# ======================
fnames = c()
fba_files = c()
pfba_files = c()
lmoma_files = c()

# Get a list with the names of the files containing the data
for (file in list.files('Results/Datasets/')) {
  if (!grepl('reactions', file)) {
    if (grepl('_fba', file) & !grepl('d7_fba', file)) {
      fba_files = c(fba_files, strsplit(file, split = '.csv')[[1]])
    } else if (grepl('pfba', file) & !grepl('d7_pfba', file)) {
      pfba_files = c(pfba_files, strsplit(file, split = '.csv')[[1]])
    } else if (grepl('lmoma', file) & !grepl('d7_lmoma', file)) {
      lmoma_files = c(lmoma_files, strsplit(file, split = '.csv')[[1]])
    }
    fnames = c(fnames, strsplit(file, split = '.csv')[[1]])
  }
}

# Create a list with datasets from files containing the data
createListOfDatasets <- function(fname_list) {
  res = list()
  for (file in fname_list) {
    df = read.csv(paste0('Results/Datasets/', file, '.csv'))
    res[[file]] = df
  }
  res
}


fba_datasets = createListOfDatasets(fba_files)
pfba_datasets = createListOfDatasets(pfba_files)
lmoma_datasets = createListOfDatasets(lmoma_files)

# Remove fluxes that are too high (1000)
removeOutliers <- function(dataset_list, threshold) {
  res = list()
  for (name in names(dataset_list)) {
    df = dataset_list[[name]]
    sim_col = grep('sim|Sim|fluxes', colnames(df))
    df = subset(df, df[, sim_col] < threshold)
    res[[name]] = df
  }
  res
}


# ==========================
#   Create plots (ggplot2)
# ==========================
library(ggplot2)

# Create a ggplot from a list of datasets
createGGPlotFromListOfDatasets <- function(dataset_list, title = '', nrow_legend = 25, xy_line_size = 0.5) {
  dataset_name = deparse(substitute(dataset_list))
  p <- ggplot() + 
    ggtitle(title) +
    scale_x_continuous(name = 'Simulation flux') +
    scale_y_continuous(name = 'Experimental flux') +
    geom_abline(intercept = 0, color = 'tomato3', size = xy_line_size, linetype = 'dashed') +
    labs(x = '', colour = 'Experiment') +
    guides(col = guide_legend(nrow = nrow_legend))

  for (name in names(dataset_list)) {
    print(name)
    
    df = dataset_list[[name]]
    exp_col = grep('exp|real|Exp|Real', colnames(df))
    sim_col = grep('sim|Sim|fluxes', colnames(df))
  
    exp_data = paste0('geom_point(data = ', dataset_name, '$', name, ', aes(x = ', dataset_name, '$', name, '[,', sim_col, '], y = ', dataset_name, '$', name, '[,', exp_col, '], color = "', name, '"))')
    
    p = p + eval(parse(text = exp_data))
  }
  p
}

# GLOBAL DATASETS

# FBA datasets
fba_datasets_proc = removeOutliers(fba_datasets, 500)
createGGPlotFromListOfDatasets(fba_datasets_proc, xy_line_size = 1.2)

# pFBA datasets
pfba_datasets_proc = removeOutliers(pfba_datasets, 500)
createGGPlotFromListOfDatasets(pfba_datasets_proc, xy_line_size = 1.2)

# LMOMA datasets
lmoma_datasets_proc = removeOutliers(lmoma_datasets, 500)
createGGPlotFromListOfDatasets(lmoma_datasets_proc, nrow_legend = 17, xy_line_size = 1.2)



# =======================
#     Global Metrics
# =======================

# Root Mean Squared Error
rmse <- function(actual, predicted) {
  error = actual - predicted
  sqrt(mean(error ^ 2))
}

# Mean Absolute Error
mae <- function(actual, predicted) {
  error = actual - predicted
  mean(abs(error))
}

# R Squared
rsquared <- function (actual, predicted) {
  cor(actual, predicted) ^ 2
}


# Create Dataframe w/ R2 and RMSE for each experiment
names = c()
for (name in names(fba_datasets_proc)) {
  print(name)
  names = c(names, strsplit(name, split = '_fba')[[1]])
}
df = data.frame(matrix(ncol = 6, nrow = length(names)))
rownames(df) = names
colnames(df) = c('R2_fba', 'RMSE_fba', 'R2_pfba', 'RMSE_pfba', 'R2_lmoma', 'RMSE_lmoma')

# Add FBA metrics
for (exp in rownames(df)) {
  dframe = fba_datasets_proc[[paste0(exp, '_fba')]]
  exp_col = grep('exp|real|Exp|Real', colnames(dframe))
  sim_col = grep('sim|Sim|fluxes', colnames(dframe))
  df[exp, 'R2_fba'] = rsquared(dframe[, exp_col], dframe[, sim_col])
  df[exp, 'RMSE_fba'] = rmse(dframe[, exp_col], dframe[, sim_col])
}

# Add pFBA metrics
for (exp in rownames(df)) {
  dframe = pfba_datasets_proc[[paste0(exp, '_pfba')]]
  exp_col = grep('exp|real|Exp|Real', colnames(dframe))
  sim_col = grep('sim|Sim|fluxes', colnames(dframe))
  df[exp, 'R2_pfba'] = rsquared(dframe[, exp_col], dframe[, sim_col])
  df[exp, 'RMSE_pfba'] = rmse(dframe[, exp_col], dframe[, sim_col])
}

# Add LMOMA metrics
for (exp in rownames(df)) {
  if (!is.null(lmoma_datasets_proc[[paste0(exp, '_lmoma')]])) {
    dframe = lmoma_datasets_proc[[paste0(exp, '_lmoma')]]
    exp_col = grep('exp|real|Exp|Real', colnames(dframe))
    sim_col = grep('sim|Sim|fluxes', colnames(dframe))
    df[exp, 'R2_lmoma'] = rsquared(dframe[, exp_col], dframe[, sim_col])
    df[exp, 'RMSE_lmoma'] = rmse(dframe[, exp_col], dframe[, sim_col]) 
  }
}

# Save results
write.csv(df, 'metrics_table.csv')



names(pfba_datasets_proc)
pfba_datasets_proc[['d13_ferm_pfba']]$X

reactions


