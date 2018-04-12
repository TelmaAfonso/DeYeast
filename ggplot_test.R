
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
createGGPlotFromListOfDatasets <- function(dataset_list, title = '', nrow_legend = 25, xy_line_size = 0.5, y_lim = c(), x_lim = c()) {
  dataset_name = deparse(substitute(dataset_list))
  p <- ggplot() + 
    ggtitle(title) +
    scale_x_continuous(name = 'Simulation flux') +
    scale_y_continuous(name = 'Experimental flux') +
    geom_abline(intercept = 0, color = 'tomato3', size = xy_line_size, linetype = 'dashed') +
    labs(x = '', colour = 'Experiment') +
    guides(col = guide_legend(nrow = nrow_legend))
  
  if (!length(y_lim) == 0) {
    p = p + ylim(y_lim)
  }
  if (!length(x_lim) == 0) {
    p = p + xlim(x_lim)
  }

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
metrics_table = df


# Get list of all reactions
reactions = c()
for (name in names(pfba_datasets_proc)) {
  reactions = c(reactions, as.character(pfba_datasets_proc[[name]]$X))
}

length(unique(reactions)) #45

python_cmd = paste(unique(reactions), collapse = "', '")

r_list = read.csv('Results/reactions_list.csv', header = F)

reactions_list = list()
for (i in 1:length(r_list$V1)) {
  if (as.character(r_list$V1[i]) != 'r_4041') {
    print(paste(as.character(r_list$V1[i]), as.character(r_list$V2[i]), sep = ': ')) 
  }
  reactions_list[[as.character(r_list$V1[i])]] = as.character(r_list$V2[i])
}

length(reactions_list) #45


# ==================================
# Best method is PFBA
# Prepare data to plot by pathway
# ==================================

# Reactions ids by pathway
# Glycolysis
gly_reactions = c('r_0486', 'r_0886', 'r_0962', 'r_0534', 'r_1054', 'r_0893', 'r_0450', 'r_0366', 'r_0892')

# TCA
tca_reactions = c('r_1239', 'r_0300', 'r_0713', 'r_0718', 'r_0467', 'r_0454', 'r_0714', 'r_0832', 'r_0452', 'r_0961', 
                  'r_2034', 'r_1022', 'r_0302', 'r_0958', 'r_2131')

# PPP
ppp_reactions = c('r_1050', 'r_0982', 'r_1048', 'r_1049', 'r_0466', 'r_0889', 'r_0091', 'r_0984')

# Other (transport, cs degradation, ethanol degradation, glyoxylate cycle(2), ...)
other_reactions = c('r_0112', 'r_0153', 'r_2116', 'r_2115', 'r_2034', 'r_1254', 'r_0959', 'r_0884', 'r_0716', 'r_0662', 
                    'r_0458', 'r_0723', 'r_0164')

length(gly_reactions) + length(tca_reactions) + length(ppp_reactions) + length(other_reactions) #45


# Subset datasets by reactions
subsetByReactions <- function(datasets_list, reactions_list) {
  res = list()
  for (name in names(datasets_list)) {
    res[[name]] = subset(datasets_list[[name]], datasets_list[[name]][['X']] %in% reactions_list)
  }
  res
}

# Glycolysis
gly_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, gly_reactions)
createGGPlotFromListOfDatasets(gly_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1)

# TCA
tca_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, tca_reactions)
createGGPlotFromListOfDatasets(tca_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1, x_lim = c(NA, 27))

# PPP
ppp_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, ppp_reactions)
createGGPlotFromListOfDatasets(ppp_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1)

# Other (NOT WORKING)
other_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, other_reactions)
createGGPlotFromListOfDatasets(other_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1)


# Create a data frame with R2 and RMSE for each pathway
createMetricsDF <- function(gly_dataset_list, tca_dataset_list, ppp_dataset_list) {
  pathway_metrics = data.frame(matrix(ncol = 6, nrow = length(names(gly_dataset_list))))
  rownames(pathway_metrics) = names(gly_dataset_list)
  colnames(pathway_metrics) = c('R2_gly', 'RMSE_gly', 'R2_tca', 'RMSE_tca', 'R2_ppp', 'RMSE_ppp')
  
  for (exp in rownames(pathway_metrics)) {
    gly_dframe = gly_dataset_list[[exp]]
    tca_dframe = tca_dataset_list[[exp]]
    ppp_dframe = ppp_dataset_list[[exp]]
    exp_col = grep('exp|real|Exp|Real', colnames(gly_dframe))
    sim_col = grep('sim|Sim|fluxes', colnames(gly_dframe))
    pathway_metrics[exp, 'R2_gly'] = rsquared(gly_dframe[, exp_col], gly_dframe[, sim_col])
    pathway_metrics[exp, 'RMSE_gly'] = rmse(gly_dframe[, exp_col], gly_dframe[, sim_col])
    pathway_metrics[exp, 'R2_tca'] = rsquared(tca_dframe[, exp_col], tca_dframe[, sim_col])
    pathway_metrics[exp, 'RMSE_tca'] = rmse(tca_dframe[, exp_col], tca_dframe[, sim_col])
    pathway_metrics[exp, 'R2_ppp'] = rsquared(ppp_dframe[, exp_col], ppp_dframe[, sim_col])
    pathway_metrics[exp, 'RMSE_ppp'] = rmse(ppp_dframe[, exp_col], ppp_dframe[, sim_col])
  }
  pathway_metrics
}

pathway_metrics = createMetricsDF(gly_pfba_datasets_proc, tca_pfba_datasets_proc, ppp_pfba_datasets_proc)
# Save results
write.csv(pathway_metrics, 'pathway_metrics_table.csv')


# ==================================
# Prepare data to plot by type
# ==================================
















