
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
createGGPlotFromListOfDatasets <- function(dataset_list, title = '', nrow_legend = 25, xy_line_size = 0.5, y_lim = c(), x_lim = c(), 
                                           legend_pos = c(1, 0), legend_alpha = 0.5, abline_color = 'tomato3') {
  dataset_name = deparse(substitute(dataset_list))
  p <- ggplot() + 
    ggtitle(title) +
    # scale_x_continuous(name = 'Simulation flux') +
    # ylab('Experimental flux') +
    geom_abline(intercept = 0, color = abline_color, size = xy_line_size, linetype = 'dashed') +
    labs(x = '', colour = 'Experiment') +
    guides(col = guide_legend(nrow = nrow_legend)) 
  
  if (!length(y_lim) == 0) {
    p = p + scale_y_continuous(name = 'Experimental flux', limits = y_lim)
  } else {
    p = p + scale_y_continuous(name = 'Experimental flux')
  }
  if (!length(x_lim) == 0) {
    p = p + scale_x_continuous(name = 'Simulation flux', limits = x_lim)
  } else {
    p = p + scale_x_continuous(name = 'Simulation flux')
  }
  if (!length(legend_pos) == 0) {
    p = p + theme(legend.justification = c(1, 0), 
                  legend.position = legend_pos,
                  legend.background = element_rect(fill = alpha('white', legend_alpha))) 
  } 

  for (name in names(dataset_list)) {
    print(name)
    if (!length(rownames(dataset_list[[name]])) == 0) {
      df = dataset_list[[name]]
      exp_col = grep('exp|real|Exp|Real', colnames(df))
      sim_col = grep('sim|Sim|fluxes', colnames(df))
      
      exp_data = paste0('geom_point(data = ', dataset_name, '$', name, ', aes(x = ', dataset_name, '$', name,
                        '[,', sim_col, '], y = ', dataset_name, '$', name, '[,', exp_col, '], color = "', name, '"))')
      exp_label = paste0('geom_text(data = ', dataset_name, '$', name, ', aes(x = ', dataset_name, '$', name,
                        '[,', sim_col, '], y = ', dataset_name, '$', name, '[,', exp_col, '], label = "', name, '"), 
                        size = 3, check_overlap = TRUE, hjust = 0, nudge_x = 0.05)')
      
      p = p + eval(parse(text = exp_data))
      p = p + eval(parse(text = exp_label))
    } else {
      print(paste0('Dataset ', name, ' is of length zero, not used.'))
    }
  }
  p
}

library(ggrepel)

createSingleGGPlotFromDataset <- function(dataset, geompoint_color = NULL, title = '', nrow_legend = 2, xy_line_size = 1, 
                                          y_lim = c(), x_lim = c(), legend_pos = c(1, 0), legend_alpha = 0.5, abline_color = 'lightblue4',
                                          label_size = 4, axis_label_size = 13, repel_level = 5, min_seg_length = 0.5,
                                          error_thr = NULL, abs_err = T, point_padding = 1, arrow_color = 'gray60',
                                          point_size = 2, subplot = F, subYLim = c(), subXLim = c(), alpha_legend = F, arrows = T) {
  exp_col = grep('exp|real|Exp|Real', colnames(dataset))
  sim_col = grep('sim|Sim|fluxes', colnames(dataset))
  dataset_name = deparse(substitute(dataset))
  
  # Set labels according to chosen metric
  labels = character(length(rownames(dataset)))
  if (!is.null(error_thr)) {
    if (abs_err) {
      Scale = absError(dataset, sim_col, exp_col)
      if (subplot) {
        labels[which(Scale > error_thr & dataset[, sim_col] > subXLim[1] & dataset[, sim_col] < subXLim[2], arr.ind = FALSE, useNames = TRUE)] = as.character(dataset$X)[which(Scale > error_thr & dataset[, sim_col] > subXLim[1] & dataset[, sim_col] < subXLim[2], arr.ind = FALSE, useNames = TRUE)]
      } else {
        labels[which(Scale > error_thr, arr.ind = FALSE, useNames = TRUE)] = as.character(dataset$X)[which(Scale > error_thr, arr.ind = FALSE, useNames = TRUE)]
      }
    } else {
      Scale = relError(dataset, sim_col, exp_col)
      if (subplot) {
        labels[which(Scale > error_thr & dataset[, sim_col] > subXLim[1] & dataset[, sim_col] < subXLim[2], arr.ind = FALSE, useNames = TRUE)] = as.character(dataset$X)[which(Scale > error_thr & dataset[, sim_col] > subXLim[1] & dataset[, sim_col] < subXLim[2], arr.ind = FALSE, useNames = TRUE)]
      } else {
        labels[which(Scale > error_thr, arr.ind = FALSE, useNames = TRUE)] = as.character(dataset$X)[which(Scale > error_thr, arr.ind = FALSE, useNames = TRUE)]
      }
    }
  }
  
  # Start plot
  p <- ggplot() + 
    ggtitle(title) +
    geom_abline(intercept = 0, color = abline_color, size = xy_line_size, linetype = 'dashed') +
    labs(x = '', colour = 'Experiment') +
    guides(col = guide_legend(nrow = nrow_legend))
    
  # Show alpha legend
  if (!alpha_legend) p = p + scale_alpha_continuous(guide = F)
  
  # If it's a subplot
  if (subplot) {
    p <- p + 
      theme_void(base_size = axis_label_size, base_family = "") +
      coord_cartesian(xlim = subXLim, ylim = subYLim)
  } else {
    p <- p + theme_classic(base_size = axis_label_size, base_family = "")
  }
  
  # Level of alpha to plot
  if (is.null(error_thr)) {
    p <- p + geom_point(data = dataset, size = point_size, aes(x = dataset[, sim_col], y = dataset[, exp_col], color = dataset_name, alpha = 0.5))
  } else {
    p <- p + geom_point(data = dataset, size = point_size, aes(x = dataset[, sim_col], y = dataset[, exp_col], color = dataset_name, alpha = Scale))
  }
  
  # Add labels
  if (arrows) {
    p <- p + geom_text_repel(data = dataset, size = label_size, arrow = arrow(length = unit(0.01, 'npc')),
                             point.padding = unit(point_padding, 'lines'), segment.color = arrow_color, force = repel_level, min.segment.length = min_seg_length,
                             aes(x = dataset[, sim_col], y = dataset[, exp_col], label = labels, fontface = "bold"))
  } else {
    p <- p + geom_text_repel(data = dataset, size = label_size, segment.color = NA,
                             point.padding = unit(point_padding, 'lines'),  force = repel_level,
                             aes(x = dataset[, sim_col], y = dataset[, exp_col], label = labels, fontface = "bold"))
  }
  
  # Custom axes limits
  if (!length(y_lim) == 0) {
    p = p + scale_y_continuous(name = 'Experimental flux', limits = y_lim)
  } else {
    p = p + scale_y_continuous(name = 'Experimental flux')
  }
  if (!length(x_lim) == 0) {
    p = p + scale_x_continuous(name = 'Simulation flux', limits = x_lim)
  } else {
    p = p + scale_x_continuous(name = 'Simulation flux')
  }
  
  # Legend customization
  if (!length(legend_pos) == 0) {
    p = p + theme(axis.line = element_line(colour = 'black', size = 1),
                  axis.ticks = element_line(colour = 'black', size = 1),
                  legend.justification = c(1, 0), 
                  legend.position = legend_pos,
                  legend.background = element_rect(fill = alpha('white', legend_alpha))) 
  } 
  
  # Set point colors
  if (!is.null(geompoint_color)) {
    color_expr = paste0('scale_color_manual(values = c(', dataset_name, ' = geompoint_color))')
    p = p + eval(parse(text = color_expr))
  }
  p
}

#d13_ferm_fba = fba_datasets_proc[['d13_ferm_fba']]
createSingleGGPlotFromDataset(d13_ferm_fba, abs_err = T, error_thr = 1, geompoint_color = 'deepskyblue4')
createSingleGGPlotFromDataset(d13_ferm_fba, abs_err = T, error_thr = 1, geompoint_color = 'deepskyblue4', subplot = T, 
                              subXLim = c(-0.2, 1), subYLim = c(-2, 2.1), legend_pos = c(), arrows = F, point_padding = 0.05, label_size = 3.5, alpha_legend = T)


df = fba_datasets_proc[['d13_ferm_fba']]
exp_col = grep('exp|real|Exp|Real', colnames(df))
sim_col = grep('sim|Sim|fluxes', colnames(df))
labels = paste0(as.character(df$X), ' (d13_ferm_fba)')

library(ggrepel)
ggplot() + 
  ggtitle('') +
  scale_color_manual(values = c('d13_ferm_fba' = 'deepskyblue4')) +
  theme_classic(base_size = 13, base_family = "") +
  # theme_void(base_size = 13, base_family = "") +
  scale_y_continuous(name = 'Experimental flux') +
  scale_x_continuous(name = 'Simulation flux') + 
  geom_abline(intercept = 0, color = 'lightblue4', size = 1, linetype = 'dashed') +
  labs(x = '', colour = 'Experiment') +
  guides(col = guide_legend(nrow = 2))  +
  theme(legend.justification = c(1, 0), 
        # panel.border = element_rect(colour = 'black', fill = NA, size = 1.5),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(1, 0),
        legend.background = element_rect(fill = alpha('white', 0.5))) +
  # coord_cartesian(xlim = c(-0.5, 1), ylim = c(-2.5, 5)) +
  scale_alpha_continuous(guide = F) +
  geom_point(data = fba_datasets_proc$d13_ferm_fba, size = 2, show.legend = T, #alpha = 0.5,
             aes(x = fba_datasets_proc$d13_ferm_fba[, 3], y = fba_datasets_proc$d13_ferm_fba[, 2], color = "d13_ferm_fba",
                 alpha = df2)) +
  geom_text_repel(data = fba_datasets_proc$d13_ferm_fba, size = 4, arrow = arrow(length = unit(0.01, 'npc')), #nudge_x = 1, nudge_x = 1,#angle = 0,
                  point.padding = unit(1, 'lines'), segment.color = 'gray60', force = 5, min.segment.length = 0.5, #direction = 'x',
                  aes(x = fba_datasets_proc$d13_ferm_fba[, 3], y = fba_datasets_proc$d13_ferm_fba[, 2], label = x, fontface = "bold"))

xlims = c(-0.5, 1)
df2 = absError(df, sim_col, exp_col)
df3 = relError(df, sim_col, exp_col)
x = character(length(df2))
alpha_vector = rep(0.25, length(df2))
alpha_vector[which(df2 > 1.5, arr.ind = FALSE, useNames = TRUE)] = 0.5
x[which(df2 > 1.5, arr.ind = FALSE, useNames = TRUE)] = as.character(df$X)[which(df2 > 1.5, arr.ind = FALSE, useNames = TRUE)]
x[which(abs(df3) > 2, arr.ind = FALSE, useNames = TRUE)] = as.character(df$X)[which(abs(df3) > 2, arr.ind = FALSE, useNames = TRUE)]
x

# labels between limits
x = character(length(df2))
x[which(df2 > 1.5 & df$Sim.Flux > xlims[1] & df$Sim.Flux < xlims[2], arr.ind = FALSE, useNames = TRUE)] = as.character(df$X)[which(df2 > 1.5 & df$Sim.Flux > xlims[1] & df$Sim.Flux < xlims[2], arr.ind = FALSE, useNames = TRUE)]
x
# labels outside limits
x = character(length(df2))
x[which(df2 > 1.5 & (df$Sim.Flux < xlims[1] | df$Sim.Flux > xlims[2]), arr.ind = FALSE, useNames = TRUE)] = as.character(df$X)[which(df2 > 1.5 & (df$Sim.Flux < xlims[1] | df$Sim.Flux > xlims[2]), arr.ind = FALSE, useNames = TRUE)]
x

#Any old plot
a_plot

#A viewport taking up a fraction of the plot area
# library(grid)
# vp <- viewport(width = 0.3, height = 0.3, x = 0.26, y = 0.82)
# print(a_plot)
# print(b_plot, vp = vp)

metrics_table


absError <- function (dataset, sim_col, exp_col) {
  res = abs(dataset[, sim_col] - dataset[, exp_col])
}

relError <- function (dataset, sim_col, exp_col) {
  res = abs(dataset[, sim_col] - dataset[, exp_col]) / dataset[, sim_col]
}

ggplot(nba, aes(x= MIN, y= PTS, colour="green", label=Name))+
  geom_point() +
  geom_text(aes(label=ifelse(PTS>24,as.character(Name),'')),hjust=0,vjust=0)


# GLOBAL DATASETS

# FBA datasets
fba_datasets_proc = removeOutliers(fba_datasets, 500)
createGGPlotFromListOfDatasets(fba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(fba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4', legend_pos = c())
# createGGPlotFromListOfDatasets(fba_datasets_proc, xy_line_size = 1.2)

# pFBA datasets
pfba_datasets_proc = removeOutliers(pfba_datasets, 500)
createGGPlotFromListOfDatasets(pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(pfba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4', legend_pos = c())
# createGGPlotFromListOfDatasets(pfba_datasets_proc, xy_line_size = 1.2)

# LMOMA datasets
lmoma_datasets_proc = removeOutliers(lmoma_datasets, 500)
createGGPlotFromListOfDatasets(lmoma_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(lmoma_datasets_proc))/3),
                               legend_alpha = 0, abline_color = 'lightblue4')
# createGGPlotFromListOfDatasets(lmoma_datasets_proc, nrow_legend = 17, xy_line_size = 1.2)



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
createGGPlotFromListOfDatasets(gly_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(gly_pfba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4', legend_pos = c())
# createGGPlotFromListOfDatasets(gly_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1)

# TCA
tca_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, tca_reactions)
createGGPlotFromListOfDatasets(tca_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(tca_pfba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4', legend_pos = c(), x_lim = c(NA, 27))
# createGGPlotFromListOfDatasets(tca_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1, x_lim = c(NA, 27))

# PPP
ppp_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, ppp_reactions)
createGGPlotFromListOfDatasets(ppp_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(ppp_pfba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4', legend_pos = c())
# createGGPlotFromListOfDatasets(ppp_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1)

# Other (datasets d6_mae1_pfba and d6_wt_pfba don't have any of these reactions)
other_pfba_datasets_proc = subsetByReactions(pfba_datasets_proc, other_reactions)
createGGPlotFromListOfDatasets(other_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(other_pfba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4', legend_pos = c(), x_lim = c(-3, NA))
# createGGPlotFromListOfDatasets(other_pfba_datasets_proc, title = '', nrow_legend = 25, xy_line_size = 1)


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

# ==================================
# Experiments by carbon source
# ==================================
glucose = c('d13_oxi_pfba', 'd13_rferm_pfba', 'd13_ferm_pfba', 'd14_wt_pfba', 'd3_glu_pfba', 'd5_glu_pfba', 'd6_mae1_pfba',
            'd6_wt_pfba', 'd8_glu_pfba', 'd9_batch_pfba', 'd9_chem_pfba', 'd9_hxk2_pfba')

mannose = c('d3_man_pfba')
galactose = c('d3_gal_pfba', 'd8_gal_pfba')
pyruvate = c('d3_pyr_pfba')
ethanol = c('d5_etoh_pfba', 'd8_etoh_pfba')
glycerol = c('d8_gly_pfba')

# Add case 7 experiments to glucose cases
d7 = names(gly_pfba_datasets_proc)[grepl('d7', names(gly_pfba_datasets_proc))]
for (exp in d7) {
  glucose = c(glucose, exp)
}

# Subset datasets by experimets
subsetByExperiments <- function(datasets_list, experiments_list) {
  res = list()
  for (name in experiments_list) {
    res[[name]] = datasets_list[[name]]
  }
  res
}

# Glucose
glucose_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, glucose)
createGGPlotFromListOfDatasets(glucose_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(glucose_pfba_datasets_proc))/3),
                               legend_alpha = 0, abline_color = 'lightblue4')
# Mannose
mannose_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, mannose)
createGGPlotFromListOfDatasets(mannose_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = 1, abline_color = 'lightblue4',
                               legend_alpha = 0)
# Galactose
galactose_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, galactose)
createGGPlotFromListOfDatasets(galactose_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = 1, abline_color = 'lightblue4',
                               legend_alpha = 0)
# Pyruvate
pyruvate_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, pyruvate)
createGGPlotFromListOfDatasets(pyruvate_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = 1, abline_color = 'lightblue4',
                               legend_alpha = 0, y_lim = c(-1.4, NA))
# Ethanol
ethanol_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, ethanol)
createGGPlotFromListOfDatasets(ethanol_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = 1, abline_color = 'lightblue4',
                               legend_alpha = 0, x_lim = c(-5, NA))
# Glycerol
glycerol_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, glycerol)
createGGPlotFromListOfDatasets(glycerol_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = 1, abline_color = 'lightblue4',
                               legend_alpha = 0, y_lim = c(-2.3, NA))


# ==================================
# Experiments by mode of operation
# ==================================
batch = c('d3_glu_pfba', 'd3_man_pfba', 'd3_gal_pfba', 'd3_pyr_pfba', 'd8_glu_pfba', 'd8_gal_pfba', 'd8_gly_pfba', 'd8_etoh_pfba',
          'd9_batch_pfba', 'd9_chem_pfba', 'd9_hxk2_pfba', 'd5_glu_pfba', 'd5_etoh_pfba')

chemostat = c('d6_mae1_pfba', 'd6_wt_pfba', 'd9_chem_pfba', 'd13_oxi_pfba', 'd13_rferm_pfba', 'd13_ferm_pfba')

# Add case 7 experiments to chemostat cases
for (exp in d7) {
  chemostat = c(chemostat, exp)
}

# Batch
batch_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, batch)
createGGPlotFromListOfDatasets(batch_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(batch_pfba_datasets_proc))/3),
                               legend_alpha = 0, abline_color = 'lightblue4')

# Chemostat
chemostat_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, chemostat)
createGGPlotFromListOfDatasets(chemostat_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(chemostat_pfba_datasets_proc))/3 + 1),
                               legend_alpha = 0, abline_color = 'lightblue4')




# ==================================
# Experiments by knockout type
# ==================================

case7_genes = c('SOL4', 'PGM1', 'PCK1', 'ZWF1', 'IDP1', 'MAE1', 'ALD6', 'ICL1', 'SOL3', 'OAC1', 'MDH1', 'GCV2', 
                'SFC1', 'CTP1', 'ADH3', 'COX5A', 'LSC1', 'GPD1', 'ALD5', 'IDP2', 'RPE1', 'SER33', 'PDA1', 'SDH1', 
                'DAL7', 'GND2', 'GLY1', 'MDH3', 'TAL1', 'PGM2', 'FUM1')

length(case7_genes) #31

ko_tca = c('d7_idp1_pfba', 'd7_mae1_pfba', 'd7_mdh1_pfba', 'd7_lsc1_pfba', 'd7_idp2_pfba', 'd7_pda1_pfba',
           'd7_sdh1_pfba', 'd7_fum1_pfba', 'd6_mae1_pfba')
ko_gly = c('d14_hxk2_pfba')
ko_ppp = c('d7_sol4_pfba', 'd7_zwf1_pfba', 'd7_sol3_pfba', 'd7_rpe1_pfba', 'd7_gnd2_pfba', 'd7_tal1_pfba')

ko_other = c('d7_pck1_pfba', 'd7_ald6_pfba', 'd7_icl1_pfba', 'd7_oac1_pfba', 'd7_pgm1_pfba', 'd7_gcv2_pfba',
             'd7_sfc1_pfba', 'd7_ctp1_pfba', 'd7_adh3_pfba', 'd7_cox5a_pfba', 'd7_gpd1_pfba', 'd7_ald5_pfba',
             'd7_ser33_pfba', 'd7_dal7_pfba', 'd7_gly1_pfba', 'd7_mdh3_pfba', 'd7_pgm2_pfba')
# d7_pck1_pfba - gluconeogenesis, d7_icl1_pfba - glyoxylate, d7_oac1_pfba - oxaloacetate transporter
# d7_pgm1_pfba - glycogen biosynthesis, d7_gcv2_pfba - glycine cleavage complex, d7_sfc1_pfba - succinate-fumarate transporter
# d7_ctp1_pfba - citrate transporter, d7_adh3_pfba - ethanol degradation, d7_cox5a_pfba - electron transport chain
# d7_gpd1_pfba - glycerol biosynthesis, d7_ald5_pfba - mitochondrial aldehyde dehydrogenase, d7_ser33_pfba - serine biosynthesis
# d7_dal7_pfba - glyoxylate cycle, d7_gly1_pfba - glycine biosynthesis, d7_mdh3_pfba - glyoxylate cycle
# d7_pgm2_pfba - glycogen biosynthesis

length(ko_tca) + length(ko_gly) + length(ko_ppp) + length(ko_other) # 33

# TCA Knockouts
tca_ko_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, ko_tca)
createGGPlotFromListOfDatasets(tca_ko_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(tca_ko_pfba_datasets_proc))/2 + 0.5),
                               legend_alpha = 0, abline_color = 'lightblue4')
# Glycolysis Knockouts
gly_ko_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, ko_gly)
createGGPlotFromListOfDatasets(gly_ko_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = 1,
                               legend_alpha = 0, abline_color = 'lightblue4')
# PPP Knockouts
ppp_ko_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, ko_ppp)
createGGPlotFromListOfDatasets(ppp_ko_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(ppp_ko_pfba_datasets_proc))/2),
                               legend_alpha = 0, abline_color = 'lightblue4')
# Other Knockouts
other_ko_pfba_datasets_proc = subsetByExperiments(pfba_datasets_proc, ko_other)
createGGPlotFromListOfDatasets(other_ko_pfba_datasets_proc, title = '', xy_line_size = 1,
                               nrow_legend = round(length(names(other_ko_pfba_datasets_proc))/2 + 0.5),
                               legend_alpha = 0, abline_color = 'lightblue4')


# Convert decimal places in csv files
write.csv(format(metrics_table, digits = 4), 'C:/Users/Telma/Desktop/DeYeast Paper/Files/global_metrics.csv')
write.csv(format(pathway_metrics, digits = 4), 'C:/Users/Telma/Desktop/DeYeast Paper/Files/pathway_metrics.csv')

format(metrics_table, digits = 4)



