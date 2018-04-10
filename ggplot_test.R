
library(ggplot2)

df = read.csv("C:/Users/Telma/Desktop/test.csv")
head(df)
df_zwf1 = df[c('ZWF1_real_flux', 'ZWF1_sim_flux')]
# names(df_zwf1) = c('x','y')
df_adh3 = df[c('ADH3_real_flux', 'ADH3_sim_flux')]
# names(df_adh3) = c('x','y')


ggplot(data = df_zwf1, aes(x, y)) + 
  geom_point(color = ) +
  geom_point(data = df_adh3, colour='red') +

ggplot() + 
geom_point(data = df_adh3, aes(x = ADH3_sim_flux, y = ADH3_real_flux, color = 'ADH3')) +
geom_point(data = df_zwf1, aes(x = ZWF1_sim_flux, y = ZWF1_real_flux, color = 'ZWF1')) +
xlab('Simulated flux') +
ylab('Real flux') +
geom_abline(intercept = 0, color = 'slategrey', size = 0.5, linetype = 'dashed')


geom_line(aes(y = TempMax, colour = "TempMax"))



