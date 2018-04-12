
'''
Functions to create summary graphs with the data

Author: Telma Afonso
'''

from yeastpack_test import PhenomenalySim
import pickle
import pandas as pd
import numpy as np
from ggplot import *

class createGraphs (PhenomenalySim):

    def __int__(self, cobra_model = None):
        super(createGraphs, self).__int__(self, cobra_model)
        # self.cases = [3, 5, 6, 8, 9, 13, 14]

    def loadResults(self, case_list):
        for i in case_list:
            res = self.loadObjectFromFile('Results/case{}_all_res.sav'.format(str(i)))
            globals()['case' + str(i)] = res

    def printDictsKeys(self, case_list):
        for i in case_list:
            assert 'case' + str(i) in globals(), 'case' + str(i) + " variable doesn't exist"
            print('\n' + '='*10 + str(list(globals()['case' + str(i)].keys())[0]).split(sep = '_')[0].upper() + '='*10)
            print(globals()['case' + str(i)].keys())

    def makeCase7DictKeysVariables(self, case7_dict, addToDict = True):
        df_fba = case7['d7_fba']
        df_pfba = case7['d7_pfba']
        df_fva = case7['d7_fva']
        df_lmoma = case7['d7_lmoma']

        for g in case7_dict['d7_genes']:
            if addToDict:
                case7_dict['d7_' + g.lower() + '_fba'] = df_fba.filter(regex = g)
                case7_dict['d7_' + g.lower() + '_pfba'] = df_pfba.filter(regex = g)
                case7_dict['d7_' + g.lower() + '_fva'] = df_fva.filter(regex = g)
                case7_dict['d7_' + g.lower() + '_lmoma'] = df_lmoma.filter(regex = g)
            else:
                # Make global
                globals()['d7_' + g.lower() + '_fba'] = df_fba.filter(regex = g)
                globals()['d7_' + g.lower() + '_pfba'] = df_pfba.filter(regex = g)
                globals()['d7_' + g.lower() + '_fva'] = df_fva.filter(regex = g)
                globals()['d7_' + g.lower() + '_lmoma'] = df_lmoma.filter(regex = g)

    def saveDatasetsToCSVFile(self, case_list, output_folder):
        for i in case_list:
            assert 'case' + str(i) in globals(), 'case' + str(i) + " variable doesn't exist"
            dict = globals()['case' + str(i)]
            for k in dict.keys():
                print(k)
                if isinstance(dict[k], pd.core.frame.DataFrame):
                    dict[k].to_csv(path_or_buf  = output_folder + k + '.csv')
                elif isinstance(dict[k], pd.core.series.Series):
                    dict[k].to_csv(path  = output_folder + k + '.csv')
                else:
                    print(k + ': Dictionary entry will not be exported.')


sorted(globals().keys())

cases = [3, 5, 6, 7, 8, 9, 13, 14]
cg = createGraphs()
cg.loadResults(cases) # Creates variables caseX (dicts) which contain information of each case
cg.printDictsKeys(cases)
cg.makeCase7DictKeysVariables(case7)
cg.saveDatasetsToCSVFile(cases, output_folder = 'Results/Datasets/')

case7.keys()
case7['d7_genes'].keys()
case7['d7_adh3_pfba'].columns
case7['d7_genes'].__class__

case3.keys()
case3['d3_reactions'].to_csv(path  = "C:/Users/Telma/Desktop/test2.csv")

isinstance(case3['d3_reactions'], pd.core.series.Series)
isinstance(case3['d3_gal_pfba'], pd.core.frame.DataFrame)


# Save file with the list of reactions
all_reactions = ['r_0534', 'r_0962', 'r_0961', 'r_2131', 'r_0713', 'r_0300', 'r_1022', 'r_0454', 'r_1054',
                 'r_0486', 'r_0450', 'r_0452', 'r_0302', 'r_0892', 'r_0893', 'r_0889', 'r_1049', 'r_1048',
                 'r_1050', 'r_0466', 'r_0467', 'r_0886', 'r_0959', 'r_2116', 'r_0112', 'r_0366', 'r_0832',
                 'r_0091', 'r_0718', 'r_0958', 'r_2034', 'r_1239', 'r_4041', 'r_2115', 'r_0458', 'r_0723',
                 'r_1254', 'r_0163', 'r_0716', 'r_0714', 'r_0662', 'r_0884', 'r_0984', 'r_0982', 'r_0164']

from case_8 import Case8
test = Case8()
test.model = test.loadObjectFromFile('model_yeast_76.sav')
test.checkReaction('r_1239')
test.getReactionInfo('r_0884')

reactions = {}
for r in all_reactions:
    reactions[r] = test.checkReaction(r)

import csv

with open('Results/reactions_list.csv', 'w') as f:
    w = csv.writer(f)
    w.writerows(reactions.items())




# GGPLOT TESTS
# p = ggplot(aes(x = 'ADH3_sim_flux', y = 'ADH3_real_flux'), data = case7['d7_adh3_pfba']) \
#     + geom_point(color = 'steelblue') \
#     + ggtitle('Beautiful graph') \
#     + xlab('Simulated flux') \
#     + ylab('Experimental flux') \
#     + theme_gray() \
#     + geom_abline(intercept = 0, color = 'slategrey', size = 1, linetype = 'dashed') # x = y line
# p
#
#
# p = ggplot(aes(x = 'ADH3_sim_flux', y = 'ADH3_real_flux'), data = case7['d7_pfba']) + geom_point(color = 'steelblue') \
#     + geom_point(aes(x = 'ZWF1_sim_flux', y = 'ZWF1_real_flux'), data = case7['d7_pfba'], color = 'tomato') \
#     + ggtitle('Beautiful graph') \
#     + xlab('Simulated flux') \
#     + ylab('Experimental flux') \
#     + theme_gray() \
#     + geom_abline(intercept = 0, color = 'slategrey', size = 1, linetype = 'dashed') # x = y line
# p
# # p.save("C:/Users/Telma/Desktop/example.png")
#
# p = ggplot(aes(x = 'ADH3_sim_flux', y = 'ADH3_real_flux', color = 'ADH3'), data = case7['d7_adh3_pfba']) \
# + geom_point(aes(x = 'ZWF1_sim_flux', y = 'ADH3_real_flux', color = 'ZWF1'), data = case7['d7_zwf1_pfba']) \
# + xlab('Simulated flux') \
# + ylab('Real flux')
# p
#
# df = case7['d7_pfba']
# df.__class__
# df.to_csv(path_or_buf  = "C:/Users/Telma/Desktop/test.csv")
#
# # + theme_bw() \
#
