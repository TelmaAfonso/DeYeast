
'''
Jupyter Notebook for Case 6

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook() #create a new notebook object

text = """\
# Case 6 Report
This report contains the results with case 6 simulations.
"""

code = """\
import warnings
from case_6 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case6 = Case6()
case6.model = case6.loadObjectFromFile('model_yeast_76.sav')
case6.model.solver = 'optlang-cplex'
case6.setMedium('MINIMAL')

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
exp_dataset, reactions = case6.loadExperimentalRes('Results/Case 6/case6_experimental_fluxes.csv')
pd.DataFrame(reactions)
# NO EXPERIMENTAL ETHANOL FLUXES TO ADJUST O2 LB
"""

wt_text = """\
# Wild Type
"""

wt_fba_text = """\
## Flux Balance Analysis (FBA) Simulation
"""

wt_fba_datasets = """\
wt_fba_res, wt_fba_exp_sim, wt_fba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'fba', res_exists = True)
pd.concat([reactions, wt_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

wt_pfba_text = """\
## Parsimonious Flux Balance Analysis (pFBA) Simulation
"""

wt_pfba_datasets = """\
wt_pfba_res, wt_pfba_exp_sim, wt_pfba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'pfba', res_exists = True)
pd.concat([reactions, wt_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

wt_fva_text = """\
## Flux Variability Analysis (FVA) Simulation
"""

wt_fva_datasets = """\
wt_fva_res, wt_fva_exp_sim, _ = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'fva', res_exists = True)
pd.concat([reactions, wt_fva_exp_sim], axis = 1, join = 'inner')
"""


mae1_text = """\
# MAE1 Deletion
"""

mae1_fba_datasets = """\
mae1_fba_res, mae1_fba_exp_sim, mae1_fba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'fba', res_exists = True)
pd.concat([reactions, mae1_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

mae1_pfba_datasets = """\
mae1_pfba_res, mae1_pfba_exp_sim, mae1_pfba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'pfba', res_exists = True)
pd.concat([reactions, mae1_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

mae1_fva_datasets = """\
mae1_fva_res, mae1_fva_exp_sim, _ = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'fva', res_exists = True)
pd.concat([reactions, mae1_fva_exp_sim], axis = 1, join = 'inner')
"""

mae1_lmoma_text = """\
## Linear Minimization of Metabolic Adjustment (LMOMA) Simulation
"""

mae1_lmoma_datasets = """\
mae1_lmoma_res, mae1_lmoma_exp_sim, mae1_lmoma_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'lmoma', res_exists = True)
pd.concat([reactions, mae1_lmoma_exp_sim_errors], axis = 1, join = 'inner')
"""


#Generate cells with plots
html_figs = {}
l = ['wt_fba', 'wt_pfba', 'mae1_fba', 'mae1_pfba', 'c']
for name in l:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 6/'+ name +'_exp_sim_plot.png", width = 100%></p>'

html_figs.keys()




nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code),

               nbf.v4.new_markdown_cell(wt_text),
               nbf.v4.new_markdown_cell(wt_fba_text),
               nbf.v4.new_code_cell(wt_fba_datasets),
               nbf.v4.new_markdown_cell(wt_fba_plot),

               nbf.v4.new_markdown_cell(wt_pfba_text),
               nbf.v4.new_code_cell(wt_pfba_datasets),
               nbf.v4.new_markdown_cell(wt_pfba_plot),

               nbf.v4.new_markdown_cell(wt_fva_text),
               nbf.v4.new_code_cell(wt_fva_datasets),

               nbf.v4.new_markdown_cell(mae1_text),
               nbf.v4.new_markdown_cell(wt_fba_text),
               nbf.v4.new_code_cell(mae1_fba_datasets),
               nbf.v4.new_markdown_cell(mae1_fba_plot),

               nbf.v4.new_markdown_cell(wt_pfba_text),
               nbf.v4.new_code_cell(mae1_pfba_datasets),
               nbf.v4.new_markdown_cell(mae1_pfba_plot),

               nbf.v4.new_markdown_cell(wt_fva_text),
               nbf.v4.new_code_cell(mae1_fva_datasets),

               nbf.v4.new_markdown_cell(mae1_lmoma_text),
               nbf.v4.new_code_cell(mae1_lmoma_datasets),
               nbf.v4.new_markdown_cell(mae1_lmoma_plot)
               ]


with open('case6.ipynb', 'w') as f:
    nbf.write(nb, f)



