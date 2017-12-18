
'''
Jupyter Notebook for Case 9

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 9 Report
This report contains the results with case 9 simulations.
"""

code = """\
import warnings
from case_9 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case9 = Case9()
case9.model = case9.loadObjectFromFile('model_yeast_76.sav')
case9.model.solver = 'optlang-cplex'
case9.setMedium('MINIMAL')
case9.dictsForCase9()

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#General datasets
exp_dataset, reactions = case9.loadExperimentalRes('Results/Case 9/case9_experimental_fluxes.csv')
pd.DataFrame(reactions)
"""

fba_text = """\
## Flux Balance Analysis (FBA) Simulation
"""

pfba_text = """\
## Parsimonious Flux Balance Analysis (pFBA) Simulation
"""

fva_text = """\
## Flux Variability Analysis (FVA) Simulation
"""


# ===== BATCH =====

b_text = """\
# Batch culture
"""

b_fba_datasets = """\
b_exp_df = case9.getColumnWithoutNAs(exp_dataset, 0, 'X')
# O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

b_fba_res, b_fba_exp_sim, b_fba_exp_sim_errors = case9.simulationPipeline(b_exp_df, cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 9/res_fba_batch_case9.sav')
pd.concat([reactions, b_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

b_pfba_datasets = """\
b_pfba_res, b_pfba_exp_sim, b_pfba_exp_sim_errors = case9.simulationPipeline(b_exp_df, cs = 'glucose',type = 'pfba', res_exists = True, fname = 'Results/Case 9/res_pfba_batch_case9.sav')
pd.concat([reactions, b_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

b_fva_datasets = """\
b_fva_res, b_fva_exp_sim, _ = case9.simulationPipeline(b_exp_df, cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 9/res_fva_batch_case9.sav')
pd.concat([reactions, b_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== CHEMOSTAT =====

c_text = """\
# Chemostat culture
"""

c_fba_datasets = """\
c_exp_df = case9.getColumnWithoutNAs(exp_dataset, 1, 'X')
# O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

c_fba_res, c_fba_exp_sim, c_fba_exp_sim_errors = case9.simulationPipeline(c_exp_df, cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 9/res_fba_chemostat_case9.sav')
pd.concat([reactions, c_fba_exp_sim_errors], axis = 1, join = 'inner') #Plot not showing r_0302
"""

c_pfba_datasets = """\
c_pfba_res, c_pfba_exp_sim, c_pfba_exp_sim_errors = case9.simulationPipeline(c_exp_df, cs = 'glucose',type = 'pfba', res_exists = True, fname = 'Results/Case 9/res_pfba_chemostat_case9.sav')
pd.concat([reactions, c_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

c_fva_datasets = """\
c_fva_res, c_fva_exp_sim, _ = case9.simulationPipeline(c_exp_df, cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 9/res_fva_chemostat_case9.sav')
pd.concat([reactions, c_fva_exp_sim], axis = 1, join = 'inner')
"""


#Generate cells with plots
x = sum([['b' + i, 'c' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 9/'+ name +'_exp_sim_plot.png", width = 100%></p>'

#List with nbformat expressions
cultures = ['b', 'c']
nbcells = [['nbf.v4.new_markdown_cell(' + cult + '_text)',
            'nbf.v4.new_markdown_cell(fba_text)',
            'nbf.v4.new_code_cell(' + cult + '_fba_datasets)',
            'nbf.v4.new_markdown_cell(' + cult + '_fba_plot)',
            'nbf.v4.new_markdown_cell(pfba_text)',
            'nbf.v4.new_code_cell(' + cult + '_pfba_datasets)',
            'nbf.v4.new_markdown_cell(' + cult + '_pfba_plot)',
            'nbf.v4.new_markdown_cell(fva_text)',
            'nbf.v4.new_code_cell(' + cult + '_fva_datasets)'] for cult in cultures]

nbcells = [item for sublist in nbcells for item in sublist]


nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code)] + [eval(exp) for exp in nbcells]



with open('case9.ipynb', 'w') as f:
    nbf.write(nb, f)


