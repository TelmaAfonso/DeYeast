
'''
Jupyter Notebook for Case 8

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 8 Report
This report contains the results with case 8 simulations.
"""

code = """\
import warnings
from case_8 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case8 = Case8()
case8.model = case8.loadObjectFromFile('model_yeast_76.sav')
case8.model.solver = 'optlang-cplex'
case8.setMedium('MINIMAL')
case8.dictsForCase8()

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#S. Cerevisiae BY4741 deltas
genes = ['HIS3', 'LEU2', 'MET17', 'URA3']
genes = list(case8.convertStdToSyst(genes).values())

exp_dataset, reactions = case8.loadExperimentalRes('Results/Case 8/case8_experimental_fluxes.csv')
"""

fba_text = """\
## Flux Balance Analysis (FBA) Simulation
"""

pfba_text = """\
## Parsimonious Flux Balance Analysis (pFBA) Simulation
"""

lmoma_text = """\
## Linear Minimization of Metabolic Adjustment (LMOMA) Simulation
"""

fva_text = """\
## Flux Variability Analysis (FVA) Simulation
"""


# ===== GLUCOSE =====

g_text = """\
# Glucose carbon source
"""

g_fba_datasets = """\
g_exp_df = case8.getColumnWithoutNAs(exp_dataset, 0, 'X')
# O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case8.simulationPipeline(g_exp_df, cs = 'glucose', geneko = genes, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_glucose_case5.sav')
pd.concat([reactions, g_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_pfba_datasets = """\
g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case8.simulationPipeline(g_exp_df, cs = 'glucose', geneko = genes, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_glucose_case5.sav')
pd.concat([reactions, g_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_lmoma_datasets = """\
g_lmoma_res, g_lmoma_exp_sim, g_lmoma_exp_sim_errors = case8.simulationPipeline(g_exp_df, cs = 'glucose', geneko = genes, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_glucose_case5.sav')
pd.concat([reactions, g_lmoma_exp_sim_errors], axis = 1, join = 'inner')
"""

g_fva_datasets = """\
g_fva_res, g_fva_exp_sim, _ = case8.simulationPipeline(g_exp_df, cs = 'glucose', geneko = genes, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_glucose_case5.sav')
pd.concat([reactions, g_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== GALACTOSE =====

gal_text = """\
# Galactose carbon source
"""

gal_fba_datasets = """\
gal_exp_df = case8.getColumnWithoutNAs(exp_dataset, 1, 'X')
# O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

gal_fba_res, gal_fba_exp_sim, gal_fba_exp_sim_errors = case8.simulationPipeline(gal_exp_df, cs = 'galactose', geneko = genes, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_galactose_case5.sav')
pd.concat([reactions, gal_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

gal_pfba_datasets = """\
gal_pfba_res, gal_pfba_exp_sim, gal_pfba_exp_sim_errors = case8.simulationPipeline(gal_exp_df, cs = 'galactose', geneko = genes, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_galactose_case5.sav')
pd.concat([reactions, gal_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

gal_lmoma_datasets = """\
gal_lmoma_res, gal_lmoma_exp_sim, gal_lmoma_exp_sim_errors = case8.simulationPipeline(gal_exp_df, cs = 'galactose', geneko = genes, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_galactose_case5.sav')
pd.concat([reactions, gal_lmoma_exp_sim_errors], axis = 1, join = 'inner')
"""

gal_fva_datasets = """\
gal_fva_res, gal_fva_exp_sim, _ = case8.simulationPipeline(gal_exp_df, cs = 'galactose', geneko = genes, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_galactose_case5.sav')
pd.concat([reactions, gal_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== GLYCEROL =====

gly_text = """\
# Glycerol carbon source
"""

gly_fba_datasets = """\
gly_exp_df = case8.getColumnWithoutNAs(exp_dataset, 2, 'X')
# O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

gly_fba_res, gly_fba_exp_sim, gly_fba_exp_sim_errors = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', geneko = genes, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_glycerol_case5.sav')
pd.DataFrame(reactions).join(gly_fba_exp_sim_errors, how = 'inner')
"""

gly_pfba_datasets = """\
gly_pfba_res, gly_pfba_exp_sim, gly_pfba_exp_sim_errors = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', geneko = genes, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_glycerol_case5.sav')
pd.DataFrame(reactions).join(gly_pfba_exp_sim_errors, how = 'inner')
"""

gly_lmoma_datasets = """\
gly_lmoma_res, gly_lmoma_exp_sim, gly_lmoma_exp_sim_errors = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', geneko = genes, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_glycerol_case5.sav')
pd.DataFrame(reactions).join(gly_lmoma_exp_sim_errors, how = 'inner')
"""

gly_fva_datasets = """\
gly_fva_res, gly_fva_exp_sim, _ = case8.simulationPipeline(gly_exp_df, cs = 'glycerol', geneko = genes, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_glycerol_case5.sav')
pd.DataFrame(reactions).join(gly_fva_exp_sim, how = 'inner')
"""


# ===== ETHANOL =====

e_text = """\
# Ethanol carbon source
"""

e_fba_datasets = """\
e_exp_df = case8.getColumnWithoutNAs(exp_dataset, 3, 'X')

e_fba_res, e_fba_exp_sim, e_fba_exp_sim_errors = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = genes, type = 'fba', res_exists = True, fname = 'Results/Case 8/res_fba_ethanol_case5.sav')
pd.DataFrame(reactions).join(e_fba_exp_sim_errors, how = 'inner')
"""

e_pfba_datasets = """\
e_pfba_res, e_pfba_exp_sim, e_pfba_exp_sim_errors = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = genes, type = 'pfba', res_exists = True, fname = 'Results/Case 8/res_pfba_ethanol_case5.sav')
pd.DataFrame(reactions).join(e_pfba_exp_sim_errors, how = 'inner')
"""

e_lmoma_datasets = """\
e_lmoma_res, e_lmoma_exp_sim, e_lmoma_exp_sim_errors = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = genes, type = 'lmoma', res_exists = True, fname = 'Results/Case 8/res_lmoma_ethanol_case5.sav')
pd.DataFrame(reactions).join(e_lmoma_exp_sim_errors, how = 'inner')
"""

e_fva_datasets = """\
e_fva_res, e_fva_exp_sim, _ = case8.simulationPipeline(e_exp_df, cs = 'ethanol', geneko = genes, type = 'fva', res_exists = True, fname = 'Results/Case 8/res_fva_ethanol_case5.sav')
pd.DataFrame(reactions).join(e_fva_exp_sim, how = 'inner')
"""



#Generate cells with plots
x = sum([['g' + i, 'gal' + i, 'gly' + i, 'e' + i] for i in ['_fba', '_pfba', '_lmoma']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 8/'+ name +'_exp_sim_plot.png", width = 100%></p>'

#List with nbformat expressions
cs = ['g', 'gal', 'gly', 'e']
nbcells = [['nbf.v4.new_markdown_cell(' + s + '_text)',
            'nbf.v4.new_markdown_cell(fba_text)',
            'nbf.v4.new_code_cell(' + s + '_fba_datasets)',
            'nbf.v4.new_markdown_cell(' + s + '_fba_plot)',
            'nbf.v4.new_markdown_cell(pfba_text)',
            'nbf.v4.new_code_cell(' + s + '_pfba_datasets)',
            'nbf.v4.new_markdown_cell(' + s + '_pfba_plot)',
            'nbf.v4.new_markdown_cell(lmoma_text)',
            'nbf.v4.new_code_cell(' + s + '_lmoma_datasets)',
            'nbf.v4.new_markdown_cell(' + s + '_lmoma_plot)',
            'nbf.v4.new_markdown_cell(fva_text)',
            'nbf.v4.new_code_cell(' + s + '_fva_datasets)'] for s in cs]

nbcells = [item for sublist in nbcells for item in sublist]


nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code)] + [eval(exp) for exp in nbcells]



with open('case8.ipynb', 'w') as f:
    nbf.write(nb, f)


