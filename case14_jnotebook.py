
'''
Jupyter Notebook for Case 14

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 14 Report (Raghevendran et al, 2004)

This report contains the results with case 14 simulations.

Paper: [*Phenotypic characterization of glucose repression mutants of Saccharomyces cerevisiae using experiments with 13 C-labelled glucose*](http://onlinelibrary.wiley.com/doi/10.1002/yea.1136/abstract)

**Abstract**
In the field of metabolic engineering and functional genomics, methods for analysis of
metabolic fluxes in the cell are attractive as they give an overview of the phenotypic
response of the cells at the level of the active metabolic network. This is unlike several
other high-throughput experimental techniques, which do not provide information
about the integrated response a specific genetic modification has on the cellular
function. In this study we have performed phenotypic characterization of several
mutants of the yeast Saccharomyces cerevisiae through the use of experiments with
13
C-labelled glucose. Through GC–MS analysis of the 13 C incorporated into the amino
acids of cellular proteins, it was possible to obtain quantitative information on the
function of the central carbon metabolism in the different mutants. Traditionally,
such labelling data have been used to quantify metabolic fluxes through the use of a
suitable mathematical model, but here we show that the raw labelling data may also
be used directly for phenotypic characterization of different mutant strains. Different
glucose derepressed strains investigated employed are the disruption mutants reg1,
hxk2, grr1, mig1 and mig1mig2 and the reference strain CEN.PK113-7D. Principal
components analysis of the summed fractional labelling data show that deleting the
genes HXK2 and GRR1 results in similar phenotype at the fluxome level, with a
partial alleviation of glucose repression on the respiratory metabolism. Furthermore,
deletion of the genes MIG1, MIG1/MIG2 and REG1 did not result in a significant
change in the phenotype at the fluxome level. 

**NOTES**
- From all the knockouts in this study ('REG1', 'MIG1', 'MIG2', 'GRR1', 'HXK2') only the latter had a gene correspondence in the model.

"""

code = """\
import warnings
from case_14 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case14 = Case14()
case14.model = case14.loadObjectFromFile('model_yeast_76.sav')
case14.model.solver = 'optlang-cplex'
case14.setMedium('MINIMAL')

genes = ['REG1', 'MIG1', 'MIG2', 'GRR1', 'HXK2'] #Knockouts in this study
HXK2 = case14.convertStdToSyst(genes)['HXK2'] # Gene match only for HXK2 gene

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#General datasets
exp_dataset, reactions = case14.loadExperimentalRes('Results/Case 14/case14_experimental_fluxes.csv')
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

o2_text = """\
# Oxygen flux lb estimation
"""

# ===== WILD TYPE =====

wt_text = """\
# Wild Type
"""

wt_o2 = """\
wt_etOH = case14.loadObjectFromFile('Results/Case 14/wt_dict_etOH_O2_fluxes.sav')
wt_o2_lb = case14.plotO2vsEtOH(wt_etOH, real_EtOH_flux = 2.1080, fname = 'Results/Case 14/wt_etOH_plot.png')
"""

wt_fba_datasets = """\
wt_fba_res, wt_fba_exp_sim, wt_fba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,0], o2_lb = wt_o2_lb, type = 'fba', res_exists = True, fname = 'Results/Case 14/res_fba_wt_case14.sav')
pd.concat([reactions, wt_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

wt_pfba_datasets = """\
wt_pfba_res, wt_pfba_exp_sim, wt_pfba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,0], o2_lb = wt_o2_lb, type = 'pfba', res_exists = True, fname = 'Results/Case 14/res_pfba_wt_case14.sav')
pd.concat([reactions, wt_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

wt_fva_datasets = """\
wt_fva_res, wt_fva_exp_sim, _ = case14.simulationPipeline(exp_dataset.ix[:,0], o2_lb = wt_o2_lb, type = 'fva', res_exists = True, fname = 'Results/Case 14/res_fva_wt_case14.sav')
pd.concat([reactions, wt_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== HXK2 DELETION =====

hxk2_text = """\
# HXK2 Knockout
"""

hxk2_o2 = """\
hxk2_etOH = case14.loadObjectFromFile('Results/Case 14/hxk2_dict_etOH_O2_fluxes.sav')
hxk2_02_lb = case14.plotO2vsEtOH(hxk2_etOH, real_EtOH_flux = 1.4221, fname = 'Results/Case 14/hxk2_etOH_plot.png')
"""

hxk2_fba_datasets = """\
hxk2_fba_res, hxk2_fba_exp_sim, hxk2_fba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, type = 'fba', res_exists = True, fname = 'Results/Case 14/res_fba_hxk2_case14.sav')
pd.concat([reactions, hxk2_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

hxk2_pfba_datasets = """\
hxk2_pfba_res, hxk2_pfba_exp_sim, hxk2_pfba_exp_sim_errors = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, type = 'pfba', res_exists = True, fname = 'Results/Case 14/res_pfba_hxk2_case14.sav')
pd.concat([reactions, hxk2_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

hxk2_fva_datasets = """\
hxk2_fva_res, hxk2_fva_exp_sim, _ = case14.simulationPipeline(exp_dataset.ix[:,5], o2_lb = hxk2_02_lb, type = 'fva', res_exists = True, fname = 'Results/Case 14/res_fva_hxk2_case14.sav')
pd.concat([reactions, hxk2_fva_exp_sim], axis = 1, join = 'inner')
"""



#Generate cells with plots
x = sum([['wt' + i, 'hxk2' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 14/'+ name +'_exp_sim_plot.png", width = 100%></p>'

#List with nbformat expressions
cs = ['wt', 'hxk2']
nbcells = [['nbf.v4.new_markdown_cell(' + s + '_text)',
            'nbf.v4.new_markdown_cell(o2_text)',
            'nbf.v4.new_code_cell(' + s + '_o2)',
            'nbf.v4.new_markdown_cell(fba_text)',
            'nbf.v4.new_code_cell(' + s + '_fba_datasets)',
            'nbf.v4.new_markdown_cell(' + s + '_fba_plot)',
            'nbf.v4.new_markdown_cell(pfba_text)',
            'nbf.v4.new_code_cell(' + s + '_pfba_datasets)',
            'nbf.v4.new_markdown_cell(' + s + '_pfba_plot)',
            'nbf.v4.new_markdown_cell(fva_text)',
            'nbf.v4.new_code_cell(' + s + '_fva_datasets)'] for s in cs]

nbcells = [item for sublist in nbcells for item in sublist]


nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code)] + [eval(exp) for exp in nbcells]



with open('case14.ipynb', 'w') as f:
    nbf.write(nb, f)


