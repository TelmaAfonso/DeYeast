
'''
Jupyter Notebook for Case 10

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 10 Report (Daran-Lapujade et al, 2004)

This report contains the results with case 10 simulations.

Paper: [*Role of Transcriptional Regulation in Controlling Fluxes in Central Carbon Metabolism of Saccharomyces cerevisiae*](http://www.jbc.org/content/279/10/9125.long)

**Abstract**
In contrast to batch cultivation, chemostat cultivation allows the identification of carbon source responses
without interference by carbon-catabolite repression, accumulation of toxic products, and differences in specific 
growth rate. This study focuses on the yeast Saccharomyces cerevisiae, grown in aerobic, carbon-limited
chemostat cultures. Genome-wide transcript levels and in vivo fluxes were compared for growth on two sugars,
glucose and maltose, and for two C2-compounds, ethanol and acetate. In contrast to previous reports on batch
cultures, few genes (180 genes) responded to changes of the carbon source by a changed transcript level. Very
few transcript levels were changed when glucose as the growth-limiting nutrient was compared with maltose
(33 transcripts), or when acetate was compared with ethanol (16 transcripts). Although metabolic flux analy-
sis using a stoichiometric model revealed major changes in the central carbon metabolism, only 117 genes exhib-
ited a significantly different transcript level when sugars and C2-compounds were provided as the growth-
limiting nutrient. Despite the extensive knowledge on carbon source regulation in yeast, many of the carbon
source-responsive genes encoded proteins with unknown or incompletely characterized biological functions. 
In silico promoter analysis of carbon source-responsive genes confirmed the involvement of several
known transcriptional regulators and suggested the involvement of additional regulators. Transcripts 
involved in the glyoxylate cycle and gluconeogenesis showed a good correlation with in vivo fluxes. This 
correlation was, however, not observed for other important pathways, including the pentose-phosphate pathway,
tricarboxylic acid cycle, and, in particular, glycolysis. These results indicate that in vivo fluxes in the central
carbon metabolism of S. cerevisiae grown in steady-state, carbon-limited chemostat cultures are controlled
to a large extent via post-transcriptional mechanisms.

**NOTES**
- Wild-Type S.Cerevisiae CEN.PK 113-7D (MATa MAL2-8c SUC2) used in this study is very similar to that used to build the model (S288C)
- Maltose only participates in two reactions (transport and exchange) in the model (maltose metabolism not incorporated into Yeast 7?)

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 10/daran_lapujade_2004_fig1", width = 80%></p>

"""

code = """\
import warnings
from case_10 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case10 = Case10()
case10.model = case10.loadObjectFromFile('model_yeast_76.sav')
case10.model.solver = 'optlang-cplex'
case10.setMedium('MINIMAL')
case10.dictsForCase10()

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#General datasets
exp_dataset, reactions = case10.loadExperimentalRes('Results/Case 10/case10_experimental_fluxes.csv')
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


# ===== GLUCOSE =====

g_text = """\
# Glucose carbon source
"""

g_fba_datasets = """\
g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,0], cs = 'glucose', type = 'fba',res_exists = True, fname = 'Results/Case 10/res_fba_glucose_case10.sav')
pd.concat([reactions, g_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_pfba_datasets = """\
g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,0], cs = 'glucose', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_glucose_case10.sav')
pd.concat([reactions, g_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_fva_datasets = """\
g_fva_res, g_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,0], cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_glucose_case10.sav')
pd.concat([reactions, g_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== MALTOSE =====

m_text = """\
# Maltose carbon source
"""

m_fba_datasets = """\
m_fba_res, m_fba_exp_sim, m_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,1], cs = 'maltose', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_maltose_case10.sav')
pd.concat([reactions, m_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

m_pfba_datasets = """\
m_pfba_res, m_pfba_exp_sim, m_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,1], cs = 'maltose', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_maltose_case10.sav')
pd.concat([reactions, m_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

m_fva_datasets = """\
m_fva_res, m_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,1], cs = 'maltose', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_maltose_case10.sav')
pd.concat([reactions, m_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== ETHANOL =====

e_text = """\
# Ethanol carbon source
"""

e_fba_datasets = """\
e_fba_res, e_fba_exp_sim, e_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,2], cs = 'ethanol', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_ethanol_case10.sav')
pd.concat([reactions, e_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

e_pfba_datasets = """\
e_pfba_res, e_pfba_exp_sim, e_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,2], cs = 'ethanol', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_ethanol_case10.sav')
pd.concat([reactions, e_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

e_fva_datasets = """\
e_fva_res, e_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,2], cs = 'ethanol', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_ethanol_case10.sav')
pd.concat([reactions, e_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== ACETATE =====

a_text = """\
# Acetate carbon source
"""

a_fba_datasets = """\
a_fba_res, a_fba_exp_sim, a_fba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,3], cs = 'acetate', type = 'fba', res_exists = True, fname = 'Results/Case 10/res_fba_acetate_case10.sav')
pd.concat([reactions, a_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

a_pfba_datasets = """\
a_pfba_res, a_pfba_exp_sim, a_pfba_exp_sim_errors = case10.simulationPipeline(exp_dataset.ix[:,3], cs = 'acetate', type = 'pfba', res_exists = True, fname = 'Results/Case 10/res_pfba_acetate_case10.sav')
pd.concat([reactions, a_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

a_fva_datasets = """\
a_fva_res, a_fva_exp_sim, _ = case10.simulationPipeline(exp_dataset.ix[:,3], cs = 'acetate', type = 'fva', res_exists = True, fname = 'Results/Case 10/res_fva_acetate_case10.sav')
pd.concat([reactions, a_fva_exp_sim], axis = 1, join = 'inner')
"""


#Generate cells with plots
x = sum([['g' + i, 'm' + i, 'e' + i, 'a' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 10/'+ name +'_exp_sim_plot.png", width = 100%></p>'

#List with nbformat expressions
cs = ['g', 'm', 'e', 'a']
nbcells = [['nbf.v4.new_markdown_cell(' + s + '_text)',
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



with open('case10.ipynb', 'w') as f:
    nbf.write(nb, f)


