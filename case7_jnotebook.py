
'''
Jupyter Notebook for Case 7

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook() #create a new notebook object

text = """\
# Case 7 Report
This report contains the results with case 7 simulations.
"""

code = """\
import warnings
from phenomenaly.io import load_yeast_76
from case_7 import *

warnings.filterwarnings('ignore')

#Initialization
case7 = Case7()
case7.model = case7.loadObjectFromFile('model_yeast_76.sav')
case7.model.solver = 'optlang-cplex'
case7.setMedium('MINIMAL_CASE7')
case7.dictsForCase7();
"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2 = case7.case7Pipeline(type = 'pfba', res_exists = True)
pd.set_option('display.max_colwidth', -1)
pd.DataFrame(reactions)
"""


fba_text = """\
## Flux Balance Analysis (FBA) Simulations
"""

fba_datasets = """\
res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors = case7.fbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, res_exists = True)
df_fba_exp_sim_errors
"""

pfba_text = """\
## Parsimonious Flux Balance Analysis (pFBA) Simulations
"""

pfba_datasets = """\
res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors = case7.pfbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, saveGenesPlot = False, plotReacts = False, plotGenes = False, res_exists = True)
df_pfba_exp_sim_errors
"""

fva_text = """\
## Flux Variability Analysis (FVA) Simulations
"""

fva_datasets = """\
res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim = case7.fvaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, res_exists = True)
df_fva_exp_sim
"""

lmoma_text = """\
## Linear Minimization of Metabolic Adjustment (LMOMA) Simulations
"""

lmoma_datasets = """\
res_lmoma, res_lmoma_df, wt_lmoma_df, df_lmoma_exp_sim, df_lmoma_exp_sim_errors = case7.lmomaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, reference_dict = res_pfba, plotGenes = False, plotReacts = False, saveGenesPlot = False, res_exists = True)
df_lmoma_exp_sim_errors
"""

# FBA Figures
fba_r = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/FBA_figs/' + str(i) + '_reacts.png", width = 100%></p>' for i in range(0,11)]
fba_reactions = 'Below are the simulated vs experimental values plotted for each reaction across the different gene knockout experiments.' + ''.join(fba_r)

# pFBA Figures
pfba_r = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/pFBA_figs/' + str(i) + '_reacts.png", width = 100%></p>' for i in range(0,11)]
pfba_reactions = 'Below are the simulated vs experimental values plotted for each reaction across the different gene knockout experiments.' + ''.join(pfba_r)

# LMOMA Figures
lmoma_r = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/LMOMA_figs/' + str(i) + '_reacts.png", width = 100%></p>' for i in range(0,11)]
lmoma_reactions = 'Below are the simulated vs experimental values plotted for each reaction across the different gene knockout experiments.' + ''.join(lmoma_r)


# Case specific results
genes_text = """\
# Analysis by Gene Knockout
"""

wt_dataset = """\
wt_res = case7.createResultsDataframeWT(reactions, wt_fba_df, wt_pfba_df, wt_lmoma_df, wt_fva_df)
wt_res
"""

genes_dataset = """\
genes_res = case7.createResultsDictByGene(df_fba_exp_sim_errors, df_pfba_exp_sim_errors, df_lmoma_exp_sim_errors, df_fva_exp_sim)
"""

genes = sorted(case7.l.keys())

#Generate cells with gene headers
for gene in genes:
    vars()[gene + '_text'] = '## ' + gene # local variable dictionary vars()


#Generate cells with gene datasets
for gene in genes:
    vars()[gene + '_dataset'] = 'genes_res["' + gene + '"]'


#Generate cells with gene images
genes_html_figs = {}
for gene in genes:
    html = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/FBA_figs/'+ gene +'_genes.png", width = 100%>Simulated vs Experiental values for FBA</p>',
            '<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/pFBA_figs/'+ gene +'_genes.png", width = 100%>Simulated vs Experiental values for pFBA</p>',
            '<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/LMOMA_figs/'+ gene +'_genes.png", width = 100%>Simulated vs Experiental values for LMOMA</p>',
            '<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/EtOH_figs/'+ gene +'_etOH_pFBA.png", width = 100%>Ethanol fluxes vs Oxygen fluxes (used to estimate Oxygen flux input)</p>']

    genes_html_figs[gene] = ''.join([html[0], html[1], html[2], html[3]])

for gene in genes:
    vars()[gene + '_images'] = genes_html_figs[gene]


#List with nbformat expressions
nbcells = [['nbf.v4.new_markdown_cell(' + gene + '_text)', 'nbf.v4.new_code_cell(' + gene + '_dataset)', 'nbf.v4.new_markdown_cell(' + gene + '_images)'] for gene in genes]
nbcells = [item for sublist in nbcells for item in sublist]


nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code),
               nbf.v4.new_markdown_cell(fba_text),
               nbf.v4.new_code_cell(fba_datasets),
               nbf.v4.new_markdown_cell(fba_reactions),
               nbf.v4.new_markdown_cell(pfba_text),
               nbf.v4.new_code_cell(pfba_datasets),
               nbf.v4.new_markdown_cell(pfba_reactions),
               nbf.v4.new_markdown_cell(fva_text),
               nbf.v4.new_code_cell(fva_datasets),
               nbf.v4.new_markdown_cell(lmoma_text),
               nbf.v4.new_code_cell(lmoma_datasets),
               nbf.v4.new_markdown_cell(lmoma_reactions),
               nbf.v4.new_markdown_cell(genes_text),
               nbf.v4.new_code_cell(wt_dataset),
               nbf.v4.new_code_cell(genes_dataset)
               ] + [eval(exp) for exp in nbcells]

with open('case7.ipynb', 'w') as f:
    nbf.write(nb, f)

