
'''
Jupyter Notebook for Case 7

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook() #create a new notebook object

text = """\
# Case 7 Report (Blank et al, 2005)

This report contains the results with case 7 simulations.

Paper: [*Large-scale 13 C-flux analysis reveals mechanistic principles of metabolic network robustness to null mutations in yeast*](http://genomebiology.com/content/6/6/R49)

**Background:** Quantification of intracellular metabolite fluxes by 13 C-tracer experiments is
maturing into a routine higher-throughput analysis. The question now arises as to which mutants
should be analyzed. Here we identify key experiments in a systems biology approach with a
genome-scale model of Saccharomyces cerevisiae metabolism, thereby reducing the workload for
experimental network analyses and functional genomics.

**Results:** Genome-scale 13 C flux analysis revealed that about half of the 745 biochemical reactions
were active during growth on glucose, but that alternative pathways exist for only 51 gene-encoded
reactions with significant flux. These flexible reactions identified in silico are key targets for
experimental flux analysis, and we present the first large-scale metabolic flux data for yeast,
covering half of these mutants during growth on glucose. The metabolic lesions were often
counteracted by flux rerouting, but knockout of cofactor-dependent reactions, as in the adh1, ald6,
cox5A, fum1, mdh1, pda1, and zwf1 mutations, caused flux responses in more distant parts of the
network. By integrating computational analyses, flux data, and physiological phenotypes of all
mutants in active reactions, we quantified the relative importance of 'genetic buffering' through
alternative pathways and network redundancy through duplicate genes for genetic robustness of
the network.

**Conclusions:** The apparent dispensability of knockout mutants with metabolic function is
explained by gene inactivity under a particular condition in about half of the cases. For the remaining
207 viable mutants of active reactions, network redundancy through duplicate genes was the major
(75%) and alternative pathways the minor (25%) molecular mechanism of genetic network
robustness in S. cerevisiae.

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 7/blank_2004_fig2.png", width = 80%></p>
<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 7/blank_2004_fig2_legend.png", width = 80%></p>

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
case7.setMedium('MINIMAL')
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
pd.concat([reactions, df_fba_exp_sim_errors], axis = 1, join = "inner")

"""

pfba_text = """\
## Parsimonious Flux Balance Analysis (pFBA) Simulations
"""

pfba_datasets = """\
res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors = case7.pfbaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, saveGenesPlot = False, plotReacts = False, plotGenes = False, res_exists = True)
pd.concat([reactions, df_pfba_exp_sim_errors], axis = 1, join = "inner")
"""

fva_text = """\
## Flux Variability Analysis (FVA) Simulations
"""

fva_datasets = """\
res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim = case7.fvaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, res_exists = True)
pd.concat([reactions, df_fva_exp_sim], axis = 1, join = "inner")
"""

lmoma_text = """\
## Linear Minimization of Metabolic Adjustment (LMOMA) Simulations
"""

lmoma_datasets = """\
res_lmoma, res_lmoma_df, df_lmoma_exp_sim, df_lmoma_exp_sim_errors = case7.lmomaPipeline(fluxes_O2 = fluxes_O2, exp_dataset = exp_dataset, reference_dict = res_pfba, plotGenes = False, plotReacts = False, saveGenesPlot = False, res_exists = True)
pd.concat([reactions, df_lmoma_exp_sim_errors], axis = 1, join = "inner")
"""

# Case specific results
wt_text = """\
## Wild Type
"""

genes_text = """\
# Analysis by Gene Knockout
"""

wt_dataset = """\
wt_res = case7.createResultsDataframeWT(reactions, wt_fba_df, wt_pfba_df, wt_fva_df)
wt_res
"""

genes_dataset = """\
genes_res = case7.createResultsDictByGene(df_fba_exp_sim_errors, df_pfba_exp_sim_errors, df_lmoma_exp_sim_errors, df_fva_exp_sim)
"""

from phenomenaly.io import load_yeast_76
from case_7 import *
import warnings

warnings.filterwarnings('ignore')

#Initialization
case7 = Case7()
case7.model = case7.loadObjectFromFile('model_yeast_76.sav')
case7.dictsForCase7();

genes = sorted(case7.l.keys())

#Generate cells with gene headers
for gene in genes:
    vars()[gene + '_text'] = '## ' + gene # local variable dictionary vars()


#Generate cells with gene datasets
for gene in genes:
    vars()[gene + '_dataset'] = 'pd.concat([reactions, genes_res["' + gene + '"]], axis = 1, join = "inner")'


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


#Generate cells with reactions images
exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2 = case7.case7Pipeline(plot = False, makeFigs = False, type = 'pfba', res_exists = True)

reactions = list(exp_dataset.index)

# FBA Figures
fba_r = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/FBA_figs/' + str(i) + '_reaction.png", width = 100%></p>' for i in reactions]
fba_reactions = 'Below are the simulated vs experimental values plotted for each reaction across the different gene knockout experiments.' + ''.join(fba_r)

# pFBA Figures
pfba_r = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/pFBA_figs/' + str(i) + '_reaction.png", width = 100%></p>' for i in reactions]
pfba_reactions = 'Below are the simulated vs experimental values plotted for each reaction across the different gene knockout experiments.' + ''.join(pfba_r)

# LMOMA Figures
lmoma_r = ['<p style="float: left; font-size: 9pt; text-align: center; width: 50%;"><img src = "Results/Case 7/LMOMA_figs/' + str(i) + '_reaction.png", width = 100%></p>' for i in reactions]
lmoma_reactions = 'Below are the simulated vs experimental values plotted for each reaction across the different gene knockout experiments.' + ''.join(lmoma_r)



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
               nbf.v4.new_markdown_cell(wt_text),
               nbf.v4.new_code_cell(wt_dataset),
               nbf.v4.new_code_cell(genes_dataset)
               ] + [eval(exp) for exp in nbcells]

with open('case7.ipynb', 'w') as f:
    nbf.write(nb, f)

