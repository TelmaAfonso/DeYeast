
'''
Jupyter Notebook for Case 13

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 13 Report (Frick et al, 2005)

This report contains the results with case 13 simulations.

Paper: [*Characterization of the metabolic shift between oxidative and fermentative growth in Saccharomyces cerevisiae by comparative 13C flux analysis*](http://www.microbialcellfactories.com/content/4/1/30)

**Background:** One of the most fascinating properties of the biotechnologically important organism
Saccharomyces cerevisiae is its ability to perform simultaneous respiration and fermentation at high
growth rate even under fully aerobic conditions. In the present work, this Crabtree effect called
phenomenon was investigated in detail by comparative 13 C metabolic flux analysis of S. cerevisiae
growing under purely oxidative, respiro-fermentative and predominantly fermentative conditions.

**Results:** The metabolic shift from oxidative to fermentative growth was accompanied by complex
changes of carbon flux throughout the whole central metabolism. This involved a flux redirection
from the pentose phosphate pathway (PPP) towards glycolysis, an increased flux through pyruvate
carboxylase, the fermentative pathways and malic enzyme, a flux decrease through the TCA cycle,
and a partial relocation of alanine biosynthesis from the mitochondrion to the cytosol. S. cerevisiae
exhibited a by-pass of pyruvate dehydrogenase in all physiological regimes. During oxidative growth
this by-pass was mainly provided via pyruvate decarboxylase, acetaldehyde dehydrogenase, acetyl-
CoA synthase and transport of acetyl-CoA into the mitochondrion. During fermentative growth
this route, however, was saturated due to limited enzyme capacity. Under these conditions the
cells exhibited high carbon flux through a chain of reactions involving pyruvate carboxylase, the
oxaloacetate transporter and malic enzyme. During purely oxidative growth the PPP alone was
sufficient to completely supply NADPH for anabolism. During fermentation, it provided only 60 %
of the required NADPH.

**Conclusion:** We conclude that, in order to overcome the limited capacity of pyruvate
dehydrogenase, S. cerevisiae possesses different metabolic by-passes to channel carbon into the
mitochondrion. This involves the conversion of cytosolic pyruvate either into acetyl CoA or
oxaloacetate followed by intercompartmental transport of these metabolites. During oxidative
growth mainly the NAD specific isoforms of acetaldehyde dehydrogenase and isocitrate
dehydrogenase catalyze the corresponding reactions in S. cerevisiae, whereas NADPH supply under
fermentative conditions involves significant contribution of sources other than the PPP such as e.
g. NADPH specific acetaldehyde dehydrogenase or isocitrate dehydrogenase.


**NOTES**
- No ethanol fluxes available for O2 flux estimation

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 13/frick_2005_fig3", width = 80%></p>


"""

code = """\
import warnings
from case_13 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case13 = Case13()
case13.model = case13.loadObjectFromFile('model_yeast_76.sav')
case13.model.solver = 'optlang-cplex'
case13.setMedium('MINERAL')
case13.dictsForCase13()

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#General datasets
exp_dataset, reactions = case13.loadExperimentalRes('Results/Case 13/case13_experimental_fluxes.csv')
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


# ===== OXIDATIVE GROWTH =====

o_text = """\
# Oxidative Growth
"""

o_fba_datasets = """\
o_exp_df = case13.getColumnWithoutNAs(exp_dataset, 0, 'X')
# NO EtOH fluxes available for O2 flux estimation

o_fba_res, o_fba_exp_sim, o_fba_exp_sim_errors = case13.simulationPipeline(o_exp_df, cs = 'g_oxidative', type = 'fba', res_exists = True, fname = 'Results/Case 13/res_fba_oxidative_case13.sav')
pd.DataFrame(reactions).join(o_fba_exp_sim_errors, how = 'inner')
"""

o_pfba_datasets = """\
o_pfba_res, o_pfba_exp_sim, o_pfba_exp_sim_errors = case13.simulationPipeline(o_exp_df, cs = 'g_oxidative',type = 'pfba', res_exists = True, fname = 'Results/Case 13/res_pfba_oxidative_case13.sav')
pd.DataFrame(reactions).join(o_pfba_exp_sim_errors, how = 'inner')
"""

o_fva_datasets = """\
o_fva_res, o_fva_exp_sim, _ = case13.simulationPipeline(o_exp_df, cs = 'g_oxidative', type = 'fva', res_exists = True, fname = 'Results/Case 13/res_fva_oxidative_case13.sav')
pd.DataFrame(reactions).join(o_fva_exp_sim, how = 'inner')
"""


# ===== RESPIRO-FERMENTATIVE GROWTH =====

rf_text = """\
# Respiro-fermentative growth
"""

rf_fba_datasets = """\
rf_exp_df = case13.getColumnWithoutNAs(exp_dataset, 1, 'X')
# NO EtOH fluxes available for O2 flux estimation

rf_fba_res, rf_fba_exp_sim, rf_fba_exp_sim_errors = case13.simulationPipeline(rf_exp_df, cs = 'g_resp_fermentative', type = 'fba', res_exists = True, fname = 'Results/Case 13/res_fba_resp_fermentative_case13.sav')
pd.DataFrame(reactions).join(rf_fba_exp_sim_errors, how = 'inner')
"""

rf_pfba_datasets = """\
rf_pfba_res, rf_pfba_exp_sim, rf_pfba_exp_sim_errors = case13.simulationPipeline(rf_exp_df, cs = 'g_resp_fermentative',type = 'pfba', res_exists = True, fname = 'Results/Case 13/res_pfba_resp_fermentative_case13.sav')
pd.DataFrame(reactions).join(rf_pfba_exp_sim_errors, how = 'inner')
"""

rf_fva_datasets = """\
rf_fva_res, rf_fva_exp_sim, _ = case13.simulationPipeline(rf_exp_df, cs = 'g_resp_fermentative', type = 'fva', res_exists = True, fname = 'Results/Case 13/res_fva_resp_fermentative_case13.sav')
pd.DataFrame(reactions).join(rf_fva_exp_sim, how = 'inner')
"""


# ===== FERMENTATIVE GROWTH =====

f_text = """\
# Fermentative growth
"""

f_fba_datasets = """\
f_exp_df = case13.getColumnWithoutNAs(exp_dataset, 2, 'X')
# NO EtOH fluxes available for O2 flux estimation

f_fba_res, f_fba_exp_sim, f_fba_exp_sim_errors = case13.simulationPipeline(f_exp_df, cs = 'g_fermentative', type = 'fba', res_exists = True, fname = 'Results/Case 13/res_fba_fermentative_case13.sav')
pd.DataFrame(reactions).join(f_fba_exp_sim_errors, how = 'inner')
"""

f_pfba_datasets = """\
f_pfba_res, f_pfba_exp_sim, f_pfba_exp_sim_errors = case13.simulationPipeline(f_exp_df, cs = 'g_fermentative',type = 'pfba', res_exists = True, fname = 'Results/Case 13/res_pfba_fermentative_case13.sav')
pd.DataFrame(reactions).join(f_pfba_exp_sim_errors, how = 'inner')
"""

f_fva_datasets = """\
f_fva_res, f_fva_exp_sim, _ = case13.simulationPipeline(f_exp_df, cs = 'g_fermentative', type = 'fva', res_exists = True, fname = 'Results/Case 13/res_fva_fermentative_case13.sav')
pd.DataFrame(reactions).join(f_fva_exp_sim, how = 'inner')
"""



#Generate cells with plots
x = sum([['o' + i, 'rf' + i, 'f' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 13/'+ name +'_exp_sim_plot.png", width = 100%></p>'

#List with nbformat expressions
growth = ['o', 'rf', 'f']
nbcells = [['nbf.v4.new_markdown_cell(' + g + '_text)',
            'nbf.v4.new_markdown_cell(fba_text)',
            'nbf.v4.new_code_cell(' + g + '_fba_datasets)',
            'nbf.v4.new_markdown_cell(' + g + '_fba_plot)',
            'nbf.v4.new_markdown_cell(pfba_text)',
            'nbf.v4.new_code_cell(' + g + '_pfba_datasets)',
            'nbf.v4.new_markdown_cell(' + g + '_pfba_plot)',
            'nbf.v4.new_markdown_cell(fva_text)',
            'nbf.v4.new_code_cell(' + g + '_fva_datasets)'] for g in growth]

nbcells = [item for sublist in nbcells for item in sublist]


nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code)] + [eval(exp) for exp in nbcells]



with open('case13.ipynb', 'w') as f:
    nbf.write(nb, f)


