# -*- coding: utf-8 -*-

'''
Jupyter Notebook for Case 9

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 9 Report (Gombert et al, 2001)

This report contains the results with case 9 simulations.

Paper: [*Network Identification and Flux Quantification in the Central Metabolism of Saccharomyces cerevisiae under Different Conditions of Glucose Repression*](http://jb.asm.org/content/183/4/1441.long)

**Abstract**
The network structure and the metabolic fluxes in central carbon metabolism were characterized in aerobically 
grown cells of Saccharomyces cerevisiae. The cells were grown under both high and low glucose concentrations, 
i.e., either in a chemostat at steady state with a specific growth rate of 0.1 h ⴚ1 or in a batch culture with
a specific growth rate of 0.37 h ⴚ1 . Experiments were carried out using [1- 13 C]glucose as the limiting substrate,
and the resulting summed fractional labelings of intracellular metabolites were measured by gas chromatography 
coupled to mass spectrometry. The data were used as inputs to a flux estimation routine that involved
appropriate mathematical modelling of the central carbon metabolism of S. cerevisiae. The results showed that
the analysis is very robust, and it was possible to quantify the fluxes in the central carbon metabolism under
both growth conditions. In the batch culture, 16.2 of every 100 molecules of glucose consumed by the cells entered 
the pentose-phosphate pathway, whereas the same relative flux was 44.2 per 100 molecules in the chemostat. 
The tricarboxylic acid cycle does not operate as a cycle in batch-growing cells, in contrast to the
chemostat condition. Quantitative evidence was also found for threonine aldolase and malic enzyme activities,
in accordance with published data. Disruption of the MIG1 gene did not cause changes in the metabolic network 
structure or in the flux pattern.

**NOTES**
- Wild-Type S.Cerevisiae CEN.PK 113-7D (MATa MAL2-8c SUC2) used in this study is very similar to that used to build the model (S288C)
- O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 9/gombert_2001_fig2", width = 80%></p>

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


# ===== SUMMARY =====

b_summary = """\

Brief summary of the results shown below:

- Glycolysis simulated fluxes are lower than those in experimental data (72% of the D-Glucose-6-phosphate is converted into D-Fructose-6-phosphate; 27% is converted into D-Glucose-1-phosphate);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 1% of D-Glucose-6-phosphate is converted into 6-phosphono-D-glucono-1,5-lactone);
- TCA cycle activity is higher in simulated fluxes. This could be due to the fact that experimentally there is a considerable production of ethanol, whereas simulated fluxes show that
only 3% of produced piruvate is used for fermentative purposes;
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0454 	Succinate-mit <==> Fumarate
    - r_0302 	Citrate <==> Isocitrate
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate
    - r_0163 	Ethanol <==> Acetaldehyde

"""

c_summary = """\

Brief summary of the results shown below:

- Overall, glycolysis simulated fluxes are similar to those in experimental data (72% of the D-Glucose-6-phosphate is converted into D-Fructose-6-phosphate; 27% is converted into D-Glucose-1-phosphate);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 1% of D-Glucose-6-phosphate is converted into 6-phosphono-D-glucono-1,5-lactone);
- TCA cycle activity is lower in simulated fluxes. This could be due to the fact that only about 51% of the produced pyruvate enters the mitochondrion and from this only 57% is converted into Acetyl-CoA
(37% is converted into 2-acetyllactic acid);
- Although produced in low amounts, simulated fluxes show that acetate is being produced as observed experimentally;
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0454 	Succinate-mit <==> Fumarate
    - r_0302 	Citrate <==> Isocitrate
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate
    - r_0163 	Ethanol <==> Acetaldehyde

"""


#Generate cells with plots
x = sum([['b' + i, 'c' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 80%;"><img src = "Results/Case 9/'+ name +'_exp_sim_plot.png", width = 100%></p>'

#List with nbformat expressions
cultures = ['b', 'c']
nbcells = [['nbf.v4.new_markdown_cell(' + cult + '_text)',
            'nbf.v4.new_markdown_cell(' + cult + '_summary)',
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



with open('case9.ipynb', 'w', encoding = "utf-8") as f:
    nbf.write(nb, f)


