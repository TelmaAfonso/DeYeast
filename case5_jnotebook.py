
'''
Jupyter Notebook for Case 5

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 5 Report (Çakir et al, 2004)
This report contains the results with case 5 simulations.

Paper: [*Metabolic Pathway Analysis of Yeast Strengthens the Bridge Between Transcriptomics and Metabolic Networks*](http://onlinelibrary.wiley.com/doi/10.1002/bit.20020/abstract)

**Abstract**
Central carbon metabolism of the yeast Saccharomyces cerevisiae was analyzed using metabolic pathway analysis tools. Elementary flux modes for three substrates
(glucose, galactose, and ethanol) were determined using the catabolic reactions occurring in yeast. Resultant elementary modes were used for gene deletion phenotype
analysis and for the analysis of robustness of the central metabolism and network functionality. Control-effective fluxes, determined by calculating the efficiency of each
mode, were used for the prediction of transcript ratios of metabolic genes in different growth media (glucose – ethanol and galactose – ethanol). A high correlation was
obtained between the theoretical and experimental expression levels of 38 genes when ethanol and glucose media were considered. Such analysis was shown to be a
bridge between transcriptomics and fluxomics. Control-effective flux distribution was found to be promising in *in silico* predictions by incorporating functionality and 
regulation into the metabolic network structure. 

**NOTES**
- O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)
- Authors did not provide specific rate values (used Sophia's rates instead)

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 5/Çakir_2004_fig4", width = 50%></p>

"""

code = """\
import warnings
from case_5 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case5 = Case5()
case5.model = case5.loadObjectFromFile('model_yeast_76.sav')
case5.model.solver = 'optlang-cplex'
case5.setMedium('MINIMAL')
case5.dictsForCase5()

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#General datasets
exp_dataset, reactions = case5.loadExperimentalRes('Results/Case 5/case5_experimental_fluxes.csv')
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
g_exp_df = case5.getColumnWithoutNAs(exp_dataset, 0)
# O2 flux estimation not possible (ethanol flux of 0 independently of O2 flux)

g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case5.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 5/res_fba_glucose_case5.sav')
pd.concat([reactions, g_fba_exp_sim_errors], axis = 1, join = 'inner') #Plot not showing r_0302
"""

g_pfba_datasets = """\
g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case5.simulationPipeline(g_exp_df, cs = 'glucose', type = 'pfba', res_exists = True, fname = 'Results/Case 5/res_pfba_glucose_case5.sav')
pd.concat([reactions, g_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_fva_datasets = """\
g_fva_res, g_fva_exp_sim, _ = case5.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 5/res_fva_glucose_case5.sav')
pd.concat([reactions, g_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== ETHANOL =====

e_text = """\
# Ethanol carbon source
"""

e_fba_datasets = """\
e_exp_df = case5.getColumnWithoutNAs(exp_dataset, 1)

e_fba_res, e_fba_exp_sim, e_fba_exp_sim_errors = case5.simulationPipeline(e_exp_df, cs = 'ethanol', type = 'fba', res_exists = False, fname = 'Results/Case 5/res_fba_ethanol_case5.sav')
pd.concat([reactions, e_fba_exp_sim_errors], axis = 1, join = 'inner') #Plot not showing r_0302
"""

e_pfba_datasets = """\
e_pfba_res, e_pfba_exp_sim, e_pfba_exp_sim_errors = case5.simulationPipeline(e_exp_df, cs = 'ethanol', type = 'pfba', res_exists = False, fname = 'Results/Case 5/res_pfba_ethanol_case5.sav')
pd.concat([reactions, e_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

e_fva_datasets = """\
e_fva_res, e_fva_exp_sim, _ = case5.simulationPipeline(e_exp_df, cs = 'ethanol', type = 'fva', res_exists = False, fname = 'Results/Case 5/res_fva_ethanol_case5.sav')
pd.concat([reactions, e_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== SUMMARY =====

g_summary = """\
Brief summary of results shown below:

- Overall glycolysis simulated fluxes are similar to the experimental data (72% of the D-Glucose-6-phosphate is converted into D-Fructose-6-phosphate; 27% is converted into D-Glucose-1-phosphate);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 1% of D-Glucose-6-phosphate is converted into 6-phosphono-D-glucono-1,5-lactone);
- TCA cycle overall has similar simulated and experimental fluxes, except for the fluxes of succinate and fumarate production, which are lower. This happens because malate 
is being produced on the citoplasm and enters the mitochondrion by antiport with 2-oxoglutarate, which becomes unavailable in the mitochondrion to produce succinyl-CoA;
- About 50% of the produced pyruvate enters the mitochondrion, with the remaining 38% being converted into oxaloacetate and 3% converted into acetaldehyde;
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0454 	Succinate-mit <==> Fumarate 
    - r_0302 	Citrate <==> Isocitrate 		
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate
    - r_0163 	Ethanol <==> Acetaldehyde 
    - r_0714 	(S)-Malate <==> Oxaloacetate    
    - r_2034 	Pyruvate <==> Pyruvate-mit

"""

e_summary = """\

Brief summary of results shown below:

- Glycolysis flux values are very discrepant between simulated and experimental fluxes. Experimentally the cell performs glycolysis, 
whereas simulated fluxes show the cell is performing gluconeogenesis (this could be due to excessive ethanol being given as carbon source(?) - specific rate values not provided by authors);
- Pentose phosphate pathway activity reduced in simulated fluxes;
- Overall experimental and simulated TCA cycle flux values are similar, with most discrepancy observed in the beginning of the cycle (most fluxes close to zero). This could be due to he fact that 
not much pyruvate is produced (being 100% produced from malate, with a flux of 0.296) and about 64% of it being used or the production of 2-acetyllactic acid;
- FVA shows most reactions are fixed, except:
    - r_2131 	Isocitrate <==> 2-Oxoglutarate + CO2-mit
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0454 	Succinate-mit <==> Fumarate 
    - r_0302 	Citrate <==> Isocitrate 		
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0163 	Ethanol <==> Acetaldehyde 
    - r_0714 	(S)-Malate <==> Oxaloacetate    
    - r_0718 	(S)-Malate-mit <==> Pyruvate + CO2
    
"""


#Generate cells with plots
x = sum([['g' + i, 'e' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 5/'+ name +'_exp_sim_plot.png", width = 80%></p>'

#List with nbformat expressions
cs = ['g', 'e']
nbcells = [['nbf.v4.new_markdown_cell(' + s + '_text)',
            'nbf.v4.new_markdown_cell(' + s + '_summary)',
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



with open('case5.ipynb', 'w') as f:
    nbf.write(nb, f)


