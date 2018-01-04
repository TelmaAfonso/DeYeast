
'''
Jupyter Notebook for Case 6

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook() #create a new notebook object

text = """\
# Case 6 Report (dos Santos et al, 2003)
This report contains the results with case 6 simulations.

Paper: [*Identification of In Vivo Enzyme Activities in the Cometabolism of Glucose and Acetate by Saccharomyces cerevisiae by Using 13C-Labeled Substrates*](http://ec.asm.org/content/2/3/599.long)

**Abstract**
A detailed characterization of the central metabolic network of Saccharomyces cerevisiae CEN.PK 113-7D was
carried out during cometabolism of different mixtures of glucose and acetate, using aerobic C-limited chemostats 
in which one of these two substrates was labeled with 13 C. To confirm the role of malic enzyme, an isogenic
strain with the corresponding gene deleted was grown under the same conditions. The labeling patterns of proteinogenic 
amino acids were analyzed and used to estimate metabolic fluxes and/or make inferences about the in vivo
activities of enzymes of the central carbon metabolism and amino acid biosynthesis. Malic enzyme flux increased 
linearly with increasing acetate fraction. During growth on a very-high-acetate fraction, the activity of
malic enzyme satisfied the biosynthetic needs of pyruvate in the mitochondria, while in the cytosol pyruvate was
supplied via pyruvate kinase. In several cases enzyme activities were unexpectedly detected, e.g., the glyoxylate
shunt for a very-low-acetate fraction, phosphoenolpyruvate carboxykinase for an acetate fraction of 0.46 C-mol
of acetate/C-mol of substrate, and glucose catabolism to CO 2 via the tricarboxylic acid cycle for a very-high-acetate 
fraction. Cytoplasmic alanine aminotransferase activity was detected, and evidence was found that alpha-iso-propylmalate 
synthase has two active forms in vivo, one mitochondrial and the other a short cytoplasmic form.

**NOTES**
- Wild-Type S.Cerevisiae CEN.PK 113-7D (MATaMAL2-8c SUC2) used in this study is very similar to that used to build the model (S288C)
- No experimental ethanol fluxes available to adjust O2 lb
- Authors did not provide specific rate values (used Sophia's rates instead)

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 6/dos_Santos_2003_fig3", width = 60%></p>

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


fba_text = """\
## Flux Balance Analysis (FBA) Simulation
"""

pfba_text = """\
## Parsimonious Flux Balance Analysis (pFBA) Simulation
"""

fva_text = """\
## Flux Variability Analysis (FVA) Simulation
"""

lmoma_text = """\
## Linear Minimization of Metabolic Adjustment (LMOMA) Simulation
"""

# ===== WILD TYPE =====

wt_text = """\
# Wild Type
"""

wt_fba_datasets = """\
wt_fba_res, wt_fba_exp_sim, wt_fba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'fba', res_exists = True)
pd.concat([reactions, wt_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

wt_pfba_datasets = """\
wt_pfba_res, wt_pfba_exp_sim, wt_pfba_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'pfba', res_exists = True)
pd.concat([reactions, wt_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

wt_fva_datasets = """\
wt_fva_res, wt_fva_exp_sim, _ = case6.simulationPipeline(exp_dataset.ix[:,1], type = 'fva', res_exists = True)
pd.concat([reactions, wt_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== MAE1 DELETION =====

mae1_text = """\
# MAE1 Deletion
"""

mae1_fba_datasets = """\
mae1 = case6.convertStdToSyst(['MAE1'])['MAE1'] #Get corresponding gene ID in the model

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

mae1_lmoma_datasets = """\
mae1_lmoma_res, mae1_lmoma_exp_sim, mae1_lmoma_exp_sim_errors = case6.simulationPipeline(exp_dataset.ix[:,0], geneko = mae1, type = 'lmoma', res_exists = True)
pd.concat([reactions, mae1_lmoma_exp_sim_errors], axis = 1, join = 'inner')
"""


# ===== SUMMARY =====

wt_summary = """\

Brief summary of the results shown below:

- Overall glycolysis simulated fluxes are similar to the experimental data (72% of the D-Glucose-6-phosphate is converted into D-Fructose-6-phosphate; 27% is converted into D-Glucose-1-phosphate);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 1% of D-Glucose-6-phosphate is converted into 6-phosphono-D-glucono-1,5-lactone);
- TCA cycle activity is lower in simulated fluxes. This could be due to the fact that only about 57% of the pyruvate that enters the mitochondrion is converted into Acetyl-CoA 
(37% is converted into 2-acetyllactic acid);
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0454 	Succinate-mit <==> Fumarate 
    - r_0302 	Citrate <==> Isocitrate 		
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate

"""

mae1_summary = """\

Brief summary of the results shown below:

- Overall glycolysis simulated fluxes are similar to the experimental data (71% of the D-Glucose-6-phosphate is converted into D-Fructose-6-phosphate; 26% is converted into D-Glucose-1-phosphate);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 3% of D-Glucose-6-phosphate is converted into 6-phosphono-D-glucono-1,5-lactone);
- TCA cycle activity is lower in simulated fluxes. This could be due to the fact that only about 54% of the produced pyruvate enters the mitochondrion and from this only 57% is converted into Acetyl-CoA 
(37% is converted into 2-acetyllactic acid);
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0454 	Succinate-mit <==> Fumarate 
    - r_0302 	Citrate <==> Isocitrate 		
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate

"""


#Generate cells with plots
x = sum([['wt' + i, 'mae1' + i] for i in ['_fba', '_pfba', '_lmoma']], [])
x.remove('wt_lmoma')

for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 6/'+ name +'_exp_sim_plot.png", width = 80%></p>'

#List with nbformat expressions
cases = ['wt', 'mae1']
nbcells = [['nbf.v4.new_markdown_cell(' + c + '_text)',
            'nbf.v4.new_markdown_cell(' + c + '_summary)',
            'nbf.v4.new_markdown_cell(fba_text)',
            'nbf.v4.new_code_cell(' + c + '_fba_datasets)',
            'nbf.v4.new_markdown_cell(' + c + '_fba_plot)',
            'nbf.v4.new_markdown_cell(pfba_text)',
            'nbf.v4.new_code_cell(' + c + '_pfba_datasets)',
            'nbf.v4.new_markdown_cell(' + c + '_pfba_plot)',
            'nbf.v4.new_markdown_cell(fva_text)',
            'nbf.v4.new_code_cell(' + c + '_fva_datasets)'] for c in cases]

mae1_cells = ['nbf.v4.new_markdown_cell(lmoma_text)',
              'nbf.v4.new_code_cell(mae1_lmoma_datasets)',
              'nbf.v4.new_markdown_cell(mae1_lmoma_plot)']

nbcells[1][7:7] = mae1_cells

nbcells = [item for sublist in nbcells for item in sublist]


nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code),
               nbf.v4.new_markdown_cell(datasets_text),
               nbf.v4.new_code_cell(datasets_code)] + [eval(exp) for exp in nbcells]


with open('case6.ipynb', 'w') as f:
    nbf.write(nb, f)



