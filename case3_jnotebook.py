
'''
Jupyter Notebook for Case 3

Author: Telma Afonso
'''

import nbformat as nbf

nb = nbf.v4.new_notebook()

text = """\
# Case 3 Report (Fendt et al, 2010)
This report contains the results with case 3 simulations.

Paper: [*Transcriptional regulation of respiration in yeast metabolizing differently repressive carbon substrates*](http://www.biomedcentral.com/1752-0509/4/12)

**Background:** Depending on the carbon source, Saccharomyces cerevisiae displays various degrees of respiration.
These range from complete respiration as in the case of ethanol, to almost complete fermentation, and thus very
low degrees of respiration on glucose. While many key regulators are known for these extreme cases, we focus
here on regulators that are relevant at intermediate levels of respiration.

**Results:** We address this question by linking the functional degree of respiration to transcriptional regulation via
enzyme abundances. Specifically, we investigated aerobic batch cultures with the differently repressive carbon
sources glucose, mannose, galactose and pyruvate. Based on 13 C flux analysis, we found that the respiratory
contribution to cellular energy production was largely absent on glucose and mannose, intermediate on galactose
and highest on pyruvate. In vivo abundances of 40 respiratory enzymes were quantified by GFP-fusions under each
condition. During growth on the partly and fully respired substrates galactose and pyruvate, several TCA cycle and
respiratory chain enzymes were significantly up-regulated. From these enzyme levels and the known regulatory
network structure, we determined the probability for a given transcription factor to cause the coordinated
expression changes. The most probable transcription factors to regulate the different degrees of respiration were
Gcr1p, Cat8p, the Rtg-proteins and the Hap-complex. For the latter three ones we confirmed their importance for
respiration by quantifying the degree of respiration and biomass yields in the corresponding deletion strains.

**Conclusions:** Cat8p is required for wild-type like respiration, independent of its known activation of gluconeogenic
genes. The Rtg-proteins and the Hap-complex are essential for wild-type like respiration under partially respiratory
conditions. Under fully respiratory conditions, the Hap-complex, but not the Rtg-proteins are essential for
respiration.

**NOTES**
- No ethanol fluxes available to estimate O2 lb

<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 3/Fendt_2010_fig1", width = 60%></p>

"""

code = """\
import warnings
from case_3 import *

pd.set_option('display.max_colwidth', -1)
warnings.filterwarnings('ignore')

#Initialization
case3 = Case3()
case3.model = case3.loadObjectFromFile('model_yeast_76.sav')
case3.model.solver = 'optlang-cplex'
case3.setMedium('MINIMAL')
case3.dictsForCase3()

"""

datasets_text = """\
## General datasets
"""

datasets_code = """\
#General datasets
exp_dataset, reactions = case3.loadExperimentalRes('Results/Case 3/case3_experimental_fluxes.csv')
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
g_exp_df = case3.getColumnWithoutNAs(exp_dataset, 0)

g_fba_res, g_fba_exp_sim, g_fba_exp_sim_errors = case3.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_glucose_case3.sav')
pd.concat([reactions, g_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_pfba_datasets = """\
g_pfba_res, g_pfba_exp_sim, g_pfba_exp_sim_errors = case3.simulationPipeline(g_exp_df, cs = 'glucose', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_glucose_case3.sav')
pd.concat([reactions, g_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

g_fva_datasets = """\
g_fva_res, g_fva_exp_sim, _ = case3.simulationPipeline(g_exp_df, cs = 'glucose', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_glucose_case3.sav')
pd.concat([reactions, g_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== MANNOSE =====

m_text = """\
# Mannose carbon source
"""

m_fba_datasets = """\
m_exp_df = case3.getColumnWithoutNAs(exp_dataset, 1)

m_fba_res, m_fba_exp_sim, m_fba_exp_sim_errors = case3.simulationPipeline(m_exp_df, cs = 'mannose', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_mannose_case10.sav')
pd.concat([reactions, m_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

m_pfba_datasets = """\
m_pfba_res, m_pfba_exp_sim, m_pfba_exp_sim_errors = case3.simulationPipeline(m_exp_df, cs = 'mannose', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_mannose_case10.sav')
pd.concat([reactions, m_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

m_fva_datasets = """\
m_fva_res, m_fva_exp_sim, _ = case3.simulationPipeline(m_exp_df, cs = 'mannose', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_mannose_case10.sav')
pd.concat([reactions, m_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== GALACTOSE =====

gal_text = """\
# Galactose carbon source
"""

gal_fba_datasets = """\
gal_exp_df = case3.getColumnWithoutNAs(exp_dataset, 2)

gal_fba_res, gal_fba_exp_sim, gal_fba_exp_sim_errors = case3.simulationPipeline(gal_exp_df, cs = 'galactose', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_galactose_case3.sav')
pd.concat([reactions, gal_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

gal_pfba_datasets = """\
gal_pfba_res, gal_pfba_exp_sim, gal_pfba_exp_sim_errors = case3.simulationPipeline(gal_exp_df, cs = 'galactose', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_galactose_case3.sav')
pd.concat([reactions, gal_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

gal_fva_datasets = """\
gal_fva_res, gal_fva_exp_sim, _ = case3.simulationPipeline(gal_exp_df, cs = 'galactose', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_galactose_case3.sav')
pd.concat([reactions, gal_fva_exp_sim], axis = 1, join = 'inner')
"""


# ===== PYRUVATE =====

p_text = """\
# Pyruvate carbon source
"""

p_fba_datasets = """\
p_exp_df = case3.getColumnWithoutNAs(exp_dataset, 3)

p_fba_res, p_fba_exp_sim, p_fba_exp_sim_errors = case3.simulationPipeline(p_exp_df, cs = 'pyruvate', type = 'fba', res_exists = True, fname = 'Results/Case 3/res_fba_pyruvate_case3.sav')
pd.concat([reactions, p_fba_exp_sim_errors], axis = 1, join = 'inner')
"""

p_pfba_datasets = """\
p_pfba_res, p_pfba_exp_sim, p_pfba_exp_sim_errors = case3.simulationPipeline(p_exp_df, cs = 'pyruvate', type = 'pfba', res_exists = True, fname = 'Results/Case 3/res_pfba_pyruvate_case3.sav')
pd.concat([reactions, p_pfba_exp_sim_errors], axis = 1, join = 'inner')
"""

p_fva_datasets = """\
p_fva_res, p_fva_exp_sim, _ = case3.simulationPipeline(p_exp_df, cs = 'pyruvate', type = 'fva', res_exists = True, fname = 'Results/Case 3/res_fva_pyruvate_case3.sav')
pd.concat([reactions, p_fva_exp_sim], axis = 1, join = 'inner')
"""

# ===== SUMMARY =====

g_summary = """\
Brief summary of results shown below:

- TCA cycle more active in simulated fluxes (more pyruvate enters the cell, maybe because it is not being converted to acetate/ethanol; pyruvate FVA flux range: 0 - 4.45);
- Glycolysis 1/4 less active compared to experimental fluxes (27% of the D-Glucose-6-phosphate is being converted to D-Glucose-1-phosphate with a flux of 4.34, not shown here);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 1% of D-Glucose-6-phosphate is converted into 6-phosphono-D-glucono-1,5-lactone) ;
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate
    - r_2034 	Pyruvate <==> Pyruvate-mit

"""

m_summary = """\

Brief summary of results shown below:

- TCA cycle more active in simulated fluxes (more pyruvate enters the cell, maybe because it is not being converted to acetate/ethanol; pyruvate FVA flux range: 0 - 3.50);
- Glycolysis less active in simualted fluxes (31% of the D-Fructose-6-phosphate is converted in D-Glucose-6-phosphate which about 94% is being converted to D-Glucose-1-phosphate with a flux of 3.67, not shown here);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 2% and 5% of the produced D-Fructose-6-phosphate and D-Glucose-6-phosphate, respectively, are directed to this path);
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate
    - r_2034 	Pyruvate <==> Pyruvate-mit
    
"""

gal_summary = """\

Brief summary of results shown below:

- TCA cycle activity is similar in simulated fluxes (although more pyruvate enters the cell, maybe because it is not being converted to acetate/ethanol; pyruvate FVA flux range: 0 - 1.23);
- Overall, glycolysis activity is similar in simualted fluxes (D-galactose is practically all converted into D-Glucose-6-phosphate of which 98% is converted into D-Fructose-6-phosphate);
- Pentose phosphate pathway activity reduced in simulated fluxes (only 2% of D-Fructose-6-phosphate is converted into D-Glucono-1,5-lactone-6-phosphate);
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate 	
    - r_1048 	Sedoheptulose-7-phosphate + D-Glyceraldehyde-3-phosphate <==> D-Erythrose-4-phosphate + D-Fructose-6-phosphate
    - r_0886 	D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate
    - r_2034 	Pyruvate <==> Pyruvate-mit
    
"""

p_summary = """\

Brief summary of results shown below:

- TCA cycle more active in simulated fluxes (practically all pyruvate enters the mitochondrion (96%); no acetate is being produced, however, ethanol is being produced with a flux of 0.737, not shown here);
- In TCA cycle the conversion of Isocitrate <==> 2-Oxoglutarate + CO2-mit is low because isocitrate is being converted to succinate in the cytoplasm (fluxes not shown here);
- Reaction D-Fructose-6-phosphate <==> D-Fructose-1,6-bisphosphate (r_0886) experimental values show it is an inverse reaction, therefore the correct comparison should be made with reaction id r_0449 which has a flux value of 0.34 (not shown here);
- Pentose phosphate pathway practically inactive comparing to experimental fluxes;
- FVA shows most reactions are fixed, except:
    - r_0713 	(S)-Malate-mit <==> Oxaloacetate-mit
    - r_1022 	Succinate-mit + CoA <==> Succinyl-CoA 	
    - r_0450 	D-Fructose-1,6-bisphosphate <==> Glycerone-phosphate + D-Glyceraldehyde-3-phosphate
    
"""

#Generate cells with plots
x = sum([['g' + i, 'm' + i, 'gal' + i, 'p' + i] for i in ['_fba', '_pfba']], [])
for name in x:
    vars()[name + '_plot'] = '<p style="float: center; font-size: 9pt; text-align: center; width: 100%;"><img src = "Results/Case 3/'+ name +'_exp_sim_plot.png", width = 80%></p>'

#List with nbformat expressions
cs = ['g', 'm', 'gal', 'p']
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



with open('case3.ipynb', 'w') as f:
    nbf.write(nb, f)


