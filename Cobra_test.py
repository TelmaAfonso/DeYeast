
'''
Script just for testing cobraby functionalities

Author: Telma Afonso

'''

# from __future__ import print_function
import cobra
import cobra.test

def modelInfo (model):
    # The reactions, metabolites, and genes attributes of the cobrapy
    # model are a special type of list called a cobra.DictList
    print('Number of reactions:', len(model.reactions))
    print('Number of metabolites:', len(model.metabolites))
    print('Number of genes:', len(model.genes))
    print('Gene number 34:', model.genes[34]) #retrieve item by index
    print('Cytosolic atp metabolite:', model.metabolites.get_by_id('atp_c'))
    print('Glucose exchange reaction bounds:', model.reactions.EX_glc__D_e.bounds)

def reactions (model):
    # Glucose 6-phosphate isomerase, which interconverts glucose 6-phosphate and fructose 6-phosphate.
    pgi = model.reactions.get_by_id('PGI')
    print('Reaction name:', pgi.name)
    print('Reaction:', pgi.reaction)
    # Because the pgi.lower_bound < 0, and pgi.upper_bound > 0, pgi is reversible.
    print('Reaction bounds:', pgi.lower_bound, '< pgi <', pgi.upper_bound)
    print('Is the reaction reversible?', pgi.reversibility)
    print('Is the reaction mass balanced? (If empty yes)', pgi.check_mass_balance())
    print('Added a metabolite to the reaction')
    pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
    print('Is the reaction mass balanced? (If empty yes)', pgi.check_mass_balance())
    pgi.subtract_metabolites({model.metabolites.get_by_id("h_c"): -1}) # Removed metabolite
    print('Create a new reaction from string: "g6p_c --> f6p_c + h_c + green_eggs + ham"')
    pgi.reaction = 'g6p_c --> f6p_c + h_c + green_eggs + ham'
    print('New reaction:', pgi.reaction)
    pgi.reaction = "g6p_c <=> f6p_c"

def metabolites (model):
    atp = model.metabolites.get_by_id('atp_c')
    print('Metabolite name:', atp.name)
    print('Metabolite compartment:', atp.compartment)
    print('Metabolite charge:', atp.charge)
    print('Metabolite chemical formula', atp.formula)
    print('How many reactions use this metabolite:', len(atp.reactions))

def genes (model):
    pgi = model.reactions.get_by_id('PGI')
    gpr = pgi.gene_reaction_rule
    print('Gene requirements for this reaction to be active:', gpr)
    #print('Gene requirements for this reaction to be active:', pgi.genes)
    pgi_gene = model.genes.get_by_id("b4025")
    pgi_gene.reactions # Each gene keeps track of the reactions it catalyzes
    pgi.gene_reaction_rule = "(spam or eggs)" #Altering the gene_reaction_rule will create new gene objects if necessary and update all relationships.
    pgi_gene.reactions
    model.genes.get_by_id("spam") # Newly created genes are also added to the model

    #The delete_model_genes function will evaluate the GPR and set the upper and lower bounds
    # to 0 if the reaction is knocked out. This function can preserve existing deletions or
    # reset them using the cumulative_deletions flag.

    cobra.manipulation.delete_model_genes(model, ["spam"], cumulative_deletions = True)
    print("after 1 KO: %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))
    cobra.manipulation.delete_model_genes(model, ["eggs"], cumulative_deletions = True)
    print("after 2 KO: %4d < flux_PGI < %4d" % (pgi.lower_bound, pgi.upper_bound))
    # The undelete_model_genes can be used to reset a gene deletion
    cobra.manipulation.undelete_model_genes(model)
    print(pgi.lower_bound, "< pgi <", pgi.upper_bound)

def makingChangesReversibly (model):
    for reaction in model.reactions[:5]: #First five reactions
        with model as model:
            reaction.knock_out()
            model.optimize()
            print('%s blocked (bounds: %s), new growth rate %f' %
                 (reaction.id, str(reaction.bounds), model.objective.value))

    print('All bounds reverted: \n', [reaction.bounds for reaction in model.reactions[:5]])

    print('original objective: ', model.objective.expression)
    with model:
        model.objective = 'ATPM'
        print('print objective in first context:', model.objective.expression)
        with model:
            model.objective = 'ACALD'
            print('print objective in second context:', model.objective.expression)
        print('objective after exiting second context:', model.objective.expression)
    print('back to original objective:', model.objective.expression)

def buildModel():
    from cobra import Model, Reaction, Metabolite
    # Best practise: SBML compliant IDs
    model = Model('example_model')
    reaction = Reaction('3OAS140')
    reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
    reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0. # This is the default
    reaction.upper_bound = 1000. # This is the default

    # Create metabolites
    ACP_c = Metabolite('ACP_c', formula = 'C11H21N2O7PRS', name = 'acyl-carrier-protein',
                        compartment = 'c')
    omrsACP_c = Metabolite('3omrsACP_c', formula = 'C25H45N2O9PRS',
                           name = '3-Oxotetradecanoyl-acyl-carrier-protein', compartment = 'c')
    co2_c = Metabolite('co2_c', formula = 'CO2', name = 'CO2', compartment = 'c')
    malACP_c = Metabolite('malACP_c', formula = 'C14H22N2O10PRS',
                          name = 'Malonyl-acyl-carrier-protein', compartment = 'c')
    h_c = Metabolite('h_c', formula = 'H', name = 'H', compartment = 'c')
    ddcaACP_c = Metabolite('ddcaACP_c', formula = 'C23H43N2O8PRS',
                            name = 'Dodecanoyl-ACP-n-C120ACP', compartment = 'c')
    # Adding metabolites to a reaction requires using a dictionary of the metabolites and their stoichiometric coefficients
    reaction.add_metabolites({
        malACP_c: -1.0,
        h_c: -1.0,
        ddcaACP_c: -1.0,
        co2_c: 1.0,
        ACP_c: 1.0,
        omrsACP_c: 1.0
    })
    print(reaction.reaction) # This gives a string representation of the reaction

    reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
    print(reaction.genes) # Gene rules

    #Model info (still empty)
    print('%i reactions initially' % len(model.reactions))
    print('%i metabolites initially' % len(model.metabolites))
    print('%i genes initially' % len(model.genes))

    model.add_reactions([reaction])
    # Now there are things in the model
    print('%i reaction' % len(model.reactions))
    print('%i metabolites' % len(model.metabolites))
    print('%i genes' % len(model.genes))

    # Iterate through the the objects in the model
    print("Reactions")
    print("---------")
    for x in model.reactions:
        print("%s : %s" % (x.id, x.reaction))
    print("")
    print("Metabolites")
    print("-----------")
    for x in model.metabolites:
        print('%9s : %s' % (x.id, x.formula))
    print("")
    print("Genes")
    print("-----")
    for x in model.genes:
        associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" % (x.id, "{" + ", ".join(associated_ids) + "}"))

    # Last we need to set the objective of the model. Here, we just want this to be the
    # maximization of the flux in the single reaction we added and we do this by assigning
    # the reaction’s identifier to the objective property of the model.
    model.objective = '3OAS140'
    print(model.objective.expression)
    print(model.objective.direction) # shows that the solver will maximize the flux in the forward direction.


def readingModels ():
    import os
    from os.path import join
    data_dir = cobra.test.data_dir
    print("mini test files: ")
    print(", ".join(i for i in os.listdir(data_dir) if i.startswith("mini")))
    textbook_model = cobra.test.create_test_model("textbook")
    ecoli_model = cobra.test.create_test_model("ecoli")
    salmonella_model = cobra.test.create_test_model("salmonella")

    # SBML models
    cobra.io.read_sbml_model(join(data_dir, "mini_fbc2.xml"))
    cobra.io.write_sbml_model(textbook_model, "test_fbc2.xml")

def simulatingFBA (model):
    solution = model.optimize()
    print(solution.objective_value)
    print(solution.status)
    print(solution.fluxes)
    print(solution.shadow_prices)

    #print(%%time)
    print(model.slim_optimize()) #faster

    model.summary() #model summary
    model.metabolites.nadh_c.summary() #examine the overall redox balance of the model
    model.metabolites.atp_c.summary() #main energy production and consumption reactions

    biomass_rxn = model.reactions.get_by_id("Biomass_Ecoli_core")

    from cobra.util.solver import linear_reaction_coefficients
    print(linear_reaction_coefficients(model))

    # The objective function can be changed by assigning Model.objective,
    # which can be a reaction object (or just it’s name), or a dict of
    # {Reaction: objective_coefficient}.

    model.objective = "ATPM" # change the objective to ATPM

    # The upper bound should be 1000, so that we get the actual optimal value
    model.reactions.get_by_id("ATPM").upper_bound = 1000.
    print(linear_reaction_coefficients(model))

    print(model.optimize().objective_value)


def simulatingFVA (model):
    # FBA will not give always give unique solution, because multiple flux
    # states can achieve the same optimum. FVA (or flux variability analysis)
    # finds the ranges of each metabolic flux at the optimum.

    from cobra.flux_analysis import flux_variability_analysis
    print(flux_variability_analysis(model, model.reactions[:10]))
    # Setting parameter fraction_of_optimium=0.90 would give the flux ranges for reactions at 90% optimality.
    print(cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:10],
        fraction_of_optimum = 0.9, loopless = True))

    model.optimize()
    model.summary(fva = 0.95) # Flux variability analysis can also be embedded in calls to summary methods.
    model.metabolites.pyr_c.summary(fva = 0.95)

def simulatingpFBA (model):
    # Parsimonious FBA (often written pFBA) finds a flux distribution which
    # gives the optimal growth rate, but minimizes the total sum of flux.
    # This involves solving two sequential linear programs, but is handled transparently by cobrapy.
    model.objective = 'Biomass_Ecoli_core'
    fba_solution = model.optimize()
    pfba_solution = cobra.flux_analysis.pfba(model)
    print(abs(fba_solution.fluxes["Biomass_Ecoli_core"] - pfba_solution.fluxes["Biomass_Ecoli_core"]))

def simulatingDels ():
    import pandas
    from time import time
    import cobra.test
    from cobra.flux_analysis import (single_gene_deletion, single_reaction_deletion,
                                     double_gene_deletion, double_reaction_deletion)
    cobra_model = cobra.test.create_test_model("textbook")
    ecoli_model = cobra.test.create_test_model("ecoli")

    print('complete model: ', cobra_model.optimize())
    with cobra_model:
        cobra_model.reactions.PFK.knock_out()
        print('pfk knocked out: ', cobra_model.optimize())

    print('complete model: ', cobra_model.optimize())
    with cobra_model:
        cobra_model.genes.b1723.knock_out()
        print('pfkA knocked out: ', cobra_model.optimize())
        cobra_model.genes.b3916.knock_out()
        print('pfkB knocked out: ', cobra_model.optimize())

    deletion_results = single_gene_deletion(cobra_model)

    print(single_gene_deletion(cobra_model, cobra_model.genes[:20])) #subset
    print(single_reaction_deletion(cobra_model, cobra_model.reactions[:20]))
    print('Hello world!')
    #Double gene deletion (NOT WORKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
    # print(double_gene_deletion(cobra_model, cobra_model.genes[-5:], return_frame = True).round(4))
    # start = time() # start timer()
    # double_gene_deletion(ecoli_model, ecoli_model.genes[:300], number_of_processes = 2)
    # t1 = time() - start
    # print("Double gene deletions for 200 genes completed in "
    #       "%.2f sec with 2 cores" % t1)
    # start = time() # start timer()
    # double_gene_deletion(ecoli_model, ecoli_model.genes[:300], number_of_processes=1)
    # t2 = time() - start
    # print("Double gene deletions for 200 genes completed in "
    #       "%.2f sec with 1 core" % t2)
    # print("Speedup of %.2fx" % (t2 / t1))

    # print(double_reaction_deletion(cobra_model, cobra_model.reactions[2:7], return_frame = True).round(4))

def productionEnvelopes (): #aka phenotype phase planes
    # will show distinct phases of optimal growth with different use of two different substrates.
    import cobra.test
    from cobra.flux_analysis import production_envelope
    model = cobra.test.create_test_model("textbook")

    prod_env = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
    prod_env.head()
    prod_env = production_envelope(model, ["EX_o2_e"], objective = "EX_ac_e", c_source = "EX_glc__D_e")
    prod_env.head()
    prod_env[prod_env.direction == 'maximum'].plot(kind = 'line', x = 'EX_o2_e', y = 'carbon_yield')

def fluxSampling (model):
    from cobra.test import create_test_model
    from cobra.flux_analysis import sample
    model = create_test_model("textbook")
    s = sample(model, 100) #number of samples to generate
    s.head()
    s = sample(model, 1000)
    s

    #The sampling process can be controlled on a lower level by using the sampler classes directly.
    from cobra.flux_analysis.sampling import OptGPSampler, ACHRSampler
    achr = ACHRSampler(model, thinning = 10) #“Thinning” means only recording samples every n iterations. A higher thinning factor means less correlated samples but also larger computation times.
    optgp = OptGPSampler(model, processes=4)  #an additional processes argument specifying how many processes are used to create parallel sampling chains.

    #For OptGPSampler the number of samples should be a multiple of the number of
    # processes, otherwise it will be increased to the nearest multiple automatically.
    s1 = achr.sample(100)
    s2 = optgp.sample(100)

    # Sampling and validation
    import numpy as np
    bad = np.random.uniform(-1000, 1000, size=len(model.reactions))
    achr.validate(np.atleast_2d(bad)) #should not be feasible
    achr.validate(s1)

    # Batch sampling
    counts = [np.mean(s.Biomass_Ecoli_core > 0.1) for s in optgp.batch(100, 10)]
    print("Usually {:.2f}% +- {:.2f}% grow...".format(
    np.mean(counts) * 100.0, np.std(counts) * 100.0))

    # Adding constraints
    co = model.problem.Constraint(model.reactions.Biomass_Ecoli_core.flux_expression, lb=0.1)
    model.add_cons_vars([co])

    # Note that this is only for demonstration purposes. usually you could set
    # the lower bound of the reaction directly instead of creating a new constraint.
    s = sample(model, 10)
    print(s.Biomass_Ecoli_core)

def solvers (model):
    import cobra.test
    model = cobra.test.create_test_model('textbook')
    model.solver = 'glpk'
    # or if you have cplex installed
    # model.solver = 'cplex'


def deYeast ():
    model = cobra.io.read_sbml_model('yeast_6.11_rui.xml')
    print('Number of reactions:', len(model.reactions))
    print('Number of metabolites:', len(model.metabolites))
    print('Number of genes:', len(model.genes))
    print('Gene number 34:', model.genes[34]) #retrieve item by index

    #FBA
    model.objective = "r_4041" # yeast 8 biomass pseudoreaction
    print(model.objective.expression)
    solution = model.optimize()
    print(solution.objective_value)
    print(solution.status)
    print(solution.fluxes)
    print(solution.shadow_prices)


if __name__ == '__main__':
    import cobra.test
    model = cobra.test.create_test_model("textbook")
    modelInfo(model)
    # reactions(model)
    # metabolites(model)
    # genes(model)
    # makingChangesReversibly(model)
    # buildModel()
    # readingModels()
    # simulatingFBA(model)
    # simulatingFVA(model)
    # simulatingpFBA(model)
    # simulatingDels()
    # productionEnvelopes(model)
    # solvers(model)







