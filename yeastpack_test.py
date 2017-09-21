#!/usr/bin/env python

'''
File for testing yeastpack features

Author: Telma Afonso
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file
from yeastpack.io import load_yeast_76
from yeastpack.simulation import fba, fva, pfba, lmoma, simulate_auxotrophy, simulate_essentiality
from types import *
import pickle

def checkModelInfo (model, n = 5):
    ''' Prints model information for a quick overview.

        Parameters
        ----------
        model (cobra.core.Model object) : model from cobra.Model
        n (integer) : determines amount of information shown
    '''
    assert isinstance(n, int), 'n must be an integer'

    print('Model ID: ', model.id)
    print('Number of reactions:', len(model.reactions))
    print('Number of metabolites:', len(model.metabolites))
    print('Number of genes:', len(model.genes))
    try:
        print('First %i reactions:' % n)
        [print('\t', x.id, x.name) for x in model.reactions[:n]]
    except IndexError:
        print('First %i reactions:' % len(model.reactions))
        [print('\t', x.id, x.name) for x in model.reactions[:len(model.reactions)]]
    try:
        print('First %i metabolites:' % n)
        [print('\t', x.id, x.name) for x in model.metabolites[:n]]
    except IndexError:
        print('First %i metabolites:' % len(model.metabolites))
        [print('\t', x.id, x.name) for x in model.metabolites[:len(model.metabolites)]]
    try:
        print('First %i genes:' % n)
        [print('\t', x.id, x.name, end = " ") for x in model.genes[:n]]
    except IndexError:
        print('First %i genes:' % len(model.genes))
        [print('\t', x.id, x.name, end = " ") for x in model.genes[:len(model.genes)]]
    print('Compartments:')
    [print('\t', ID, c) for ID, c in model.compartments.items()]
    print('Yeastpack model description:')
    [print('\t', key, value) for key, value in model.describe().items()]

def saveModelToFile (model, filename):
    ''' Saves a model object to a file for later use

        Parameters
        ----------
        model (cobra.core.Model object) : model from cobra.Model
        filename (str) : the name of the file to store the model
    '''
    assert isinstance(filename, str), 'filename must be a string'

    with open(filename, 'wb') as f:
        pickle.dump(model, f)

def loadModelFromFile (filename):
    ''' Loads a model object from a file

        Parameters
        ----------
        filename (str) : the name of the file where the model is contained

        Returns
        ----------
        model (cobra.core.Model object) : model from cobra.Model
    '''
    assert isinstance(filename, str), 'filename must be a string'

    with open(filename, 'rb') as f:
        model = pickle.load(f)

    return model


def writeMetReactToFile (filename, type = 'metabolites'):
    ''' Writes the list of metabolites/reactions to a file

    Parameters
    ----------
    filename (str) : the name of the file to be created
    type (str): whether to write the reactions or metabolites to the file
    '''
    assert isinstance(filename, str), 'filename must be a string'
    assert type in ['metabolites', 'reactions'], "type must be either 'metabolites' or 'reactions'"

    dict = {}
    if type == 'metabolites':
        [dict.update({x.id: x.name}) for x in model.metabolites[:len(model.metabolites)]]
    else:
        [dict.update({x.id: x.name}) for x in model.reactions[:len(model.reactions)]]

    with open(filename, 'w') as f:
        for key, val in dict.items():
            f.write(key + '\t' + val + '\n')



if __name__ == '__main__':
    # dir(a) #check object attributes

    #Importing directly with cobra
    model_cobra = create_cobra_model_from_sbml_file('.\DATA\Models\Yeast_7\yeast_7.6_original\yeast_7.6_cobra.xml')

    #Importing using yeastpack
    model = load_yeast_76()
    checkModelInfo(model, n = 10)

    # Save model to file for later use
    saveModelToFile(model, 'model_yeast_76.sav')

    # Load the model from file
    model = loadModelFromFile('model_yeast_76.sav')
    model.summary()

    model.metabolites.s_1203.summary() #examine the overall redox balance of the model
    model.metabolites.s_0434.summary() #main energy production and consumption reactions

    #SIMULATIONS
    #FBA
    res_fba = fba(model)
    res_fba.objective_value     # The (optimal) value for the objective function.
    res_fba.status              # The solver status related to the solution.
    res_fba.fluxes              # Contains the reaction fluxes (primal values of variables).
    res_fba.reduced_costs       # Contains reaction reduced costs (dual values of variables).
    res_fba.shadow_prices       # Contains metabolite shadow prices (dual values of constraints).

    #PFBA
    res_pfba = pfba(model)
    res_pfba.f                  # The objective value
    res_pfba.solver             # A string indicating which solver package was used.
    res_pfba.x                  # List or Array of the fluxes (primal values).
    res_pfba.x_dict             # A dictionary of reaction IDs that maps to the respective primal values.
    res_pfba.y                  # List or Array of the dual values.
    res_pfba.y_dict             # A dictionary of reaction IDs that maps to the respective dual values.

    #FVA
    res_fva = fva(model, reaction_list = model.reactions) #pandas dataframe
    res_fva.to_csv('fva_res.csv', sep = '\t')

    
    writeMetReactToFile('reactions.csv', type = 'reactions')






