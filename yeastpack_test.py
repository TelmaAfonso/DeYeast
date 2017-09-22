#!/usr/bin/env python

'''
File for testing yeastpack features

Author: Telma Afonso
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file
from yeastpack.io import load_yeast_76
from yeastpack.simulation import fba, fva, pfba, lmoma, simulate_auxotrophy, simulate_essentiality
from yeastpack.data import Media
from types import *
import pickle
import pandas as pd

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


    def translate_medium_modified (medium, model):
        ''' Modification of translate_medium function from yeastpack.utils

        Parameters
        ----------
        medium (str) : the name of the medium to use
        model (cobra.core.Model object) : model from cobra.Model

        Returns
        ----------
        df (pandas dataframe) : a dataframe with the exchange reactions and respective lower/upper bounds
        '''
        media = [m for m in dir(Media) if m[0] is not '_']
        assert medium in media, 'Please provide one of the following as medium:' + str(media)
        medium_filename = getattr(Media, medium)

        df = pd.read_csv(medium_filename, sep=";")
        index = []
        # Translate each exchange in medium. If not present in model, remove from medium
        for i, row in df.iterrows():
            try:
                df.set_value(
                    i, "Exchange", model.convert_exchange_name(Mod_name=row["Exchange"]))
            except:
                index.append(i)

        df.drop(df.index[index], inplace=True)

        return df



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
    res_pfba.objective_value
    res_pfba.status
    res_pfba.fluxes
    res_pfba.reduced_costs
    res_pfba.shadow_prices

    #FVA
    res_fva = fva(model, reaction_list = model.reactions) #pandas dataframe
    res_fva.to_csv('fva_res.csv', sep = '\t')


    #RANDOM TESTS
    model.get_model_id()
    model.get_biomass_id()
    model.get_biomass()
    model.get_exchanges_ids()
    model.get_genes_ids()
    model.get_minimum_lower_bound()
    model.get_maximum_upper_bound()

    writeMetReactToFile('reactions.csv', type = 'reactions')

    model.set_medium(translate_medium_modified('MINIMAL_CASE7', model))

    l = ['r_1654', 'r_2060', 'r_2020', 'r_2005', 'r_1861', 'r_2049', 'r_2020', 'r_1671', 'r_1967', 'r_1915',  'r_1947']
    for r in l:
        print(model.reactions.get_by_id(r).name)
        print(model.reactions.get_by_id(r).lower_bound)





