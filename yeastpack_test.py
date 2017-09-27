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

def saveObjectToFile (object, filename):
    ''' Saves an object to a file for later use

        Parameters
        ----------
        object : a python object to save to file
        filename (str) : the name of the file to store the model
    '''
    assert isinstance(filename, str), 'filename must be a string'

    with open(filename, 'wb') as f:
        pickle.dump(object, f)

def loadObjectFromFile (filename):
    ''' Loads an object object from a file

        Parameters
        ----------
        filename (str) : the name of the file where the object is contained

        Returns
        ----------
        object : the python object contained in the file
    '''
    assert isinstance(filename, str), 'filename must be a string'

    with open(filename, 'rb') as f:
        object = pickle.load(f)

    return object


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
        for key, val in sorted(dict.items()):
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

def convertStdToSyst (gene_name_list):
    gs = pd.read_csv('genes_summary.csv', delim_whitespace = True)
    d = dict(zip(gs['symbol'], gs['input']))
    res = {}

    for g in gene_name_list:
        if g in d.keys():
            res[g] = d[g]
        else:
            print('No correspondence in this model for gene ' + g)
    return res

def createResultsDataset (res_dict):
    df = pd.DataFrame()

    for gene, r in sorted(res_dict.items()):
        df = pd.concat([df, pd.DataFrame(r.fluxes)], axis = 1)

    df.columns = sorted(res_dict.keys())

    return df



if __name__ == '__main__':
    # dir(a) #check object attributes

    #Importing directly with cobra
    model_cobra = create_cobra_model_from_sbml_file('.\DATA\Models\Yeast_7\yeast_7.6_original\yeast_7.6_cobra.xml')

    #Importing using yeastpack
    model = load_yeast_76()
    checkModelInfo(model, n = 10)

    # Save model to file for later use
    saveObjectToFile(model, 'model_yeast_76.sav')

    # Load the model from file
    model = loadObjectFromFile('model_yeast_76.sav')
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


    # =======================================
    #                  CASE 7
    # =======================================

    d = {'ADH3': -14.4, 'ALD5': -16.2, 'ALD6': -4.7, 'COX5A': -13.3, 'CTP1': -13.9, 'DAL7': -14.3,
         'FUM1': -9.1, 'GND2': -15.5, 'GCV2': -16, 'GLY1': -13.1, 'GPD1': -16.4, 'ICL1': -15.7,
         'IDP1': -15, 'IDP2': -13.4, 'LSC1': -17.8, 'MAE1': -12.5, 'MDH1': -10, 'MHD2': -13,
         'MDH3': -14.4, 'MSL1': -17.5, 'OAC1': -11.7, 'PCK1': -13.4, 'PDA1': -8, 'PGM1': -14.7,
         'PGM2': -16.1, 'RPE1': -4.6, 'SDH1': -13.4, 'SER33': -14.5, 'SFC1': -12.9, 'SOL1': -14.4,
         'SOL2': -15.6, 'SOL3': -12.2, 'SOL4': -15.9, 'TAL1': -12.2, 'YGR043C': -16, 'ZWF1': -6.5}

    l = convertStdToSyst(d.keys())              #Dict with gene : yeastpack gene
    l_inv = {v: k for k, v in l.items()}        #Dict with yeastpack gene : gene
    g_lb = {v: d[k] for k, v in l.items()}      #Dict with yeastpack gene : lb

    model = loadObjectFromFile('model_yeast_76.sav')
    model.set_medium(translate_medium_modified('MINIMAL_CASE7', model))

    # ======================
    #           FBA
    # ======================
    res = {}
    for g, lb in g_lb.items():
        print('======= ' + l_inv[g] + ' (' + g + ')' ' =======')
        with model as model:
            model.set_carbon_source('r_1714', lb = lb)
            model.set_environmental_conditions(gene_knockout = g)
            key = ''.join(l_inv[g] + ' (%s)' % g)
            res[key] = fba(model)
            print('Objective value:', res[key].objective_value, '\n')

    #Save results to file
    saveObjectToFile(res, 'res_fba_case7.sav')

    #Load results from file
    result = loadObjectFromFile('res_fba_case7.sav')

    #Create dataset with results
    res_df = createResultsDataset(res)

    #Save dataset to file
    res_df.to_csv('res_fba_case7.csv', sep = ';')


    # ======================
    #          pFBA
    # ======================
    res_pfba = {}
    for g, lb in g_lb.items():
        print('======= ' + l_inv[g] + ' (' + g + ')' ' =======')
        with model as model:
            model.set_carbon_source('r_1714', lb = lb)
            model.set_environmental_conditions(gene_knockout = g)
            key = ''.join(l_inv[g] + ' (%s)' % g)
            res_pfba[key] = pfba(model)
            print('Objective value:', res_pfba[key].objective_value, '\n')

    #Save results to file
    saveObjectToFile(res_pfba, 'res_pfba_case7.sav')

    #Load results from file
    result_pfba = loadObjectFromFile('res_pfba_case7.sav')

    #Create dataset with results
    res_pfba__df = createResultsDataset(res_pfba)

    #Save dataset to file
    res_pfba__df.to_csv('res_pfba_case7.csv', sep = ';')

    #res['ADH3 (YMR083W)'].fluxes - result_pfba['ADH3 (YMR083W)'].fluxes


    # ======================
    #          FVA
    # ======================
    res_fva = {}
    for g, lb in g_lb.items():
        print('======= ' + l_inv[g] + ' (' + g + ')' ' =======')
        with model as model:
            try:
                model.set_carbon_source('r_1714', lb = lb)
                model.set_environmental_conditions(gene_knockout = g)
                key = ''.join(l_inv[g] + ' (%s)' % g)
                res_fva[key] = fva(model, reaction_list = model.reactions, fix_biomass = False)
                print('Done!', '\n')
            except:
                print('FVA solution status infeasible for case ' + key + '\n')



















