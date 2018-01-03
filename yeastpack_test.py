#!/usr/bin/env python

'''
File for testing phenomenaly features

Author: Telma Afonso
'''

import cplex
from cobra.io.sbml import create_cobra_model_from_sbml_file
from phenomenaly.io import load_yeast_76
from phenomenaly.simulation import fba, fva, pfba, lmoma, simulate_auxotrophy, simulate_essentiality
from phenomenaly.variables import Media
from types import *
import pickle
import pandas as pd
import numpy as np
import sys, io
import re


class PhenomenalySim (object):

    def __init__(self, cobra_model = None):
        self.model = cobra_model
        self.medium = None

    def checkModelInfo (self, n = 5):
        ''' Prints model information for a quick overview.

            Parameters
            ----------
            model (cobra.core.Model object) : model from cobra.Model
            n (integer) : determines amount of information shown
        '''
        assert isinstance(n, int), 'n must be an integer'

        print('Model ID: ', self.model.id)
        print('Number of reactions:', len(self.model.reactions))
        print('Number of metabolites:', len(self.model.metabolites))
        print('Number of genes:', len(self.model.genes))
        try:
            print('First %i reactions:' % n)
            [print('\t', x.id, x.name) for x in self.model.reactions[:n]]
        except IndexError:
            print('First %i reactions:' % len(self.model.reactions))
            [print('\t', x.id, x.name) for x in self.model.reactions[:len(self.model.reactions)]]
        try:
            print('First %i metabolites:' % n)
            [print('\t', x.id, x.name) for x in self.model.metabolites[:n]]
        except IndexError:
            print('First %i metabolites:' % len(self.model.metabolites))
            [print('\t', x.id, x.name) for x in self.model.metabolites[:len(self.model.metabolites)]]
        try:
            print('First %i genes:' % n)
            [print('\t', x.id, x.name, end = " ") for x in self.model.genes[:n]]
        except IndexError:
            print('First %i genes:' % len(self.model.genes))
            [print('\t', x.id, x.name, end = " ") for x in self.model.genes[:len(self.model.genes)]]
        print('Compartments:')
        [print('\t', ID, c) for ID, c in self.model.compartments.items()]
        print('Yeastpack model description:')
        [print('\t', key, value) for key, value in self.model.describe().items()]

    def saveObjectToFile (self, object, filename):
        ''' Saves an object to a file for later use

            Parameters
            ----------
            object : a python object to save to file
            filename (str) : the name of the file to store the model
        '''
        assert isinstance(filename, str), 'filename must be a string'

        with open(filename, 'wb') as f:
            pickle.dump(object, f)

    def loadObjectFromFile (self, filename):
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

    def writeMetReactToFile (self, filename, type = 'metabolites'):
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
            [dict.update({x.id: x.name}) for x in self.model.metabolites[:len(self.model.metabolites)]]
        else:
            [dict.update({x.id: x.name}) for x in self.model.reactions[:len(self.model.reactions)]]

        with open(filename, 'w') as f:
            for key, val in sorted(dict.items()):
                f.write(key + '\t' + val + '\n')

    def translate_medium_modified (self, medium):
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
                    i, "Exchange", self.model.convert_exchange_name(Mod_name=row["Exchange"]))
            except:
                index.append(i)

        df.drop(df.index[index], inplace=True)

        return df

    def convertStdToSyst (self, gene_name_list):
        gs = pd.read_csv('genes_summary.csv', delim_whitespace = True)
        d = dict(zip(gs['symbol'], gs['input']))
        res = {}

        for g in gene_name_list:
            if g in d.keys():
                res[g] = d[g]
            else:
                print('No correspondence in this model for gene ' + g)
        return res

    def createResultsDataset (self, res_dict):
        res_dict = res_dict.copy()
        df = pd.DataFrame()
        wt = pd.DataFrame()

        for gene, r in sorted(res_dict.items()):
            if hasattr(r, 'x_dict'): #if legacy solution
                if gene != 'WildType':
                    df = pd.concat([df, pd.DataFrame(pd.DataFrame(list(r.x_dict.items())).set_index(0))], axis = 1)
                else:
                    wt = pd.concat([wt, pd.DataFrame(pd.DataFrame(list(r.x_dict.items())).set_index(0))], axis = 1)
                    wt.columns = [gene]
                    del res_dict[gene]
            else:
                if gene != 'WildType':
                    df = pd.concat([df, pd.DataFrame(r.fluxes)], axis = 1)
                else:
                    wt = pd.concat([wt, pd.DataFrame(r.fluxes)], axis = 1)
                    wt.columns = [gene]
                    del res_dict[gene]

        df.columns = sorted(res_dict.keys())

        return df, wt

    def createResultsDatasetLMOMA (self, res_dict):
        #Because it is of type LegacySolution (deprecated)
        res_dict = res_dict.copy()
        df = pd.DataFrame()

        for gene, r in sorted(res_dict.items()):
            res = pd.DataFrame(list(r.x_dict.items())).set_index([0])
            df = pd.concat([df, pd.DataFrame(res)], axis = 1)

        df.columns = sorted(res_dict.keys())

        return df

    def createResultsDatasetFVA (self, res_dict):
        res_dict = res_dict.copy()
        df = pd.DataFrame()
        wt = pd.DataFrame()

        for gene, r in sorted(res_dict.items()):
            if gene != 'WildType':
                df = pd.concat([df, r], axis = 1)
            else:
                wt = pd.concat([wt, r], axis = 1)
                wt.columns = [gene + '_maximum', gene + '_minimum']
                del res_dict[gene]

        keys = sorted(res_dict.keys())
        col_names = [[key + '_maximum', key + '_minimum'] for key in keys]
        df.columns = [item for sublist in col_names for item in sublist]

        return df, wt

    def convertKeggID (self, keggID):
        id_df = pd.read_csv('kegg_yeast_ids.csv', sep = '\t', index_col = False, header = None)
        id_df.columns = ['kegg', 'yeast', 'bool']

        id_dict = dict(zip(id_df.kegg, id_df.yeast))
        if keggID in id_dict.keys():
            return id_dict[keggID]
        else:
            print('No correspondence in this model for KEGG ID ' + keggID)

    def relativeError (self, actual, simulated):
        return (self.absoluteError(actual, simulated)/np.absolute(actual))*100

    def absoluteError (self, actual, simulated):
        return np.absolute(actual - simulated)

    def divide_list(self, l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def checkReaction (self, reaction_id):
        print('Reation (' + reaction_id + '): ' + self.model.reactions.get_by_id(reaction_id).name)
        r = self.model.reactions.get_by_id(reaction_id).reaction.split()
        r2 = [(self.model.metabolites.get_by_id(met).name if 's' in met else (met)) for met in r]
        print('\n', ' '.join(r2))

        # For plots
        s = ' '.join(r2)
        s = s.replace('cytoplasm', 'c')
        s = s.replace('mitochondrion', 'm')
        s = s.replace('nucleus', 'n')
        s = s.replace('extracellular', 'e')
        s = s.replace('peroxisome', 'p')

        return s

    def setMedium (self, medium):
        self.medium = medium
        self.model.set_medium(self.translate_medium_modified(medium))

    def printDict (self, dict):
        for key, val in dict.items():
            print(key, '\t', val)

    def correctReversedReactions (self, dataset, reactions = None):
        if reactions is None:
            reactions = ['r_0962', 'r_0300', 'r_1022', 'r_0454', 'r_1054', 'r_0452', 'r_0892', 'r_0893', 'r_1049']

        dat = dataset.copy()
        dat.update(dat.loc[reactions] * -1)

        return dat

    def checkReactionLB (self, reaction_id):
        print('Reation (' + reaction_id + '): ' + self.model.reactions.get_by_id(reaction_id).name)
        lb = self.model.reactions.get_by_id(reaction_id).lower_bound
        print('Lower bound: {}'.format(lb))

    def checkReactionUB (self, reaction_id):
        print('Reation (' + reaction_id + '): ' + self.model.reactions.get_by_id(reaction_id).name)
        ub = self.model.reactions.get_by_id(reaction_id).upper_bound
        print('Upper bound: {}'.format(ub))

    def singleSimulation(self, carbon_source = 'r_1714', cs_lb = -1.5, geneko = None, o2_lb = None, type = 'fba'):
        with self.model as m:
            m.set_carbon_source(carbon_source, lb = cs_lb)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb)
            if geneko is not None:
                if type == 'lmoma':
                    r = pfba(m)
                    m.set_environmental_conditions(gene_knockout = geneko)
                    res = lmoma(m, r)
                else:
                    m.set_environmental_conditions(gene_knockout = geneko)

            if type == 'fba':
                res = fba(m)
            elif type == 'pfba':
                res = pfba(m)
            elif type == 'fva':
                res = fva(m, reaction_list = m.reactions, fix_biomass = True)

        return res

    def getReactionInfo(self, reaction_id):
        r = self.model.reactions.get_by_id(reaction_id)
        for k,v in vars(r).items():
            print(k, '\t', v)

    def getSimFluxesRange (self, sim_res, min = None, max = None, dropReacts = []):
        if min is not None:
            if max is not None:
                df = sim_res.fluxes[(sim_res.fluxes > min) & (sim_res.fluxes < max)]
            else:
                df = sim_res.fluxes[(sim_res.fluxes > min)]
        elif max is not None:
            df = sim_res.fluxes[(sim_res.fluxes < max)]

        df = df.drop([r for r in dropReacts if r in df.index])
        row_names = [self.checkReaction(r) + '(' + str(r) + ')' for r in df.index]
        df.index = row_names

        return df.sort_values()

    def getMetaboliteSummary (self, metabolite_id, analysis_result):
        stdout = sys.stdout
        sys.stdout = io.StringIO()

        self.model.metabolites.get_by_id(metabolite_id).summary(analysis_result)

        # get output and restore sys.stdout
        output = sys.stdout.getvalue()
        sys.stdout = stdout

        return output

    def getMetaboliteName (self, metabolite_id, regex = True):
        if regex:
            metabolite_id = metabolite_id.group()
        return self.model.metabolites.get_by_id(metabolite_id).name

    def replaceMetaboliteIdsInString (self, string):
        return re.sub(r's_\w+', self.getMetaboliteName, string)

    def getMetaboliteSummaryWithNames (self, metabolite_id, analysis_result):
        output = self.replaceMetaboliteIdsInString(self.getMetaboliteSummary(metabolite_id, analysis_result))
        print(output)

    def getListOfMetabolitesSummary (self, analysis_result, metabolite_dict = None):
        if metabolite_dict is None:
            metabolite_dict = {'D-Glucose-6-phosphate (c)': 's_0568',
                               'Phosphoenolpyruvate (c)': 's_1360',
                               'Phosphoenolpyruvate (m)': 's_1361',
                               'Pyruvate (c)': 's_1399',
                               'Pyruvate (m)': 's_1401',
                               'Oxaloacetate (c)': 's_1271',
                               'Oxaloacetate (m)': 's_1273',
                               # 'Acetaldehyde (c)': 's_0359',
                               'Acetaldehyde (m)': 's_0361',
                               'Succinate (c)': 's_1458',
                               'Succinate (m)': 's_1460',
                               'Ethanol (c)': 's_0680',
                               'Ethanol (m)': 's_0682',
                               'Acetate (c)': 's_0362',
                               'Acetate (c)': 's_0365'
                               }

        for key, value in sorted(metabolite_dict.items()):
            print('='*20 + key.upper() + '='*20)
            output = self.replaceMetaboliteIdsInString(self.getMetaboliteSummary(value, analysis_result))
            print(output)


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

    #Reaction Bounds
    for r in model.reactions:
        print(r.id, '\t', r.lower_bound, '\t', r.upper_bound)


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

    model.set_medium(translate_medium_modified('MINIMAL', model))

    l = ['r_1654', 'r_2060', 'r_2020', 'r_2005', 'r_1861', 'r_2049', 'r_2020', 'r_1671', 'r_1967', 'r_1915',  'r_1947']
    for r in model.genes.get_by_id('YMR083W').reactions:
        print(r)
        #print(model.reactions.get_by_id('r_0187').reaction)
        #print(model.reactions.get_by_id(r).lower_bound)


    model.reactions.get_by_id('r_0713').reaction
    model.metabolites.get_by_id('s_1273').name

    model.genes.get_by_id('YMR083W').reactions


    # USING NEW CLASS
    ySim = YeastpackSim(loadObjectFromFile('model_yeast_76.sav'), 'MINIMAL')
    ySim.checkModelInfo()

