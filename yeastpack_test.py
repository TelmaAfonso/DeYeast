#!/usr/bin/env python

'''
File for testing yeastpack features

Author: Telma Afonso
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file
from yeastpack.io import load_yeast_76
from types import *
import pickle

def checkModelInfo (model, n = 5):
    ''' Prints model information for a quick overview.

        Parameters
        ----------
        model (cobra.core.Model object) : model from cobra.Model
        n (integer) : determines amount of information shown
    '''
    assert type(n) == int, "%r is not an integer" % n

    print('Model ID: ', model.id)
    print('Number of reactions:', len(model.reactions))
    print('Number of metabolites:', len(model.metabolites))
    print('Number of genes:', len(model.genes))
    try:
        print('First %i reactions:' % n)
        [print('\t', x.id, x.reaction) for x in model.reactions[:n]]
    except IndexError:
        print('First %i reactions:' % len(model.reactions))
        [print('\t', x.id, x.reaction) for x in model.reactions[:len(model.reactions)]]
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


def saveModelToFile (model, filename):
    ''' Saves a model object to a file for later use

        Parameters
        ----------
        model (cobra.core.Model object) : model from cobra.Model
        filename (str) : the name of the file to store the model
    '''
    assert type(filename) == str, "%r is not a string" % filename

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
    assert type(filename) == str, "%r is not a string" % filename

    with open(filename, 'rb') as f:
        model = pickle.load(f)

    return model



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




