#!/usr/bin/env python

'''
File for testing yeastpack features

Author: Telma Afonso
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file
from yeastpack.io import load_yeast_76

def checkModelInfo (model, n = 5):
    '''
        Simple function that prints model information for a quick overview.

        model (cobra.core.Model object) : model from cobra.Model
        n (integer) : determines amount of information shown

    '''
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


if __name__ == '__main__':
    # dir(a) #check object attributes

    #Importing directly with cobra
    model_cobra = create_cobra_model_from_sbml_file('D:\Dropbox\Trabalhos\\16.17\DeYeastLib Project\yeastpack\src\yeastpack\DATA\Models\Yeast_7\yeast_7.6_original\yeast_7.6_cobra.xml')

    #Importing using yeastpack
    model = load_yeast_76()
    checkModelInfo(model, n = 10)




