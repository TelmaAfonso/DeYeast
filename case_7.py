#!/usr/bin/env python

'''
File case 7 simulations

Author: Telma Afonso
'''

from yeastpack.simulation import fba, fva, pfba, lmoma
from yeastpack.data import Media
from types import *
import pickle
import string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.backends.backend_pdf import PdfPages

from yeastpack_test import YeastpackSim


class Case7 (YeastpackSim):

    def __int__(self, cobra_model = None):
        super(Case7, self).__int__(self, cobra_model)
        self.fba_res = None
        self.pfba_res = None
        self.fva_res = None
        self.dictsForCase7() #Not working for some reason

    def dictsForCase7 (self):
        # Glucose lbs
        d = {'ADH3': -14.4, 'ALD5': -16.2, 'ALD6': -4.7, 'COX5A': -13.3, 'CTP1': -13.9, 'DAL7': -14.3,
             'FUM1': -9.1, 'GND2': -15.5, 'GCV2': -16, 'GLY1': -13.1, 'GPD1': -16.4, 'ICL1': -15.7,
             'IDP1': -15, 'IDP2': -13.4, 'LSC1': -17.8, 'MAE1': -12.5, 'MDH1': -10, 'MHD2': -13,
             'MDH3': -14.4, 'MSL1': -17.5, 'OAC1': -11.7, 'PCK1': -13.4, 'PDA1': -8, 'PGM1': -14.7,
             'PGM2': -16.1, 'RPE1': -4.6, 'SDH1': -13.4, 'SER33': -14.5, 'SFC1': -12.9, 'SOL1': -14.4,
             'SOL2': -15.6, 'SOL3': -12.2, 'SOL4': -15.9, 'TAL1': -12.2, 'YGR043C': -16, 'ZWF1': -6.5}

        exp_biomass = {'ADH3':0.3716, 'ALD5': 0.4126, 'ALD6': 0.1387, 'COX5A': 0.256, 'CTP1': 0.37, 'DAL7': 0.35,
                       'FUM1': 0.2097, 'GND2': 0.338, 'GCV2': 0.3745, 'GLY1': 0.318, 'GPD1': 0.34, 'ICL1': 0.386,
                       'IDP1': 0.372, 'IDP2': 0.356, 'LSC1': 0.423, 'MAE1': 0.411, 'MDH1': 0.29, 'MHD2': 0.36,
                       'MDH3': 0.405, 'MSL1': 0.4, 'OAC1': 0.29, 'PCK1': 0.4, 'PDA1': 0.17, 'PGM1': 0.37,
                       'PGM2': 0.41, 'RPE1': 0.1316, 'SDH1': 0.293, 'SER33': 0.33, 'SFC1': 0.339, 'SOL1': 0.37,
                       'SOL2': 0.402, 'SOL3': 0.3557, 'SOL4': 0.3843, 'TAL1': 0.361, 'YGR043C': 0.3745, 'ZWF1': 0.16}

        l = self.convertStdToSyst(d.keys())         #Dict with gene : yeastpack gene
        l_inv = {v: k for k, v in l.items()}        #Dict with yeastpack gene : gene
        g_lb = {v: d[k] for k, v in l.items()}      #Dict with yeastpack gene : lb
        self.d, self.l, self.l_inv, self.g_lb, self.exp_biomass = d, l, l_inv, g_lb, exp_biomass

    def case7fba (self, o2_lb = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res = {}
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                m.set_carbon_source('r_1714', lb = lb)
                m.set_environmental_conditions(gene_knockout = g)
                if o2_lb is not None:
                    m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                key = ''.join(self.l_inv[g] + ' (%s)' % g)
                res[key] = fba(m)
                print('Objective value:', res[key].objective_value, '\n')
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7) # GLC WT = -16.7
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res['WildType'] = fba(m)
            print('Objective value:', res['WildType'].objective_value, '\n')
        return res

    def case7pfba (self, o2_lb = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res_pfba = {}
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                m.set_carbon_source('r_1714', lb = lb)
                m.set_environmental_conditions(gene_knockout = g)
                if o2_lb is not None:
                    m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                key = ''.join(self.l_inv[g] + ' (%s)' % g)
                res_pfba[key] = pfba(m)
                print('Objective value:', res_pfba[key].objective_value, '\n')
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res_pfba['WildType'] = pfba(m)
            print('Objective value:', res_pfba['WildType'].objective_value, '\n')
        return res_pfba

    def case7lmoma (self, o2_lb = None, reference_dict = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res_lmoma = {}
        i = 1
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                m.set_carbon_source('r_1714', lb = lb)
                m.set_environmental_conditions(gene_knockout = g)
                if o2_lb is not None:
                    m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                key = ''.join(self.l_inv[g] + ' (%s)' % g)
                res_lmoma[key] = lmoma(m, reference_dict[key])
                print('Objective value:', res_lmoma[key].f, '\n')
                print('(' + (i/len(self.g_lb.keys())) * 100 + '% Complete)')
                i += 1
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res_lmoma['WildType'] = lmoma(m, reference_dict['WildType'])
            print('Objective value:', res_lmoma['WildType'].f, '\n')
        return res_lmoma

    def case7fva (self, o2_lb = None):
        # g_lb - dict with yeastpack gene : lb
        # l_inv - dict with yeastpack gene : gene
        res_fva = {}
        for g, lb in self.g_lb.items():
            print('======= ' + self.l_inv[g] + ' (' + g + ')' ' =======')
            with self.model as m:
                try:
                    m.set_carbon_source('r_1714', lb = lb)
                    m.set_environmental_conditions(gene_knockout = g)
                    if o2_lb is not None:
                        m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb[self.l_inv[g]]) # oxygen exchange
                    key = ''.join(self.l_inv[g] + ' (%s)' % g)
                    res_fva[key] = fva(m, reaction_list = m.reactions, fix_biomass = True)
                    print('Done!', '\n')
                except:
                    print('FVA solution status infeasible for case ' + key + '\n')
        print('======= ' + 'Wild Type' + ' =======')
        with self.model as m:
            m.set_carbon_source('r_1714', lb = -16.7)
            if o2_lb is not None:
                m.reactions.get_by_id('r_1992').lower_bound = float(o2_lb['WildType'])
            res_fva['WildType'] = fva(m, reaction_list = m.reactions, fix_biomass = True)
            print('Done!', '\n')
        return res_fva

    def readExpFluxes (self, filename, sep = ';'):
        exp_res = pd.read_csv(filename, sep = sep)
        col_names = list(exp_res.columns)

        genes = [gene for gene in col_names if len(gene) <= 5]

        for i in range(0, len(col_names)):
            if 'Real' in col_names[i]:
                col_names[i] = col_names[i-1] + '_real_flux'

        exp_res.columns = col_names

        return exp_res, genes

    def loadExperimentalRes (self, filename_exp):
        exp_fluxes, genes = self.readExpFluxes(filename_exp)
        exp_fluxes.set_index('Yeast7_ID', inplace = True)
        reactions = exp_fluxes.iloc[:,0]
        sub_exp_fluxes = exp_fluxes.drop([col for col in list(exp_fluxes.columns) if '_real_flux' not in col], 1)
        sub_exp_fluxes = sub_exp_fluxes.apply(lambda x: x.str.replace(',', '.')).astype('float')

        return sub_exp_fluxes, reactions

    def createDatasetExpVsSimul (self, dataset_exp, dataset_sim):
        # For FBA and pFBA
        dsim = dataset_sim.copy()
        dsim.columns = [n.split()[0] + '_sim_flux' for n in list(dsim.columns)]
        df = pd.concat([dataset_exp, dsim], axis = 1, join = 'inner')
        df = df[list(sum(zip(dataset_exp.columns, dsim.columns), ()))]

        return df

    def createDatasetExpVsSimulFVA (self, dataset_exp, dataset_sim):
        dsim = dataset_sim.copy()
        dsim.columns = [(i.split()[0] + '_minimum' if count % 2 == 1 else (i.split()[0] + '_maximum')) for count, i in enumerate(list(dsim.columns))]
        df = pd.concat([dataset_exp, dsim], axis = 1, join = 'inner')
        df = df.reindex_axis(sorted(df.columns), axis = 1)

        return df

    def createDatasetWithAbsRelError (self, dataset_exp_vs_sim):
        df = dataset_exp_vs_sim.copy()
        for ind in range(len(list(dataset_exp_vs_sim.columns)), 0,  -2):
            ae = self.absoluteError(df.ix[:, ind - 2], df.ix[:, ind - 1])
            colname_ae = df.ix[:, ind - 2].name.split('_')[0] + '_abs_error'
            df.insert(loc = ind, column = colname_ae, value = ae)

            re = self.relativeError(df.ix[:, ind - 2], df.ix[:, ind - 1])
            colname_re = df.ix[:, ind - 2].name.split('_')[0] + '_rel_error'
            df.insert(loc = ind + 1, column = colname_re, value = re)

        return df

    def getBiomassObj (self, res_dict):
        biom = {}

        for gene, res in res_dict.items():
            biom[gene] = res.objective_value

        return  biom

    def getEthanolFux (self, dataset_exp, reaction_id = 'r_2115'):
        df = dataset_exp.copy()
        # sub_df = df.filter(regex = '_real_flux')
        df.columns = [col.split('_')[0] for col in list(df.columns)]
        # df.columns = list(sorted(l_inv.keys()))
        return dict(df.ix[reaction_id, :])

    def testO2EthanolProd (self, do_fba = True, g_knockout = None, gluc_lb = -10, range_o2 = list(np.arange(-20, 0, 2))):
        res = {}
        for i in range_o2:
            with self.model as m:
                m.set_carbon_source('r_1714', lb = gluc_lb)                   # glucose
                m.reactions.get_by_id('r_1992').lower_bound = float(i)        # oxygen exchange
                m.set_environmental_conditions(gene_knockout = g_knockout)
                if do_fba: r = fba(m)
                else: r = pfba(m)
                res[str(i)] = r.fluxes['r_2115']
        for key, val in sorted(res.items()): print(key, '\t', val)
        return res

    def testAllO2EthanolProd (self, range_o2 = list(np.arange(-20, 0, 2)), do_fba = False):
        res_EtOH_fuxes = {}
        genes = sorted(list(self.l_inv.values()))
        for gene in genes:
            print('Gene ' + gene + ':')
            res_EtOH_fuxes[gene] = self.testO2EthanolProd(g_knockout = self.l[gene], gluc_lb = self.d[gene], range_o2 = range_o2, do_fba = do_fba)
            print('Done!')
        print('Wild Type:')
        res_EtOH_fuxes['WildType'] = self.testO2EthanolProd(g_knockout = None, gluc_lb = -16.7, range_o2 = range_o2, do_fba = do_fba)
        print('Done!')
        return res_EtOH_fuxes

    def fixEtO2FluxesForPlotting (self, dict_EtOH_fuxes):
        # Remove O2 fluxes that yield no EtOH (for plot slope reasons)
        res = dict_EtOH_fuxes.copy()
        for key in list(dict_EtOH_fuxes.keys()):
            for k in list(dict_EtOH_fuxes[key].keys()):
                if res[key][k] == 0:
                    del res[key][k]
        return res

    def getO2flux (self, dict_EtOH_fluxes, dict_real_EtOH_fluxes):
        res = {}
        for gene in dict_EtOH_fluxes.keys():
            x = sorted([int(x) for x in dict_EtOH_fluxes[gene].keys()])
            y = [dict_EtOH_fluxes[gene][str(key)] for key in x]
            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            real_O2 = lambda x0: (y0 - intercept) / slope
            y0 = dict_real_EtOH_fluxes[gene]
            res[gene] = format(real_O2(y0), '.2f')

        return res

    def correctReversedReactions (self, dataset, reactions = None):
        if reactions is None:
            reactions = ['r_0962', 'r_0300', 'r_1022', 'r_1054', 'r_0452', 'r_0892', 'r_0893', 'r_1049', 'r_1048']

        dat = dataset.copy()
        dat.update(dat.loc[reactions] * -1)

        return dat

    def convertPandasDFToExcel (self, reactions, dataframe,  title = 'FBA Results Case 7', filename = 'Results', imgFolder = None, type = 'fba'):
        df = pd.concat([reactions, dataframe], axis = 1)
        writer = pd.ExcelWriter(filename, engine = 'xlsxwriter') # Create a Pandas Excel writer using XlsxWriter as the engine.
        df.to_excel(writer, sheet_name = 'Results', startrow = 2) # Convert the dataframe to an XlsxWriter Excel object.
        workbook  = writer.book # Get the xlsxwriter objects from the dataframe writer object.
        worksheet = writer.sheets['Results'] # Equivalent to xlsxwriter.Workbook('filename.xlsx').add_worksheet()

        # Formatting
        header_format = workbook.add_format({'bold': True, 'text_wrap': False, 'valign': 'vcenter', 'fg_color': '#d2d7bb', 'border': 0, 'align': 'center'})
        header_format2 = workbook.add_format({'bold': True, 'text_wrap': False, 'valign': 'vcenter', 'fg_color': '#4e7792', 'border': 0, 'font_color': '#FFFFFF', 'align': 'center'})
        num_format = workbook.add_format({'num_format': '#,##0.0000', 'valign': 'center', 'font_size': 10, 'fg_color': '#FFFFFF'})
        title_format = workbook.add_format({'bold': True, 'text_wrap': False, 'valign': 'vcenter', 'font_size': 12})
        reactions_format = workbook.add_format({'text_wrap': False, 'valign': 'right', 'font_size': 8, 'fg_color': '#FFFFFF'})
        fva_format = workbook.add_format({'num_format': '#,##0.0000', 'valign': 'vcenter', 'font_size': 10, 'fg_color': '#fbf9f1'})

        worksheet.write(0, 1, title, title_format) #Title
        worksheet.write(2, 0, 'Yeast7_ID', header_format2)
        for col_num, value in enumerate(list(df.columns.values)):
            worksheet.write(2, col_num + 1, value, header_format2) # Write the column headers with the defined format.
        for row_num, value in enumerate(df.index):
            worksheet.write(row_num + 3, 0, value, header_format)
        worksheet.set_column(2, df.shape[1], width = 15, cell_format = num_format) # Set column width
        worksheet.set_row(2, height = 20)

        # Conditional Formatting
        inds = [i for i in range(4, df.shape[1], 4)] #Columns to colour
        for i in inds:
            worksheet.conditional_format(first_row = 3, last_row = df.shape[0] + 3, first_col = i, last_col = i + 1,
                                         options = {'type': '3_color_scale', 'min_color': '#63BE7B',  'mid_color': '#FFEB84', 'max_color': '#F8696B'})

        worksheet.set_column(1, 1, width = 40, cell_format = reactions_format)
        worksheet.set_column(0, 0, width = 11)

        #Add Images
        worksheet.write(35, 1, 'Experimental vs Simulated Fluxes (Reactions)', title_format)
        worksheet.insert_image(36, 0, imgFolder + '/0_reacts.png', {'x_scale': 0.7, 'y_scale': 0.7})
        col = 5
        for i in range(1, 11):
            worksheet.insert_image(36, col, imgFolder + '/' + str(i) + '_reacts.png', {'x_scale': 0.7, 'y_scale': 0.7})
            col += 6

        # ============== Add case specific sheets ==============
        l = [value for col_num, value in enumerate(list(dataframe.columns.values))]
        cases = [vals for vals in iter(self.divide_list(l, 4))]

        for ind, gene_l in enumerate(cases):
            gene = gene_l[0].split('_')[0]
            sub_df = pd.concat([reactions, dataframe[gene_l]], axis = 1)
            sub_df.to_excel(writer, sheet_name = gene, startrow = 2)
            setattr(self, 'worksheet' + str(ind), writer.sheets[gene])
            getattr(self, 'worksheet' + str(ind)).write(0, 1, title + ' [' + gene + ']', title_format)


            getattr(self, 'worksheet' + str(ind)).write(2, 0, 'Yeast7_ID', header_format2)
            for col_num, value in enumerate(list(sub_df.columns.values)):
                getattr(self, 'worksheet' + str(ind)).write(2, col_num + 1, value, header_format2) # Write the column headers with the defined format.
            for row_num, value in enumerate(sub_df.index):
                getattr(self, 'worksheet' + str(ind)).write(row_num + 3, 0, value, header_format)
            getattr(self, 'worksheet' + str(ind)).set_column(2, sub_df.shape[1], width = 15, cell_format = num_format) # Set column width

            # Conditional Formatting
            inds = [i for i in range(4, sub_df.shape[1])] #Columns to colour
            for i in inds:
                getattr(self, 'worksheet' + str(ind)).conditional_format(first_row = 3, last_row = sub_df.shape[0] + 3, first_col = i, last_col = i + 1,
                                                                         options = {'type': '3_color_scale', 'min_color': '#63BE7B',  'mid_color': '#FFEB84', 'max_color': '#F8696B'})
            #Add Images
            if type == 'pfba' or type == 'lmoma':
                getattr(self, 'worksheet' + str(ind)).insert_image('G3', 'Results/Case 7/EtOH_figs/' + gene + '_etOH_pFBA.png', {'x_scale': 0.7, 'y_scale': 0.7})
            elif type == 'fba':
                getattr(self, 'worksheet' + str(ind)).insert_image('G3', 'Results/Case 7/EtOH_figs/' + gene + '_etOH.png', {'x_scale': 0.7, 'y_scale': 0.7})

            getattr(self, 'worksheet' + str(ind)).insert_image('G17', imgFolder + '/' + gene + '_genes.png', {'x_scale': 0.7, 'y_scale': 0.7})

            getattr(self, 'worksheet' + str(ind)).set_column(1, 1, width = 40, cell_format = reactions_format)
            getattr(self, 'worksheet' + str(ind)).set_column(0, 0, width = 11)
            getattr(self, 'worksheet' + str(ind)).set_row(2, height = 20)

        workbook.close()

    def convertPandasDFToExcelFVA (self, reactions, dataframe, title = 'FVA Results Case 7', filename = 'Results', imgFolder = None):
        df = pd.concat([reactions, dataframe], axis = 1)
        writer = pd.ExcelWriter(filename, engine = 'xlsxwriter') # Create a Pandas Excel writer using XlsxWriter as the engine.
        df.to_excel(writer, sheet_name = 'Results', startrow = 2) # Convert the dataframe to an XlsxWriter Excel object.
        workbook  = writer.book # Get the xlsxwriter objects from the dataframe writer object.
        worksheet = writer.sheets['Results'] # Equivalent to xlsxwriter.Workbook('filename.xlsx').add_worksheet()

        # Formatting
        header_format = workbook.add_format({'bold': True, 'text_wrap': False, 'valign': 'vcenter', 'fg_color': '#d2d7bb', 'border': 0, 'align': 'center'})
        header_format2 = workbook.add_format({'bold': True, 'text_wrap': False, 'valign': 'vcenter', 'fg_color': '#4e7792', 'border': 0, 'font_color': '#FFFFFF', 'align': 'center'})
        num_format = workbook.add_format({'num_format': '#,##0.0000', 'valign': 'center', 'font_size': 10, 'fg_color': '#FFFFFF'})
        title_format = workbook.add_format({'bold': True, 'text_wrap': False, 'valign': 'vcenter', 'font_size': 12})
        reactions_format = workbook.add_format({'text_wrap': False, 'valign': 'right', 'font_size': 8, 'fg_color': '#FFFFFF'})
        fva_format = workbook.add_format({'num_format': '#,##0.0000', 'valign': 'vcenter', 'font_size': 10, 'fg_color': '#fbf9f1', 'align': 'center'})

        worksheet.write(0, 1, title, title_format) #Title
        worksheet.write(2, 0, 'Yeast7_ID', header_format2)
        for col_num, value in enumerate(list(df.columns.values)):
            worksheet.write(2, col_num + 1, value, header_format2) # Write the column headers with the defined format.
        for row_num, value in enumerate(df.index):
            worksheet.write(row_num + 3, 0, value, header_format)
        worksheet.set_column(2, df.shape[1], width = 15, cell_format = num_format) # Set column width
        worksheet.set_row(2, height = 20)

        # Conditional Formatting
        inds = [i for i in range(4, df.shape[1], 3)]
        for i in inds:
            worksheet.set_column(i, i, width = 15, cell_format = fva_format)

        worksheet.set_column(1, 1, width = 50, cell_format = reactions_format)
        worksheet.set_column(0, 0, width = 11)

        # ============== Add case specific sheets ==============
        l = [value for col_num, value in enumerate(list(dataframe.columns.values))]
        cases = [vals for vals in iter(self.divide_list(l, 3))]

        for ind, gene_l in enumerate(cases):
            gene = gene_l[0].split('_')[0]
            sub_df = pd.concat([reactions, dataframe[gene_l]], axis = 1)
            sub_df.to_excel(writer, sheet_name = gene, startrow = 2)
            setattr(self, 'worksheet' + str(ind), writer.sheets[gene])
            getattr(self, 'worksheet' + str(ind)).write(0, 1, title + ' [' + gene + ']', title_format)


            getattr(self, 'worksheet' + str(ind)).write(2, 0, 'Yeast7_ID', header_format2)
            for col_num, value in enumerate(list(sub_df.columns.values)):
                getattr(self, 'worksheet' + str(ind)).write(2, col_num + 1, value, header_format2) # Write the column headers with the defined format.
            for row_num, value in enumerate(sub_df.index):
                getattr(self, 'worksheet' + str(ind)).write(row_num + 3, 0, value, header_format)
            getattr(self, 'worksheet' + str(ind)).set_column(2, sub_df.shape[1], width = 15, cell_format = num_format) # Set column width

            # Conditional Formatting
            getattr(self, 'worksheet' + str(ind)).set_column(4, 4, width = 15, cell_format = fva_format)

            # getattr(self, 'worksheet' + str(ind)).conditional_format('E4:E33', {'type':     'formula',
            #                                                                     'criteria': '=AND(E4 >= D4, E4 <= C4)',
            #                                                                     'format':   fva_format})
            getattr(self, 'worksheet' + str(ind)).set_column(1, 1, width = 60, cell_format = reactions_format)
            getattr(self, 'worksheet' + str(ind)).set_column(0, 0, width = 11)
            getattr(self, 'worksheet' + str(ind)).set_row(2, height = 20)

        workbook.close()

    # PLOTS

    def plotO2vsEtOH (self, dict_gene_EtOH, dict_real_EtOH_fulxes, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', gene = 'ADH3'):
        plt.figure(figsize = (10, 10))
        x = sorted([int(x) for x in dict_gene_EtOH.keys()])
        y = [dict_gene_EtOH[str(key)] for key in x]
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        plt.plot(x, y, 'o', x, line)
        plt.axhline(y = dict_real_EtOH_fulxes[gene])
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.legend([gene, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f' % dict_real_EtOH_fulxes[gene]])
        plt.show()

    def multiplePlotO2vsEtOH (self, dict_EtOH_fluxes, dict_real_EtOH_fluxes, n = 3, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', pdf_filename = None):
        plt.rcParams["figure.figsize"] = (20,3)
        g_list = [genes for genes in iter(self.divide_list(sorted(list(dict_EtOH_fluxes.keys())), n))]

        for g_ind in range(0, len(g_list)):
            plt.figure(g_ind, figsize = (10, 10))
            i = 1
            for gene in g_list[g_ind]:
                x = sorted([int(x) for x in dict_EtOH_fluxes[gene].keys()])
                y = [dict_EtOH_fluxes[gene][str(key)] for key in x]
                slope, intercept, r_value, p_value, std_err = linregress(x, y)
                line = [slope * x + intercept for x in x]
                real_O2 = lambda x0: (y0 - intercept) / slope
                y0 = dict_real_EtOH_fluxes[gene]
                plt.subplot(n, 1, i)
                plt.plot(x, y, 'o', x, line)
                plt.axhline(y = dict_real_EtOH_fluxes[gene], ls = 'dashed')
                plt.ylabel(ylab)
                # plt.title(title)
                plt.legend([gene, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f (O2 flux of %.2f)' % (dict_real_EtOH_fluxes[gene], real_O2(y0))])

                i += 1
            plt.xlabel(xlab)

        if pdf_filename is not None:
            pdf = PdfPages(pdf_filename)
            for i in plt.get_fignums():
                pdf.savefig(plt.figure(i))
            pdf.close()

        return [plt.figure(i) for i in plt.get_fignums()]

    def multiplePlotO2vsEtOHSaveFigs (self, dict_EtOH_fluxes, dict_real_EtOH_fluxes, xlab = 'O2 Flux', ylab = 'EtOH Flux', title = 'Ethanol production with O2 flux', folder = None, type = 'fba'):
        plt.rcParams["figure.figsize"] = (10, 4)

        for i, gene in enumerate(sorted(list(dict_EtOH_fluxes.keys()))):
            plt.figure(i)
            x = sorted([int(x) for x in dict_EtOH_fluxes[gene].keys()])
            y = [dict_EtOH_fluxes[gene][str(key)] for key in x]
            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            line = [slope * x + intercept for x in x]
            real_O2 = lambda x0: (y0 - intercept) / slope
            y0 = dict_real_EtOH_fluxes[gene]
            plt.plot(x, y, 'o', x, line)
            plt.axhline(y = dict_real_EtOH_fluxes[gene], ls = 'dashed')
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)
            plt.legend([gene, 'R2: %.4f' % r_value**2, 'Real EtOH flux: %.2f (O2 flux of %.2f)' % (dict_real_EtOH_fluxes[gene], real_O2(y0))])
            if type == 'pfba' or type == 'lmoma':
                plt.savefig(folder + '/' + gene + '_etOH_pFBA.png')
            elif type == 'fba':
                plt.savefig(folder + '/' + gene + '_etOH_FBA.png')


    def plotGeneExpVsSim (self, absRelErrorDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'ADH3', pdf_filename = None, gene = 'ADH3'):
        plt.rcParams["figure.figsize"] = (10,5)

        x = absRelErrorDataset[gene + '_real_flux']
        y = absRelErrorDataset[gene + '_sim_flux']
        gene_IDs = list(absRelErrorDataset.index)
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        meanRelErr = absRelErrorDataset[gene + '_rel_error'].mean()
        corr = absRelErrorDataset[gene + '_real_flux'].corr(absRelErrorDataset[gene + '_sim_flux'])

        plt.plot(x, y, 'o', x, line)
        for ind, gene_ID in enumerate(gene_IDs):
            plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.plot([], [], ' ') # To show correlation in legend
        plt.plot([], [], ' ') # To show mean relative error in legend
        plt.legend([gene, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

    def multipleGenesPlotExpVsSim (self, absRelErrorDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', pdf_filename = None):
        plt.rcParams["figure.figsize"] = (20,3)

        genes = sorted(set([name.split('_')[0] for name in list(absRelErrorDataset.columns)]))
        g_list = [genes for genes in iter(self.divide_list(genes, n))]

        for g_ind in range(0, len(g_list)):
            plt.figure(g_ind, figsize = (10, 10))
            i = 1
            for gene in g_list[g_ind]:
                x = absRelErrorDataset[gene + '_real_flux']
                y = absRelErrorDataset[gene + '_sim_flux']
                gene_IDs = list(absRelErrorDataset.index)
                slope, intercept, r_value, p_value, std_err = linregress(x, y)
                line = [slope * x + intercept for x in x]
                meanRelErr = absRelErrorDataset[gene + '_rel_error'].mean()
                corr = absRelErrorDataset[gene + '_real_flux'].corr(absRelErrorDataset[gene + '_sim_flux'])

                plt.subplot(n, 1, i)
                plt.plot(x, y, 'o', x, line)
                for ind, gene_ID in enumerate(gene_IDs):
                    plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))

                plt.ylabel(ylab)
                plt.plot([], [], ' ') # To show correlation in legend
                plt.plot([], [], ' ') # To show mean relative error in legend
                plt.legend([gene, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

                i += 1
            plt.xlabel(xlab)

        if pdf_filename is not None:
            pdf = PdfPages(pdf_filename)
            for i in plt.get_fignums():
                pdf.savefig(plt.figure(i))
            pdf.close()

        return [plt.figure(i) for i in plt.get_fignums()]

    def multipleGenesPlotExpVsSimSaveFigs (self, absRelErrorDataset, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'Experimental vs Simulated Fluxes', folder = None):
        plt.rcParams["figure.figsize"] = (10, 4)
        genes = sorted(set([name.split('_')[0] for name in list(absRelErrorDataset.columns)]))

        for i, gene in enumerate(genes):
            plt.figure(i)
            x = absRelErrorDataset[gene + '_real_flux']
            y = absRelErrorDataset[gene + '_sim_flux']
            gene_IDs = list(absRelErrorDataset.index)
            slope, intercept, r_value, p_value, std_err = linregress(x, y)
            line = [slope * x + intercept for x in x]
            meanRelErr = absRelErrorDataset[gene + '_rel_error'].mean()
            corr = absRelErrorDataset[gene + '_real_flux'].corr(absRelErrorDataset[gene + '_sim_flux'])

            plt.plot(x, y, 'o', x, line)
            for ind, gene_ID in enumerate(gene_IDs):
                plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))
            plt.xlabel(xlab)
            plt.ylabel(ylab)
            plt.title(title)
            plt.plot([], [], ' ') # To show correlation in legend
            plt.plot([], [], ' ') # To show mean relative error in legend
            plt.legend([gene, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])
            plt.savefig(folder + '/' + gene + '_genes.png')

    def plotReactExpVsSim (self, simVsExpDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'r_0534', reaction = 'r_0534'):
        plt.rcParams["figure.figsize"] = (10,5)
        # Prepare dataset
        df = simVsExpDataset.copy()
        df_sim = df[list(df.columns[1::2])].transpose()
        df_exp = df[list(df.columns[0::2])].transpose()

        x = df_exp[reaction]
        y = df_sim[reaction]
        genes = [gene.split('_')[0] for gene in list(df_exp.index)]
        x.index = y.index = genes
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        line = [slope * x + intercept for x in x]
        meanRelErr = self.relativeError(x, y).mean()
        corr = x.corr(y)

        plt.plot(x, y, 'o', x, line)
        for ind, gene_ID in enumerate(genes):
            plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind], y[ind]))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.plot([], [], ' ') # To show correlation in legend
        plt.plot([], [], ' ') # To show mean relative error in legend
        plt.legend([reaction, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

    def multipleReactsPlotExpVsSim (self, simVsExpDataset, n = 3, xlab = 'Experimental Flux', ylab = 'Simulated Flux', pdf_filename = None, folder = None):
        plt.rcParams["figure.figsize"] = (10, 10)
        # Prepare dataset
        df = simVsExpDataset.copy()
        df_sim = df[list(df.columns[1::2])].transpose()
        df_exp = df[list(df.columns[0::2])].transpose()
        r_list = [r for r in iter(self.divide_list(list(df_exp.columns), n))]
        genes = [gene.split('_')[0] for gene in list(df_exp.index)]

        for r_ind in range(0, len(r_list)):
            plt.figure(r_ind)
            i = 1
            for reaction in r_list[r_ind]:
                x = df_exp[reaction]
                y = df_sim[reaction]
                x.index = y.index = genes
                slope, intercept, r_value, p_value, std_err = linregress(x, y)
                line = [slope * x + intercept for x in x]
                meanRelErr = self.relativeError(x, y).mean()
                corr = x.corr(y)
                plt.subplot(n, 1, i)
                plt.plot(x, y, 'o', x, line)
                for ind, gene_ID in enumerate(genes):
                    plt.annotate(gene_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind], y[ind]))

                plt.ylabel(ylab)
                plt.plot([], [], ' ') # To show correlation in legend
                plt.plot([], [], ' ') # To show mean relative error in legend
                plt.legend([reaction, 'R2: %.4f' % r_value**2, 'Pearson correlation: %.4f' % corr, 'Mean relative error: %.4f' % meanRelErr])

                i += 1

            plt.xlabel(xlab)
            if folder is not None:
                plt.savefig(folder + '/' + str(r_ind) + '_reacts.png')

        if pdf_filename is not None:
            pdf = PdfPages(pdf_filename)
            for i in plt.get_fignums():
                pdf.savefig(plt.figure(i))
            pdf.close()

        return [plt.figure(i) for i in plt.get_fignums()]

    def plotGeneExpVsSimFVA (self, simVsExpDataset, xlab = 'Experimental Flux', ylab = 'Simulated Flux', title = 'ADH3', gene = 'ADH3'):
        # NOT FINISHED
        plt.rcParams["figure.figsize"] = (10,5)

        x = simVsExpDataset[gene + '_real_flux']
        ymin = simVsExpDataset[gene + '_minimum']
        ymax = simVsExpDataset[gene + '_maximum']
        y = [np.mean(y) for y in list(zip(ymin, ymax))]

        plt.errorbar(x, y, yerr = [ymin, ymax], fmt = 'o')
        #plt.plot(x, 'o')
        #lt.fill_between(range(len(ymin)), ymin, ymax, alpha=.1)

        # react_IDs = list(simVsExpDataset.index)
        #
        # for ind, react_ID in enumerate(react_IDs):
        #     plt.annotate(react_ID, (x[ind], y[ind]), fontsize = 8, xytext = (x[ind] + 0.25, y[ind] + 0.5))

        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.show()



# Pipelines

def case7Pipeline (plot = False, makeFigs = False, type = 'fba'):
    # Experimental Fluxes
    exp_dataset, reactions = case7.loadExperimentalRes('Results/Case 7/case7_experimental_fluxes.csv')

    #Add Biomass Values (r_4041)
    reactions['r_4041'] = 'Biomass'
    genes = [gene.split('_')[0] for gene in list(exp_dataset.columns)]
    exp_dataset.loc['r_4041'] = [case7.exp_biomass[gene] for gene in genes]

    # Get real EtOH fluxes
    real_EtOH_fluxes = case7.getEthanolFux(exp_dataset, 'r_2115')
    real_EtOH_fluxes['WildType'] = 23.6318184606019 #From authors file

    if type == 'pfba' or type == 'lmoma':
        # Testing EtOH fluxes with O2 consumption
        # sim_EtOH_O2_fluxes = case7.testAllO2EthanolProd(do_fba = ispFBA)
        # case7.saveObjectToFile(sim_EtOH_O2_fluxes, 'Results/Case 7/dict_etOH_O2_fluxes_pfba.sav')
        sim_EtOH_O2_fluxes = case7.loadObjectFromFile('Results/Case 7/dict_etOH_O2_fluxes_pfba.sav')
    elif type == 'fba':
        # Testing EtOH fluxes with O2 consumption
        # sim_EtOH_O2_fluxes = case7.testAllO2EthanolProd()
        # case7.saveObjectToFile(sim_EtOH_O2_fluxes, 'Results/Case 7/dict_etOH_O2_fluxes.sav')
        sim_EtOH_O2_fluxes = case7.loadObjectFromFile('Results/Case 7/dict_etOH_O2_fluxes.sav')

    # Fix EtOH fluxes with O2 consumption for plotting
    sim_EtOH_O2_fluxes_fixed = case7.fixEtO2FluxesForPlotting(sim_EtOH_O2_fluxes)
    if plot:
        if type == 'pfba' or type == 'lmoma':
            case7.multiplePlotO2vsEtOH(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes, pdf_filename = 'Results/Case 7/etOH_O2_fluxes_plot_pfba.pdf')
        elif type == 'fba':
            case7.multiplePlotO2vsEtOH(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes, pdf_filename = 'Results/Case 7/etOH_O2_fluxes_plot.pdf')
    if makeFigs:
        if type == 'pfba' or type == 'lmoma':
            case7.multiplePlotO2vsEtOHSaveFigs(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes, folder = 'Results/Case 7/EtOH_figs', type = 'pfba')
        elif type == 'fba':
            case7.multiplePlotO2vsEtOHSaveFigs(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes, folder = 'Results/Case 7/EtOH_figs', type = 'fba')

    # Get correct O2 flux according to exp EtOH flux
    fluxes_O2 = case7.getO2flux(sim_EtOH_O2_fluxes_fixed, real_EtOH_fluxes)
    #case7.printDict(fluxes_O2)

    return exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2

def fbaPipeline (plotGenes = False, plotReacts = False, saveGenesPlot = False):
    # res_fba = case7.case7fba(fluxes_O2)
    # case7.saveObjectToFile(res_fba, 'Results/Case 7/res_fba_case7.sav')
    res_fba = case7.loadObjectFromFile('Results/Case 7/res_fba_case7.sav')
    res_fba_df, wt_fba_df = case7.createResultsDataset(res_fba)
    res_fba_df = case7.correctReversedReactions(res_fba_df) # Reactions fixed
    wt_fba_df = case7.correctReversedReactions(wt_fba_df)

    # Dataset with experimental vs simulated fluxes
    df_fba_exp_sim = case7.createDatasetExpVsSimul(exp_dataset, res_fba_df)

    # Dataset with absolute and realtive errors
    df_fba_exp_sim_errors = case7.createDatasetWithAbsRelError(df_fba_exp_sim)

    if plotGenes:
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_fba_exp_sim_errors, pdf_filename = 'Results/Case 7/fba_genes_exp_vs_sim_plots.pdf')
    if plotReacts:
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_fba_exp_sim, pdf_filename = 'Results/Case 7/fba_reacts_exp_vs_sim_plots.pdf', folder = 'Results/Case 7/FBA_figs')
    if saveGenesPlot:
        case7.multipleGenesPlotExpVsSimSaveFigs(df_fba_exp_sim_errors, folder = 'Results/Case 7/FBA_figs')

    return res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors

def pfbaPipeline (plotGenes = False, plotReacts = False, saveGenesPlot = False):
    # res_pfba = case7.case7pfba(fluxes_O2)
    # case7.saveObjectToFile(res_pfba, 'Results/Case 7/res_pfba_case7.sav')
    res_pfba = case7.loadObjectFromFile('Results/Case 7/res_pfba_case7.sav')
    res_pfba_df, wt_pfba_df = case7.createResultsDataset(res_pfba)
    res_pfba_df = case7.correctReversedReactions(res_pfba_df) # Reactions fixed
    wt_pfba_df = case7.correctReversedReactions(wt_pfba_df)

    # Dataset with experimental vs simulated fluxes
    df_pfba_exp_sim = case7.createDatasetExpVsSimul(exp_dataset, res_pfba_df)

    # Dataset with absolute and realtive errors
    df_pfba_exp_sim_errors = case7.createDatasetWithAbsRelError(df_pfba_exp_sim)

    if plotGenes:
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_pfba_exp_sim_errors, pdf_filename = 'Results/Case 7/pfba_genes_exp_vs_sim_plots.pdf')
    if plotReacts:
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_pfba_exp_sim, pdf_filename = 'Results/Case 7/pfba_reacts_exp_vs_sim_plots.pdf', folder = 'Results/Case 7/pFBA_figs')
    if saveGenesPlot:
        case7.multipleGenesPlotExpVsSimSaveFigs(df_pfba_exp_sim_errors, folder = 'Results/Case 7/pFBA_figs')

    return res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors

def fvaPipeline ():
    # res_fva = case7.case7fva(fluxes_O2)
    # case7.saveObjectToFile(res_fva, 'Results/Case 7/res_fva_case7.sav')
    res_fva = case7.loadObjectFromFile('Results/Case 7/res_fva_case7.sav')
    res_fva_df, wt_fva_df = case7.createResultsDatasetFVA(res_fva)
    res_fva_df = case7.correctReversedReactions(res_fva_df) # Reactions fixed
    wt_fva_df = case7.correctReversedReactions(wt_fva_df)

    # Dataset with experimental vs simulated fluxes
    df_fva_exp_sim = case7.createDatasetExpVsSimulFVA(exp_dataset, res_fva_df)

    #NOT FINISHED - PLOTS NEXT (?)

    # if plot == 'genes':
    #     #Plots w/ errors for genes
    #     # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
    #     case7.multipleGenesPlotExpVsSim(df_pfba_exp_sim_errors, pdf_filename = 'Results/Case 7/pfba_genes_exp_vs_sim_plots.pdf')
    # elif plot == 'reactions':
    #     #Plots w/ errors for reactions
    #     # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
    #     case7.multipleReactsPlotExpVsSim(df_fba_exp_sim, pdf_filename = 'Results/Case 7/pfba_reacts_exp_vs_sim_plots.pdf')

    return res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim

def lmomaPipeline (reference_dict = None, plotGenes = False, plotReacts = False, saveGenesPlot = False):
    #res_lmoma = case7.case7lmoma(fluxes_O2, reference_dict)
    #case7.saveObjectToFile(res_pfba, 'Results/Case 7/res_lmoma_case7.sav')
    res_lmoma = case7.loadObjectFromFile('Results/Case 7/res_lmoma_case7.sav')
    res_lmoma_df, wt_lmoma_df = case7.createResultsDataset(res_lmoma)
    res_lmoma_df = case7.correctReversedReactions(res_lmoma_df) # Reactions fixed
    wt_lmoma_df = case7.correctReversedReactions(wt_lmoma_df)

    # Dataset with experimental vs simulated fluxes
    df_lmoma_exp_sim = case7.createDatasetExpVsSimul(exp_dataset, res_lmoma_df)

    # Dataset with absolute and realtive errors
    df_lmoma_exp_sim_errors = case7.createDatasetWithAbsRelError(df_lmoma_exp_sim)

    if plotGenes:
        #Plots w/ errors for genes
        # case7.plotGeneExpVsSim(df_fba_exp_sim_errors, gene = 'ADH3')
        case7.multipleGenesPlotExpVsSim(df_lmoma_exp_sim_errors, pdf_filename = 'Results/Case 7/lmoma_genes_exp_vs_sim_plots.pdf')
    if plotReacts:
        #Plots w/ errors for reactions
        # case7.plotReactExpVsSim(df_fba_exp_sim, reaction = 'r_1054')
        case7.multipleReactsPlotExpVsSim(df_lmoma_exp_sim, pdf_filename = 'Results/Case 7/lmoma_reacts_exp_vs_sim_plots.pdf', folder = 'Results/Case 7/LMOMA_figs')
    if saveGenesPlot:
        case7.multipleGenesPlotExpVsSimSaveFigs(df_lmoma_exp_sim_errors, folder = 'Results/Case 7/LMOMA_figs')

    return res_lmoma, res_lmoma_df, wt_lmoma_df, df_lmoma_exp_sim, df_lmoma_exp_sim_errors

if __name__ == '__main__':

    #Initialization
    case7 = Case7()
    case7.model = case7.loadObjectFromFile('model_yeast_76.sav')
    case7.setMedium('MINIMAL_CASE7')
    case7.dictsForCase7()

    # General datasets
    exp_dataset, reactions, real_EtOH_fluxes, sim_EtOH_O2_fluxes_fixed, fluxes_O2 = case7Pipeline(plot = False, makeFigs = False, type = 'lmoma')

    # FBA
    res_fba, res_fba_df, wt_fba_df, df_fba_exp_sim, df_fba_exp_sim_errors = fbaPipeline(saveGenesPlot = False, plotReacts = False, plotGenes = False)

    # pFBA
    res_pfba, res_pfba_df, wt_pfba_df, df_pfba_exp_sim, df_pfba_exp_sim_errors = pfbaPipeline(saveGenesPlot = False, plotReacts = False, plotGenes = False)

    # FVA
    res_fva, res_fva_df, wt_fva_df, df_fva_exp_sim = fvaPipeline()

    # LMOMA
    res_lmoma, res_lmoma_df, wt_lmoma_df, df_lmoma_exp_sim, df_lmoma_exp_sim_errors = lmomaPipeline(reference_dict = res_pfba, plotGenes = False, plotReacts = False, saveGenesPlot = False)

    # Create xlsx with results
    # case7.convertPandasDFToExcel(reactions, df_fba_exp_sim_errors, filename = 'Results/Case 7/fba_results_case7_test.xlsx', imgFolder = 'Results/Case 7/FBA_figs')
    # case7.convertPandasDFToExcel(reactions, df_pfba_exp_sim_errors, title = 'pFBA Results Case 7', filename = 'Results/Case 7/pfba_results_case7_test.xlsx', imgFolder = 'Results/Case 7/pFBA_figs', type = 'pfba')
    # case7.convertPandasDFToExcelFVA(reactions, df_fva_exp_sim, title = 'FVA Results Case 7', filename = 'Results/Case 7/fva_results_case7.xlsx')
    case7.convertPandasDFToExcel(reactions, df_lmoma_exp_sim_errors, title = 'LMOMA Results Case 7', filename = 'Results/Case 7/lmoma_results_case7.xlsx', imgFolder = 'Results/Case 7/LMOMA_figs', type = 'lmoma')


    # Compare biomasses
    # biom = case7.getBiomassObj(res_fba)
    # case7.printDict(biom)
    # case7.printDict(case7.exp_biomass)

    # Check reactions
    # checkReaction('r_4041') #Biomass
    # checkReaction('r_2115') #EtOH

    # xlsxwriter
    # http://xlsxwriter.readthedocs.io/working_with_pandas.html
    # http://xlsxwriter.readthedocs.io/worksheet.html

