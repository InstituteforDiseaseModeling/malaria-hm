import itertools
import os, re
import sqlite3
from collections import OrderedDict

import pandas as pd
import numpy as np
import seaborn as sns
#Plotting
import matplotlib.pyplot as plt
from scipy.special import gammaln

from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from calibtool.analyzers.Helpers import season_channel_age_density_csv_to_pandas
from calibtool.LL_calculators import dirichlet_multinomial_pandas

from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from dtk.utils.parsers.malaria_summary import summary_channel_to_pandas

from calibtool import LL_calculators
from calibtool.analyzers.Helpers import \
    convert_to_counts, age_from_birth_cohort, season_from_time, aggregate_on_index



def grouped_df(df, pfprdict, index, column_keep, column_del):
    """
    Recut dataframe to recategorize data into desired age and parasitemia bins

    Args:
        df: Dataframe to be rebinned
        pfprdict: Dictionary mapping postive counts per slide view (http://garkiproject.nd.edu/demographic-parasitological-surveys.html)
                to density of parasites/gametocytes per uL
        index: Multi index into which 'df' is rebinned
        column_keep: Column (e.g. parasitemia) to keep
        column_del: Column (e.g. gametocytemia) to delete
    """
    dftemp = df.copy()
    del dftemp[column_del]

    dftemp['PfPR Bin'] = df[column_keep]
    dftemp = aggregate_on_index(dftemp, index)

    dfGrouped = dftemp.groupby(['Season', 'Age Bin', 'PfPR Bin'])

    dftemp = dfGrouped[column_keep].count()
    dftemp = dftemp.unstack().fillna(0).stack()
    dftemp = dftemp.rename(column_keep).reset_index()
    dftemp['PfPR Bin'] = [pfprdict[p] for p in dftemp['PfPR Bin']]

    dftemp = dftemp.set_index(['Season', 'Age Bin', 'PfPR Bin'])

    return dftemp

def get_reference_data(self):
    dir_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ),'..', 'reference data'))
    dfFileName = os.path.join(dir_path, 'Garki_df.csv')
    df = pd.read_csv(dfFileName)
    df = df.loc[df['Village']==self.metadata['village']]
    pfprBinsDensity = self.metadata['density_bins']
    uL_per_field = 0.5 / 200.0  # from Garki PDF - page 111 - 0.5 uL per 200 views
    pfprBins = 1 - np.exp(-np.asarray(pfprBinsDensity) * uL_per_field)
    seasons = self.metadata['seasons']
    pfprdict = dict(zip(pfprBins, pfprBinsDensity))

    bins = OrderedDict([
        ('Season', self.metadata['seasons']),
        ('Age Bin', self.metadata['age_bins']),
        ('PfPR Bin', pfprBins)
    ])
    bin_tuples = list(itertools.product(*bins.values()))
    index = pd.MultiIndex.from_tuples(bin_tuples, names=bins.keys())

    df = df.loc[df['Seasons'].isin(seasons)]
    df = df.rename(columns={'Seasons': 'Season', 'Age': 'Age Bin'})

    df2 = grouped_df(df, pfprdict, index, 'Parasitemia', 'Gametocytemia')
    df3 = grouped_df(df, pfprdict, index, 'Gametocytemia', 'Parasitemia')
    dfJoined = df2.join(df3).fillna(0)
    dfJoined = pd.concat([dfJoined['Gametocytemia'], dfJoined['Parasitemia']])
    dfJoined.name = 'Counts'
    dftemp = dfJoined.reset_index()
    dftemp['Channel'] = 'PfPR by Gametocytemia and Age Bin'
    dftemp.loc[len(dftemp) / 2:, 'Channel'] = 'PfPR by Parasitemia and Age Bin'
    dftemp = dftemp.rename(columns={'Seasons': 'Season', 'PfPR Bins': 'PfPR Bin', 'Age Bins': 'Age Bin'})
    dftemp = dftemp.set_index(['Channel', 'Season', 'Age Bin', 'PfPR Bin'])


    return dftemp

class Summary_Prevalence_Analyzer(BaseAnalyzer):

    def __init__(self, **kwargs):
        super().__init__(filenames=['output/MalariaSummaryReport_Monthly_Report.json'])

        self.metadata =  {
        'density_bins': [0, 50, 200, 500, np.inf],  # (, 0] (0, 50] ... (50000, ]
        'density_bin_edges':['0', '50', '200', '500'],
        'age_bins': [1, 4, 8, 18, 28, 43, np.inf],  # (, 5] (5, 15] (15, ],
        'age_bin_labels':['<1', '1-4', '4-8', '8-18', '18-28', '28-43', '>43'],
        'seasons': ['DC2', 'DH2', 'W2'],
        'seasons_by_month': {
            'Apr': 'DH2',
            'June/Aug': 'W2',
            'Dec/Jan': 'DC2'
        },
        'village': 'Sugungum'
    }
        self.reference = get_reference_data(self)
        # Get channels to extract from 'Channel' level of reference MultiIndex
        ref_ix = self.reference.index
        channels_ix = ref_ix.names.index('Channel')
        self.channels = ref_ix.levels[channels_ix].values

        self.seasons = kwargs.get('seasons')

    def select_simulation_data(self, data, simulation):
        # The multi-index that was used in the reference, should match the bins that were used in the call to the summary report in sim commissioning
        #build the multindex
        #slice the df object by just those months corresponding to seasons of interest
        #Add that dimension to the multindex and d the groupby and count steps
        #return as metadata and sim_df
        ref_df_multiindex_levels = self.reference.index.levels
        summary_data_by_month = data[self.filenames[0]]['DataByTimeAndPfPRBinsAndAgeBins']['PfPR by Parasitemia and Age Bin']
        #the months of sim time corresponding to ref seasons are (0-indexed) 3-Apr, 5-Jun, and 12, for subsequent Jan
        sim_df_season_subset = [summary_data_by_month[x] for x in [3,6,11]]
        sim_df = pd.DataFrame(sim_df_season_subset).T
        #columns should now be season subset
        sim_df.columns = ('DH2', 'W2', 'DC2')
        #index should be age_bins level of multindex
        sim_df.index = ref_df_multiindex_levels[3]

        sim_df.sample = simulation.tags.get("__sample_index__")
        sim_df.sim_id = simulation.tags.get("sim_id")

        return sim_df

    def finalize(self, all_data):
        # compare the sim_data and ref objects by doing a diff and norm
        # append result to a results csv with sampleid, simid, and value
        #plot all results on a facet grid


        channels = self.reference.index.levels[0].values
        seasons = self.reference.index.levels[1].values
        age_bins = self.reference.index.levels[2].values
        density_bins = self.reference.index.levels[3].values



        results = pd.DataFrame()
        for sim, data in all_data.items():
            exp_id = sim.experiment_id


            for season in range(len(seasons)):
                for density_bin in density_bins:
                    for age_bin in range(len(age_bins)):
                        sim_value = data.loc[density_bin, seasons[season]][age_bin]




                        sim_id = sim.id
                        sim.tags.get("__sample_index__")
                        sample = sim.tags.get("__sample_index__")
                        sub_results = pd.DataFrame({'sample': sample,
                                                    'sim_id': sim_id,
                                                    'season': seasons[season],
                                                    'age_bin': age_bins[age_bin],
                                                    'density_bin': density_bin,
                                                    'value': [sim_value]})
                        results = pd.concat([results, sub_results])

        #define exp_id

        results.to_csv(os.path.join('..', 'iter0', f'{exp_id}','analyzer_results.csv'))


        ######


if __name__ == '__main__':

    am = AnalyzeManager('d8ccb6e9-1e04-ea11-a2c3-c4346bcb1551',analyzers=Summary_Prevalence_Analyzer())
    am.analyze()
