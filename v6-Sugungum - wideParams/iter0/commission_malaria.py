import os, copy, re
import math
import random
import pandas as pd
from pyDOE import lhs
import numpy as np

from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.utils.builders.TemplateHelper import TemplateHelper
from dtk.utils.builders.ConfigTemplate import ConfigTemplate
from dtk.utils.builders.TaggedTemplate import CampaignTemplate, DemographicsTemplate
from history_matching import quick_read
from malaria.reports.MalariaReport import add_patient_report, add_survey_report, add_summary_report
from dtk.interventions.input_EIR import add_InputEIR

iteration = int(re.search(r'iter(\d+)', os.getcwd()).group(1))
N_rep_per_sample = 1
N_samples = 512

params = quick_read( os.path.join('..', 'Params.xlsx'), 'Params').set_index('Name')

plugin_files_dir = os.path.join('..', 'InputFiles')
cfg = ConfigTemplate.from_file( os.path.join(plugin_files_dir, 'config.json') )
cpn = CampaignTemplate.from_file( os.path.join(plugin_files_dir, 'campaign.json') )
# demog = DemographicsTemplate.from_file( os.path.join(plugin_files_dir, 'Malariatherapy_demographics.json') )

templates = TemplateHelper()
table_base = {
    'ACTIVE_TEMPLATES': [cfg, cpn],
    'TAGS': {'BayesianHistoryMatching':None, 'Iteration': iteration}
}

# Standard DTKConfigBuilder
config_builder = DTKConfigBuilder.from_files(
    os.path.join(plugin_files_dir, 'config.json'),
    os.path.join(plugin_files_dir, 'campaign.json')
)
years = 3

# config_builder.update_params({
#         'Base_Population_Scale_Factor': 0.1,
#         'Use_Fixed_Scale_Factor': 0,
#         'Enable_Vital_Dynamics': 1,
#         'Enable_Disease_Mortality':0,
#         'Demographics_Filenames' : ['Namawala_single_node_demographics.json'],
#         'Max_Individual_Infections': 30,
#         'Use_Fixed_Wave_Period': 0
#         # 'Report_Detection_Threshold_Blood_Smear_Gametocytes': 40
#     })

add_summary_report(config_builder,
                   start = 7300,
                   description='Monthly_Report',
                   interval = 365/12,
                   nreports = years*12,
                   age_bins = [1, 4, 8, 18, 28, 43, 125],
                   parasitemia_bins = [0, 50, 200, 500, 2000000]
                   )
# add_patient_report(config_builder)

# "Sugungum":     [2, 1, 0.1, 0.1, 0.1, 5, 9.625, 10.375, 20, 23.25, 4.125, 4],
#     "Matsari":      [4, 1, 0.1, 0.1, 1, 2.5, 7.75, 11.625, 56.25, 27.125, 12, 6],
#     'Siaya':        [55*m for m in [0.07, 0.08, 0.04, 0.02, 0.03, 0.04, 0.22, 0.13, 0.18, 0.09, 0.07, 0.05]],   # EIR =  56
#     'Laye':         [1, 1, 1, 1, 1, 7, 8, 9, 5, 4, 6, 1],                                                       # EIR =  45
#     'Kilifi_South': [3.6, 1.4, 0.0, 2.8, 0.4, 10.3, 10.6, 3.1, 0.4, 0.0, 0.7, 0.8],                             # EIR =  34
#     "Ndiop":        [0.39, 0.19, 0.77, 0, 0, 0, 6.4, 2.2, 4.7, 3.9, 0.87, 0.58],                                # EIR =  20
#     "Rafin_Marke":  [1, 1, 0.5, 1, 1, 2, 3.875, 7.75, 15.0, 3.875, 1, 1],
scale_factor = 1
monthly_EIR = [1, 1, 0.5, 1, 1, 2, 3.875, 7.75, 15.0, 3.875, 1, 1]
add_InputEIR(config_builder, [x * scale_factor for x in monthly_EIR])

def map_sample_to_model_input(config_builder, sample_idx, replicate_idx, sample):
    table = copy.deepcopy(table_base)
    table['TAGS'].update({'__sample_index__': sample_idx, '__replicate_index__': replicate_idx})
    table['Run_Number'] = random.randint(0, 1e6)
    table['TAGS'].update({'[SAMPLE] %s' % k: v for k, v in sample.items()})

    for param_name,p in params.iterrows():
        if param_name in sample and 'MapTo' in p:
            if isinstance(p['MapTo'],float) and math.isnan(p['MapTo']):
                continue
            table[p['MapTo']] = sample.pop( param_name )

    for name,value in sample.items():
        print('UNUSED PARAMETER:', name)
    assert( len(sample) == 0 ) # All params used

    return templates.mod_dynamic_parameters(config_builder, table)


def choose_and_scale_samples_unconstrained(num_samples):
    N_dim = params.shape[0]
    samples = pd.DataFrame(lhs(N_dim, samples = num_samples), columns=params.index.tolist())

    for param_name in samples.columns.values:
        pmin,pmax = (params.loc[param_name,'Min'], params.loc[param_name,'Max'])
        samples[param_name] = pmin + samples[param_name]*(pmax-pmin)
    samples.index.name = 'Sample'
    return samples


def choose_and_scale_samples(num_samples):
    if iteration == 0:
        samples = choose_and_scale_samples_unconstrained(num_samples)
        # samples_unconstrained['Days'] = samples_unconstrained[['Env Ramp Up', 'Env Ramp Down', 'Env Cutoff']].sum(axis=1)
        # samples = samples_unconstrained.loc[ samples_unconstrained['Days'] < 365, :]

        remaining = num_samples - samples.shape[0]
        while remaining > 0:
            samples_unconstrained = choose_and_scale_samples_unconstrained(num_samples)
            samples_unconstrained['Days'] = samples_unconstrained[['Env Ramp Up', 'Env Ramp Down', 'Env Cutoff']].sum(axis=1)
            samples = pd.concat([samples, samples_unconstrained.loc[ samples_unconstrained['Days'] < 365, :] ], ignore_index=True)
            remaining = num_samples - samples.shape[0]

        samples.index.name = 'Sample'
        return samples.iloc[:num_samples,:]#.drop('Days', axis=1)
    else:
        #samples = pd.read_excel(os.path.join('..', 'iter%d'%(iteration-1), 'Candidates_for_iter%d.xlsx'%iteration), sheetname='Values')
        samples = pd.read_hdf(os.path.join('..', 'iter%d'%(iteration-1), 'Candidates_for_iter%d.hd5'%iteration), key='values')

        samples.index.name = 'Sample'
        # assert( pd.Series.all( samples[['Env Ramp Up', 'Env Ramp Down', 'Env Cutoff']].sum(axis=1) < 365 ))
        return samples.iloc[:num_samples,:]


try:
    xlsx = pd.ExcelFile('Samples.xlsx')
    samples = pd.read_excel(xlsx, 'Samples')
    samples.set_index('Sample', inplace=True)
    # input('NOTE: Using previously compted samples from Samples.xlsx.  Enter to continue ...')
except IOError:
    samples = choose_and_scale_samples(N_samples)

    writer = pd.ExcelWriter('Samples.xlsx')
    samples.to_excel(writer, sheet_name='Samples')
    params.to_excel(writer, sheet_name='Params')
    writer.save()

exp_builder = ModBuilder.from_combos(
    [
        ModFn(map_sample_to_model_input,
            sample[0],  # <-- sample index
            rep,        # <-- replicate index
            {k:v for k,v in zip(samples.columns.values, sample[1:])})
        for sample in samples.itertuples() for rep in range(N_rep_per_sample)
    ])

run_sim_args =  {'config_builder': config_builder,
                 'exp_builder': exp_builder,
                 'exp_name': 'Rafin Marke_20yrBurn_Iter%d'%iteration}

if __name__ == "__main__":

    if not SetupParser.initialized:
        SetupParser.init('HPC')


    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())