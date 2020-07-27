import os, copy, re
import math
import random
import pandas as pd
import json
from pyDOE import lhs
import numpy as np
from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.utils.builders.TemplateHelper import TemplateHelper
from dtk.utils.builders.ConfigTemplate import ConfigTemplate
from dtk.utils.builders.TaggedTemplate import CampaignTemplate, DemographicsTemplate
from dtk.utils.Campaign.utils.CampaignEncoder import CampaignEncoder
from history_matching import quick_read
from malaria.reports.MalariaReport import add_patient_report, add_survey_report, add_summary_report
from dtk.interventions.input_EIR import add_InputEIR
import malaria.site.input_EIR_by_site

##################### Calibration Setup #########################

iteration = int(re.search(r'iter(\d+)', os.getcwd()).group(1))
N_rep_per_sample = 1
N_samples = 200
params = quick_read( os.path.join('..', 'Params.xlsx'), 'Params').set_index('Name')
exp_name = 'sensitiveparams_singlesite_sugungum_Iter%d'%iteration


##################### Simulation Setup #########################

site = 'Sugungum'
mEIR_profile = malaria.site.input_EIR_by_site.study_site_monthly_EIRs[site]
years_burnin = 20
years_start_to_report = 1
years_to_report = 2

simulation_duration = (years_burnin + years_start_to_report+years_to_report+2)*365 #in days

input_files_dir = os.path.join('..', 'InputFiles')
default_config = DTKConfigBuilder.from_files(os.path.join(input_files_dir, 'config.json'),
                                             os.path.join(input_files_dir, 'campaign.json'))
#I want to make experiment specific changes to that config json
#MALARIA_FUNCTIONAL_MODEL is the model type flag to trigger Malaria2.0 style infection and immune model
default_config.update_params({"Malaria_Model": "MALARIA_FUNCTIONAL_MODEL"})

config_params_to_update = {
#Generic sim config updates
"Base_Population_Scale_Factor": 1,
"Enable_Vital_Dynamics": 1,
"Enable_Disease_Mortality":0,
"Simulation_Duration": simulation_duration,
#Pyrogenic Threshold params
"Base_PT": 15000.0,
"ActiveInfectionsScaleFactorPT": 1.0,
"AgeScaleFactorPT": 1.0,
"Age_PT_Increment": 1.0,
"ClearedInfectionsScaleFactorPT": 1.0,
"ExposureScaleFactorPT": 0.0001,
"FreeIntervalScaleFactorPT": 1.0,
"ExposureThreshold": 100.0,
"Individual_Variation_Magnitude_PT": 1.0,
"NumWavesScaleFactorPT": 1.0,
#Immune Modifier Inclusion Params
"Immune_Modifier_Include_Age": 1,
"Immune_Modifier_Include_Cumulative_Exposure": 1,
"Immune_Modifier_Include_Recent_Exposure": 0,
"Immune_Modifier_Include_Strain_Diversity": 0,
"Immune_Modifier_Include_Malaria_Free_Interval": 0,

#PPP Params
"First_Peak_Mean": 4.7474,
"First_Peak_Std": 0.42609,
"Biological_Age_Immune_Coefficient_PPP": 3,
"Cumulative_Exposure_Immune_Coefficient_PPP": 1.0,
"Malaria_Free_Interval_Immune_Coefficient_PPP": 0.0,
"Recent_Exposure_Immune_Coefficient_PPP": 0.0,
"Strain_Diversity_Immune_Coefficient_PPP": 0.0,
#TM Params
"Biological_Age_Immune_Coefficient_TM": 1,
"Cumulative_Exposure_Immune_Coefficient_TM": 1.0,
"Malaria_Free_Interval_Immune_Coefficient_TM": 0.0,
"Recent_Exposure_Immune_Coefficient_TM": 0.0,
"Strain_Diversity_Immune_Coefficient_TM": 0.0,
"Parasite_Density_Bin_Edges": [
         0,
         100,
         1000,
         10000,
         100000,
         1000000
      ],
"Parasite_Peak_Density_Probabilities": [
         [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
         ],
         [
            0.373134328,
            0.358208955,
            0.179104478,
            0.074626866,
            0.014925373,
            0.0
         ],
         [
            0.277108434,
            0.289156627,
            0.301204819,
            0.13253012,
            0.0,
            0.0
         ],
         [
            0.042253521,
            0.126760563,
            0.232394366,
            0.542253521,
            0.056338028,
            0.0
         ],
         [
            0.168421053,
            0.010526316,
            0.126315789,
            0.410526316,
            0.284210526,
            0.0
         ],
         [
            0.176470588,
            0.058823529,
            0.0,
            0.235294118,
            0.529411765,
            0.0
         ]
      ],
#Immunity Params
"Concurrent_Infections_a": 4.0,
"Concurrent_Infections_b": 0.5,
"Concurrent_Infections_c": 7.5,
"Diversity_From_Marginal_Prevalence": 40.0,
"Parasite_Density_Wave_Sigma": 5,
"Prevalence_Threshold_For_Diversity": 0.1,
"Scale_Factor_Age_a": 4.0,
"Scale_Factor_Age_b": 0.2,
"Scale_Factor_Cum_Exp_shape": 2.0,
"Scale_Factor_Free_Interval_a": 0.01,
"Scale_Factor_Rec_Exp_midpoint": 4.0,
"Scale_Factor_Rec_Exp_steepness": 2.0,
"Scale_Factor_Strain_Diversity_n": 2.0,
"Wave_Vs_Infection_Relative_Weight": 1.0,
"Use_Fixed_Scale_Factor": 0,
"Use_Fixed_Wave_Period": 0,
"use_PRISM_age_modifiers": 0,
#Sexual Stage Params
"Gametocyte_Conversion_Draw_Mean": 1.8,
"Gametocyte_Conversion_Draw_Std_Dev": 1.0,
#Genetics Params
"Genome_Markers": [],
"Number_Basestrains": 1,
"Number_Substrains": 1,
#Climate and site demographic Params
"Geography": site,
"Default_Geography_Initial_Node_Population": 100,
"Demographics_Filenames": [f"Namawala_single_node_demographics.json"],
"Air_Temperature_Filename": f"Namawala_single_node_air_temperature_daily.bin",
"Land_Temperature_Filename": f"Namawala_single_node_land_temperature_daily.bin",
"Rainfall_Filename": f"Namawala_single_node_rainfall_daily.bin",
"Relative_Humidity_Filename": f"Namawala_single_node_relative_humidity_daily.bin",

# "Demographics_Filenames": [f"{site}_single_node_demographics.json"],
# "Air_Temperature_Filename": f"{site}_single_node_air_temperature_daily.bin",
# "Land_Temperature_Filename": f"{site}_single_node_land_temperature_daily.bin",
# "Rainfall_Filename": f"{site}_single_node_rainfall_daily.bin",
# "Relative_Humidity_Filename": f"{site}_single_node_relative_humidity_daily.bin",
#Vector Params
"Vector_Species_Names": [],
"Vector_Sampling_Type": "VECTOR_COMPARTMENTS_NUMBER",
"x_Temporary_Larval_Habitat": 0.01

}
default_config.update_params(config_params_to_update)

### Reset the mEIR profile to be used for this site
default_config.campaign.Events[0].Event_Coordinator_Config.Intervention_Config.Monthly_EIR = mEIR_profile
default_config.campaign.save_to_file(filename = os.path.join('..', 'InputFiles','template_campaign' ))

def save_to_file(self, filename=None):
    if filename is None:
        filename = 'output_class'
    content = self.to_json()
    f = open('{}.json'.format(filename), 'w')
    f.write(content)
    f.close()

#I want to save that json locally
def write_file(working_directory, name, content):
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
    filename = os.path.join(working_directory, '%s' % name)
    with open(filename, 'w') as f:
        json.dump(content,f)

write_file(os.path.join('..', 'InputFiles'),'template_config.json', default_config.config)


plugin_files_dir = os.path.join('..', 'InputFiles')
#I want to read from that json to make a template for subsequent exp_builder modification.
cfg = ConfigTemplate.from_file( os.path.join(input_files_dir, 'template_config.json') )
cpn = CampaignTemplate.from_file( os.path.join(input_files_dir, 'template_campaign.json') )

config_params = {
        'Base_Population_Scale_Factor': 0.1,
        'Use_Fixed_Scale_Factor': 0,
        'Enable_Vital_Dynamics': 1,
        'Enable_Disease_Mortality':0,
        'Demographics_Filenames' : ['Namawala_single_node_demographics.json'],
        'Max_Individual_Infections': 30,
        'Use_Fixed_Wave_Period': 0
        # 'Report_Detection_Threshold_Blood_Smear_Gametocytes': 40
    }


templates = TemplateHelper()

table_base = {
    'ACTIVE_TEMPLATES': [cfg, cpn],
    'TAGS': {'BayesianHistoryMatching':None, 'Iteration': iteration}
}
# Standard DTKConfigBuilder
config_builder = DTKConfigBuilder.from_files(
    os.path.join(plugin_files_dir, 'template_config.json'),
    os.path.join(plugin_files_dir, 'template_campaign.json')
)


add_summary_report(config_builder,
                   start = years_burnin+years_start_to_report*365,
                   description='Monthly_Report',
                   interval = 365/12,
                   nreports = years_to_report*12,
                   age_bins = [1, 4, 8, 18, 28, 43, 125],
                   parasitemia_bins = [0, 50, 200, 500, 2000000]
                   )

##################### Sampling Functions ##########################

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
        if "Include" in param_name:
            samples[param_name] = round(pmin + samples[param_name]*(pmax-pmin))
        else:
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


##################### Sim Commissioning #########################
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
                 'exp_name': exp_name}

if __name__ == "__main__":

    if not SetupParser.initialized:
        SetupParser.init('HPC')


    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())