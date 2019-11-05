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
from malaria.reports.MalariaReport import add_patient_report


iteration = int(re.search(r'iter(\d+)', os.getcwd()).group(1))


params = quick_read( os.path.join('..', 'Params.xlsx'), 'Params').set_index('Name')

plugin_files_dir = os.path.join('..', 'InputFiles')

# Standard DTKConfigBuilder
config_builder = DTKConfigBuilder.from_files(
    os.path.join(plugin_files_dir, 'burnin_config.json'),
    os.path.join(plugin_files_dir, 'burnin_campaign.json')
)

config_builder.update_params({
    "Simulation_Duration": 365*41,
    "Serialization_Time_Steps": [
        365*40
    ],
    "Serialization_Type": "TIMESTEP"
})
# add_patient_report(config_builder)


run_sim_args =  {'config_builder': config_builder,
                 'exp_name': 'Malaria Challenge burnin for v2'}

if __name__ == "__main__":

    if not SetupParser.initialized:
        SetupParser.init('HPC')


    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())