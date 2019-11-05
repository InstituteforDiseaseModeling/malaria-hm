#title: GhanaKNDcohortsite
#
#description: A site for developing the history matching worflow with Ghana KND cohort infection durations by age as an initial target
#
#author: Jon Russell
#
#date: 7/19/2019
#
#notes and dependencies:
#
#Institute for Disease Modeling, Bellevue, WA


from abc import ABCMeta

from calibtool.CalibSite import CalibSite
from malaria.study_sites.site_setup_functions import *
import os
import pandas as pd

from malaria.reports.MalariaReport import add_patient_report


class ChallengeInfectionCalibSite(CalibSite):
    __metaclass__ = ABCMeta
    def __init__(self):
        super(ChallengeInfectionCalibSite, self).__init__(name='Navrongo')


    reference_dict = {
        "Age bin": [],
        "MeanDuration": [],
        "CI_lower": [],
        "CI_upper": []
    }

    def get_setup_functions(self):

        return [add_challenge_trial_fn(start_day=0)]

    def get_reference_data(self,reference_type):

        return ref_df

    def get_analyzers(self):
        from analyzers.CompareMozambiqueInfectionsAnalyzer import CompareMozambiqueInfectionsAnalyzer
        return [CompareMozambiqueInfectionsAnalyzer(site=self)]
