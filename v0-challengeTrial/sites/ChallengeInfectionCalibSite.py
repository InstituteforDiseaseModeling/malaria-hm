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



class NavronogoCalibSite(CalibSite):
    __metaclass__ = ABCMeta
    def __init__(self):
        super(NavronogoCalibSite, self).__init__(name='Navrongo')


    reference_dict = {

        "Age bin": [],
        "MeanDuration": [],
        "MeanDuration": [],
    }


    def get_setup_functions(self):


    def get_reference_data(self,reference_type):

        def collect_reference_trajectories(df):

            infection_trajectories = {}
            infection_trajectories['microscopy'] = pd.DataFrame()
            df_name = df[df['diagnostic'] == 'microscopy'].reset_index()
            for id, data in df_name.groupby('individual_ID'):
                data = data.reset_index()
                patient_parasites = pd.DataFrame()

                for infection_day in range(data.shape[0]):
                    patient_parasites[data.day[infection_day]] = [data.parasitemia[infection_day]]

                infection_trajectories['microscopy'] = pd.concat(
                    [infection_trajectories['microscopy'], patient_parasites], axis=0)

            # nonan_df = pd.DataFrame()
            # for ind, patient in infection_trajectories['microscopy'].iterrows():
            #     parasitemia = pd.Series(patient).dropna()
            #     nonan_df = pd.concat([nonan_df,parasitemia],axis = 1)

            return infection_trajectories['microscopy']


        ref_path = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'data', 'Mozambique',
                                'parasitemia',
                                'asymptomatic_men_CID_2018', 'galatas_CID_2018_supp_parasitemias.csv')
        moz_df = pd.read_csv(ref_path)
        ref_df = collect_reference_trajectories(df=moz_df)

        return ref_df

    def get_analyzers(self):
        from analyzers.CompareMozambiqueInfectionsAnalyzer import CompareMozambiqueInfectionsAnalyzer
        return [CompareMozambiqueInfectionsAnalyzer(site=self)]
