import os, re
import logging
import pandas as pd
import numpy as np      # for finfo tiny


# Plotting
import matplotlib.pyplot as plt
import seaborn as sns


from calibtool import LL_calculators
from simtools.Analysis.BaseAnalyzers import  BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager


class InfectionDurationsByAgeAnalyzer(BaseAnalyzer):
    def __init__(self,reference_path):
        super().__init__(filenames=['output/MalariaPatientReport.json'])
        self.channels = ['initial_age', 'true_asexual_parasites', 'true_gametocytes','temps']
        self.reference_path = reference_path

    def select_simulation_data(self, data, simulation):
        patient_data = data[self.filenames[0]]['patient_array']
        patient_dict = {'patient %d' % x['id']: {channel : x[channel] for channel in self.channels} for x in patient_data}
        patient_df = pd.DataFrame.from_dict(patient_dict,orient='index')
        patient_df.reset_index(drop=True)
        #add age as a feature to patient_df
        #patient id .. age .. asexual list .. gametocyte list .. fever list

        return patient_df

    def finalize(self, all_data):
        #first extract the age bins from the reference data structure
        ref_df = pd.read_excel(self.reference_path)
        age_bins = ref_df['Age bin'].values * 365


 ######
        results = pd.DataFrame()

        for sim, data in all_data.items():
            #for each sim I want to group by each age bin and then calculate a mean infection duration and then append to a dataframe a row that has the columns:
            # sample. simid. category. and value from sim

            sample = sim.tags.get("__sample_index__")
            sim_id = sim.id
            for age in range(len(age_bins[:-1])):
                sim_subset = data[(data['initial_age'] >= age_bins[age])&(data['initial_age'] < age_bins[age+1])]
                simulation_durations = []
                for iid, patient in sim_subset.iterrows():
                    (positive_day_indices,) = np.where(np.array(patient['true_asexual_parasites'][0]) > 0)
                    try:
                        simulation_durations.append(max(positive_day_indices))
                    except:
                        print('')

                age_cat = f'{age_bins[age]/365}-{age_bins[age+1]/365}'
                value = np.mean(simulation_durations)
                sub_results = pd.DataFrame({'sample' : sample,
                                            'sim_id': sim_id,
                                            'age_cat': age_cat,
                                            'value': [value]})
                #Could log the individual id, ind age and ind duration.

                results = pd.concat([results,sub_results])
            #loop through all the age bins
            #evaluate their mean duration for each age bin
        results.to_csv(os.path.join('..', 'iter0', 'analyzer_results.csv'))





if __name__ == '__main__':
    ref_data_path = os.path.join('..', 'reference data', 'GhanaDurations.xlsx')
    am = AnalyzeManager('8435ba1a-12b8-e911-a2c1-c4346bcb1555', analyzers=InfectionDurationsByAgeAnalyzer(reference_path= ref_data_path))
    am.analyze()

#move samples and results into a exp_id dir in iter0 after every iter



