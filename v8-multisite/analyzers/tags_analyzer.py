from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
import pandas as pd


class TagsAnalyzer(BaseAnalyzer):

    def __init__(self, **kwargs):
        self.filter_tags = kwargs['filter_tags'] if 'filter_tags' in kwargs else {}
        super().__init__(**kwargs)
        print('WD:', self.working_dir)


    def filter(self, sim_metadata):
        if '__replicate_index__' in sim_metadata.tags:
            return sim_metadata.tags['__replicate_index__'] == 0
        return True


    def select_simulation_data(self, data, simulation):
        print('SELECT:', simulation)
        print(dir(simulation))
        print(type(simulation.tags))

        print(simulation.tags)
        tags = {t[9:]:v for t,v in simulation.tags.items() if 'SAMPLE' in t}
        tags['Sim_Id'] = simulation.id
        tags['Run_Number'] = simulation.tags['Run_Number']
        tags['__sample_index__'] = simulation.tags['__sample_index__']
        print('TAGS:', tags)
        return tags

    def finalize(self, all_data):
        print('All Data:\n', all_data)
        print('All Data:\n', all_data.keys())
        all_results = pd.DataFrame(all_data).transpose()
        all_results.reset_index(drop=True)
        all_results.set_index('__sample_index__', inplace=True)
        all_results.index.name = 'Sample'

        fn = 'Results_%s.csv'%(self.__class__.__name__)
        print('--> Writing %s' % fn)
        all_results.to_csv(fn)

if __name__ == '__main__':

    am = AnalyzeManager('9be99541-56ff-e911-a2c3-c4346bcb1551',analyzers=TagsAnalyzer())
    am.analyze()

