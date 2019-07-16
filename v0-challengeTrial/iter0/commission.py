import os, sys, copy, re
import math
import random
import pandas as pd
from pyDOE import lhs
import numpy as np
from scipy import ndimage
from scipy.misc import imresize
from scipy.stats import gamma

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.utils.builders.TemplateHelper import TemplateHelper
from dtk.utils.builders.ConfigTemplate import ConfigTemplate
from dtk.utils.builders.TaggedTemplate import CampaignTemplate, DemographicsTemplate

from history_matching import quick_read

sys.path.insert(0, os.path.join('..', '..', '..', 'Scripts')) # <-- Ew
from utils import make_campaign_template

iteration = int(re.search(r'iter(\d+)', os.getcwd()).group(1))

N_samples = 5000
N_rep_per_sample = 1

params = quick_read( os.path.join('..', 'Params.xlsx'), 'Params').set_index('Name')

plugin_files_dir = os.path.join('..', '..', '..', 'InputFiles', 'Zimbabwe', 'v3')
cfg = ConfigTemplate.from_file( os.path.join(plugin_files_dir, 'config.json') )
cpn = make_campaign_template( plugin_files_dir,
    with_risk_reduction_events = True,
    with_medical_male_circumcision_events = True,
    with_sti_coinfection_events = True,
    with_traditional_male_circumcision_events = True,
    with_commercial_sex_events = True,
    with_PrEP = False
)
demog = DemographicsTemplate.from_file( os.path.join(plugin_files_dir, 'Demographics.json') )
demog_pfa = DemographicsTemplate.from_file( os.path.join(plugin_files_dir, 'PFA_Overlay.json') )
demog_acc = DemographicsTemplate.from_file( os.path.join(plugin_files_dir, 'Accessibility_and_Risk3_IP_Overlay.json') )
demog_asrt = DemographicsTemplate.from_file( os.path.join(plugin_files_dir, 'Risk_Assortivity_Overlay.json') )

print 'SETTING STATIC PARAMETERS:'
static_params = {}

base_year = 1900
static_params.update({
    'Society__KP_Defaults.COMMERCIAL.Pair_Formation_Parameters.Formation_Rate_Constant': 0.4,
    'Base_Year': base_year,
    'PFA_Burnin_Duration_In_Days': (1975-base_year)*365.,
    'Base_Population_Scale_Factor' : 1/400.0,
    'Simulation_Duration': math.floor( (2017.1 - base_year)*365 ),
    'Simulation_Timestep': 365/12.0,
    'Report_HIV_ByAgeAndGender_Collect_Age_Bins_Data': [ 0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75 ],
})
#'Report_HIV_ByAgeAndGender_Collect_Age_Bins_Data': [ 0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 ],
# 'Report_HIV_ByAgeAndGender_Start_Year': 1900
# 'Demographic_Coverage__KP_Seeding': 0.05,

# Set static parameters
cfg.set_params( static_params )
cpn.set_params( static_params )
demog.set_params( static_params )
demog_pfa.set_params( static_params )
demog_acc.set_params( static_params )
demog_asrt.set_params( static_params )


# Prepare templates
templates = TemplateHelper()
table_base = {
    'ACTIVE_TEMPLATES': [cfg, cpn, demog, demog_pfa, demog_acc, demog_asrt],
    'TAGS': {'Scenario':'StatusQuo_Baseline'}
}

# Standard DTKConfigBuilder
config_builder = DTKConfigBuilder.from_files(
    os.path.join(plugin_files_dir, 'config.json'),
    os.path.join(plugin_files_dir, 'campaign_StatusQuo.json')
)

s2 = np.sqrt(2)/2.
rot = np.array([[s2, -s2], [s2, s2]]) # Rotation matrix

q = np.arange(8, 68, 2)
#q = np.linspace(0, 99, 100)
x, y = np.meshgrid(q, q)
pos = np.dstack((x, y))
posr = np.dot(pos, rot) # Rotated

# Marital
lower = 6
upper = 54
fn = os.path.join('..', 'ZimbabweAge1stMarriage_age6thru53_padded_Mcol_Frow_DJK.csv')
inc_1st_marriages_raw = pd.read_csv(fn, skipinitialspace=True, index_col=0)

# Pad to 0-99
inc_1st_marriages_raw = np.pad(inc_1st_marriages_raw, pad_width=[(lower,100-upper), (lower,100-upper)], mode='constant')

s = 2*inc_1st_marriages_raw.shape[0]
c = s/2

def map_sample_to_model_input(config_builder, sample_idx, replicate_idx, sample):
    table = copy.deepcopy(table_base)
    table['Run_Number'] = random.randint(0, 1e6)
    table['TAGS'].update({'__sample_index__': sample_idx, '__replicate_index__': replicate_idx, 'Run_Number':table['Run_Number']})
    table['TAGS'].update(sample)

    if 'Commercial Condom Min' in sample:
        value = sample.pop('Commercial Condom Min')
        table['Society__KP_Defaults.COMMERCIAL.Relationship_Parameters.Condom_Usage_Probability.Min'] = value

    if 'Commercial Condom Max' in sample:
        value = sample.pop('Commercial Condom Max')
        table['Society__KP_Defaults.COMMERCIAL.Relationship_Parameters.Condom_Usage_Probability.Max'] = value

    if 'Commercial Condom Mid' in sample:
        value = sample.pop('Commercial Condom Mid')
        table['Society__KP_Defaults.COMMERCIAL.Relationship_Parameters.Condom_Usage_Probability.Mid'] = value

    if 'Seeding Year' in sample:
        value = sample.pop('Seeding Year')
        table['Start_Year__KP_Seeding_Year'] = value

    if 'STI Per Act LOW MED' in sample:
        value = sample.pop('STI Per Act LOW MED')
        table['Choices__KP_AcquireSTI_LOW_MEDIUM.AcquireSTI'] = value
        table['Choices__KP_AcquireSTI_LOW_MEDIUM.NoTrigger'] = 1-value

    if 'STI Per Act HIGH' in sample:
        value = sample.pop('STI Per Act HIGH')
        table['Choices__KP_AcquireSTI_HIGH.AcquireSTI'] = value
        table['Choices__KP_AcquireSTI_HIGH.NoTrigger'] = 1-value

    if 'STI Frac Temporary' in sample:
        value = sample.pop('STI Frac Temporary')
        table['Choices__KP_Prob_STI_Temporary.Need_Delay_To_STI_Cure'] = value
        table['Choices__KP_Prob_STI_Temporary.Lifelong_STI'] = 1-value

    if 'LOG Base Infectivity' in sample:
        value = sample.pop('LOG Base Infectivity')
        table['Base_Infectivity'] = np.exp(value)

    if 'HIGH Size Mult' in sample:
        value = sample.pop('HIGH Size Mult')

        before = cpn.get_param('Time_Value_Map__KP_FSW.Values')[1][0]
        table['Time_Value_Map__KP_FSW.Values'] = [value*p for p in before]

        before = cpn.get_param('Time_Value_Map__KP_Clients.Values')[1][0]
        table['Time_Value_Map__KP_Clients.Values'] = [value*p for p in before]

    if 'Trans Form Period' in sample:
        value = sample.pop('Trans Form Period')
        table['Society__KP_Defaults.TRANSITORY.Pair_Formation_Parameters.Formation_Rate_Constant'] = 1./float(value)

    if 'Informal Form Period' in sample:
        value = sample.pop('Informal Form Period')
        table['Society__KP_Defaults.INFORMAL.Pair_Formation_Parameters.Formation_Rate_Constant'] = 1./float(value)

    if 'Marital Form Period' in sample:
        value = sample.pop('Marital Form Period')
        table['Society__KP_Defaults.MARITAL.Pair_Formation_Parameters.Formation_Rate_Constant'] = 1./float(value)

    if 'Male To Female Young' and 'Male To Female Old' in sample:
        young = sample.pop('Male To Female Young')
        old = sample.pop('Male To Female Old')
        table['Male_To_Female_Relative_Infectivity_Multipliers'] = [young, young, old]

    risk_reduction_fraction = sample.pop( 'Risk Ramp Annual Max' ) if 'Risk Ramp Annual Max' in sample else float('NaN')
    risk_ramp_rate = sample.pop( 'Risk Ramp Rate' ) if 'Risk Ramp Rate' in sample else float('NaN')
    risk_ramp_midyear = sample.pop( 'Risk Ramp Mid' ) if 'Risk Ramp Mid' in sample else float('NaN')

    for province in ['Bulawayo', 'Harare', 'Manicaland', 'Mashonaland', 'Masvingo', 'Matabeleland', 'Midlands']:
        key = '%s: LOW Risk' % province
        if key in sample:
            value = sample.pop(key)
            param = 'Initial_Distribution__KP_Risk_%s' % province
            table[param] = [value, 1-value, 0]

        key = '%s: TI Condoms Max' % province
        if key in sample:
            value = sample.pop(key)
            table['Society__KP_%s.TRANSITORY.Relationship_Parameters.Condom_Usage_Probability.Max' % province] = value
            table['Society__KP_%s.INFORMAL.Relationship_Parameters.Condom_Usage_Probability.Max' % province] = value

        key = '%s: PreStaging Loss' % province
        if key in sample:
            value = sample.pop(key)
            table['Choices__KP_PreStagingLoss_%s.HCTUptakePostDebut3' % province] = value
            table['Choices__KP_PreStagingLoss_%s.ARTStaging2' % province] = 1-value

        if not math.isnan(risk_reduction_fraction):
            param = 'Actual_IndividualIntervention_Config__KP_Medium_Risk_%s.Ramp_Max' % province
            table[param] = risk_reduction_fraction

        if not math.isnan(risk_ramp_rate):
            param = 'Actual_IndividualIntervention_Config__KP_Medium_Risk_%s.Ramp_Rate' % province
            table[param] = risk_ramp_rate

        if not math.isnan(risk_ramp_midyear):
            param = 'Actual_IndividualIntervention_Config__KP_Medium_Risk_%s.Ramp_MidYear' % province
            table[param] = risk_ramp_midyear

    if 'Rate Ratio' in sample:
        v = sample.pop('Rate Ratio')
        table['Pair_Formation_Parameters__KP_All.Extra_Relational_Rate_Ratio_Male'] = v
        table['Pair_Formation_Parameters__KP_All.Extra_Relational_Rate_Ratio_Female'] = v

    if 'Risk Assortivity' in sample:
        v = sample.pop('Risk Assortivity')
        table['Weighting_Matrix_RowMale_ColumnFemale__KP_RiskAssortivity'] = [
            [ v, 1-v, 0 ],
            [ 1-v, v, v ],
            [ 0, v, 1-v ] ]

    if set(['Informal PFA Diag A', 'Informal PFA Diag Loc', 'Informal PFA Diag Scale', 'Informal PFA Cross A', 'Informal PFA Cross Loc', 'Informal PFA Cross Scale']).issubset(set(sample.keys())):
        # Make the informal PFA
        diag_a = sample.pop('Informal PFA Diag A')
        diag_loc = sample.pop('Informal PFA Diag Loc')
        diag_scale = sample.pop('Informal PFA Diag Scale')
        cross_a = sample.pop('Informal PFA Cross A')
        cross_loc = sample.pop('Informal PFA Cross Loc')
        cross_scale = sample.pop('Informal PFA Cross Scale')

        diag_rv = gamma(a=diag_a, loc=diag_loc, scale=diag_scale)
        cross_rv = gamma(a=cross_a, loc=cross_loc, scale=cross_scale)

        X = diag_rv.pdf(x=posr[:,:,0])
        Y = cross_rv.pdf(x=-posr[:,:,1]) # <-- Note negative!

        Z = np.multiply(X, Y)
        L = Z.transpose().tolist() # <-- Transpose for input

        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Number_Age_Bins_Male'] = Z.shape[1]
        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Number_Age_Bins_Female'] = Z.shape[0]
        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Age_of_First_Bin_Edge_Male'] = q[0]
        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Age_of_First_Bin_Edge_Female'] = q[0]
        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Years_Between_Bin_Edges_Male'] = q[1] - q[0]
        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Years_Between_Bin_Edges_Female'] = q[1] - q[0]
        table['Society__KP_Defaults_All_Nodes.INFORMAL.Pair_Formation_Parameters.Joint_Probabilities'] = L

        '''
        print 'Informal PFA Diag A', diag_a
        print 'Informal PFA Diag Loc', diag_loc
        print 'Informal PFA Diag Scale', diag_scale
        print 'Informal PFA Cross A', cross_a
        print 'Informal PFA Cross Loc', cross_loc
        print 'Informal PFA Cross Scale', cross_scale
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        ax.contourf(x, y, Z)
        ax.plot(q,q,'r-', lw=2)
        plt.show()
        '''

    if set(['Marital Center Weight', 'Marital Diag InvRate', 'Marital Diag A', 'Marital Cross A', 'Martial Cross Loc', 'Marital Cross Scale']).issubset(set(sample.keys())):
        # Make the marital PFA
        center_weight = sample.pop('Marital Center Weight')
        diag_invrate = sample.pop('Marital Diag InvRate')
        diag_a = sample.pop('Marital Diag A')
        cross_a = sample.pop('Marital Cross A')
        cross_loc = sample.pop('Martial Cross Loc')
        cross_scale = sample.pop('Marital Cross Scale')

        diag_rv = gamma(a=diag_a, loc=0, scale=diag_invrate)
        cross_rv = gamma(a=cross_a, loc=cross_loc, scale=cross_scale)

        # Build kernel
        K = np.zeros((s,s))
        for i in range(s):
            for j in range(s-i-1,s):
                diag_dist = np.sqrt( 2*((i+j)/2. - c)**2 )
                diag_f = diag_rv.pdf(x=diag_dist)

                lat_dist = np.sqrt( (i-(i+j)/2.)**2 + (j-(i+j)/2.)**2 )
                theta = np.arctan2(lat_dist, diag_dist)
                if i > j:
                    theta *= -1

                K[i,j] = diag_f * cross_rv.pdf(x=theta)

        K[c,c] += center_weight*np.amax(K)

        convolved = ndimage.convolve(inc_1st_marriages_raw, K, mode='constant', cval=0.0)
        years_between_bins = 2 # Other values not tested, probably won't work
        smaller = imresize(convolved, 1/float(years_between_bins), interp='bicubic')
        first_bin = 6
        # Chop below first bin
        Z = smaller[first_bin/years_between_bins:, first_bin/years_between_bins:]
        L = Z.transpose().tolist() # <-- Transpose for input

        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Number_Age_Bins_Male'] = Z.shape[1]
        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Number_Age_Bins_Female'] = Z.shape[0]
        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Age_of_First_Bin_Edge_Male'] = first_bin
        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Age_of_First_Bin_Edge_Female'] = first_bin
        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Years_Between_Bin_Edges_Male'] = years_between_bins
        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Years_Between_Bin_Edges_Female'] = years_between_bins
        table['Society__KP_Defaults_All_Nodes.MARITAL.Pair_Formation_Parameters.Joint_Probabilities'] = L

        if False:
            print 'Marital Center Weight', center_weight
            print 'Marital PFA Diag A', diag_a
            print 'Marital PFA Diag InvRate', diag_invrate
            print 'Marital PFA Cross A', cross_a
            print 'Marital PFA Cross Loc', cross_loc
            print 'Marital PFA Cross Scale', cross_scale

            import matplotlib.pyplot as plt
            fig,ax = plt.subplots(1,1)
            #ax.contourf(x, y, Z)
            ax.imshow(Z)
            ax.plot([0,Z.shape[0]], [0,Z.shape[1]], 'r-', lw=2)

            plt.show()


    for param_name,p in params.iterrows():
        if param_name in sample and 'MapTo' in p:
            if isinstance(p['MapTo'],float) and math.isnan(p['MapTo']):
                continue
            table[p['MapTo']] = sample.pop( param_name )

    for name,value in sample.iteritems():
        print 'UNUSED PARAMETER:', name

    assert( len(sample) == 0 ) # All params used


    ###
    mod_tags = templates.mod_dynamic_parameters(config_builder, table)
    #print '\n'.join(mod_tags.keys())
    #table['TAGS'].update(sample.to_dict())

    return table['TAGS']

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
        samples_unconstrained = choose_and_scale_samples_unconstrained(num_samples)
        # Constraints:
        #  --> 'Male To Female Young' >= 'Male To Female Old'
        samples = samples_unconstrained.loc[ samples_unconstrained['Male To Female Young'] >= samples_unconstrained['Male To Female Old'], :]
        #samples = samples_unconstrained

        remaining = num_samples - samples.shape[0]
        while remaining > 0:
            samples_unconstrained = choose_and_scale_samples_unconstrained(num_samples)
            samples = pd.concat([samples, samples_unconstrained.loc[ samples_unconstrained['Male To Female Young'] >= samples_unconstrained['Male To Female Old'], :] ], ignore_index=True)
            remaining = num_samples - samples.shape[0]

        samples.index.name = 'Sample'
        return samples.iloc[:num_samples,:]
    else:
        #samples = pd.read_excel(os.path.join('..', 'iter%d'%(iteration-1), 'Candidates_for_iter%d.xlsx'%iteration), sheetname='Values')
        samples = pd.read_hdf(os.path.join('..', 'iter%d'%(iteration-1), 'Candidates_NS_for_iter%d.hd5'%iteration), key='values')
        print 'CANDIDATES:\n', samples.head()
        samples.index.name = 'Sample'
        assert( pd.Series.all( samples['Male To Female Young'] >= samples['Male To Female Old'] ))

        return samples.iloc[:num_samples,:]


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
                 'exp_name': 'HIV-Zim BHMv3 Iter%d'%iteration}




# Notes
# * Contraceptive mult baked into M-->F mult
# * Neglecting pregnancy
# * MTCT not right due to low maternal prevalence
# * STI for life

    # CSW POPULATION SIZE!!!
    # STI PARAMS
    # Rate ratio?

    # RISK reduction - really?


'''
Circumcision Effectiveness	0.5	0.6	Circumcision_Reduced_Acquire	Boily: Siegfried Cochrane Review (2013), Mehta AIDS (2013)
ART Effectiveness	0.64	0.92	ART_Viral_Suppression_Multiplier	Boily
Condom Effectiveness	0.6	0.95	Condom_Transmission_Blocking_Probability	Giannou, Cochrane, Liu
    {
        'Name': 'Max Trns M&F MED',
        #'MapTo': [  'Society__KP_Defaults_All_Nodes.TRANSITORY.Concurrency_Parameters.MEDIUM.Max_Simultaneous_Relationships_Male',
                    'Society__KP_Defaults_All_Nodes.TRANSITORY.Concurrency_Parameters.MEDIUM.Max_Simultaneous_Relationships_Female' ],
        'Min': 1.0,
        'Max': 3.0
    },
    {
        'Name': 'Max Infmrl M&F MED',
        #'MapTo': [  'Society__KP_Defaults_All_Nodes.INFORMAL.Concurrency_Parameters.MEDIUM.Max_Simultaneous_Relationships_Male',
                    'Society__KP_Defaults_All_Nodes.INFORMAL.Concurrency_Parameters.MEDIUM.Max_Simultaneous_Relationships_Female' ],
        'Min': 1.0,
        'Max': 3.0
    },
    {
        'Name': 'Acute Duration',
        'MapTo': 'Acute_Duration_In_Months',
        'Min': 0.55,
        'Max': 6.8,
        'Source': 'Boily/Bellan'
    },
    {
        'Name': 'Acute RR',
        'MapTo': 'Acute_Stage_Infectivity_Multiplier',
        'Min': 1,
        'Max': 30,
        'Source': 'Boily, Bellan'
    },
    ]
'''