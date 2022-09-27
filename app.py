import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

def parse_table(df):
    data = pd.read_csv(df, header = 4)
    final_col = data.columns[-1]
    if 'Unnamed' in final_col:
        data = data.drop(final_col, axis = 1)
    if data.Content[0][:4] == 'Time':
        data = data.drop('Content', axis = 1)
        times = list(data.iloc[0])[1:]
        data = data.drop(0)
    else:
        data = data.drop('Content', axis = 1)
        times = 0
    rows = data.apply(lambda x: x['Well'][0], axis = 1)
    data['Row'] = rows
    columns = data.apply(lambda x: int(x['Well'][1:]), axis = 1)
    data['Column'] = columns
    data = data.set_index(['Row', 'Column']).drop('Well', axis = 1)
    data.columns = [[x.split('.')[0] for x in data.columns], times]
    return data, set(times), rows, columns

concentration_units = {
    'molar': 1,
    'millimolar': 10**-3,
    'micromolar': 10**-6,
    'nanomolar': 10**-9
}

st.title('TR-FRET Data Wrangler')
'''
I've made a thing! It takes the .csv file from the PheraStar and separates out the wells you've used as you specify.
At the minute, it only does this for competition curves run with replicates in separate columns. I'll add functionality for other experiments one day.
'''
uploaded_file = st.file_uploader('Upload the data for your experiment in .csv format')
time_delay = st.number_input('First time point', min_value=0, value = 0)
time_units = st.selectbox('Time Units:', ['seconds', 'minutes', 'hours'])
test_compound_string = st.text_input('Test Compounds (separate with commas)')

if uploaded_file is not None:
    tr_fret_data, tr_fret_times, available_rows, available_columns = parse_table(uploaded_file)
    timepoints = len(tr_fret_times)
    ordered_times = [x + time_delay for x in sorted(list(tr_fret_times))]
    test_compounds = test_compound_string.split(', ')
    compound_dict = {}

    st.text(f'{timepoints} time points:')
    st.text(ordered_times)
    st.header('Options')
    '''
    Once you've specified which wells your experiment was run in, it'll display your data below. This can be copied and pasted into a spreadsheet.

    You have to specify concentrations in the same order you specify the rows used for that concentration. This is done separately for each compound tested.

    Don't worry if the concentrations appear to be zero, if you paste them elsewhere, they'll show properly. I tried to get this to display scientific notation but it just wouldn't work.
    '''
    for compound in test_compounds:
        compound_dict[compound] = {}
        st.subheader(f'{compound} options')
        compound_dict[compound]['units'] = st.selectbox(f'{compound} Concentration Units', ['nanomolar', 'micromolar', 'millimolar', 'molar'])
        compound_dict[compound]['raw_concs'] = st.text_input(f'{compound} Compound concentrations (separate with commas)')
        if compound_dict[compound]['raw_concs'] is not '':
            if '.' in compound_dict[compound]['raw_concs']:
                conc_numbers = [float(x) for x in compound_dict[compound]['raw_concs'].split(', ')]
            else:
                conc_numbers = [int(x) for x in compound_dict[compound]['raw_concs'].split(', ')]
            compound_dict[compound]['concentrations'] = np.array(conc_numbers)*concentration_units[compound_dict[compound]['units']]
        compound_dict[compound]['replicate cols'] = st.multiselect(f'Columns containing {compound} replicates', available_columns.unique())
        compound_dict[compound]['conc rows'] = st.multiselect(f'Rows containing {compound} concentrations', available_rows.unique())
       
        if 'concentrations' in compound_dict[compound].keys():
            if len(compound_dict[compound]['concentrations']) > len(compound_dict[compound]['conc rows']):
                st.text(f'Select more rows for {compound} concentrations')
            elif len(compound_dict[compound]['concentrations']) < len(compound_dict[compound]['conc rows']):
                st.text(f"You've selected too many rows for {compound} concentrations")
            else:
                conc_df = pd.DataFrame({
                    'Row': compound_dict[compound]['conc rows'],
                    'Concentration': compound_dict[compound]['concentrations']
                })
                reps_df = pd.DataFrame({
                    'Column': compound_dict[compound]['replicate cols'],
                    'Replicate': [f'Replicate {str(x+1)}' for x in range(len(compound_dict[compound]['replicate cols']))]
                })
                time_dfs = dict([(f'{t+time_delay} {time_units}', tr_fret_data.xs(t, axis = 1, level = 1)) for t in sorted(list(tr_fret_times))])
                compound_dict[compound]['display timepoint'] = st.selectbox(f'Choose a time-point to display for {compound}', time_dfs.keys())
                
                compound_dict[compound]['timepoints'] = dict([
                    (t, df.reset_index()
                    .merge(conc_df, on = 'Row')
                    .merge(reps_df, on = 'Column')
                    .drop(['Row', 'Column'], axis = 1)
                    .pivot(index = 'Concentration', columns = 'Replicate'))
                    for t, df in time_dfs.items()
                ])
                display_table = compound_dict[compound]['timepoints'][compound_dict[compound]['display timepoint']]
                fret_ratio = display_table.groupby(['Replicate'], axis = 1).apply(lambda x: x['Raw Data (337/665 A)']/x['Raw Data (337/620 B)'])
                mean_f_r = fret_ratio.apply(np.mean, axis = 1)
                sd_f_r = fret_ratio.apply(lambda x: np.std(x)/sqrt(len(x)), axis = 1)
                fig, ax = plt.subplots()

                ax.errorbar(mean_f_r.index, mean_f_r, yerr = sd_f_r, color = 'black', fmt = 'o', capsize = 3, ls = 'None')

                ax.set_xscale('log')
                ax.set_xlabel(f'[{compound}] (M)')
                ax.set_ylabel('FRET ratio')

                st.dataframe(display_table)
                st.pyplot(fig)
