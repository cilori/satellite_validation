#!/home/claysa/anaconda3/bin/python3

import numpy as np
import pandas as pd


#==============================================================================
# DESCRIPTION
#==============================================================================
#
# 27 Mar 2017
# Stephanie Clay
#
# Put in situ data from .txt into a data frame and write it to an Excel file.
# Input data must start with the following column names, in this order:
#    cruise_id ID DEPTH LATDEC LONGDEC YEAR DOY
# followed by columns containing data for other requested variables.
# 
# Format notes:
#    Separate columns by spaces.
#    The number of sample measurement columns does not matter.
#    Longitude must be a negative number (for the AZMP region).
#    Year must be a 4-digit integer.
#    NA values should be marked NA, NaN, or -NaN.


#==============================================================================
# VARIABLES TO CHANGE
#==============================================================================

# Path to the this python script and the input file.
inpath = 'C:/Users/ClaySA/Desktop/claysa/satellite_validation/01_in_situ_data/'

# Input filename
filename = 'chl_in_situ_gosl_1997-2019.txt'

# List of variables of datatype "string".
string_vars = ['method', 'datetime']

outfile = ''.join([inpath, 'in_situ_extNA_', filename, '_formatted.csv'])

# WARNING:
#   See 2 short sections of code below marked "POTENTIAL ISSUE #1" - might have
#   to uncomment them to make it work properly, depending on input file format.



#==============================================================================
# MAIN CODE
#==============================================================================

if (filename.endswith('.csv')):
    
    df = pd.read_csv(''.join([inpath, filename]))
    newcols = ['cruise_id','id','depth','latitude','longitude','year',
               'day_of_year']
    newcols.extend([x.lower() for x in df.columns[7:len(df.columns)]])
    df.columns = newcols
    
else:
    
    with open(''.join([inpath, filename]), 'r') as f:
        
        # Create column names.
        cols = f.readline().rstrip().lower().split(',')
        
        
        # POTENTIAL ISSUE #1 (if you uncomment this, you also need to uncomment the lines below)
        #cols = cols[0].split(" ")
        #cols = [i.strip('"') for i in cols]
        
        
        print(cols)

        cols[3] = 'latitude'
        cols[4] = 'longitude'
        cols[6] = 'day_of_year'
        
        print('Completed up to line...')
        # Loop through input .txt or .csv file, and add lines to the frame as rows.
        df = pd.DataFrame(np.nan, index=[0], columns=cols)
        for count,line in enumerate(f):
            if count%1000==0: print(count)
            
            line = line.rstrip().split(',')
            
            
            # POTENTIAL ISSUE #1 (note if you uncomment the lines above, you also need to uncomment this)
            #line = line[0].split(" ")
            #line = [i.strip('"') for i in line]
            


            df.loc[count,:] = line

# Convert columns to appropriate data types (ints must be converted
# individually later - could be NaN values coercing the data type to float).
exclude_strs = ['cruise_id','id']
exclude_strs.extend(string_vars)
for col in df.columns.difference(exclude_strs):
    df.loc[:,col] = pd.to_numeric(df.loc[:,col], errors='coerce')

# Sort rows by values in columns year, day of year, latitude, longitude,
# and depth, then reset index.
df.sort_values(['year', 'day_of_year', 'latitude', 'longitude', 'depth'], 
                 ascending=[True, True, True, True, True], inplace=True)
df.reset_index(drop=True, inplace=True)

# Change order of columns.
col_order = ['year', 'day_of_year', 'latitude', 'longitude', 'cruise_id', 
             'id', 'depth']
col_order.extend(cols[7:])
df = df[col_order]

# Add new columns to the front of the HPLC frame.
df.insert(0,'extracted',False)
df.insert(1,'num_of_extractions',0)
# Add a column of empty lists. If an HPLC record has multiple extractions,
# the extraction numbers can be appended to the list for that particular
# record/row of the HPLC dataframe.
df.insert(2,'extraction_numbers',pd.Series(np.empty((len(df), 0)).tolist()))

# Write resulting HPLC data frame to csv file.
df.to_csv(outfile)
