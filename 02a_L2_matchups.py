#!/home/claysa/anaconda3/bin/python3

import numpy as np
import netCDF4 as nc4
import datetime as dt
import pandas as pd
import sys
import os
import glob
import pyproj
import operator as op
import time

#==============================================================================
# FUTURE ADDONS/FIXES
#==============================================================================
#
# If a matchup netCDF contains matchups from multiple in situ files, they won't necessarily be connected to the right record numbers so you have no way of knowing which input file they came from


#==============================================================================
# DESCRIPTION
#==============================================================================
#
# 6 Jul 2017
# Stephanie Clay
# 
#
# DESCRIPTION:
# Input in_situ_data file and find matches in NASA's satellite files (.nc format).
# WORKS FOR SEAWIFS, MODIS-AQUA, AND VIIRS
# 
#
# INPUT IN SITU DATA:
# in_situ_data must be in .txt format, starting with the following columns in this
# order:
# cruise_id, ID, DEPTH, LATDEC, LONGDEC, YEAR, DOY
# followed by any sample measurement columns, such as HPLC or absorption. The
# number of columns does not matter (do not add extra 'dummy' columns).
# Format notes:
#    Longitude must be a negative number (for the AZMP region).
#    Year must be a 4-digit integer.
#    NA values should be marked NA, NaN, -NaN.
# IMPORTANT: VALUES MUST BE SEPARATED BY COMMAS
# 
#
# INPUT SATELLITE DATA:
#   Level 2 netCDF files downloaded from OBPG.
# 
#
# SCRIPT INSTRUCTIONS:
# First line of this script must contain #![location of python]
# example: #!/usr/bin/python
#       or #!/usr/bin/python3
#       or #!/home/claysa/anaconda3/bin/python3
# NOTE: This program was built with Python3, NOT Python2.
# Change variables as needed in the "VARIABLES TO CHANGE" section below.
# To run a program, type in Cygwin: python directoryname/programname.py
# (If that does not work, try python3)
# 
# 
# NetCDF file notes:
    # Shape of numpy array of 2D variables given by:
    #    (line,pixel) = (y,x) = (row,column)
    # In Seadas, the grid might be flipped depending on whether the satellite
    #    pass was ascending or descending (info found in global attributes).
    #    If viewing the original .nc file of a match in Seadas, might need to
    #    check the coordinate given by:
    #    (number_of_lines - 1 - line, pixels_per_line - 1 - pixel)
    #    The index of the grid begins at 0, which is why you must subtract 1
    #    from the total size of each dimension.


#==============================================================================
# VARIABLES TO CHANGE
#==============================================================================

# Location of this python script, the input data file, and the output .nc file.
prog_directory = '/home/claysa/'

# Type of in_situ data to be matched with satellite files. If it's not something
# specific (e.g. hplc, absorption, turner), you can just use "in_situ"
in_situ_dataname = 'in_situ'

# Enter a list of column names after cruise_id, ID, DEPTH, LATDEC, LONGDEC, YEAR, and DOY.
# NOTE: THESE MUST BE LOWER-CASE, even if they are upper-case in the in situ file.
in_situ_cols = ["chl", "method", "nbsamples", "datetime"]
# List of in_situ_cols variables of datatype "string". (If none, leave blank - i.e. [])
# Again, keep them lower-case.
string_vars = ["method", "datetime"]

# Set of rows to check in the (SORTED) in_situ_data dataframe (first row = 0).
# To view sorted version, run .txt or .csv file through in_situ_data_to_excel.py,
# which is the in_situ_data_toframe() method within its own .py file.
# NOTE: This is used in a Pandas DataFrame, where endpoints are inclusive,
# as opposed to regular Python slicing where they are exclusive.
# ALSO KEEP IN MIND: Python uses zero-indexing
startrow = 0
endrow = 7102

# For first match: distance from data point to "exact" match (in metres).
exact_match_dist = 500

# Maximum radius (in metres) around data point to search for matches.
max_dist = 20000

# How to deal with a situation where an output file of the same name
# already exists.
# If False: program quits (does not append to already existing file).
# If True: program checks if already existing file contains the same variables;
#          if it does, it appends to it. If not, it quits.
append = True

# Possible sensors: modis, viirs, seawifs.
sensor = 'modis'

# Window of pixels extracted around matching point.
# window_size (one side of the square window) must be an odd integer > 0.
window_size = 5 # default = 3
# NOTES: Number of good points must be <= (window_size)^2, and low number of
# good points must be <= high number of good points.
# For second match: lowest number of good (non-NaN) pixels in window.
low_good_points = 3 # default = 3
# For third match: second lowest number of good (non-Nan) pixels in window.
high_good_points = 9 # default = 9

# Enter a list of flags in all capitals. To select no flags, type {}
# Example using the default l2 flags:
# l2_flags = ['ATMFAIL','BOWTIEDEL','LAND','HIGLINT','CLDICE','HILT','HISOLZEN','HISATZEN']
l2_flags = ['ATMFAIL','BOWTIEDEL','LAND','HIGLINT','CLDICE','HILT','HISOLZEN','HISATZEN']

# in_situ_data_filename must include the extension (.txt or .csv), and must be
# located in the same folder as the program.
# If the file was already formatted using in_situ_data_to_excel.py, make sure
# the filename contains the string "formatted".
in_situ_data_filename = 'chl_in_situ_gosl_1997-2019.txt'

# Output NetCDF file. Filename can only include spaces, alphanumeric characters,
# or the following: ( ) \' - _ . ,'
netcdf_filename = ''.join([in_situ_data_filename,'_extNA_',sensor,'_matches'])

# Location of input netCDF L2 satellite files, which MUST BE SORTED INTO YEAR SUBFOLDERS
input_dir = '/Volumes/Data-panix/modis-aqua/NASA/R2018.0/L2_LAC_OC.extNA'

# WARNING:
# The in_situ_data_toframe() method has 2 sections of code marked "POTENTIAL ISSUE #1"
# which might need to be uncommented depending on the format of  the input file.
# Try formatting it with in_situ_data_to_excel.py first to see how it behaves, and
# when you figure out if it needs to be uncommented or not, do the same to this script.



#==============================================================================
# DEFAULTS AND OTHER VARIABLES
#==============================================================================

# List of 2D variables for the output .nc file (MUST have at least one).
# Example: output_vars = ['rrs_469','rrs_531']
# These are the available variables for each sensor:
if sensor=="modis":
    output_vars = ['Rrs_412','Rrs_443','Rrs_469','Rrs_488','Rrs_531','Rrs_547', 
                   'Rrs_555','Rrs_645','Rrs_667','Rrs_678','chlor_a','Kd_490',
                   'latitude','longitude']
elif sensor=="seawifs":
    output_vars = ['Rrs_412','Rrs_443','Rrs_490','Rrs_510','Rrs_555','Rrs_670', 
                   'chlor_a','Kd_490','latitude','longitude']
elif sensor=="viirs":
    output_vars = ['Rrs_410','Rrs_443','Rrs_486','Rrs_551','Rrs_671','chlor_a',
                   'Kd_490','latitude','longitude']

# Full collection of l2 flags for SeaWiFS, MODIS, and VIIRS.
l2_flags_dict = {'PRODWARN': 4, 'HIPOL': 536870912, 'LOWLW': 16384, 
        'SEAICE': 16777216, 'ATMFAIL': 1, 'BOWTIEDEL': 268435456, 'LAND': 2, 
        'HIGLINT': 8, 'NAVFAIL': 33554432, 'FILTER': 67108864, 
        'MODGLINT': 1048576, 'PRODFAIL': 1073741824, 'CLDICE': 512, 
        'HILT': 16, 'CHLFAIL': 32768, 'NAVWARN': 65536, 'ABSAER': 131072, 
        'HISOLZEN': 4096, 'COCCOLITH': 1024, 'SPARE': 2147483648, 
        'HISATZEN': 32, 'COASTZ': 64, 'ATMWARN': 4194304, 'TURBIDW': 2048, 
        'CHLWARN': 2097152, 'MAXAERITER': 524288, 'STRAYLIGHT': 256}

# 1D variables for the output .nc file.
output = {'/netcdf4_data/': ['pass_id'],
          '/processing_data/': ['extraction_number', 'matched_pass',
                                'window_center_line_zero_based',
                                'window_center_pixel_zero_based',
                                'window_center_original_line_zero_based',
                                'window_center_original_pixel_zero_based',
                                'distance', 'nctime', 'time_str']}
in_situ_path = ''.join(['/', in_situ_dataname, '_data/', in_situ_dataname, '_'])
base_cols = ['record_number', 'year', 'day_of_year', 'latitude', 'longitude',
             'cruise_id', 'id', 'depth']
base_cols.extend(in_situ_cols)
output[in_situ_path] = base_cols



#==============================================================================
# Check for errors in input.

if not os.path.isfile(prog_directory + in_situ_data_filename):
    sys.exit(''.join([in_situ_data_filename,' does not exist in the current directory']))

if not netcdf_filename.endswith('.nc'):
    netcdf_filename = ''.join([netcdf_filename, '.nc'])

sensor = sensor.lower()
if not sensor in ['modis', 'viirs', 'seawifs']:
    sys.exit('Error: Unknown sensor (enter modis, seawifs, or viirs)')

if (window_size%2 == 0) or (window_size < 1):
    sys.exit('Error: Window size must be an odd number > 0')

if not (0 < low_good_points <= window_size**2):
    sys.exit('Enter a low number of good points between 0 and window_size^2')

if not (low_good_points <= high_good_points <= window_size**2):
    sys.exit('Enter a high number of good points between low_good_points and window_size^2')

if not output_vars:
    sys.exit('Enter at least one output variable.')

valid_sensor_vars = {
    'modis': ['Rrs_412','Rrs_443','Rrs_469','Rrs_488','Rrs_531','Rrs_547', 
               'Rrs_555','Rrs_645','Rrs_667','Rrs_678','chlor_a','Kd_490',
               'latitude','longitude'],
    'seawifs': ['Rrs_412','Rrs_443','Rrs_490','Rrs_510','Rrs_555','Rrs_670', 
                'chlor_a','Kd_490','latitude','longitude'], 
    'viirs': ['Rrs_410','Rrs_443','Rrs_486','Rrs_551','Rrs_671','chlor_a', 
              'Kd_490','latitude','longitude']}
valid_sensor_vars = valid_sensor_vars[sensor]

for var in output_vars:
    if not var in valid_sensor_vars:
        sys.exit(''.join(['Input error: ', var, ' not a recognized '
                          'variable for given sensor']))

# Add latitude and longitude to the 2D output variables for the sensor.
if not 'latitude' in output_vars:
    output_vars.append('latitude')
if not 'longitude' in output_vars:
    output_vars.append('longitude')

if (len(l2_flags) > 0):
    flags = l2_flags
    l2_flags = {}
    for flag in flags:
        # If too many commas used, ignore blank space between them.
        if not flag or flag.isspace(): continue
        flag = flag.strip()
        # Ignore duplicate flags.
        if flag in l2_flags.keys(): continue
        # Check if input is a valid l2 flag.
        if not flag in l2_flags_dict.keys():
            sys.exit(''.join(['Input error: unrecognized flag: ', flag]))
        l2_flags[flag] = l2_flags_dict[flag]

# Check if output netCDF file already exists and contains the same variables.
file_exists = False
if os.path.isfile(prog_directory + netcdf_filename):
    if not append:
        sys.exit(''.join(['Output filename ', netcdf_filename, ' already exists.']))
    else:
        with nc4.Dataset(netcdf_filename,'r') as nc:
            nc_sensor = nc.getncattr('sensor').lower()
            nc_output_vars = nc.getncattr('output_vars').split(', ')
            nc_window_size = np.int32(nc.getncattr('window_size')[0])
            nc_low_good_points = np.int32(nc.getncattr('low_min_good_points'))
            nc_high_good_points = np.int32(nc.getncattr('high_min_good_points'))
            nc_l2_flag_meanings = nc.getncattr('l2_flag_meanings').split(', ')
        # Check if the variables from the existing file match with the
        # current user-selected variables.
        if ((sensor==nc_sensor) &
            (set(output_vars)==set(nc_output_vars)) &
            (window_size==nc_window_size) &
            (low_good_points==nc_low_good_points) &
            (high_good_points==nc_high_good_points) &
            (set(l2_flags.keys())==set(nc_l2_flag_meanings))):
            
            file_exists = True
            
            with nc4.Dataset(netcdf_filename,'a') as nc:
                # If appending to a file, check if it was created with the same
                # data file, and if not, add this data file's name to the list.
                d = nc.getncattr('in_situ_data_filename')
                if d.find(in_situ_data_filename) == -1:
                    newname = ', '.join([d, in_situ_data_filename])
                    nc.setncattr('in_situ_data_filename', newname)
        else:
            sys.exit(''.join(['Output filename ', netcdf_filename, ' already ',
                              'exists and contains different variables.']))


#==============================================================================
# Put in_situ_data into a Pandas DataFrame. Convert data to appropriate data types,
# and sort columns and rows to make it easier to extract and compare data to
# .nc files.
# Input data must start with the following column names, in this order:
#    cruise_id ID DEPTH LATDEC LONGDEC YEAR DOY
# followed by columns containing data for other requested variables.
# Format notes:
#    Longitude must be a negative number (for the AZMP region).
#    Year must be a 4-digit integer.
#    NA values should be marked NA, NaN, -NaN, or -999.
# NOTE: this method exists separately in its own program too: in_situ_data_to_excel.py

def in_situ_data_toframe():
    
    # Get a list of variables that should be strings.
    table_strs = ['cruise_id','id']
    table_strs.extend(string_vars)
    
    if 'formatted' in in_situ_data_filename:
        
        # Get the table.
        df = pd.read_csv(''.join([prog_directory + in_situ_data_filename]),index_col=0)
        
        # Make sure 'year' and day of year columns are int.
        df.loc[:,'year'] = pd.to_numeric(df.loc[:,'year'], errors='coerce')
        df.loc[:,'day_of_year'] = pd.to_numeric(df.loc[:,'day_of_year'], errors='coerce')
        
        # Make sure cruise_id, id, and string_vars columns are string/object dtype.
        for col in table_strs:
            df[col] = df[col].astype(str)
        
        # Replace "extraction_nuumbers" column with a column of empty lists (reading it
        # from a formatted file will result in a column of strings instead).
        df = df.drop(labels='extraction_numbers',axis=1)
        df.insert(2,'extraction_numbers',pd.Series(np.empty((len(df),0)).tolist()))
        
    else:
    
        with open(''.join([prog_directory + in_situ_data_filename]), 'r') as f:
        
            # Create column names.
            cols = f.readline().rstrip().lower().split(',')
            
            
            # POTENTIAL ISSUE #1 (if you uncomment this, you also need to uncomment the lines below)
            # # DEPENDING ON INPUT FILE FORMATTING, YOU MIGHT NOT NEED THESE LINES (same with lines in loop just below)
            # cols = cols[0].split(" ")
            # cols = [i.strip('"') for i in cols]
            
            
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
                # # DEPENDING ON INPUT FILE FORMATTING, YOU MIGHT NOT NEED THESE LINES
                # line = line[0].split(" ")
                # line = [i.strip('"') for i in line]
                
                
                df.loc[count,:] = line
        
        # Convert columns to appropriate data types (ints must be converted
        # individually later - could be NaN values coercing the data type to float).
        for col in df.columns.difference(table_strs):
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
        
        # Add new columns to the front of the data frame.
        df.insert(0,'extracted',False)
        df.insert(1,'num_of_extractions',0)
        # Add a column of empty lists. If a in_situ_data record has multiple extractions
        # the extraction numbers can be appended to the list for that particular
        # record/row of the in_situ_data dataframe.
        df.insert(2,'extraction_numbers',pd.Series(np.empty((len(df),0)).tolist()))
        
    return df


#==============================================================================
# Create the output .nc file that will contain the data from in_situ_data rows that
# match with latitude/longitude points in certain netCDF files, and info from
# the corresponding netCDF files, including a pixel window around the matching
# point with data from all selected variables.

def createnetcdf():
    
    with nc4.Dataset(netcdf_filename, 'w', format='NETCDF4') as n:
        
        # Add global attributes to the output netcdf file.
        n.in_situ_data_filename = in_situ_data_filename
        n.sensor = sname
        n.platform = platform
        n.output_vars = ', '.join(output_vars)
        l2_flag_meanings = []
        l2_flag_masks = []
        # Sort l2 flags by bit, then add the l2 flag meanings to a string, and
        # the corresponding bits to a numpy array.
        for meaning,mask in sorted(l2_flags.items(), key=op.itemgetter(1)):
            l2_flag_meanings.append(meaning)
            l2_flag_masks.append(mask)
        n.l2_flag_meanings = ', '.join(l2_flag_meanings)
        n.l2_flag_masks = np.asarray(l2_flag_masks)
        n.window_size = ''.join([str(window_size), 'x', str(window_size)])
        n.window_length = window_size**2
        n.low_min_good_points = low_good_points
        n.high_min_good_points = high_good_points
        
        n.createDimension('length', window_length)
        n.createDimension('records', None)
        
        n.createGroup('netcdf4_data')
        n.createGroup(''.join([in_situ_dataname, '_data']))
        n.createGroup('processing_data')
        
        # 1D variables.
        
        # Create lists of 1D output variables by type.
        str_vars = ['pass_id', 'cruise_id', 'id', 'matched_pass', 'time_str']
        str_vars.extend(string_vars)
        num_vars = ['distance', 'nctime', 'latitude', 'longitude'] + [x for x in output[in_situ_path][7:] if x not in string_vars]
        
        for path,varlist in output.items():
            for var in varlist:
                varname = path + var
                # String variables in output.
                if var in str_vars:
                    n.createVariable(varname, str, ('records',))
                # Float variables in output, and their units.
                elif (var in num_vars):
                    x = n.createVariable(varname, 'f4', ('records',))
                    if var == 'latitude':
                        x.units = 'degree_north'
                    elif var == 'longitude':
                        x.units = 'degree_east'
                    elif var == 'distance':
                        x.units = 'metres'
                    elif var == 'nctime':
                        x.units = 'seconds since 1970-01-01 00:00:00.0'
                # Integer variables in output.
                else:
                    n.createVariable(varname, 'i4', ('records',))
        
        # 2D variables, and their units.
        for ovar in output_vars:
            varname = '/processing_data/' + ovar
            x = n.createVariable(varname, 'f4', ('records', 'length'),
                                 fill_value=np.nan)
            if ovar == 'latitude':
                x.units = 'degree_north'
            elif ovar == 'longitude':
                x.units = 'degree_east'
            elif ovar == 'chlor_a':
                x.units = 'mg m^-3'
            elif ovar == 'Kd_490':
                x.units = 'm^-1'
            else:
                x.units = 'sr^-1'


#==============================================================================
# If match(es) have been found, append them to the output .nc file with data
# from the corresponding in_situ_data row/record.

def addtoOutput(in_situ_dataindex, in_situ_datarow, filename, matches):
    
    f, n = filename, netcdf_filename
    wl = window_length
    
    with nc4.Dataset(f, 'r') as inputnc, nc4.Dataset(n, 'a') as outputnc:
        
        # Get dimensions of .nc file containing the pixel that matches the
        # current in_situ_data row.
        number_of_lines = inputnc.dimensions['number_of_lines'].size
        pixels_per_line = inputnc.dimensions['pixels_per_line'].size
        
        # Get current number of records of output .nc file.
        nextrow = outputnc.dimensions['records'].size
        
        # Extract pass_id, time variables, and their units.
        pass_id = inputnc.getncattr('id')
        centre_line = number_of_lines//2
        path = '/scan_line_attributes'
        year = int(inputnc[path].variables['year'][centre_line])
        doy = int(inputnc[path].variables['day'][centre_line])
        msec = int(inputnc[path].variables['msec'][centre_line])
        time_units = outputnc['/processing_data'].variables['nctime'].units
        # Process time variables.
        month = int((dt.date(year,1,1) + dt.timedelta(doy-1)).strftime('%m'))
        day = int((dt.date(year,1,1) + dt.timedelta(doy-1)).strftime('%d'))
        time_str = (dt.datetime(year, month, day)
                    + dt.timedelta(milliseconds=msec))
        time_str = time_str.strftime('%Y-%m-%d %H:%M:%S')
        nctime = dt.datetime(year,month,day) + dt.timedelta(milliseconds=msec)
        nctime = nc4.date2num(nctime, units=time_units)
        
        # If netCDF has global attributes that give information on the original
        # pixel/line values (before extraction from NASA files), add them to
        # the output netCDF.
        try:
            l2_line = inputnc.getncattr('l2extract_start_line_zero_based')
            l2_pix = inputnc.getncattr('l2extract_start_pixel_zero_based')
            l2_exists = True
        except:
            l2_exists = False
        
        # Get 2D variables from input .nc file.
        for count,match in enumerate(matches, start=nextrow):
        
            line, pixel = match[0], match[1]
            distance = np.float32(match[2])
            maskframe = match[3]
            # Window borders, top/bottom.
            v1, v2 = match[5], match[6]
            # Window borders, left/right.
            h1, h2 = match[7], match[8]
            # If window is off edge of grid, edge(s) will be dropped and window
            # will start or end at a different line or pixel.
            # Window line start, line end.
            l1, l2 = match[9], match[10]
            # Window pixel start, pixel end.
            p1, p2 = match[11], match[12]
            
            # If window is off the edge of the grid, find the amount of
            # necessary padding for each side.
            if (l2-l1)*(p2-p1) < wl:
                
                top_pad, bottom_pad, left_pad, right_pad = 0, 0, 0, 0
                if v1 < 0:
                    top_pad = -v1
                if v2 > number_of_lines:
                    bottom_pad = v2 - number_of_lines
                if h1 < 0:
                    left_pad = -h1
                if h2 > pixels_per_line:
                    right_pad = h2 - pixels_per_line
                padding = ((top_pad,bottom_pad), (left_pad, right_pad))
            
            data = {}
            geopath =  '/geophysical_data'
            # Collect all variables except latitude and longitude, which don't
            # require any masking.
            for var in output_vars[:-2]:
                # Possible error: "reshape" can't be used on a Pandas series or dataframe
                # https://stackoverflow.com/questions/42240376/dataframe-object-has-no-attribute-reshape
                # "pandas.dataframe doesn't have a built-in reshape method, but you can use .values to access the underlying numpy array object and call reshape on it"
                # Note that you didn't have this problem when using "reshape" in other places because you converted the object to a numpy array first.
                mask = maskframe.loc[:,var].values.reshape(l2-l1,p2-p1)
                
                # Extract a window of pixels around the matching point.
                window = inputnc[geopath].variables[var][l1:l2, p1:p2]
                
                # Pad resulting window and mask if they're off edge of grid.
                if (l2-l1)*(p2-p1) < wl:
                    window_pad = np.pad(window, padding, 'constant',
                                        constant_values=np.nan)
                    mask_pad = np.pad(mask, padding, 'constant',
                                      constant_values=True)
                    window = window_pad
                    mask = mask_pad
                
                # Collect the data in a flattened window.
                data[var] = np.where(mask, np.nan, window).ravel()
                
            navpath = '/navigation_data'
            # Collect latitude and longitude here, padding if necessary.
            for var in ['latitude', 'longitude']:
                
                window = inputnc[navpath].variables[var][l1:l2, p1:p2]
                if (l2-l1)*(p2-p1) < wl:
                    window_pad = np.pad(window, padding, 'constant',
                                        constant_values=np.nan)
                    window = window_pad
                data[var] = window.ravel()
            
            # Calculate original line/pixel from NASA netCDF file.
            if l2_exists:
                orig_line = np.int32(line + l2_line)
                orig_pix = np.int32(pixel + l2_pix)
            else:
                orig_line, orig_pix = 0, 0
            
            # Compile collected variable data for this particular match.
            values = {'pass_id': pass_id, 'extraction_number': count,
                      'matched_pass': f, 'distance': distance,
                      'window_center_line_zero_based': np.int32(line),
                      'window_center_pixel_zero_based': np.int32(pixel),
                      'window_center_original_line_zero_based': orig_line,
                      'window_center_original_pixel_zero_based': orig_pix,
                      'nctime': nctime, 'time_str': time_str}
            
            # Append 1D variables to new .nc file.
            for path,varlist in output.items():
                for var in varlist:
                    # Write data variables to output file.
                    if path == in_situ_path:
                        if var == 'record_number':
                            d = np.int32(in_situ_dataindex)
                        else:
                            d = in_situ_datarow[var]
                        var = ''.join([in_situ_dataname, '_', var])
                        if isinstance(d, str):
                            outputnc[path[:-path_chars]].variables[var][count] = d
                        else:
                            if not np.isnan(d): outputnc[path[:-path_chars]].variables[var][count] = d
                            else: outputnc[path[:-path_chars]].variables[var][count] = np.nan
                    # Write pass_id and processing variables to output file.
                    else:
                        outputnc[path[:-1]].variables[var][count] = values[var]
                        
            # Append 2D variables to new .nc file.
            path = '/processing_data'
            for ovar in output_vars:
                outputnc[path].variables[ovar][count,:] = data[ovar]
            
            # Add new extraction number to in_situ_data record in dataframe.
            in_situ_data.loc[in_situ_dataindex, 'extraction_numbers'].append(count)


#==============================================================================
# If a duplicate in_situ_data record is found (for example, same year, day of year,
# latitude, and longitude, but different depth), re-use the data from the
# original in_situ_data record.
# NOTE: the original in_situ_data record could have multiple matches from a single
# netCDF file or multiple netCDF files, so the duplicate record will also have
# multiple matches.

def duplicatein_situ_datarow(index, row, extraction_nums):
    
    new_extraction_nums = []
    
    with nc4.Dataset(netcdf_filename, 'a') as outputnc:
        
        # Get starting point for new extraction records.
        records = outputnc.dimensions['records'].size
        
        # Copy data from original data record matches to duplicate data record.
        for nextrow,record in enumerate(extraction_nums, start=records):
            
            new_extraction_nums.append(nextrow)
            
            # 1D variables.
            for path,varlist in output.items():
                for var in varlist:
                    # Write data variables to output file.
                    if path == in_situ_path:
                        dv = ''.join([in_situ_dataname, '_', var])
                        # Copy old variables from original record.
                        if var == 'record_number':
                            d = np.int32(index)
                        elif var in varlist[1:6]:
                            d = outputnc[path[:-path_chars]].variables[dv][record]
                        # Get new variables from duplicate row.
                        else:
                            d = row[var]
                        outputnc[path[:-path_chars]].variables[dv][nextrow] = d
                    # Write pass_id and processing variables to output file.
                    else:
                        if var == 'extraction_number':
                            d = np.int32(nextrow)
                        else:
                            d = outputnc[path[:-1]].variables[var][record]
                        outputnc[path[:-1]].variables[var][nextrow] = d
            
            # 2D variables.
            path = '/processing_data'
            for ovar in output_vars:
                d = outputnc[path].variables[ovar][record,:]
                outputnc[path].variables[ovar][nextrow,:] = d
    
    # Return new extraction numbers of duplicate row to in_situ_data record.
    return new_extraction_nums


#==============================================================================
# Extract a window around each potential match in good_points, and mask pixels.
# Pixels are masked by a 'general mask', which includes l2 flags and pixels
# that have >1 negative or "fill value" Rrs pixels. Each variable/window also
# has a specific mask, which may or may not include more masked pixels (for
# example, if a pixel has a negative Rrs value for Rrs_412, but good values for
# all other variables, that pixel would not be included in the general mask,
# but it should be masked from Rrs_412 specifically).

def extractandmask(filename, lines, pixels, good_points, mask, ws):
    
    lgp = low_good_points
    output_other = [v for v in output_vars if not (v.startswith('Rrs')
                    or v in ['latitude', 'longitude'])]
    output_nav = ['latitude', 'longitude']
    output_Rrs = [v for v in output_vars if v.startswith('Rrs')]
    
    # Check the window around each potential matching pixel in good_points and
    # weed out pixels with flagged or null values, and pixels with >1 negative
    # Rrs values (suggesting the pixel has bad data).
    results = []

    with nc4.Dataset(filename,'r') as nc:
        
        path = '/geophysical_data'
        
        for line,pixel,dist in good_points:
            
            # Get pixel coordinates of window around the matching point.
            # (v = vertical/line, h = horizontal/pixel)
            v1 = np.int32(line - ((ws-1)/2))
            v2 = np.int32(line + ((ws+1)/2))
            h1 = np.int32(pixel - ((ws-1)/2))
            h2 = np.int32(pixel + ((ws+1)/2))
            
            # Trim edges of window if they're on the edge of the grid.
            # Window line start, line end.
            l1, l2 = max(0, v1), min(lines, v2)
            # Window pixel start, pixel end.
            p1, p2 = max(0, h1), min(pixels, h2)
            
            # Window length for this particular match: includes situations
            # where the matching pixel is at the edge of the grid, shrinking
            # the window (if so, padding will be added in addtoOutput method).
            wl = (l2-l1)*(p2-p1)
            
            # Create subarray of l2 mask for this window and check if it has
            # enough unflagged pixels; if not, ignore this window.
            l2_mask = mask[l1:l2, p1:p2]
            if wl - np.count_nonzero(l2_mask) < lgp: continue
            
            # Create a dataframe to collect specific masks for variables.
            # Note: False = unmasked
            cols = output_vars
            ind = list(range(wl))
            maskframe = pd.DataFrame(True, index=ind, columns=cols)
            
            # If there are non-Rrs variables (excluding lat/long), extract a
            # window around the matching point, and mask out l2 flags and fill
            # values.
            for var in output_other:
                window = nc[path].variables[var][l1:l2, p1:p2]
                l2_fill_mask = np.array(l2_mask)
                if np.ma.is_masked(window):
                    fill = window.get_fill_value()
                    l2_fill_mask = np.where(l2_mask | (window.data == fill),
                                             True, False)
                maskframe.loc[:,var] = l2_fill_mask.ravel()
            
            # If there are no Rrs variables, use non-Rrs variables to check
            # for sufficient number of good points in this window for at least
            # one of the variables.
            if len(output_Rrs) == 0:
                
                goodwindows = [True for name,col in maskframe.iteritems()
                               if col[~col].size >= lgp]
                if len(goodwindows) == 0: continue
            
            # Or if there is only one Rrs variable, find its mask by cutting
            # out fill values and negative values.
            elif len(output_Rrs) == 1:
                
                rrs = output_Rrs[0]
                window = nc[path].variables[rrs][l1:l2, p1:p2]
                # Mask l2 flags, fill values (if any), and negative Rrs.
                total_mask = np.where(l2_mask | (window < 0), True, False)
                if np.ma.is_masked(window):
                    fill = window.get_fill_value()
                    total_mask = np.where(l2_mask | (window.data == fill)
                                          | (window < 0), True, False)
                if (wl - np.count_nonzero(total_mask)) < lgp: continue
                maskframe.loc[:,rrs] = total_mask.ravel()
                # Update the masks for other output as well, since they are
                # calculated using Rrs values for that pixel.
                for var in output_other:
                    old = np.array(maskframe.loc[:,var]).reshape(l2-l1,p2-p1)
                    maskframe.loc[:,var] = np.where(old | total_mask,
                                                    True, False).ravel()
            
            # Or if >1 Rrs variables used, test for >1 negative Rrs values in
            # each pixel of the window (if true, pixel likely has bad data).
            elif len(output_Rrs) > 1:
                
                # Create dataframe that will be used to check if a pixel (row)
                # has more than one negative or "fill" Rrs (column) value.
                cols = list(output_Rrs)
                cols.append('Rrs_mask')
                ind = list(range(wl))
                bad_Rrs = pd.DataFrame(np.nan, index=ind, columns=cols)
                
                for count,var in enumerate(output_Rrs):
                    
                    # NOTE: Do not use the l2 mask in this section of code--
                    # apply it afterward, since it applies to all variables in
                    # the same way (not the individual masks).
                    window = nc[path].variables[var][l1:l2, p1:p2]
                    condition = window < 0
                    if np.ma.is_masked(window):
                        fill = window.get_fill_value()
                        fill_mask = np.where(window.data == fill, True, False)
                        condition = (window < 0) | fill_mask
                        rrs_window = np.array(window)
                        window = np.where(fill_mask, -1, rrs_window)
                    maskframe.loc[:,var] = np.where(condition,
                                                    True, False).ravel()
                    bad_Rrs.iloc[:,count] = window.ravel()
                
                # Add another column for results of the test (True = bad data).
                bad_Rrs.loc[:,'Rrs_mask'] = pd.Series(False, index=ind)
                
                for index,row in bad_Rrs.iterrows():
                    if np.count_nonzero([row < 0]) > 1:
                        bad_Rrs.loc[index,'Rrs_mask'] = True
                
                # Take the final column, and reshape it into original window to
                # get the logical array containing fill or negative Rrs pixels.
                mask_col = bad_Rrs.loc[:,'Rrs_mask']
                fill_rrs_mask = np.array(mask_col).reshape(l2-l1, p2-p1)
                
                # Combine masks from l2 flags and fill/negative Rrs test, and
                # test if there are still enough good pixels in the window.
                total_mask = np.where(l2_mask | fill_rrs_mask, True, False)
                if wl - np.count_nonzero(total_mask) < lgp: continue
                
                # Combine masks for flagged pixels, fill value pixels, and
                # pixels with >1 negative Rrs values, and update individual
                # masks for each selected variable, excluding lat/long.
                for var in output_vars[:-2]:
                    old = np.array(maskframe.loc[:,var]).reshape(l2-l1,p2-p1)
                    maskframe.loc[:,var] = np.where(old | total_mask,
                                                    True, False).ravel()
            
            # Count the number of masked points in each column/window (i.e. for
            # each variable) and get the lowest number (best window). Exclude
            # lat/long from this since they are only masked by l2 flags.
            # If this mixes you up as much as it mixed me up:
                # The window with the minimum number of bad points is the one
                # with the highest number of good points. When the results are
                # sorted in the checkncfile method, they'll be sorted by
                # increasing min number of bad points, i.e. decreasing max
                # number of good points, so your first results will be the
                # windows that have variables with the highest number of good
                # points.
            minbadpts = min([sum(col) for name,col in maskframe.iteritems()])
            
            # If window has passed all mask checks, apply l2 flags to lat/long
            # variables.
            for var in output_nav:
                maskframe.loc[:,var] = np.where(l2_mask, True, False).ravel()
            
            # Create a tuple with info for the potential match, and append it
            # to the list of potential matches.
            result = (line, pixel, dist, maskframe, minbadpts,
                      v1, v2, h1, h2, l1, l2, p1, p2)
            results.append(result)
            
    return results


#==============================================================================
# Check a given .nc file for matches near in_situ_data point, and extract a window
# around the point, given by the variable window_size. The 3 potential matches
# to keep are:
#     match 1: exact location (< exact_match_dist from in_situ_data point),
#              enough good points in window (>= low_good_points)
#     match 2: < max_dist metres from in_situ_data point,
#              enough good points in window (>= low_good_points)
#     match 3: < max_dist metres from in_situ_data point,
#              enough good points in window (>= high_good_points)
# The number of good points in the window are counted after applying the mask
# of selected l2 flags, removing fill values, and removing pixels that have
# >1 negative Rrs values (meaning the pixel likely has bad data).
# The final match(es) are then sent to the "addtoOutput" method to be
# appended to the output .nc file.

def checkncfile(in_situ_dataindex, in_situ_datarow, filename, sensor, platform):
    
    with nc4.Dataset(filename,'r') as nc:
        
        # Confirm date from .nc file metadata.
        lines = nc.dimensions['number_of_lines'].size
        pixels = nc.dimensions['pixels_per_line'].size
        centre = lines//2
        year = np.int32(nc['/scan_line_attributes'].variables['year'][centre])
        doy = np.int32(nc['/scan_line_attributes'].variables['day'][centre])
        if in_situ_datarow['year'] != year or in_situ_datarow['day_of_year'] != doy:
            return 0
        
        # Confirm sensor/instrument and platform from .nc file metadata.
        instr = nc.getncattr('instrument').lower()
        plat = nc.getncattr('platform').lower()
        # From input.
        s, p = sensor, platform.lower()
        if s != instr or p != plat: return 0
        
        lat = nc['/navigation_data'].variables['latitude'][:,:]
        long = nc['/navigation_data'].variables['longitude'][:,:]
        nc_l2_flags = nc['/geophysical_data'].variables['l2_flags'][:,:]
    
    in_situ_data_lat = np.float32(in_situ_datarow['latitude'])
    in_situ_data_long = np.float32(in_situ_datarow['longitude'])
    
    # Gnomonic projection with 2D array and in_situ_data point plotted.
    proj = pyproj.Proj(proj='gnom', ellps='WGS84', lat_0=in_situ_data_lat,
                       lon_0=in_situ_data_long, datum='WGS84', units='m')
    x, y = proj(long, lat)
    in_situ_data_x, in_situ_data_y = proj(in_situ_data_long, in_situ_data_lat)
    # Distance calculations between in_situ_data point and other coordinates.
    dist1 = np.sqrt((x - in_situ_data_x)**2 + (y - in_situ_data_y)**2)
    
    # Points within a certain distance (radius) of in_situ_data point.
    good_points = dist1 < max_dist
    
    if np.sum(good_points)==0: return 0
    
    # NaNs trigger error message from geoid.inv:
    # ValueError: undefined inverse geodesic (may be an antipodal point)
        # Check for masked cells and change them to False so they don't cause error in the distance calculations.
    if (np.ma.is_masked(good_points)):
        good_points[good_points.mask] = False
        print(good_points)
    
    # MORE ACCURATE DISTANCE CALCULATION (code from Brendan DeTracey)
    # Note: this will reduce the distance since it's affected by the Earth's curvature
    num_idx = np.count_nonzero(good_points)
    geoid = pyproj.Geod(ellps='WGS84')
    __, __, dist1 = geoid.inv(long[good_points], lat[good_points], np.tile(in_situ_data_long, num_idx), np.tile(in_situ_data_lat, num_idx))
    
    # Get pixel coordinates of good points, and add them to a tuple with their
    # actual distance to the in_situ_data point.
    #good_points = list(zip(np.where(good_points)[0], np.where(good_points)[1],
    #                       np.sqrt((x[good_points] - in_situ_data_x)**2
    #                               + (y[good_points] - in_situ_data_y)**2)))
    
    good_points = list(zip(np.where(good_points)[0], np.where(good_points)[1], dist1))
    
    # Create a logical array with flagged pixels = True.
    flagmask = flagbit
    mask = (nc_l2_flags & flagmask) != 0
    hgp = high_good_points
    ws = window_size
    
    # Extract window around potential matches in good points, and mask pixels
    # that are flagged, fill values, or have >1 negative Rrs values.
    results = extractandmask(filename, lines, pixels, good_points, mask, ws)
    
    # No match.
    if len(results) == 0: return 0
    
    matches = []
    # Sort results by distance and number of acceptable points in window to
    # find a match within a specified range of in_situ_data point, in metres (close
    # enough to call it "exact"). Sort by increasing bad points/decreasing good
    # points.
    closest = sorted(results, key=lambda x: (x[2], x[4]))[0]
    if closest[2] < exact_match_dist:
        matches.append(closest)
        results.remove(closest)
    # Sort remaining results by number of good pixels in window, then distance.
    if len(results) > 0:
        results = sorted(results, key=lambda x: (x[4], x[2]))
        # Get second match ( < max_dist, >= low_good_points).
        matches.append(results[0])
        if len(results) > 1:
            l1, l2 = results[1][9], results[1][10]
            p1, p2 = results[1][11], results[1][12]
            wl = (l2-l1)*(p2-p1)
            # Window length and number of good pixels in window of next result.
            if (wl - results[1][4]) >= hgp:
                # Get third match ( < max_dist, >= high_good_points).
                matches.append(results[1])
    
    # Add match(es) to output .nc file.
    addtoOutput(in_situ_dataindex, in_situ_datarow, filename, matches)
    
    return len(matches)


start_time = time.time()
#==============================================================================
# Start the program, calling methods to add the in_situ_data data to a frame
# and create an output .nc file for the results.
# Loop through the rows of the in_situ_data frame, extracting data to form paths and
# filenames in order to find specific files that might contain matches for
# the current in_situ_data row.
# If a file could have a potential match, call the checkncfile method to
# determine if the in_situ_data point has an adequate match in that file. If it
# does, the checkncfile method will send the results to the addtoOutput method.
# If a in_situ_data row is a duplicate of another row already added to the output
# .nc file, it will be handled by calling the duplicatein_situ_datarow method.

# Get an integer used to extract the path name alone in the output variable dictionary.
# Example: take "/hplc_data" from "/hplc_data/hplc_" by removing the last 6 characters.
path_chars = len(in_situ_dataname) + 2

# Process sensor, window size, and l2 flag mask.
if sensor=='modis':
    sletter, platform, sname = 'A', 'Aqua', 'MODIS'
elif sensor=='seawifs':
    sletter, platform, sname = 'S', 'Orbview-2', 'SeaWiFS'
elif sensor=='viirs':
    sletter, platform, sname = 'V', 'Suomi-NPP', 'VIIRS'
window_length = window_size**2
flagbit = sum(l2_flags.values())

# Put data from in_situ_data file into a Pandas DataFrame.
print('\nReading in_situ_data...')
in_situ_data = in_situ_data_toframe()

if os.path.isfile(prog_directory + netcdf_filename):
    # Get existing in_situ record numbers and corresponding extraction numbers to
    # avoid repeating them when looping through the in_situ data frame.
    with nc4.Dataset(prog_directory + netcdf_filename,'r') as nc:
        bpath = ''.join([in_situ_dataname, '_data'])
        bvar = ''.join([in_situ_dataname, '_record_number'])
        in_situ_data_recs = nc[bpath].variables[bvar][:]
        exts = nc['/processing_data'].variables['extraction_number'][:]
    print('Appending results to',(prog_directory + netcdf_filename),'...')
else:
    # Create new output .nc file.
    createnetcdf()

# Determine the path pattern to use for each loop, based on sensor.
path_pattern = os.path.join(input_dir, '%Y', sletter + '%Y%j*.nc')

print('Searching for .nc files...\n')
count, duplicates = 0, 0
for index,row in in_situ_data.loc[startrow:endrow,:].iterrows():
# in_situ_data.iterrows() to check the entire in_situ_dataset.
# in_situ_data.loc[startrow:endrow,:].iterrows() to check a specific set of rows
# in the dataset.
    
    print('in_situ_data index:', index)
    
    num_of_ext = 'num_of_extractions'
    extnums = 'extraction_numbers'
    
    
    
    
    # ADD AN ERROR CHECK HERE:
        # if adding new matchups to an existing netCDF based on a certain in situ record file, first check if that in situ file already exists in the netCDF
    
    
    
    
    # Check if in_situ_data record index has already been matched and added to
    # the existing .nc match file. If so, ignore this row.
    if file_exists and append:
        if index in list(in_situ_data_recs):
            oldextnums = exts[in_situ_data_recs == index]
            print('Matches found in existing netCDF file: extraction numbers',
                  oldextnums)
            in_situ_data.loc[index, 'extracted'] = True
            in_situ_data.loc[index, num_of_ext] = (list(in_situ_data_recs)).count(index)
            in_situ_data.loc[index, extnums].extend(oldextnums)
            continue
    
    # If year, day of year, or latitude/longitude are NaN, skip this row.
    if row.iloc[3:7].isnull().values.any(): continue
    # in_situ_data is sorted by year, day of year, latitude, and longitude, so check
    # the preceding row to see if this one is a duplicate.
    elif (row.iloc[3:7].equals(in_situ_data.iloc[index-1, 3:7])
          and in_situ_data.loc[index-1,'extracted'] == True):
        
        in_situ_data.loc[index, 'extracted'] = True
        in_situ_data.loc[index, num_of_ext] = in_situ_data.loc[index-1, num_of_ext]
        # In the output netcdf file, copy data from original in_situ_data record's
        # matches to duplicate in_situ_data record, and return extraction numbers
        # of the new rows.
        dnums = duplicatein_situ_datarow(index, row, in_situ_data.loc[index-1, extnums])
        in_situ_data.loc[index, extnums].extend(dnums)
        duplicates += len(dnums)
        continue
    
    # Extract year and day of year from this row to build paths to .nc files.
    in_situ_datayear = int(row.loc['year'])
    in_situ_datadoy = int(row.loc['day_of_year'])
    # Calculate month and day of month from year and day of year.
    in_situ_datamonth = int((dt.date(in_situ_datayear,1,1)
                         + dt.timedelta(in_situ_datadoy-1)).strftime('%m'))
    in_situ_dataday = int((dt.date(in_situ_datayear,1,1)
                       + dt.timedelta(in_situ_datadoy-1)).strftime('%d'))
    
    expanded_path = dt.date(in_situ_datayear, in_situ_datamonth,
                            in_situ_dataday).strftime(path_pattern)
    
    # Check directories and subdirectories for files matching given sensor,
    # year, and month.
    for filename in glob.iglob(expanded_path, recursive=True):
        
        match_num = checkncfile(index, row, filename, sensor, platform)
        
        # If match found, update first 3 columns of in_situ_data dataframe.
        if match_num > 0:
            in_situ_data.loc[index, 'extracted'] = True
            in_situ_data.loc[index, num_of_ext] += match_num
            print(filename,': match found')
            count += match_num
            

#==============================================================================
end_time = time.time()

# Get time of search completion from system settings and convert it to a string
# to be used in the results below.
tc = str(dt.datetime.today())[:-7].split(':')
time_completed = ''.join([tc[0], 'h', tc[1], 'm', tc[2], 's'])

# Write the relevant rows of the resulting in_situ_data dataframe to an Excel file.
excel_name = ''.join([in_situ_dataname,'_rows',str(startrow),'-',
					  str(endrow),'_',netcdf_filename[:-2],'xlsx'])
writer = pd.ExcelWriter(excel_name)
in_situ_data.iloc[startrow:endrow+1,:].to_excel(writer,'in_situ_data')
writer.save()

# Print results of the search.
print('Elapsed time was %g seconds' % (end_time - start_time))
print('Time completed:', time_completed)
print('Summary:')
print('Start row:', startrow, '\nEnd row:', endrow)
print(count + duplicates, 'records found:', count, 'unique matches and',
      duplicates, 'duplicates')
print('Input filename:', in_situ_data_filename)
print('Output filename:', netcdf_filename)
print('Sensor:', sname)
print('Window length:', window_length)
print('Low number of good points in window:', low_good_points)
print('High number of good points in window:', high_good_points)
print('L2 flags:', list(l2_flags.keys()))
print('Output variables:', output_vars)
print(''.join(['See "', excel_name, '" for resulting in_situ_data dataframe.']))
