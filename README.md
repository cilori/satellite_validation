# Satellite Validation  

This repository contains scripts to extract in situ data from the BIOCHEM database, and match it to satellite files using an R script to extract satellite L3b matchups (the easiest and fastest option), or a Python script to extract L2 matchups.  



Files in this repository:  

| Filename                               | Description                                                       |
| -------------------------------------- | ----------------------------------------------------------------- |
| 01_BC_Data_Retrievals_08Mar2017.txt    | BIOCHEM IDs for variables to extract                              |
| 01_database_functions.R                | SQL/Database functions used for extracting data                   |
| 01_in_situ_data_extraction.R           | R script to that reads SQL script and formats data                |
| 01_in_situ_data_extraction.sql         | SQL script to query database                                      |
| 02a_L2_matchups.py                     | Get in situ / satellite matchups from L2 files                    |
| 02b_L3b_matchups.R                     | Get in situ / satellite matchups from L3b files                   |
| 03a_matchup_analysis.Rmd               | Create Rmarkdown report summarizing results from multiple sensors |
| 03b_matchup_analysis_single_sensor.Rmd | Create Rmarkdown report summarize results for one sensor          |




#### TO DO list:  

- Adjust the L2 matchup script to write matchups to csv in the same format as the L3b script, for consistency and easier access.  


--------------------------------------------------------------------------------

## INSTRUCTIONS

### STEP 1: Get in situ data

**01_in_situ_data_extraction.R**  

### STEP 2: Get satellite matchups

*Note: Access to PanCanadian dataset (for L3b matchups) or Panix (for L2 matchups) required for this step*  

**02a_L2_matchups.py**  

    *NOTE: This program was scripted in Python3*.  

    [to be updated soon]  
    
    To make changes, open the program in Spyder or Notepad++.  
    Check to make sure directories and maximum searching distances are correct in the ‘Global Variables’ section, near the top of the script.  
    Move program to Panix where the netCDF files are located.  
    Run the program in the Linux terminal (or login to Linux using Cygwin on Windows) by typing:  
    `python3 directoryname/programname.py`
    Follow instructions and enter required variables.  
    
    **NOTES:**  
    
    Input HPLC database and find matches in NASA's satellite files (.nc format).  HPLC database must be in a specific format, containing columns cruise_id, ID, DEPTH, LATDEC, LONGDEC, YEAR, DOY, HPLCHLA, HPLCPHAE, ZEA, CHLB, HEX19, BUT19, ALLOX, FUCOX, PERID, CHLC3  
    
    NetCDF file notes:  
    Shape of numpy array of 2D variables given by:  
    `(line,pixel) = (y,x) = (row,column)`  
    In Seadas, the grid might be flipped depending on whether the satellite pass was ascending or descending (info found in global attributes). If viewing the original .nc file of a match in Seadas, might need to check the coordinate given by:  
    `(number_of_lines - 1 - line, pixels_per_line - 1 - pixel)`  
    The index of the grid begins at 0, which is why you must subtract 1 from the total size of each dimension.  
    
    **DEFAULTS FOR INPUT VARIABLES:**  

    default_l2_flags = {'ATMFAIL': 1, 'BOWTIEDEL': 268435456, 'LAND': 2, 'HIGLINT': 8, 'CLDICE': 512, 'HILT': 16, 'HISOLZEN': 4096, 'HISATZEN': 32}  
    
    Input HPLC database file:  
    default_hplc_filename = 'HPLC_mbeck_2016_20170130.txt'  
    
    Output .nc file:  
    default_netcdf_filename = 'HPLC_NetCDF_matches.nc'  
    
    Window of pixels extracted around matching point:  
    default_window_size = 3  
    
    For second match: lowest number of good (non-NaN) pixels in window:  
    default_low_good_points = 3  
    
    For third match: second lowest number of good (non-Nan) pixels in window:  
    default_high_good_points = 9  


**02b_L3b_matchups.py**  

