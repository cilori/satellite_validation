# Pieces of code and scripts taken from Benoit Casault's sample code to query Bedford Basin.
# Biochem_Example.zip, 2016

# NOTE 1: cruise_id is "NAME" in the BIOCHEM "MISSIONS" table
# NOTE 2: As of 2019-10-18, the QC in all query results was 0 (might want to
# specify the quality control level in the future)


library(dplyr)
library(tidyr)
library(openxlsx)
library(lubridate)
library(oceancolouR)

# source custom functions
source("~/.Rprofile")
source("01_database_functions.R")

# bounds on year, latitude, longitude (inclusive)
year_bounds <- c(1997,2020)
lat_bounds <- c(-Inf, Inf)
lon_bounds <- c(-Inf, Inf)
# NOTE: restrict these from within the sql code first, otherwise the query results will be huge and it might freeze

output_file <- "HPLC_1997-2020"


#-------------------------------------------------------------------------------

# get the data
sql_file <- "01_in_situ_data_extraction.sql"
sql <- Read_SQL(sql_file)

data <- Run_Database_Query(sql, my.env$host, my.env$port, my.env$sid, my.env$username, my.env$password)

# rename variables
sql_str <- data[[1]]
df_raw <- data[[2]]
rm(data)


#-------------------------------------------------------------------------------

# format output
df <- df_raw %>%
    dplyr::filter(startsWith(PARAMETER, "HPLC") & #PARAMETER %in% c("HPLC_CHLA", "HPLC_FUCOX") &
                  between(YEAR, year_bounds[1], year_bounds[2]) &
                  between(LATITUDE, lat_bounds[1], lat_bounds[2]) &
                  between(LONGITUDE, lon_bounds[1], lon_bounds[2])) %>%
    dplyr::mutate(DOY=yday(ymd(paste(YEAR,MONTH,DAY)))) %>%
    dplyr::rename(cruise_id=CRUISE_ID,
                  ID=SAMPLE_ID,
                  DEPTH=START_DEPTH,
                  LATDEC=LATITUDE,
                  LONGDEC=LONGITUDE) %>%
    dplyr::select(-METHOD) %>%
    tidyr::pivot_wider(.,names_from=PARAMETER, values_from=VALUE) %>%
    dplyr::arrange(YEAR,DOY,DEPTH) %>%
    dplyr::select(cruise_id,ID,DEPTH,LATDEC,LONGDEC,YEAR,DOY,TIME,starts_with("HPLC"))


#-------------------------------------------------------------------------------

# create output filename
output_file <- paste0("01_in_situ_data/", output_file, "_", format(Sys.Date(), "%Y%m%d"))

# save to rdata file, along with the sql string used in the query
save(file=paste0(output_file,".rda"), list=c("df","sql_str"))

# print to csv (comma separated values) file
write.csv(df, file=paste0(output_file,".csv"), row.names=F, na="")

