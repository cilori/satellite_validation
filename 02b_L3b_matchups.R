# Stephanie.Clay@dfo-mpo.gc.ca
# 17 Sep 2020

# Match L3b data to in situ data.

library(geosphere)  # to calculate distance between matches
library(pbapply)
library(parallel)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(grid)
library(data.table)
library(oceancolouR)

mult_num_cl <- 10 # or whatever works best, leaving at least one processing core free

# For MODIS/SeaWiFS/VIIRS-SNPP: CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR
# For OLCI-A: CHL1, CHL2, CHL-OC5
# For OLCI-B: CHL1
variable <- "CHL_POLY4"

# OLCI-A, OLCI-B, MODIS, VIIRS-SNPP, SeaWiFS
sensor <- "VIIRS-SNPP"

region <- "NWA"

sat_path <- paste0("/mnt/data3/claysa/", sensor, "/", variable, "/", region, "/")
in_situ_path <- "01_in_situ_data/"

max_dist <- 10000 # metres

max_depth <- 10 # in situ depth, metres

# maximum number of matchup pixels for each in situ record (within max_dist from in situ data)
max_matchups <- 100

olci_flag_subset <- c("INVALID", "LAND", "TURBID", "ICE", "CLOUD1", "CLOUD2")


#*******************************************************************************
# IN SITU DATA ####

# Note: Dataframe must contain the columns YEAR, DOY, LATDEC, LONGDEC, DEPTH
# (DOY = day of year, LATDEC = decimal latitude (<0 for south), LONGDEC = decimal longitude (<0 for west))

# ... from a file:
in_situ_file <- "HPLC_1997-2020_20201110.csv"
in_situ_data <- read.csv(paste0(in_situ_path, in_situ_file), header=TRUE)
in_situ_data$cruise_id <- as.character(in_situ_data$cruise_id)

# # ... or manually entered here:
# in_situ_data <- data.frame(YEAR=c(rep(2016,length(131:366)), rep(2017, 201)),
#                            DOY=c(131:366,1:201),
#                            LATDEC=56.823145,
#                            LONGDEC=-52.219467,
#                            DEPTH=0)


#*******************************************************************************
# OUTPUT FORMAT ####

output_path <- "02_matchups/"
output_name_base <- paste0(output_path, "matchups_L3b_daily_", sensor, "_", variable, "_maxdist_", max_dist, "_maxdepth_", max_depth, "_maxmatchups_", max_matchups)

# create output dataframe column names
output_colnames <- c(colnames(in_situ_data), paste0("sat_", variable), "pancan_lon", "pancan_lat", "dist_to_matchup")

# note which columns should be in what format, so that when output is created later,
# the columns have the right variable types
# NOTE: this is required to convert the output from a list to a dataframe
integer_cols <- which(colnames(in_situ_data) %in% c("YEAR", "DOY", "TIME"))
character_cols <- which(colnames(in_situ_data) %in% c("cruise_id", "ID"))
numeric_cols <- (1:ncol(in_situ_data))[!((1:ncol(in_situ_data)) %in% c(integer_cols, character_cols))]


#*******************************************************************************
# OLCI FLAGS ####

olci_flags <- NULL

if (startsWith(sensor, "OLCI")) {
    
    olci_flags <- data.frame(bit=0:15,
                             value=2^(0:15),
                             string=c("NO_MEASUREMENT", "INVALID", "OLCI_A", "LAND",
                                      "CLOUD1", "CLOUD2", "DEPTH1", "DEPTH2", "TURBID", 
                                      "ICE", "TROPHIC1", "TROPHIC2", "VIIRS_N",
                                      "SEAWIFS_or_VIIRS_J1", "MODIS", "MERIS_or_OLCI_B"))
    
    # (CLOUD2=0) + (CLOUD1=0): CF < 5%
    # (CLOUD2=0) + (CLOUD1=1): 5% <= CF < 25%
    # (CLOUD2=1) + (CLOUD1=0): 25% <= CF < 50%
    # (CLOUD2=1) + (CLOUD1=1): CF >= 50%
    # 
    # (DEPTH2=0) + (DEPTH1=0): depth < 30m
    # (DEPTH2=0) + (DEPTH1=1): 30m <= depth < 200m
    # (DEPTH2=1) + (DEPTH1=0): 200m <= depth < 1000m
    # (DEPTH2=1) + (DEPTH1=1): depth >= 1000m
    # 
    # (TROPHIC2=0) + (TROPHIC1=1): Oligotrophic water
    # (TROPHIC2=1) + (TROPHIC1=0): Mesotrophic water
    # (TROPHIC2=1) + (TROPHIC1=1): Eutrophic water
    
    # BIOLOGICAL PRODUCTIVITY SCALE:
    #       LOW --------------------------------> HIGH
    #       oligotrophic --> mesotrophic --> eutrophic
    
    # OLCI MASK - INVALID / LAND / TURBID / ICE / CLOUD > 50%
    olci_mask_value <- sum(olci_flags[olci_flags$string %in% olci_flag_subset, "value"])
    olci_mask <- as.logical(intToBits(olci_mask_value))[1:nrow(olci_flags)]
    
    output_colnames <- c(output_colnames, paste0("MASK_", olci_flags$string))
    
    # ADD THESE TO INTEGER DATATYPE COLUMN LIST
    integer_cols <- c(integer_cols, (ncol(in_situ_data)+1):(ncol(in_situ_data)+length(olci_flag_subset)))
    
}


#*******************************************************************************
# MATCHUP FUNCTION ####

find_matches <- function(i, in_situ_data, sat_path, olci_flags=NULL, pancan_lon, pancan_lat, max_dist, output_colnames, input_varname, sensor, max_matchups) {
    
    sensor_is_olci <- startsWith(sensor, "OLCI")
    
    results <- data.frame(matrix(nrow=0, ncol=length(output_colnames)), stringsAsFactors = FALSE)
    colnames(results) <- output_colnames
    
    isd_row <- in_situ_data[i,]
    year <- as.numeric(isd_row$YEAR)
    doy <- as.numeric(isd_row$DOY)
    isd_lat <- as.numeric(isd_row$LATDEC)
    isd_lon <- as.numeric(isd_row$LONGDEC)
    
    if (sensor_is_olci) {
        date <- as.Date(as.numeric(doy), origin = paste0(year-1, "-12-31"))
        month <- format(date, "%m")
        day <- format(date, "%d")
        sat_file <- list.files(paste0(sat_path, year), pattern=paste0("L3b_", year, month, day, "_"))
    } else {
        sat_file <- list.files(paste0(sat_path, year), pattern=paste0(year, doy))
    }
    
    if (length(sat_file) == 0) {return(results)}
    
    nc <- nc_open(paste0(sat_path, year, "/", sat_file))
    sat_var <- ncvar_get(nc, input_varname)
    if (sensor_is_olci) {nc_flags <- ncvar_get(nc, "flags")}
    nc_close(nc)
    
    if (sensor_is_olci) {
        
        finite_flags <- which(is.finite(nc_flags))
        nc_flags_subset <- nc_flags[finite_flags]
        # nc_mask <- rep(NaN, length(nc_flags))
        # nc_mask[finite_flags] <- sapply(1:length(nc_flags_subset), function(x) sum(as.logical(intToBits(nc_flags_subset[x])) & olci_mask) > 0)
        
        nc_flags_bits <- matrix(nrow=length(nc_flags), ncol=nrow(olci_flags))
        nc_flags_bits[finite_flags,] <- t(sapply(1:length(nc_flags_subset), function(x) as.integer(intToBits(nc_flags_subset[x]))[1:nrow(olci_flags)]))
        colnames(nc_flags_bits) <- paste0("MASK_", olci_flags$string)
        
        matchup_df <- cbind(data.frame(sat_var=sat_var,
                                       # olci_mask=nc_mask,
                                       pancan_lon=pancan_lon,
                                       pancan_lat=pancan_lat,
                                       stringsAsFactors = FALSE),
                            nc_flags_bits) %>%
               dplyr::filter(., !is.na(sat_var) & abs(pancan_lat - isd_lat) < 1 & abs(pancan_lon - isd_lon) < 1)# & !olci_mask)
        
    } else {
        
        matchup_df <- data.frame(sat_var=sat_var,
                                 pancan_lon=pancan_lon,
                                 pancan_lat=pancan_lat,
                                 stringsAsFactors = FALSE) %>%
            dplyr::filter(., !is.na(sat_var) & abs(pancan_lat - isd_lat) < 1 & abs(pancan_lon - isd_lon) < 1)
        
    }
    
    if (nrow(matchup_df)==0) {return(results)}
    
    # get accurate distance between pts, make sure it's within the range you want
    matchup_df <- matchup_df %>%
        dplyr::mutate(., dist_to_matchup = as.numeric(sapply(1:nrow(matchup_df),function(i) distGeo(as.numeric(matchup_df[i,c("pancan_lon", "pancan_lat")]), c(isd_lon, isd_lat))))) %>% # use decimal degrees
        dplyr::filter(., dist_to_matchup <= max_dist)
    
    if (nrow(matchup_df)==0) {return(results)}
    
    # reduce to max number of matchups
    last_row <- ifelse(max_matchups > nrow(matchup_df), nrow(matchup_df), max_matchups)
    matchup_df <- matchup_df %>%
        dplyr::arrange(., dist_to_matchup)
    matchup_df <- matchup_df[1:last_row,]
    
    if (sensor_is_olci) {
        # move flag columns to the end
        matchup_df <- dplyr::bind_cols(matchup_df %>% dplyr::select(!starts_with("MASK_")),
                                       matchup_df %>% dplyr::select(starts_with("MASK_")))
    }
    
    isd_output <- data.frame(matrix(isd_row, nrow=nrow(matchup_df), ncol=length(isd_row), byrow=TRUE), stringsAsFactors = FALSE)
    results <- cbind(isd_output, matchup_df)
    colnames(results) <- output_colnames
    
    return(results)
    
}


#*******************************************************************************
# MAIN CODE ####

# get 4km-resolution lats/lons in vector form
if (region=="PANCAN") {
    data("pancan_lats_4km")
    data("pancan_lons_4km")
    lats_4km <- pancan_lats_4km
    lons_4km <- pancan_lons_4km
} else if (region=="NWA") {
    data("nwa_lats_4km")
    data("nwa_lons_4km")
    lats_4km <- nwa_lats_4km
    lons_4km <- nwa_lons_4km
} else if (region=="NEP") {
    data("nep_lats_4km")
    data("nep_lons_4km")
    lats_4km <- nep_lats_4km
    lons_4km <- nep_lons_4km
}

# restrict records by in situ depth
in_situ_data <- in_situ_data %>%
    dplyr::filter(DEPTH <= 10)

if (variable=="CHL_GSM_GS") {
    input_varname <- "chl_GSM_GS"
} else if (startsWith(variable, "CHL")) {
    input_varname <- "chlor_a"
} else if (variable=="PAR") {
    input_varname <- "par"
}

if (nrow(in_situ_data)==1) {
    num_cl <- 1
} else {
    num_cl <- mult_num_cl
}

# Create clusters for parallel processing.
num_cl <- min(num_cl, detectCores()-1)
cl <- makeCluster(num_cl)
# Load necessary variables or custom functions into cluster.
clusterExport(cl,
              varlist = c("find_matches", "in_situ_data", "sat_path", "olci_flags",
                          "lons_4km", "lats_4km", "max_dist", "output_colnames",
                          "input_varname", "sensor", "max_matchups"),
              envir = environment())
# Load necessary libraries into cluster.
clusterEvalQ(cl, c(library(ncdf4), library(dplyr), library(geosphere)))
# Like lapply, but with the clusters variable cl as the first argument.
results <- pblapply(1:nrow(in_situ_data), find_matches,
                    in_situ_data=in_situ_data, sat_path=sat_path, olci_flags=olci_flags,
                    pancan_lon=lons_4km, pancan_lat=lats_4km,
                    max_dist=max_dist, output_colnames=output_colnames,
                    input_varname=input_varname, sensor=sensor, max_matchups=max_matchups,
                    cl=cl)
# Stop parallel processing and return processing power to other operations.
stopCluster(cl)

# # for testing
# results <- list()
# for (i in 1:nrow(in_situ_data)) {
#     results[[i]] <- find_matches(i, in_situ_data=in_situ_data, sat_path=sat_path, olci_flags=olci_flags,
#                                  pancan_lon=lons_4km, pancan_lat=lats_4km,
#                                  max_dist=max_dist, output_colnames=output_colnames,
#                                  input_varname=input_varname, sensor=sensor, max_matchups=max_matchups)
# }


#***************************************************************************
# FORMAT AND WRITE RESULTS TO CSV ####

# get the number of matchups
valid_match_ind <- which(sapply(1:length(results), function(i) length(results[[i]][[1]]) > 0))

if (sum(valid_match_ind)==0) {
    
    cat(paste("No matchups for", region, sensor, variable))
    
} else {
    
    # reshape into table after removing indices where there were no matchups
    matchup_df <- results[valid_match_ind]
    matchup_df <- rbindlist(matchup_df)
    matchup_df <- as.data.frame(matchup_df)
    
    # convert columns to proper variable types
    # NOTE: this is required to convert the output from a list to a dataframe
    matchup_df[,numeric_cols] <- lapply(matchup_df[,numeric_cols], as.numeric)
    matchup_df[,integer_cols] <- lapply(matchup_df[,integer_cols], as.integer)
    matchup_df[,character_cols] <- lapply(matchup_df[,character_cols], as.character)
    
    # str(matchup_df)
    
    # # see how many matchups each day gets
    # matchup_df %>%
    #     dplyr::filter(., dist_to_matchup <= max_dist) %>%
    #     dplyr::group_by(., YEAR, DOY, LATDEC, LONGDEC, DEPTH) %>%
    #     dplyr::summarize(., sat_PAR = mean(sat_PAR, na.rm=TRUE),
    #                         num_pix = n()) %>%
    #     dplyr::ungroup() %>%
    #     dplyr::arrange(., YEAR, DOY, LATDEC, LONGDEC, DEPTH)
    
    # arrange rows by YEAR, DOY, LATDEC, LONGDEC, DEPTH
    matchup_df <- matchup_df %>%
        dplyr::arrange(., YEAR, DOY, LATDEC, LONGDEC, DEPTH)
    
    write.csv(matchup_df, paste0(output_name_base, ".csv"), row.names=FALSE)
    
}
