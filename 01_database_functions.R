# Functions written by Benoit Casault to read sql data from a database.
# These will only work on a system that has Oracle Client installed, and requires
# an "~/.Rprofile" file with your database account details.

#*******************************************************************************
Read_SQL <- function(sql_file) {
    
    # Function to read a sql query from a text file.
    #
    # Input: sql_file : sql filename - string
    #
    # Output: sql : sql text string - string
    #
    # Last update: 20141015
    # Benoit.Casault@dfo-mpo.gc.ca
    
    # read sql query
    sql <- scan(file=sql_file, what=ls(), sep="\n")
    sql <- paste(sql, collapse=" ")
    
    return(sql)
}

#*******************************************************************************
# Small modification by Stephanie Clay (Nov 2020) to read SQL separately, and input string to this function.
Run_Database_Query <- function(sql, host, port, sid, username, password) {
    
    # Function to extract data from a database (e.g. Biochem).  The function reads a sql query, opens a
    # connection to a database, extracts the data and returns the extracted data and the associated
    # sql code.
    #
    # Input: sql : sql as read by Read_SQL function - string
    #        data_source : data source name for database - string
    #        username : username - string
    #        password: user password - string
    #
    # Output: data : list of two elements
    #         data[[1]] : sql code for the query - string
    #         data[[2]] : data extracted from the database - dataframe
    #
    # Last update: 20141015
    # Benoit.Casault@dfo-mpo.gc.ca
    
    # load odbc library
    library(ROracle)
    
    # declare empty list to store outputs
    data <- list()
    
    # read sql file
    data[[1]] <- sql
    
    # create an Oracle Database instance and create connection
    drv <- dbDriver("Oracle")
    
    # connect string specifications
    connect.string <- paste(
        "(DESCRIPTION=",
        "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))",
        "(CONNECT_DATA=(SID=", sid, ")))", sep = "")
    
    # use username/password authentication
    conn <- ROracle::dbConnect(drv, username=username, password=password, dbname = connect.string)
    
    # run SQL statement by creating first a resultSet object
    rs <- dbSendQuery(conn, data[[1]])
    
    # fetch records from the resultSet into a dataframe
    data[[2]] <- fetch(rs)
    
    # close database connection
    conn <- dbDisconnect(conn)
    
    return(data)
}

#*******************************************************************************
Write_Database_Table <- function(data, table_name, host, port, sid, username, password) {
    
    # Function to ...
    #
    # Input: data : data to write to table - dataframe
    #        table_name : name of table to write - string
    #        ....
    #        data_source : data source name for database - string
    #        username : username - string
    #        password: user password - string
    #
    # Last update: 20141015
    # Benoit.Casault@dfo-mpo.gc.ca
    
    # load odbc library
    library(ROracle)
    
    # create an Oracle Database instance and create connection
    drv <- dbDriver("Oracle")
    
    # connect string specifications
    connect.string <- paste(
        "(DESCRIPTION=",
        "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))",
        "(CONNECT_DATA=(SID=", sid, ")))", sep = "")
    
    # use username/password authentication
    conn <- ROracle::dbConnect(drv, username=username, password=password, dbname = connect.string)
    
    ## check if table already exixts
    tf <- dbExistsTable(conn, name=table_name, schema = NULL)
    if(tf){
        status <- dbRemoveTable(conn, name=table_name, purge=TRUE, schema=NULL)
    }
    
    ## write data to table
    status <- dbWriteTable(conn, name=table_name, value=data, row.names=FALSE, overwrite=TRUE,
                           append=FALSE, ora.number=TRUE, schema=NULL)
    
    # close database connection
    conn <- dbDisconnect(conn)
    
    return()
}

#*******************************************************************************
Remove_Database_Table <- function(table_name, host, port, sid, username, password) {
    
    # Function to ...
    #
    # Input: table_name : name of table to remove - string
    #        ....
    #        data_source : data source name for database - string
    #        username : username - string
    #        password: user password - string
    #
    # Last update: 20141015
    # Benoit.Casault@dfo-mpo.gc.ca
    
    # load odbc library
    library(ROracle)
    
    # create an Oracle Database instance and create connection
    drv <- dbDriver("Oracle")
    
    # connect string specifications
    connect.string <- paste(
        "(DESCRIPTION=",
        "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))",
        "(CONNECT_DATA=(SID=", sid, ")))", sep = "")
    
    # use username/password authentication
    conn <- ROracle::dbConnect(drv, username=username, password=password, dbname = connect.string)
    
    ## check if table exists
    tf <- dbExistsTable(conn, name=table_name, schema = NULL)
    if(tf){
        status <- dbRemoveTable(conn, name=table_name, purge=TRUE, schema=NULL)
    }
    
    # close database connection
    conn <- dbDisconnect(conn)
    
    return()
}
