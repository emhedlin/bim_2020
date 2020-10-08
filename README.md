Work for the 2020 BIM report
================


### File structure

```
└── data
│     ├── 2012-2019_report.csv  <- data used in 2019 report
│     ├── 2012-2020_report.csv  <- 2020 data added to the 2019 report data
│     ├── 2019_BIM_db.csv       <- alastair's BIM database, central location for data
│     ├── dnon.csv              <- distance to nearest neighobour, long format, generated from dnon.R
│     ├── yearly_site_occ.csv   <- for mapping. Yearly summary for each site wrt pefa/rlha occupancy
│     ├── pefa_territory.csv    <- list of pefa sites, coordinates, and territory centroid coordinates
│     ├── rlha_territory.csv    <- list of rlha sites, coordinates, and territory centroid coordinates
│     └── 2020_survey.csv       <- data from 2020
└── scripts
      ├── 01_data_cleaning.R    <- combining data from last year's sheets with this years surveys
      │                              * also verifying counts across datasheets and reports
      ├── 02_dnon.R             <- distance to nearest occupied neighbour
      ├── 03_dist_disturb.R     <- distance to disturbance
      ├── 04_territory.R        <- clustering, and solve territory conflicts
      ├── 05_occupancy.R        <- occupancy analysis
      └── 06_breeding.R         <- occupancy analysis
```

## Status
### October 4

I've improved scripts from last year that calculate distance to nearest occupied neighbour, and distance to the road. Before I delineate territories, I need to make sure the accounting makes sense across the data sheets/reports. Turns out it doesn't. Spent yesterday and today trying to figure out where the errors are coming from. 

### October 5
* Summarized the database into a table that indicates what species is occupying what territory.
* Summarized the datasheet used in last year's report in the same way
* joined tables and filtered the values to show differences

### October 6

The numbers in the report's summary table are taken from the database, but the analysis is done using the datasheet. This year, I'm hoping to use the data sheet for both. I'll check the comments in the database to see if I can find clarification on these issues, but if I can't, I'll use the database as the truth. Changes to the datasheet will be completed on the ```data/2012-2019_report.csv``` sheet. See last year's datasheet for historical values.

Discrepancies between datasheets and reports need to be solved before I can work on territory delineation.

Issues with the following sites / years - changes were made to the occupancy datasheet based on values in the database.
##### 2012
* 55 - changed to pefa on v1
* 61 - changed from RLHA to not occupied based on database
* 62 - changed to rlha on v1
* 69 - changed to unoccupied based on database
* 81 - added PEFA on v2 (seen perching on rock nearby)
* 360 - changed from unoccupied to RLHA
* 374 - changed to PEFA on v2


##### 2013
* 46 - changed from PEFA to unoccupied based on database
* 50 - changed from PEFA to unoccupied based on database
* 101 - changed from unoccupied to PEFA based on database

##### 2014
* 61 - changed from unoccupied to RLHA on v2
* 80 - changed from unoccupied to RLHA on v1
* 392 - changed from unoccupied to RLHA on v1
* 393 - changed from unoccupied to RLH Aon v1


##### 2015
* changed 71 from unoccupied to PEFA

Number of of sites with pefa and rlha line up with the report now. I can't get the number of sites "not checked" and the number of sites with "no detections" to line up. I suspect this is because I'm not account for sites that were occupied by species other than pefa and rlha. I'll investigate that tomorrow.

### October 7 

* Ignoring differences in unoccupied territories (there are lots) between data sheets and the report. The report incorporates PDA = 3, whereas the datasheets don't.
* calculated distance to nearest occupied neighbour data
* created csv for mapping yearly occupancy in pefa/rlha
* started script for territory delineation

### October 8
* completed script to programmatically search for territories that are occupied by multiple pairs in a given year, and assign new territory IDs
