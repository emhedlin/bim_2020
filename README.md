Work for the 2020 BIM report
================


### File structure

```
└── data
│     ├── 2012-2019_report.csv         <- data used in 2019 report
│     ├── 2012-2020_report.csv         <- 2020 data added to the 2019 report data
│     ├── 2019_BIM_db.csv              <- alastair's BIM database, central location for data
│     ├── dnon.csv                     <- distance to nearest neighobour, long format, generated from dnon.R
│     ├── yearly_site_occ.csv          <- for mapping. Yearly summary for each site wrt pefa/rlha occupancy
│     ├── pefa_territory.csv           <- list of pefa sites, coordinates, and territory centroid coordinates
│     ├── rlha_territory.csv           <- list of rlha sites, coordinates, and territory centroid coordinates
│     ├── pefa_territory_1km.csv       <- pefa territories using 1km tree cutoff
│     ├── rlha_territory_1km.csv       <- pefa territories using 1km tree cutoff
│     ├── unique_rlha_territory.csv    <- unique territories, made for mapping purposes only
│     ├── unique_pefa_territory.csv    <- unique territories, made for mapping purposes only
│     └── 2020_survey.csv              <- data from 2020
└── scripts -> .gitignoring for now
      ├── 01_data_cleaning.R           <- combining data from last year's sheets with this years surveys
      │                                     * also verifying counts across datasheets and reports
      ├── 02_dnon.R                    <- distance to nearest occupied neighbour
      ├── 03_dist_disturb.R            <- distance to disturbance
      ├── 04_territory.R               <- clustering, and solve territory conflicts
      ├── 05_occupancy.R               <- occupancy analysis - unmarked and Jags
      └── 06_breeding.R                <- zero-truncated counts in INLA 
```

<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/pefa_terr.png">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/rlha_terr.png">
</p>
Updated territories for peregrines (left), and rough-legs (right), and occupied sites (points) for the year 2020. Territories were delineated using cluster analysis and euclidean distance to group nest with a tree cut-off of 3500 meters. Territories occupied by two pairs of the same species within the same year were seperated into distinct territories. Maximum ndvi values from 2015 - 2020 are displayed as varying shades of green to show veg productivity in the study area.

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
* completed script to programmatically search for territories that are occupied by multiple pairs in a given year, and to and assign new territory IDs to the extra territories.
* exported csv's for territory locations to produce figures in qgis (see above maps)


### October 9
* started cleaning data for occupancy analysis, gathering covariates
* see below for diagram of additional issues with using territories.


<p align="center">
  <img width="950" src="https://github.com/emhedlin/bim_2020/blob/master/documents/territory_issue.jpg">
</p>

* based on work done in rankin, it's clear that adults sometimes utilize alternate nesting sites around a central location. This results in spatial variation in nest site locations from year to year, and in many cases, such groups of nest sites can be considered as a single territory. Establishing hard boundaries that define these territories comes with assumptions however, and it may be best to explore nest site spatial dynamics among marked birds before accepting these assumptions. 
* Maybe there's a better way than using territories. Maybe we define a grid, and use state space count models to investigate variation in nest site densities in each cell per year.

### October 10
* noticed that the report says territories were calculated using 1km tree cutoffs. Had to recalculate territories.

### October 11
