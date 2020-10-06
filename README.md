Work for the 2020 BIM report
================


### File structure

```
└── data
│     ├── 2012-2019_report.csv  <- data used in 2019 report
│     ├── 2012-2020_report.csv  <- 2020 data added to the 2019 report data
│     ├── 2019_BIM_db.csv       <- alastair's BIM database, central location for data
│     └── 2020_survey.csv       <- data from 2020
└── scripts
│     ├── 01_data_cleaning.R    <- combining data from last year's sheets with this years surveys
│     │                              * also verifying counts across datasheets and reports
│     ├── 02_dnon.R             <- distance to nearest occupied neighbour
│     ├── 03_dist_disturb.R     <- distance to disturbance

```

## Status
### October 4

I've improved scripts from last year that calculate distance to nearest occupied neighbour, and distance to the road. Before I delineate territories, I need to make sure the accounting makes sense across the data sheets/reports. Turns out it doesn't. Spent yesterday and today trying to figure out where the errors are coming from. 

### October 5
* Summarized the database into a table that indicates what species is occupying what territory.
* Summarized the datasheet used in last year's report in the same way
* joined tables and filtered the values to show differences

|	Site    | Year      | Datasheet | Database      |
|-----------|-----------|-----------|--------|
|	101		|	13		|unoccupied |	PEFA		|
|	141		|	17		|RLHA     (updated in some data sheets)  |   PEFA   |
|	374		|	12		|RLHA       |	PEFA        |
|	71		|	15		|RLHA	    |	PEFA        |

More issues with the following sites / years
##### 2012
"7"   "9"   "15"  "49"  "55"  "63"  "81"  "91"  "374"

##### 2013
"101"

##### 2015
"8"   "17"  "18"  "34"  "80"  "122" "141" "357" "360" "391" "442" "444"
"445" "447" "448" "450" "464" "465" "470" "472" "474" "476" "478"

The numbers in the report's summary table are taken from the database, but the analysis is done using the datasheet. This year, I'm hoping to use the data sheet for both. I'll check the comments in the database to see if I can find clarification on these issues, but if I can't, I'll use the database as the truth. Changes to the datasheet will be completed on the ```data/2012-2019_report.csv``` sheet. See last year's datasheet for historical values.

Discrepancies between datasheets and reports need to be solved before I can work on territory delineation.
