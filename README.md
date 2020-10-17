Work for the 2020 BIM report
================

This repo contains all of the analytical components for 2020 BIM report. High level summary of the workflow includes:
* data cleaning
* calculating the distance between all sites and the nearest occupied nest site
* distance between all sites and nearest anthropogenic disturbance
* hierarchical clustering of individual nest sites into territory clumps
* dynamic occupancy modeling - unmarked and Jags scripts included
* zero-truncated poisson spatial models for breeding productivity

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
## Territory Clustering
Updated territories for peregrines (left), and rough-legs (right), and occupied sites (points) for the year 2020 . Territories were delineated using cluster analysis and euclidean distance to group nests with a tree cut-off of 3500 meters. Territories occupied by two pairs of the same species within the same year were seperated into distinct territories. Maximum ndvi values from 2015 - 2020 are displayed as varying shades of green to show veg productivity in the study area.
<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/pefa_terr.png">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/rlha_terr.png">
</p>

* based on work done in rankin, it's clear that adults sometimes utilize alternate nesting sites around a central location. This results in spatial variation in nest site locations from year to year, and in many cases, such groups of nest sites can be considered as a single territory. Establishing hard boundaries that define these territories comes with assumptions however, and it may be best to explore nest site spatial dynamics among marked birds before accepting these assumptions. 
* Maybe there's a better way than using territories. Maybe we define a grid, and use state space count models to investigate variation in nest site densities in each cell per year.

<p align="center">
  <img width="950" src="https://github.com/emhedlin/bim_2020/blob/master/documents/territory_issue.jpg">
</p>



## Occupancy results

To investigate spatial patterns in occupancy, all nest sites were plotted for both species with a point size dependant on the number of years the site was occupied. PEFA have remained consistent since 2012, and have occupied roughly half of the sites each year. RLHAs have cycled heavily (see below), and have therefore occupied fewer sites overall.

<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/pefa_sum_occ.png">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/rlha_sum_occ.png">
</p>


<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/pefa_PAO.jpg">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/RLHA_PAO.jpg">
</p>
Occupancy dynamics corrected for detection error. Both species remain stable, with a high amount of variation in RLHA's (as to be expected).





