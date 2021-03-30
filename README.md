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
Notes for 2021.

* 2019 report indicated 76 pefa territories and 71 rlha territories.
* In 2020, code to check for occupancy among > 1 nest sites in each territory, which would split territories, was added in ```04_territory.R`` 
* this resulted in 84 pefa and 87 rlha
* seemed like a large increase from 2019 so we're manually checking.

The code is correct in separating these territories based on the rules outlined in 2019. However, we're updating the rules to be more discriminatory. According to the new rules, if two sites are occupied by the same species, are within 1 km, and produce eggs, then the territory should be split up. 

#### Territory manual checks
Territories below were flagged

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PEFA

25  - 18 - * male flushed on first survey, no detections on following survey
475 - 18 - * nest failed, no confirmation of eggs

481 - 19 

338 - 18 - * detected on 3rd survey, 2 young
59  - 18 - * no confirmation of eggs

359 - 19   * keep together
380 - 19

485 - 18 - no eggs confirmed
86  - 18 - 3 young.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RLHA

359 - 2015
460 - 2015 - no eggs at 460, don't split this site up

86 - 2015 - eggs at both sites, split up
88 - 2015

48 - 2012 - no eggs, don't split this territory up
51 - 2012

466 - 2015 - eggs at both sites, split
477

359 2015
460 2015 - no eggs at this site, don't split

21  2012 - eggs at both, split these sites
360 2012





## Breeding Success

Given a site was occupied, we modeled the chance that each site successfully reared young using mixed logistic models. We examined potential spatial patterns in breeding success via matern covariance. Three types of spatial patterns were estimated:
* the pattern remains the same across all years
* the pattern changes with each year (iid)
* the pattern changes from year to year with temporal correlation of one step (AR1)

For Peregrines, the model that performed the best included a spatial pattern that was fixed across all years (below left).

<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/mesh.png">
</p>


<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/pefa_spatial_pattern.png">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/rlha_sum_occ.png">
</p>


## Occupancy
Occupancy dynamics corrected for detection error. Both species remain stable, with a high amount of variation in RLHA's (as is expected).


<p align="center">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/pefa_PAO.jpg">
  <img width="400" src="https://github.com/emhedlin/bim_2020/blob/master/documents/RLHA_PAO.jpg">
</p>




