Work for the 2020 BIM report
================
Erik Hedlin
29/09/2020

## Current status

``` r
occ %>% group_by(SiteID, year) %>% summarize(occ = max(pefa > 0, na.rm = TRUE)) %>% ungroup() %>% 
  group_by(year) %>% summarize(sum = sum_pefa(occ, na.rm = TRUE)) 
```

Summing the occupied sites for PEFA across years (2012 - 2020) reveals
some differences between the data and what’s been reported. I compared
the data used for analysis with the historical database spreadsheet
(both in the data folder), and one of the PEFA sites in 2017 was
incorrectly written as RLHA. Some errors obviously exist and we’ll now
have to go through each year again to find where the errors are.
