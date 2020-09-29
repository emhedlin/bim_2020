library(tidyverse)


# join 2020 survey data and historical ------------------------------------

hist <- read_csv("data/2012-2020_occ.csv")
current <- read_csv("data/2020_survey.csv") 



# clean 2020 survey data --------------------------------------------------
names(current)
current <- current %>% select("SiteID" = "Site",
                              "y20v1" = "20_Sp_1",
                              "y20v2" = "20_Sp_2",
                              "y20v3" = "20_Sp_3")

full <- current %>% left_join(hist, by= "SiteID")

# filter out sites not included in 2019 analysis because they haven't been occupied since 2011
full <- full %>% filter(!SiteID %in% c(5, 23, 52, 73, 104, 117)) 



# double check numbers with last years report -----------------------------

#2019
  #pefa = 39
  #rlha = 11
#2018 
  #pefa = 48 
  #rlha = 11
occ <- full %>% select(1:4, 19:42) 
occ <- occ %>% pivot_longer(cols = 2:ncol(occ)) %>% mutate(year = str_sub(name, 2,3), survey = str_sub(name, 5,5)) %>% 
  mutate(value = ifelse(value == "not detected" | is.na(value), 0, value)) %>%
  mutate(pefa = ifelse(value == "PEFA", 1, 0), rlha = ifelse(value == "RLHA", 1, 0))

# lists the total pefa sites occupied in each year. Discrepancy with report, and this needs to be looked at year by year.
occ %>% group_by(SiteID, year) %>% summarize(occ = max(pefa > 0, na.rm = TRUE)) %>% ungroup() %>% 
  group_by(year) %>% summarize(sum = sum(occ, na.rm = TRUE)) 

View(occ %>% filter(SiteID == "2"))
  

range(grp$pefa, na.rm = TRUE)


pefa <- ifelse(occ[,2:ncol(occ)] == "PEFA", 1, 0) # convert pefa occupied to binary


?pivot_longer
pefa <- as.data.frame(pefa)

pefa <- cbind(pefa, sum = rowSums(pefa, na.rm = "TRUE")) # tally occupancy across years
pefa <- cbind(dat[,1], pefa) # add site ID's and locations back in
pefa <- subset(pefa, sum !=0) # drop sites that were never occupied 2012-2019
pefa <- inner_join(pefa, md, by = "SiteID")
pefa <- inner_join(pefa, dnon, by = "SiteID")
pefa <- pefa %>% dplyr::select(-sum)
dim(pefa) # 94 sites