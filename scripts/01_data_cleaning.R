library(tidyverse)


# join 2020 survey data and historical ------------------------------------

hist <- read_csv("data/2012-2020_occ.csv")
current <- read_csv("data/2020_survey.csv") 
db <- read_csv("data/2019_BIM_db.csv")
current %>% filter(SiteID == "52")
db_occ <- db %>% select(1, 66, 74, 82, 90, 98, 106, 114)


occ %>% filter(SiteID == "52")
unique(occ$SiteID)
# clean 2020 survey data --------------------------------------------------
names(current)
current <- current %>% select("SiteID" = "Site",
                              "y20v1" = "20_Sp_1",
                              "y20v2" = "20_Sp_2",
                              "y20v3" = "20_Sp_3")

full <- current %>% left_join(hist, by= "SiteID")
full %>% filter(SiteID == "5")
# filter out sites not included in 2019 analysis because they haven't been occupied since 2011
full <- full %>% filter(!SiteID %in% c(5, 23, 52, 73, 104, 117)) 



# tables -----------------------------

occ <- full %>% select(1:4, 6:7, 20:42, -X35) 
write.csv(occ, "data/2012-2020.csv")
occ <- occ %>% pivot_longer(cols = 2:ncol(occ)) %>% mutate(year = str_sub(name, 2,3), survey = str_sub(name, 5,5)) %>% 
  mutate(value = ifelse(value == "not detected" | is.na(value), 0, value)) %>%
  mutate(pefa = ifelse(value == "PEFA", 1, 0), rlha = ifelse(value == "RLHA", 1, 0))

# total pefa sites occupied in each year. 
occ %>% group_by(SiteID, year) %>% summarize(occ = max(pefa > 0, na.rm = TRUE)) %>% ungroup() %>% 
  group_by(year) %>% summarize(sum = sum(occ, na.rm = TRUE)) 


# total rlha sites occupied in each year.
occ %>% group_by(SiteID, year) %>% summarize(occ = max(rlha > 0, na.rm = TRUE)) %>% ungroup() %>% 
  group_by(year) %>% summarize(sum = sum(occ, na.rm = TRUE)) 



# n years occupied by site ------------------------------------------------
occ %>% group_by(SiteID, year) %>% summarize(occ = max(pefa > 0, na.rm = TRUE)) %>% ungroup() %>%
  group_by(SiteID) %>% summarize(sum = sum(occ, na.rm = TRUE)) %>% arrange(desc(sum))






# compare 2019 data with db -----------------------------------------------

db_occ[,2:ncol(db_occ)] <- ifelse(db_occ[,2:ncol(db_occ)] == "PEFA", 1, 0)


       