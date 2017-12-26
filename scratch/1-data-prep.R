library(foreign)
library(tidyverse)
library(stringr)
library(mice)
library(forcats)

mice.impute.ranger = function (y, ry, x, ntree = 10, wy = NULL, type = NULL, ...) 
{
  # had to add wy = NULL so that this conforms with
  # the updated MICE specification. It is ignored
  
  if (!requireNamespace("ranger", quietly = TRUE)) 
    stop("Package 'ranger' needed fo this function \n to work. Please install it.", 
         call. = FALSE)
  
  obs = x[ry, ]
  mis = x[!ry, ]
  obs$y = y[ry]
  do_ordered = is.numeric(y) | (is.factor(y) & length(levels(y)) == 2)
  f = ranger::ranger(dependent.variable.name = "y", data = obs, replace = FALSE, 
                     num.trees = ntree, respect.unordered.factors = do_ordered, 
                     min.node.size = 5, ...)
  o_nodes = data.frame(matrix(predict(f, data = obs, predict.all = T, 
                                      type = "terminalNodes")$predictions, ncol = ntree)) %>% 
    mutate(donor_id = 1:n()) %>% gather(tree, node, -donor_id)
  mis1 = FALSE
  if (nrow(mis) == 1) {
    mis1 = TRUE
    mis = bind_rows(mis, mis)
  }
  p = data.frame(matrix(predict(f, data = mis, predict.all = TRUE, 
                                type = "terminalNodes")$predictions, ncol = ntree)) %>% 
    mutate(id = 1:n()) %>% gather(tree, node, -id) %>% group_by(id) %>% 
    slice(sample(1:n(), size = 1)) %>% left_join(o_nodes, 
                                                 by = c("tree", "node")) %>% group_by(id) %>% slice(sample(1:n(), 
                                                                                                           size = 1))
  if (mis1) {
    p = p[1, ]
  }
  impute = obs$y[p$donor_id]
  return(impute)
}




full_data = read.spss("scratch/raw/ALL_MERGED_WEIGHTED.sav", 
                      use.value.labels = T, 
                      to.data.frame = T, 
                      trim.factor.names = T, 
                      trim_values = T, 
                      add.undeclared.levels = "no")

full_data$sample_id = str_replace(full_data$Survey_mask, "Panel ", "")
full_data$sample_id = str_replace(full_data$sample_id, "J", "MTk")

full_data$fcregion = factor(full_data$fcregion)

full_data$educcat = car::recode(as.numeric(full_data$Q0050_RECODE), "1:3='HS or Less'; 4:5='Some College'; 6:8='College Grad'; else=NA", 
                                as.factor.result=T, levels=c('HS or Less', 'Some College', 'College Grad'))

full_data$sex = car::recode(full_data$Q0042, "'Refused'=NA")
full_data$age = car::recode(full_data$Q0043, "999=NA")

# Topcode age at 85 to be consistent with CPS
full_data$age = ifelse(full_data$age > 85, 85, full_data$age)

levels(full_data$racethn) = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other", "Refused")
full_data$racethn = car::recode(full_data$racethn, "'Refused'=NA", 
                                levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other"))
# Create analytic file

np = full_data %>%
  # Filter down to nonprob samples of interest
  filter(sample_id %in% c("ABS", "ATP", "K") == FALSE) %>% 
  select(sample_id,
         rid=RID,
         # Civic Engagement Variables
         talk_neighbor = Q0002,
         trust_neighbors = Q0004,
         school_group = Q0009,
         civic_assoc = Q0010,
         recreational_assoc = Q0011,
         vote_local = Q0029,
         
         # Basic demos
         age,
         sex,
         racethn,
         educcat,
         fcregion
         
  ) %>%
  mutate(  
    # CPS Civic Engagement outcome variables
    y_talk_neighbor_weekly = talk_neighbor %in% c("Basically every day", "A few times a week"),
    y_trust_neighbors_all_most = trust_neighbors %in% c("All of the people in your neighborhood", "Most of the people in your neighborhood"),
    y_school_group_yes = school_group=="Yes",
    y_civic_assoc_yes = civic_assoc=="Yes",
    y_recreational_assoc_yes = recreational_assoc=="Yes",
    y_always_vote_local = vote_local == "Always vote in local elections"
  ) %>%
  group_by(sample_id) %>%
  mutate(np_index = 1:n()) %>%
  ungroup()

# Split nonprobs into separate samples and impute missingness for use in
# modeling
set.seed(1234)
np_imputed = np %>% split(.$sample_id) %>%
  map(~complete(mice(., m=1, method="ranger"))) %>%
  bind_rows()



saveRDS(np_imputed, "scratch/processed/cleaned_np_civic_data.RDS")



### CPS Civic Engagement
cps_civic = readRDS("scratch/raw/cps_nov13_civic_adult_benchmarks.RDS") %>% 
  # Select self responses from supplement
  filter(PWSRWGT > 0) %>%
  select(
    y_talk_neighbor_weekly = Q0002_BMCAT,
    y_trust_neighbors_all_most = Q0004_BMCAT,
    y_school_group_yes = Q0009_BMCAT,
    y_civic_assoc_yes = Q0010_BMCAT,
    y_recreational_assoc_yes = Q0011_BMCAT,
    y_always_vote_local = Q0029_BMCAT,
    age=PRTAGE,
    sex,
    racethn,
    educcat, 
    fcregion=GEREG,
    state = GESTCEN,
    pwsrwgt = PWSRWGT) %>%
  mutate(rid = 9999999 + 1:n(),
         sex = droplevels(sex),
         racethn = droplevels(racethn),
         educcat = droplevels(educcat),
         fcregion = factor(fcregion),
         state = fct_recode(factor(state), 
                            "ME" = "11", 
                            "DE" = "51", 
                            "NM" = "85",
                            "NH" = "12", 
                            "MD" = "52", 
                            "AZ" = "86",
                            "VT" = "13", 
                            "DC" = "53", 
                            "UT" = "87",
                            "MA" = "14", 
                            "VA" = "54", 
                            "NV" = "88",
                            "RI" = "15", 
                            "WV" = "55", 
                            "WA" = "91",
                            "CT" = "16", 
                            "NC" = "56", 
                            "OR" = "92",
                            "NY" = "21", 
                            "SC" = "57", 
                            "CA" = "93",
                            "NJ" = "22", 
                            "GA" = "58", 
                            "AK" = "94",
                            "PA" = "23", 
                            "FL" = "59", 
                            "HI" = "95",
                            "OH" = "31", 
                            "KY" = "61",
                            "IN" = "32", 
                            "TN" = "62",
                            "IL" = "33", 
                            "AL" = "63",
                            "MI" = "34", 
                            "MS" = "64",
                            "WI" = "35", 
                            "AR" = "71",
                            "MN" = "41", 
                            "LA" = "72",
                            "IA" = "42", 
                            "OK" = "73",
                            "MO" = "43", 
                            "TX" = "74",
                            "ND" = "44", 
                            "MT" = "81",
                            "SD" = "45", 
                            "ID" = "82",
                            "NE" = "46", 
                            "WY" = "83",
                            "KS" = "47", 
                            "CO" = "84"),
         sample_id = "CPS Civic",
         cps_index = 1:n())

saveRDS(cps_civic, "scratch/processed/cps_civic_full_edited.RDS")


