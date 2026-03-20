library(tidyverse)

# library(tidycensus)
# vars <- load_variables(2016, "acs5", cache = TRUE)
# cacs <- get_acs(
#   geography = "county",
#   year = 2016,
#   survey = "acs5",
#   variables = c(
#     
#     # Poverty
#     pov_total = "B17001_001",
#     pov_below = "B17001_002",
#     
#     # Race
#     race_total = "B03002_001",
#     white_nonhisp = "B03002_003",
#     
#     # Age total
#     age_total = "B01001_001",
#     
#     # Male 65+
#     m65_66 = "B01001_020",
#     m67_69 = "B01001_021",
#     m70_74 = "B01001_022",
#     m75_79 = "B01001_023",
#     m80_84 = "B01001_024",
#     m85p   = "B01001_025",
#     
#     # Female 65+
#     f65_66 = "B01001_044",
#     f67_69 = "B01001_045",
#     f70_74 = "B01001_046",
#     f75_79 = "B01001_047",
#     f80_84 = "B01001_048",
#     f85p   = "B01001_049"
#   )
# )
# 
# saveRDS(cacs, "../data/cacs.csv")

# Load in American Community Survey (Poverty Rate, Nonwhite Rate, Age 65+ Rate)
cacs <- readRDS("../data/cacs.csv")
cacs_wide <- cacs %>%
  select(GEOID, NAME, variable, estimate) %>%
  pivot_wider(names_from = variable, values_from = estimate) %>%
  mutate(
    poverty_rate = pov_below / pov_total,
    nonwhite_rate = 1 - white_nonhisp / race_total,
    age65plus = (
      m65_66 + m67_69 + m70_74 + m75_79 + m80_84 + m85p +
        f65_66 + f67_69 + f70_74 + f75_79 + f80_84 + f85p
    ) / age_total
  ) %>%
  select(GEOID, poverty_rate, nonwhite_rate, age65plus)
cacs_wide$STFIPS <- as.integer(substr(cacs_wide$GEOID, 1, nchar(cacs_wide$GEOID) - 3))
cacs_wide$CTYFIPS <- as.integer(substr(cacs_wide$GEOID, nchar(cacs_wide$GEOID) - 3 + 1, nchar(cacs_wide$GEOID)))

# Load in Age-adjusted Cardiovascular Disease Mortality Rate per 100,000
cvd <- read.csv('../data/cvd.csv')
cvd <- cvd[!is.na(cvd$County.Code), ]
cvd$STFIPS <- as.integer(substr(cvd$County.Code, 1, nchar(cvd$County.Code) - 3))
cvd$CTYFIPS <- as.integer(substr(cvd$County.Code, nchar(cvd$County.Code) - 3 + 1, nchar(cacs_wide$GEOID)))
cvd$CTYFIPS <- ifelse(cvd$STFIPS == 46 & cvd$CTYFIPS == 113, 102, cvd$CTYFIPS)
cvd$CTYFIPS <- ifelse(cvd$STFIPS == 2 & cvd$CTYFIPS == 270, 158, cvd$CTYFIPS)

# Load in Medicaid Expansion Status of each state
expansion <- read.csv('../data/expansion_modified.csv')
expansion <- expansion[-c(1, 53:75), ]

# Load in uninsurance rate
uninsurance <- read.csv('../data/uninsurance.csv')
uninsurance <- uninsurance %>% 
  filter(Year == 2016)
uninsurance <- uninsurance[-c(1:51), ]
uninsurance$STFIPS <- as.integer(substr(uninsurance$ID, 1, nchar(uninsurance$ID) - 3))
uninsurance$CTYFIPS <- as.integer(substr(uninsurance$ID, nchar(uninsurance$ID) - 3 + 1, nchar(cacs_wide$GEOID)))

# Load in urban/rural classification
urban <- read.csv('../data/urban.csv')

# Merge
df_list <- list(cacs_wide, cvd, uninsurance, urban)
df <- df_list %>%
  reduce(full_join, by = c("STFIPS", "CTYFIPS"))

df <- full_join(df, expansion, by = 'STFIPS') %>%
  filter(STFIPS != 72, County != "Bedford city, VA" & 
           County != "Clifton Forge city, VA" |
           is.na(County),
           CTYNAME != "Chugach Census Area" &
           CTYNAME != "Copper River Census Area" &
           CTYNAME != "Wade Hampton Census Area" &
           CTYNAME != "Prince of Wales-Outer Ket" &
           CTYNAME != "Skagway-Hoonah-Angoon Census Area" &
           CTYNAME != "Wrangell-Petersburg Census Area" &
         !grepl('Planning Region', CTYNAME) | is.na(CTYNAME),
         !(STFIPS == 46 & CTYFIPS == 113)) %>%
  select(STFIPS, CTYFIPS, Population, poverty_rate, nonwhite_rate, age65plus, Age.Adjusted.Rate, Pct, 
         CODE2013, medicaid_expansion_status) %>%
  # filter((Age.Adjusted.Rate != 'Suppressed') &  %>%
  #          !grepl('Unreliable', Age.Adjusted.Rate)) %>%
  filter(!(STFIPS == 15 & CTYFIPS == 5)) %>%
  mutate(population = as.numeric(Population),
         cvd = ifelse(Age.Adjusted.Rate == "Suppressed", "", Age.Adjusted.Rate),
         cvd = parse_number(cvd),
         uninsurance = Pct / 100,
         # rurality = factor(CODE2013, 
         #                   ordered = TRUE, 
         #                   levels = c(1, 2, 3, 4, 5, 6)),
         rurality = as.integer(CODE2013),
         medicaid = as.integer(medicaid_expansion_status)
  ) %>%
  select(-c(Population, Age.Adjusted.Rate, Pct, CODE2013, medicaid_expansion_status))




