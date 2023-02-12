library(dplyr)

allzips <- readRDS("data/superzip.rds")
allzips$latitude <- jitter(allzips$latitude)
allzips$longitude <- jitter(allzips$longitude)
allzips$zip_code <- formatC(allzips$zip_code, width=5, format="d", flag="0")
row.names(allzips) <- allzips$zip_code

cleantable <- allzips %>%
  select(
    Place=place,
    Community = community,
    Zipcode = zip_code,
    All_cancersites=all_cancer_sites,
    Tongue=tongue,
    Mouth=mouth,
    Breast=breast,
    Lungs=lungs,
    Stomach=stomach,
    Bladder=bladder,
    Sex=sex,
    Lat = latitude,
    Long = longitude
  )