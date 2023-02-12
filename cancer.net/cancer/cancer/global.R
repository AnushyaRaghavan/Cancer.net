library(dplyr)

allzips <- read.csv(file.choose(),header = T)
allzips$latitude <- jitter(allzips$latitude)
allzips$longitude <- jitter(allzips$longitude)
allzips$zipcode <- formatC(allzips$pincode, width=5, format="d", flag="0")
row.names(allzips) <- allzips$pincode

cleantable <- allzips %>%
  select(
   Zipcode = zipcode,
   Lat = latitude,
    Long = longitude
  )
