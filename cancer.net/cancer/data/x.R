install.packages("leaflet").

library("rio")
superzip<- import("superzip.csv")
export(superzip, "superzip.csv")
# 
# # convert Stata to SPSS
convert("superzip.csv", "superzip.rds")