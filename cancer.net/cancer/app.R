
library(leaflet)

########################
# Simple World Map
########################
my_map <-  leaflet() %>% 
  addTiles() 
my_map 
df <-read.csv("superzip.csv") 

prepare_map <- df %>% leaflet() %>%
  addProviderTiles(providers$MtbMap) %>%
  addProviderTiles(providers$Stamen.TonerLines,
                   options = providerTileOptions(opacity = 0.35)) %>%
  addProviderTiles(providers$Stamen.TonerLabels)


my_map <- prepare_map %>%
  addMarkers(clusterOptions = markerClusterOptions())
my_map