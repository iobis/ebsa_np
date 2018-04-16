library(robis)
library(sf)
library(dplyr)

generate_ebsa_summaries <- function() {
  poly <- sf::st_read('data/Global_EBSAs_Automated_Final_1104_2016_WGS84/Global_EBSAs_Automated_Final_1104_2016_WGS84.shp', stringsAsFactors = FALSE)
  ids <- unique(poly$GLOBAL_ID)

  for (id in ids) {
    ebsafile <- paste0('data/ebsa_', id, '.rds')
    if(!file.exists(ebsafile)) {
      ft <- poly %>% filter(GLOBAL_ID == id)
      info <- ft %>% select(NAME, Workshop, EBSA_ID, AREA_MW_KM, GLOBAL_ID) %>% distinct()
      data <- list(globalid=id, info=info)

      # get bounding box for querying
      geom <- sf::st_geometry(ft)
      geom <- sf::st_cast(geom, 'POLYGON')
      convex <- sf::st_convex_hull(geom)
      data$occ <- bind_rows(lapply(convex, function(poly) {
        wkt <- st_as_text(poly)
        occ <- robis2::occurrence(geometry = wkt)
        fields <- c('decimalLongitude', 'decimalLatitude', 'eventDate', 'depth', 'year', 'worms_id', 'valid_id', 'resource_id', 'institutionCode', 'taxonomicgroup')
        occ[!is.na(occ$worms_id),fields]
      }))
      obisid_table <- table(data$occ$valid_id)
      data$taxa <- bind_rows(lapply(unique(data$occ$valid_id), function(obisid) {
        t <- robis::taxon(obisid = obisid)
        cbind(t, ebsa_record_count = obisid_table[paste(obisid)])
      }))
      saveRDS(data, ebsafile)
    }
    outputfile <- rmarkdown::render('scripts/summaries.Rmd', output_file = paste0('ebsa_', id, '.pdf'), output_dir = 'reports', params = list(ebsafile = ebsafile))
  }
}
