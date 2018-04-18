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

      datasets <- lapply(unique(occ$resource_id), function(id) {
        r <- httr::GET(paste0('http://api.iobis.org/resource/', id))
        d <- httr::content(r)
        provider <- ifelse(is.null(d$provider$name), '', d$provider$name)
        institutes <- bind_rows(d$institutes)
        institutes$resource_id <- rep(as.integer(id), NROW(institutes))
        list(dataset=data_frame(resource_id=id, name=d$name, node=d$node$name, provider=provider, taxon_cnt=d$taxon_cnt, record_cnt=d$record_cnt),
             institutes=institutes)
      })

      datasetrecordcount <- data$occ %>%
        group_by(resource_id) %>%
        summarise(ebsa_record_count=n())
      datasettaxacount <- data$occ %>%
        select(resource_id, valid_id) %>%
        distinct() %>%
        group_by(resource_id) %>%
        summarise(ebsa_taxa_count = n())

      data$datasets <- bind_rows(lapply(datasets, function(d) d$dataset)) %>%
        left_join(datasetrecordcount) %>%
        left_join(datasettaxacount)

      data$institutes <- bind_rows(lapply(datasets, function(d) d$institutes)) %>%
        left_join(datasetrecordcount) %>%
        left_join(datasettaxacount) %>%
        group_by(id, name, acronym, parent) %>%
        summarise(datasets = n(), ebsa_record_count = sum(ebsa_record_count), ebsa_taxa_count = sum(ebsa_taxa_count))

      saveRDS(data, ebsafile)
    }
    outputfile <- rmarkdown::render('scripts/summaries.Rmd', output_file = paste0('ebsa_', id, '.pdf'), output_dir = 'reports', params = list(ebsafile = ebsafile))
  }
}
