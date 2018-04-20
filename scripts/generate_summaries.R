library(robis)
library(sf)
library(dplyr)

generate_ebsa_summaries <- function(overwrite_data=FALSE, overwrite_reports=FALSE) {
  poly <- sf::st_read('data/Global_EBSAs_Automated_Final_1104_2016_WGS84/Global_EBSAs_Automated_Final_1104_2016_WGS84.shp', stringsAsFactors = FALSE)
  ids <- unique(poly$GLOBAL_ID)

  for (id in ids) {
    print(id)
    ebsafile <- paste0('data/ebsa_', id, '.rds')
    if(!file.exists(ebsafile) || overwrite_data) {
      ft <- poly %>% filter(GLOBAL_ID == id)
      info <- ft %>% select(NAME, Workshop, EBSA_ID, AREA_MW_KM, GLOBAL_ID) %>% distinct()
      data <- list(globalid=id, info=info)

      # get bounding box for querying
      geom <- sf::st_geometry(ft)
      geom <- sf::st_cast(geom, 'POLYGON')
      data$geom <- geom
      convex <- sf::st_convex_hull(geom)
      convex <- sf::st_buffer(convex, 0.1)
      convex <- sf::st_simplify(convex, preserveTopology = TRUE, dTolerance = 0.05)
      data$occ <- bind_rows(lapply(convex, function(poly) {
        wkt <- st_as_text(poly)
        occ <- robis2::occurrence(geometry = wkt)
        fields <- c('decimalLongitude', 'decimalLatitude', 'eventDate', 'depth', 'year', 'worms_id', 'valid_id', 'resource_id', 'institutionCode', 'taxonomicgroup')
        if(NROW(occ) > 0) {
          if(!all(fields %in% colnames(occ))) {
            warning("Field not found")
            warning(fields[!fields %in% colnames(occ)])
          }
          occ <- occ[!is.na(occ$worms_id),fields]
        }
        occ
      }))

      # filter occurrences with real geometries
      occ_sp <- sf::st_as_sf(data$occ, coords=c('decimalLongitude', 'decimalLatitude'))
      st_crs(occ_sp) <- st_crs(geom)
      occ_intersects <- unlist(sf::st_intersects(geom, occ_sp, sparse=TRUE, prepared = TRUE))
      data$occ <- data$occ[occ_intersects,]

      obisid_table <- table(data$occ$valid_id)
      data$taxa <- bind_rows(lapply(unique(data$occ$valid_id), function(obisid) {
        t <- robis::taxon(obisid = obisid)
        cbind(t, ebsa_record_count = obisid_table[paste(obisid)])
      }))

      # get datasets
      datasets <- lapply(unique(data$occ$resource_id), function(id) {
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

      # get institutes
      data$institutes <- bind_rows(lapply(datasets, function(d) d$institutes)) %>%
        left_join(datasetrecordcount) %>%
        left_join(datasettaxacount) %>%
        group_by(id, name, acronym, parent) %>%
        summarise(datasets = n(), ebsa_record_count = sum(ebsa_record_count), ebsa_taxa_count = sum(ebsa_taxa_count))

      saveRDS(data, ebsafile)
    }
    if(!file.exists(paste0('reports/ebsa_', id, '.pdf')) || overwrite_reports) {
      outputfile <- rmarkdown::render('scripts/summaries.Rmd', output_file  = paste0('ebsa_', id, '.pdf'), output_dir = 'reports', params = list(id=id, ebsafile = paste0('../', ebsafile)))
    }
  }
}
# generate_ebsa_summaries()
