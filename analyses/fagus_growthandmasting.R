rm(list = ls())
wd <- '/home/victor/projects/intratradeoff'
library(ggplot2)
library(terra)
library(tidyterra)
library(rnaturalearth)
world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
suisse <- aggregate(vect(ne_states(country = c("Switzerland"), returnclass = "sf")))

# Load MastWeb data for Fagus sylvatica
years <- 1997:2018
mastdata <- lapply(years, function(y){
  read.csv(file.path(wd, 'data/mastweb/fagus_sylvatica', paste0(y,'.csv')))
})
mastdata <- as.data.frame(do.call(rbind, mastdata))
mastdata <- mastdata[,c('Species', 'Year', 'Type', 'Coordinate..WGS.84.')]
names(mastdata) <- c('vernacular_species', 'year', 'mast_type', 'coordinates')
mastdata$lat <- as.numeric(stringr::str_split_i(mastdata$coordinates, ",", 1))
mastdata$lon <- as.numeric(stringr::str_split_i(mastdata$coordinates, ",", 2))

# Load dendrometer data
years <- 1997:2018
dat <- lapply(years, function(y){
  s <- ifelse(y > 2010, 25, ifelse(y == 2010, 24, 23))
  read.table(file.path(wd, 'data/zweifel2021', paste0('Zweifel-etal_2021_TreeNet',y,'.tab')), skip = s, sep = '\t', header = TRUE)[,1:6]
})
dat <- as.data.frame(do.call(rbind, dat))
names(dat) <- c('date_time', 'id_plot', 'id_tree', 'species', 'growth', 'growth_rate')

# Sites
sites <- read.table(file.path(wd, 'data/zweifel2021', 'Zweifel-etal_2021_site-info.tab'), skip = 20, sep = '\t', header = TRUE)

# Look at Fagus growth rate
subset_sp <- na.omit(dat[dat$species == "Fagus sylvatica",])
subset_sp$date <- format(as.Date(subset_sp$date_time), "%Y-%m-%d")
subset_sp$year <- format(as.Date(subset_sp$date), "%Y")
subset_sp$month_day <- format(as.Date(subset_sp$date), "%m-%d")
subset_sp$week <- week(subset_sp$date)
count_obs <- aggregate(month_day ~ id_tree + year, FUN = function(x) length(unique(x)), data = subset_sp)
names(count_obs)[3] <- 'nobs'
subset_sp <- merge(subset_sp, count_obs)
mean_growthrate <- aggregate(growth_rate ~ week + year + id_plot, FUN = mean, data = subset_sp[subset_sp$nobs > 300,])

# Prepare dendrometer site data
fagus_sites <- unique(mean_growthrate$id_plot)
metadata_sites <- unique(sites$Site)
metadata_sites <- gsub("-","_",metadata_sites)
corrected_sites <- data.frame(original = unique(mean_growthrate$id_plot), corrected = NA)
corrected_sites[(corrected_sites$original %in% metadata_sites),'corrected'] <- 
  corrected_sites[(corrected_sites$original %in% metadata_sites),'original']
corrected_sites[!(corrected_sites$original %in% metadata_sites) & 
                  grepl('Saillon', corrected_sites$original), 'corrected'] <- 'Saillon_A780'
corrected_sites[!(corrected_sites$original %in% metadata_sites) &
                  corrected_sites$original == 'Sempach_NA', 'corrected'] <- 'Sempach_Forest'
corrected_sites[!(corrected_sites$original %in% metadata_sites) &
                  corrected_sites$original == 'Vordemwald_NA', 'corrected'] <- "Vordemwald_Forest"
corrected_sites[!(corrected_sites$original %in% metadata_sites) &
                  corrected_sites$original == 'Muri_Beech', 'corrected'] <- "Muri_Forest"
corrected_sites[!(corrected_sites$original %in% metadata_sites) &
                  corrected_sites$original == 'Lausanne_NA', 'corrected'] <- "Lausanne_Forest"
corrected_sites[!(corrected_sites$original %in% metadata_sites) &
                  corrected_sites$original == 'Birmensdorf_NA', 'corrected'] <- "Birmensdorf_Forest"
corrected_sites[!(corrected_sites$original %in% metadata_sites) &
                  corrected_sites$original == 'Laegeren_FF', 'corrected'] <- "Laegeren"
site_coordinates <- unique(sites[c('Site', 'Latitude', 'Longitude')])
site_coordinates$Latitude <- round(site_coordinates$Latitude, 2)
site_coordinates$Longitude <- round(site_coordinates$Longitude, 2)
site_coordinates <- unique(site_coordinates)
site_coordinates$Site <- gsub("-","_",site_coordinates$Site)
site_coordinates <- rbind(
  site_coordinates,
  data.frame(Site = "Laegeren", Latitude = 47.478333, Longitude = 8.364389)
)
names(site_coordinates) <- c('corrected', 'lat', 'lon')
correctedsite_coordinates <- merge(corrected_sites, site_coordinates)

# Quick map
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(suisse)
points(vect(unique(correctedsite_coordinates[c('lat', 'lon')]), geom = c('lon', 'lat')), pch = 2, cex = 2, col = 'darkgreen')
points(vect(unique(mastdata[c('lat', 'lon')]), geom = c('lon', 'lat')))

# Mast data in numeric
mastdata$numeric_mast <- NA
mastdata$numeric_mast <- ifelse(mastdata$mast_type == 'Full mast', 1,
                                ifelse(mastdata$mast_type == 'Half-mast', 0.65,
                                       ifelse(mastdata$mast_type == 'Partial mast', 0.3,
                                              ifelse(mastdata$mast_type == 'No mast', 0, NA))))


# Check buffers
dendro_sites <- vect(correctedsite_coordinates[c('lat', 'lon')], geom = c('lon', 'lat'))
crs(dendro_sites) <- 'EPSG:4326'
buffer_effect <- data.frame()
for(b in c(1e3, 5e3, 10e3, 20e3, 50e3, 75e3, 100e3)){
  dendro_buffers <- buffer(dendro_sites, b)
  mastingstate_sites <- data.frame()
  for(y in unique(mean_growthrate$year)){
    mastdata_y <- vect(mastdata[mastdata$year == y, c('lat', 'lon', 'numeric_mast')], geom = c('lon', 'lat'), crs = crs(dendro_sites))
    r <- relate(dendro_buffers, mastdata_y, "intersects")
    count <- apply(r, 1, function(i) length(mastdata_y$numeric_mast[i]))
    num <- apply(r, 1, function(i) mean(mastdata_y$numeric_mast[i]))
    mastingstate_sites <- rbind(
      mastingstate_sites,
      cbind(correctedsite_coordinates, year = y, mast = num, count = count)
    )
  }
  buffer_effect <- rbind(
    buffer_effect,
    data.frame(buffer = b/1000, unobs = sum(is.na(mastingstate_sites$mast))/(7*19), mean_count = mean(mastingstate_sites$count))
  )
}
par(mfrow = c(2,1), mar = c(4,4,1,1))
plot(unobs*100~buffer, data = buffer_effect, xlab = 'Buffer distance (km)', ylab = 'Unobserved year-site (%)', type = 'l')
plot(mean_count~buffer, data = buffer_effect, xlab = 'Buffer distance (km)', ylab = 'No. of mastweb obs. per site (mean)', type = 'l')

# Cross both data
dendro_buffers <- buffer(dendro_sites, 75e3)
mastingstate_sites <- data.frame()
for(y in unique(mean_growthrate$year)){
  mastdata_y <- vect(mastdata[mastdata$year == y, c('lat', 'lon', 'numeric_mast')], geom = c('lon', 'lat'), crs = crs(dendro_sites))
  r <- relate(dendro_buffers, mastdata_y, "intersects")
  count <- apply(r, 1, function(i) length(mastdata_y$numeric_mast[i]))
  num <- apply(r, 1, function(i) mean(mastdata_y$numeric_mast[i]))
  numsd <- apply(r, 1, function(i) sd(mastdata_y$numeric_mast[i]))
  mastingstate_sites <- rbind(
    mastingstate_sites,
    cbind(correctedsite_coordinates, year = y, mast_mean = num, mast_sd = numsd, count = count)
  )
}
mastingstate_sites <- mastingstate_sites[c('original', 'year', 'mast_mean', 'mast_sd', 'count')]
names(mastingstate_sites)[1] <- 'id_plot'
growth_masting <- merge(mean_growthrate, mastingstate_sites)

# Quick figure
dendromast_figure <- ggplot(data = growth_masting) +
  facet_grid(id_plot ~ year) +
  geom_line(aes(x = week, y = growth_rate), col = 'black', linewidth = 0.8) + 
  geom_line(aes(x = week, y = growth_rate, col = mast_mean), linewidth = 0.6) + 
  scale_color_distiller(palette = "Reds", direction = 1,  na.value = "grey80", name = 'Masting state') +
  theme_classic() +
  labs(x = 'Week', y = "Growth rate (µm/h)") +
  theme(strip.background = element_blank(), strip.text.y = element_text(size = 5),
        axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 7),
        legend.direction = "horizontal", legend.position = 'inside',
        legend.position.inside = c( 0.15, 0.05), legend.title.position = 'top',
        legend.key.height = unit(0.2, 'cm'))
ggsave(dendromast_figure, filename = file.path(wd, 'analyses/figures', 'dendromast.pdf'), width = 210, height = 297, units = "mm")
