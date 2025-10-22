library(ggplot2)
library(lubridate)
wd <- '/home/victor/projects/intrareprogrowth'

years <- 1997:2018
dat <- lapply(years, function(y){
  s <- ifelse(y > 2010, 25, ifelse(y == 2010, 24, 23))
  read.table(file.path(wd, 'data/zweifel2021', paste0('Zweifel-etal_2021_TreeNet',y,'.tab')), skip = s, sep = '\t', header = TRUE)[,1:6]
})
dat <- as.data.frame(do.call(rbind, dat))

names(dat) <- c('date_time', 'id_plot', 'id_tree', 'species', 'growth', 'growth_rate')

Ntrees <- length(unique(dat$id_tree))
Nspecies <- length(unique(dat$species))

uniq_trees_idxs <- unique(dat$id_tree)

dat$date <- format(as.Date(dat$date_time), "%Y-%m-%d")
dat$year <- format(as.Date(dat$date), "%Y")
dat$month_day <- format(as.Date(dat$date), "%m-%d")

# Look at one tree!
subset_onetree <- na.omit(dat[dat$id_tree %in% uniq_trees_idxs[1],])
subset_onetree$week <- week(subset_onetree$date)
count_obs <- aggregate(month_day ~ year, FUN = function(x) length(unique(x)), data = subset_onetree)
names(count_obs)[2] <- 'nobs'
mean_growthrate <- aggregate(growth_rate ~ week, FUN = mean, data = subset_onetree[subset_onetree$year %in% count_obs[count_obs$nobs > 340, 'year'],])
names(mean_growthrate)[1] <- 'week'
ggplot(data = mean_growthrate, aes(x = week, y = growth_rate)) +
  # geom_point() +
  geom_bar(stat = "identity")

# Look at one species!
mean_growthrate_perspecies <- data.frame()
for(sp in unique(dat$species)){
  if(sp %in% c('Sorbus aria', 'Que', 'Larix decidua')){next}
  subset_sp <- na.omit(dat[dat$species == sp,])
  subset_sp$week <- week(subset_sp$date)
  count_obs <- aggregate(month_day ~ id_tree + year, FUN = function(x) length(unique(x)), data = subset_sp)
  names(count_obs)[3] <- 'nobs'
  subset_sp <- merge(subset_sp, count_obs)
  mean_growthrate <- aggregate(growth_rate ~ week + year, FUN = mean, data = subset_sp[subset_sp$nobs > 300,])
  mean_growthrate$species <- sp
  mean_growthrate_perspecies <- rbind(mean_growthrate_perspecies, mean_growthrate)
}

mean_growthrate_perspecies$week_date <- lubridate::ymd(mean_growthrate_perspecies$year) + lubridate::weeks(mean_growthrate_perspecies$week - 1)

ggplot(data = mean_growthrate_perspecies, aes(x = week_date, y = growth_rate)) +
  facet_wrap(~species) +
  geom_bar(stat = "identity") +
  scale_x_date(date_breaks = "1 month", date_labels = "%m") +
  geom_segment(data = species_cycle, aes(x = start_fruitmat, xend = end_fruitmat, y = 0.8, yend = 0.8),
               linewidth = 2, col = 'darkorange') +
  theme_classic()

species_cycle <-
  data.frame(species = 'Fagus sylvatica',
             start_fruitmat = as.Date("2000-07-01"), end_fruitmat = as.Date("2000-09-01"))

               