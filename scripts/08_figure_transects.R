library(sf)
library(tidyverse)
library(lubridate)
library(ggpubr)


# NBAR figure -------------------------------------------------------------

prepareData <- function(transID, dem, timesat, thres) {
  transDEM <- read_sf(dem) %>%
    as_tibble() %>%
    mutate(ID = as.character(geometry)) %>%
    select(ID, DEM = AVE_DSM)
  
  files <- Sys.glob(timesat) %>% 
    str_subset(thres) %>% str_subset("M_loeFiltered")
  
  tranVIpheno <- NULL
  for (i in files) {
    VIname <- str_extract(i, "NDVI|EVI2|PPI")
    tmp <- read_csv(i) %>%
      transmute(ID, year,
                SOS = start_of_season,
                EOS = end_of_season,
                trans =  transID,
                VIname = VIname) %>% 
      filter(between(year, 2018, 2019)) %>% 
      mutate(lat = as.numeric(str_extract(ID, "[:blank:]\\d{2}\\.\\d{1,}"))) %>% 
      left_join(transDEM)
    tranVIpheno <- bind_rows(tranVIpheno, tmp)
  }
  return(tranVIpheno)
}

transScandic <- 
  prepareData("transScandic",
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic/*.csv",
              "0.25")

transSpain <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain/*.csv",
              "0.25") %>% 
  filter(between(lat, 42.15, 43.1))


transects <- bind_rows(transScandic, transSpain)

axisscale1 <- 10
tscandic <- transects %>%
  filter(year %in% c(2019), trans == "transScandic") %>%
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric",
               values_to = "pheno") %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")),
         metric = factor(metric, levels = c("SOS", "EOS"))) %>%
  
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) + 
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.15, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale1, color = "DEM"), size = 0.5, alpha = 0.1, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale1, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(0, 400),
                     sec.axis = sec_axis(trans = ~.*axisscale1,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(VIname), ncol = 3) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

axisscale2 <- 20
tspain <- transects %>%
  filter(year %in% c(2019), trans == "transSpain") %>%
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric",
               values_to = "pheno") %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")),
         metric = factor(metric, levels = c("SOS", "EOS"))) %>%
  
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) +
  scale_shape_manual("Land cover type", values = seq(0, 6)) +
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale2, color = "DEM"), size = 0.5, alpha = 0.15, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale2, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(-100, 500), breaks = seq(0, 500, 100),
                     sec.axis = sec_axis(trans = ~.*axisscale2,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(VIname), ncol = 3) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

ggarrange(tscandic, tspain, common.legend = T,
          labels = c("a", "b"), ncol = 1)

ggsave("figures/figure_transects_NBAR.pdf",
       width = 25, height = 25, units = "cm")


# TOC figure --------------------------------------------------------------

prepareData <- function(transID, dem, timesat, thres) {
  transDEM <- read_sf(dem) %>%
    as_tibble() %>%
    mutate(ID = as.character(geometry)) %>%
    select(ID, DEM = AVE_DSM)
  
  files <- Sys.glob(timesat) %>% 
    str_subset(thres) %>% str_subset("raw_loeFiltered")
  
  tranVIpheno <- NULL
  for (i in files) {
    VIname <- str_extract(i, "NDVI|EVI2|PPI")
    tmp <- read_csv(i) %>%
      transmute(ID, year,
                SOS = start_of_season,
                EOS = end_of_season,
                trans =  transID,
                VIname = VIname) %>% 
      filter(between(year, 2018, 2019)) %>% 
      mutate(lat = as.numeric(str_extract(ID, "[:blank:]\\d{2}\\.\\d{1,}"))) %>% 
      left_join(transDEM)
    tranVIpheno <- bind_rows(tranVIpheno, tmp)
  }
  return(tranVIpheno)
}

transScandic <- 
  prepareData("transScandic",
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic/*.csv",
              "0.25")

transSpain <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain/*.csv",
              "0.25") %>% 
  filter(between(lat, 42.15, 43.1))


transects <- bind_rows(transScandic, transSpain)

axisscale1 <- 10
tscandic <- transects %>%
  filter(year %in% c(2019), trans == "transScandic") %>%
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric",
               values_to = "pheno") %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")),
         metric = factor(metric, levels = c("SOS", "EOS"))) %>%
  
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) + 
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.15, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale1, color = "DEM"), size = 0.5, alpha = 0.1, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale1, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(0, 400),
                     sec.axis = sec_axis(trans = ~.*axisscale1,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(VIname), ncol = 3) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

axisscale2 <- 20
tspain <- transects %>%
  filter(year %in% c(2019), trans == "transSpain") %>%
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric",
               values_to = "pheno") %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")),
         metric = factor(metric, levels = c("SOS", "EOS"))) %>%
  
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) +
  scale_shape_manual("Land cover type", values = seq(0, 6)) +
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale2, color = "DEM"), size = 0.5, alpha = 0.15, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale2, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(-100, 500), breaks = seq(0, 500, 100),
                     sec.axis = sec_axis(trans = ~.*axisscale2,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(VIname), ncol = 3) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

ggarrange(tscandic, tspain, common.legend = T,
          labels = c("a", "b"), ncol = 1)

ggsave("figures/figure_transects_raw.pdf",
       width = 25, height = 25, units = "cm")


