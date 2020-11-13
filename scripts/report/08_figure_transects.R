library(sf)
library(tidyverse)
library(lubridate)
library(ggpubr)


# TOC figure --------------------------------------------------------------

prepareData <- function(transID, dem, timesat) {
  transDEM <- read_sf(dem) %>%
    as_tibble() %>%
    mutate(ID = as.character(geometry)) %>%
    select(ID, DEM = AVE_DSM)
  
  tranVIpheno <- read_csv(timesat) %>% 
      transmute(ID, year,
                SOS = start_of_season,
                EOS = end_of_season,
                trans =  transID) %>% 
      filter(between(year, 2018, 2019)) %>% 
      mutate(lat = as.numeric(str_extract(ID, "[:blank:]\\d{2}\\.\\d{1,}"))) %>% 
      left_join(transDEM)

  return(tranVIpheno)
}

transScandicSOS25 <- 
  prepareData("transScandic",
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic/sos_eos_0.250_PPI_raw_loeFiltered_Transect_Scandic_S2_VI_SCL45_10m_newDL.csv") %>% 
  select(-EOS)

transScandicEOS15 <- 
  prepareData("transScandic",
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic/sos_eos_0.150_PPI_raw_loeFiltered_Transect_Scandic_S2_VI_SCL45_10m_newDL.csv") %>% 
  select(-SOS)
transScandic <- left_join(transScandicSOS25, transScandicEOS15)

##########################
transSpainSOS25 <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain/sos_eos_0.250_PPI_raw_loeFiltered_Transect_Spain_S2_VI_SCL45_10m_newDL.csv") %>% 
  select(-EOS) %>% 
  filter(between(lat, 42.15, 43.1))

transSpainEOS15 <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain/sos_eos_0.150_PPI_raw_loeFiltered_Transect_Spain_S2_VI_SCL45_10m_newDL.csv") %>% 
  select(-SOS) %>% 
  filter(between(lat, 42.15, 43.1))

transSpain <- left_join(transSpainSOS25, transSpainEOS15)


transects <- bind_rows(transScandic, transSpain)

axisscale1 <- 10
transects %>%
  filter(year %in% c(2019)) %>%
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric",
               values_to = "pheno") %>%
  mutate(metric = factor(metric, levels = c("SOS", "EOS")),
         trans = factor(trans, levels = c("transScandic", "transSpain"),
                        labels = c("Northern Europe", "Southern Europe"))) %>%
  
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
  facet_wrap(vars(trans), scales = "free") +
  theme_bw(base_size = 15) +
  labs(x = "Latitude") +
  theme(legend.position = "top")

ggsave("figures/report/figure_transects_raw.pdf",
       width = 25, height = 15, units = "cm")


