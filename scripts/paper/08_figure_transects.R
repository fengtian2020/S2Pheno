library(sf)
library(tidyverse)
library(lubridate)
library(ggpubr)


# NBAR figure -------------------------------------------------------------

prepareData <- function(transID, dem, timesat) {
  transDEM <- read_sf(dem) %>%
    as_tibble() %>%
    mutate(ID = as.character(geometry)) %>%
    select(ID, Elevation = AVE_DSM)
  
  files <- Sys.glob(timesat) %>% 
    # str_subset(thres) %>% 
    str_subset("M_loeFiltered")
  
  tranVIpheno <- NULL
  for (i in files) {
    VIname <- str_extract(i, "NDVI|EVI2|PPI")
    thres <- str_extract(i, "\\d\\.\\d{2}")
    tmp <- read_csv(i) %>%
      transmute(ID, year,
                SOS = start_of_season,
                EOS = end_of_season,
                trans =  transID,
                VIname = VIname,
                thres = thres) %>% 
      filter(between(year, 2019, 2019)) %>% 
      mutate(lat = as.numeric(str_extract(ID, "[:blank:]\\d{2}\\.\\d{1,}"))) %>% 
      left_join(transDEM)
    tranVIpheno <- bind_rows(tranVIpheno, tmp)
  }
  return(tranVIpheno)
}

loessFilter <- function(tbl, col, span = 0.1) {
  col <- enquo(col)
  tbl <- tbl %>% mutate(colv = !!col)
  tbl %>% 
    mutate(losp = predict(loess(colv ~ index, data = ., span = span,
                                family = "symmetric",
                                na.action = na.exclude))) %>%
    pull(losp)
}

loessFilterSE <- function(tbl, col, span = 0.1) {
  col <- enquo(col)
  tbl <- tbl %>% mutate(colv = !!col)
  tbl %>% 
    mutate(lospse = predict(loess(colv ~ index, data = ., span = span,
                                  family = "symmetric",
                                  na.action = na.exclude),
                            se = TRUE)$se.fit) %>% 
    pull(lospse)
}

transScandicSOS <- 
  prepareData("transScandic", 
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic_all/*.csv") %>% 
  select(-c(EOS, ID, year, trans)) %>% 
  filter((VIname == "NDVI" & thres == "0.40") |
         (VIname == "EVI2" & thres == "0.40") |
         (VIname == "PPI" & thres == "0.25")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
    mutate(index = row_number()) %>% 
    group_modify(~{
      .x %>%
        mutate(SOS_loess = loessFilter(., SOS),
               SOS_loess_se = loessFilterSE(., SOS),
               elev_loess = loessFilter(., Elevation),
               # elev_loess_se = loessFilterSE(., Elevation)
               )
    }) %>% ungroup()
elevScandic <- transScandicSOS %>% group_by(lat) %>% 
  summarise(across(elev_loess, mean))

transScandicEOS <- 
  prepareData("transScandic", 
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic_all/*.csv") %>% 
  select(-SOS) %>% 
  filter((VIname == "NDVI" & thres == "0.50") |
           (VIname == "EVI2" & thres == "0.45") |
           (VIname == "PPI" & thres == "0.15")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(EOS_loess = loessFilter(., EOS),
             EOS_loess_se = loessFilterSE(., EOS))
  }) %>% ungroup()


transSpainSOS <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain_all/*.csv") %>% 
  filter(between(lat, 42.15, 43.1)) %>% 
  select(-EOS) %>% 
  filter((VIname == "NDVI" & thres == "0.40") |
           (VIname == "EVI2" & thres == "0.40") |
           (VIname == "PPI" & thres == "0.25")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(SOS_loess = loessFilter(., SOS),
             SOS_loess_se = loessFilterSE(., SOS),
             elev_loess = loessFilter(., Elevation),
             # elev_loess_se = loessFilterSE(., Elevation)
      )
  }) %>% ungroup()
elevSpain <- transSpainSOS %>% group_by(lat) %>% 
  summarise(across(elev_loess, mean))

transSpainEOS <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain_all/*.csv") %>% 
  filter(between(lat, 42.15, 43.1)) %>% 
  select(-SOS) %>% 
  filter((VIname == "NDVI" & thres == "0.50") |
           (VIname == "EVI2" & thres == "0.45") |
           (VIname == "PPI" & thres == "0.15")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(EOS_loess = loessFilter(., EOS),
             EOS_loess_se = loessFilterSE(., EOS))
  }) %>% ungroup()





axisscale1 <- 30 
axisbias1 <- 120
scandicsos <- transScandicSOS %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = SOS_loess - SOS_loess_se, 
                  ymax = SOS_loess + SOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = SOS_loess, color = VIname)) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  geom_line(data = elevScandic, aes(y = elev_loess/axisscale1 + axisbias1, linetype = "Elevation")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias1) * axisscale1,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 2000, 500),
                                         labels = seq(0, 2, 0.5))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "SOS - northern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

axisscale2 <- 30 
axisbias2 <- 240
scandiceos <- transScandicEOS %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = EOS_loess - EOS_loess_se, 
                  ymax = EOS_loess + EOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = EOS_loess, color = VIname)) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  geom_line(data = elevScandic, aes(y = elev_loess/axisscale2 + axisbias2, linetype = "Elevation")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias2) * axisscale2,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 2000, 500),
                                         labels = seq(0, 2, 0.5))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "EOS - northern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

axisscale3 <- 20 
axisbias3 <- 0
spainsos <- transSpainSOS %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = SOS_loess - SOS_loess_se, 
                  ymax = SOS_loess + SOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = SOS_loess, color = VIname)) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  geom_line(data = elevSpain, aes(y = elev_loess/axisscale3 + axisbias3, linetype = "Elevation")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias3) * axisscale3,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 5000, 1000),
                                         labels = seq(0, 5, 1))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "SOS - southern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

axisscale4 <- 20 
axisbias4 <- 150
spaineos <- transSpainEOS %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = EOS_loess - EOS_loess_se, 
                  ymax = EOS_loess + EOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = EOS_loess, color = VIname)) +
  geom_line(data = elevSpain, aes(y = elev_loess/axisscale4 + axisbias4, linetype = "Elevation")) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias4) * axisscale4,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 5000, 1000),
                                         labels = seq(0, 5, 1))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "EOS - southern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

ggarrange(scandicsos, spainsos, scandiceos, spaineos, 
          common.legend = T,
          labels = c("a", "c", "b", "d"),
          nrow = 2, ncol = 2,
          font.label = list(size = 19))

ggsave("figures/figure_transects_NBAR.pdf",
       width = 35, height = 25, units = "cm")


# TOC figure --------------------------------------------------------------

prepareData <- function(transID, dem, timesat) {
  transDEM <- read_sf(dem) %>%
    as_tibble() %>%
    mutate(ID = as.character(geometry)) %>%
    select(ID, Elevation = AVE_DSM)
  
  files <- Sys.glob(timesat) %>% 
    # str_subset(thres) %>% 
    str_subset("raw_loeFiltered")
  
  tranVIpheno <- NULL
  for (i in files) {
    VIname <- str_extract(i, "NDVI|EVI2|PPI")
    thres <- str_extract(i, "\\d\\.\\d{2}")
    tmp <- read_csv(i) %>%
      transmute(ID, year,
                SOS = start_of_season,
                EOS = end_of_season,
                trans =  transID,
                VIname = VIname,
                thres = thres) %>% 
      filter(between(year, 2019, 2019)) %>% 
      mutate(lat = as.numeric(str_extract(ID, "[:blank:]\\d{2}\\.\\d{1,}"))) %>% 
      left_join(transDEM)
    tranVIpheno <- bind_rows(tranVIpheno, tmp)
  }
  return(tranVIpheno)
}

loessFilter <- function(tbl, col, span = 0.1) {
  col <- enquo(col)
  tbl <- tbl %>% mutate(colv = !!col)
  tbl %>% 
    mutate(losp = predict(loess(colv ~ index, data = ., span = span,
                                family = "symmetric",
                                na.action = na.exclude))) %>%
    pull(losp)
}

loessFilterSE <- function(tbl, col, span = 0.1) {
  col <- enquo(col)
  tbl <- tbl %>% mutate(colv = !!col)
  tbl %>% 
    mutate(lospse = predict(loess(colv ~ index, data = ., span = span,
                                  family = "symmetric",
                                  na.action = na.exclude),
                            se = TRUE)$se.fit) %>% 
    pull(lospse)
}

transScandicSOS <- 
  prepareData("transScandic", 
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic_all/*.csv") %>% 
  select(-c(EOS, ID, year, trans)) %>% 
  filter((VIname == "NDVI" & thres == "0.40") |
           (VIname == "EVI2" & thres == "0.40") |
           (VIname == "PPI" & thres == "0.25")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(SOS_loess = loessFilter(., SOS),
             SOS_loess_se = loessFilterSE(., SOS),
             elev_loess = loessFilter(., Elevation),
             # elev_loess_se = loessFilterSE(., Elevation)
      )
  }) %>% ungroup()
elevScandic <- transScandicSOS %>% group_by(lat) %>% 
  summarise(across(elev_loess, mean))

transScandicEOS <- 
  prepareData("transScandic", 
              "data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/scandic_all/*.csv") %>% 
  select(-SOS) %>% 
  filter((VIname == "NDVI" & thres == "0.50") |
           (VIname == "EVI2" & thres == "0.45") |
           (VIname == "PPI" & thres == "0.15")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(EOS_loess = loessFilter(., EOS),
             EOS_loess_se = loessFilterSE(., EOS))
  }) %>% ungroup()


transSpainSOS <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain_all/*.csv") %>% 
  filter(between(lat, 42.15, 43.1)) %>% 
  select(-EOS) %>% 
  filter((VIname == "NDVI" & thres == "0.40") |
           (VIname == "EVI2" & thres == "0.40") |
           (VIname == "PPI" & thres == "0.25")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(SOS_loess = loessFilter(., SOS),
             SOS_loess_se = loessFilterSE(., SOS),
             elev_loess = loessFilter(., Elevation),
             # elev_loess_se = loessFilterSE(., Elevation)
      )
  }) %>% ungroup()
elevSpain <- transSpainSOS %>% group_by(lat) %>% 
  summarise(across(elev_loess, mean))

transSpainEOS <- 
  prepareData("transSpain",
              "data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp",
              "data/3_timesat_outputs/phenometrics/transects/spain_all/*.csv") %>% 
  filter(between(lat, 42.15, 43.1)) %>% 
  select(-SOS) %>% 
  filter((VIname == "NDVI" & thres == "0.50") |
           (VIname == "EVI2" & thres == "0.45") |
           (VIname == "PPI" & thres == "0.15")) %>% 
  group_by(VIname) %>% arrange(lat) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(EOS_loess = loessFilter(., EOS),
             EOS_loess_se = loessFilterSE(., EOS))
  }) %>% ungroup()





axisscale1 <- 30 
axisbias1 <- 120
scandicsos <- transScandicSOS %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = SOS_loess - SOS_loess_se, 
                  ymax = SOS_loess + SOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = SOS_loess, color = VIname)) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  geom_line(data = elevScandic, aes(y = elev_loess/axisscale1 + axisbias1, linetype = "Elevation")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias1) * axisscale1,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 2000, 500),
                                         labels = seq(0, 2, 0.5))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "SOS - northern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

axisscale2 <- 30 
axisbias2 <- 240
scandiceos <- transScandicEOS %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = EOS_loess - EOS_loess_se, 
                  ymax = EOS_loess + EOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = EOS_loess, color = VIname)) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  geom_line(data = elevScandic, aes(y = elev_loess/axisscale2 + axisbias2, linetype = "Elevation")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias2) * axisscale2,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 2000, 500),
                                         labels = seq(0, 2, 0.5))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "EOS - northern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

axisscale3 <- 20 
axisbias3 <- 0
spainsos <- transSpainSOS %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = SOS_loess - SOS_loess_se, 
                  ymax = SOS_loess + SOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = SOS_loess, color = VIname)) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  geom_line(data = elevSpain, aes(y = elev_loess/axisscale3 + axisbias3, linetype = "Elevation")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias3) * axisscale3,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 5000, 1000),
                                         labels = seq(0, 5, 1))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "SOS - southern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

axisscale4 <- 20 
axisbias4 <- 150
spaineos <- transSpainEOS %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>%
  ggplot(aes(x = lat))+
  geom_ribbon(aes(ymin = EOS_loess - EOS_loess_se, 
                  ymax = EOS_loess + EOS_loess_se,
                  fill = VIname), 
              alpha = 0.2) +
  geom_line(aes(y = EOS_loess, color = VIname)) +
  geom_line(data = elevSpain, aes(y = elev_loess/axisscale4 + axisbias4, linetype = "Elevation")) +
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  scale_fill_manual(values = c("#FE6100", "blue", "#DC267F")) +
  scale_y_continuous(name = "Day of year", 
                     sec.axis = sec_axis(trans = ~ (. - axisbias4) * axisscale4,
                                         name = "Elevation (km)",
                                         breaks = seq(0, 5000, 1000),
                                         labels = seq(0, 5, 1))) +
  theme_bw(base_size = 15) +
  labs(x = "Latitude", title = "EOS - southern transect") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2))

ggarrange(scandicsos, spainsos, scandiceos, spaineos, 
          common.legend = T,
          labels = c("a", "c", "b", "d"),
          nrow = 2, ncol = 2,
          font.label = list(size = 19))

ggsave("figures/figure_transects_raw.pdf",
       width = 35, height = 25, units = "cm")


