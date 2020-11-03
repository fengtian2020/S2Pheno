
# libraries ---------------------------------------------------------------

library(sf)
library(tidyverse)
library(readr)
library(lubridate)


# read data and calculate VIs ---------------------------------------------

fluxS2 <- read_sf("data/1_gee_extracted_data/flux/flux_drought2018_site_S2_NBAR_raw_201704_202003_100m_SCL45_mean.shp") %>%
  mutate(date = ymd(date)) %>%
  as_tibble() %>%
  select(ID, date, SZAlocal, #SCL, 
         Red_raw = B4, NIR_raw = B8, 
         NBAR_Red_M = NBAR_Red, NBAR_NIR_M = NBAR_NIR) %>%
  mutate_at(vars(ends_with("raw")), function(x) x * 0.0001) %>% 
  mutate_at(vars(contains("Red"), contains("NIR")),
            ~replace(., . > 1, NA))


# Land cover from CORINE 2018 ---------------------------------------------

level_key <- c(`313` = "Mixed forest",
               `324` = "Mixed forest",
               `312` = "Coniferous forest",
               `311` = "Broad-leaved forest",
               `231` = "Grassland",
               `321` = "Grassland",
               `412` = "Wetland",
               `411` = "Wetland",
               `211` = "Agriculture",
               `212` = "Agriculture",
               `221` = "Agriculture",
               `243` = "Agriculture",
               `244` = "Agriculture")
fluxLC <- read_sf("data/1_gee_raw_reflectance/flux/flux_drought2018_site_CORINE_LC_100m.shp") %>%
  as_tibble() %>% filter(year == 2018, !is.na(b1)) %>% 
  mutate(LC = recode(b1, !!!level_key)) %>% 
  select(ID, LC)
# lct <- fluxLC %>% group_by(LC) %>% mutate(N = n()) %>% distinct(LC, .keep_all = T)


# GPP data 2018 -----------------------------------------------------------

# from hh to day
# load("RData/fluxhh_drought2018_6_March_2020.RData")
# fluxGPP <- fluxhh %>%
#   mutate(ID = str_replace_all(ID,  "CZ-wet", "CZ-Wet")) %>% 
#   filter(year(time) %in% c(2017, 2018), ID != "IT-Cp2",
#          !(ID == "SE-Svb" & year(time) == 2017)) %>% 
#   mutate(date = date(time)) %>% 
#   select(ID, date, GPP = tower_GPP) %>% 
#   group_by(ID, date) %>% 
#   summarise(GPP = mean(GPP, na.rm = T)) %>% 
#   ungroup()

fluxGPP <- read_csv("data/0_ground_data_raw/flux/flux_GPP_Drought2018_daily_year2017_2018.csv")


# generate consistent daily time series -----------------------------------

loessFilter <- function(tbl, col, span = 0.2, thres = 1.5) {
  col <- enquo(col)
  tbl <- tbl %>% mutate(colv = !!col)
  if (length(which(!is.na(tbl$colv) == 1)) > 19) {
    tbl %>% 
      mutate(losp = predict(loess(colv ~ index, data = ., span = span,
                                  family = "symmetric",
                                  na.action = na.exclude))) %>%
      mutate(res = abs(colv - losp)) %>% 
      mutate_at(vars(colv), 
                ~replace(., res > sd(res, na.rm = T) * thres, NA)) %>% 
      select(colv) %>% pull()
  } else {
    tbl %>% select(colv) %>% pull()
  }
}


IDlist = fluxGPP %>% distinct(ID) %>% pull()
days = seq.Date(ymd("2017-4-1"), ymd("2020-3-31"), 1)
flux <- tibble(ID = rep(IDlist, each = length(days)),
               date = rep(days, length(IDlist))) %>% 
  left_join(fluxS2) %>% 
  left_join(fluxLC) %>% 
  left_join(fluxGPP) %>% 

  ## filtering and smoothing the GPP and reflectance time series
  group_by(ID) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(GPP_loess = predict(loess(GPP ~ index, data = ., span = 0.2, 
                                       family = "symmetric",
                                       na.action = na.exclude)),
             NBAR_NIR_M_loeFilter = loessFilter(., NBAR_NIR_M),
             NBAR_Red_M_loeFilter = loessFilter(., NBAR_Red_M),
             NIR_raw_loeFilter = loessFilter(., NIR_raw),
             Red_raw_loeFilter = loessFilter(., Red_raw),
             )
  }) %>% ungroup() %>% 
  filter(!is.na(LC))



# calculate VIs based on raw and smoothed reflectance ---------------------

calPPI <- function(DVI, MDVI, SZA) {
  PPI <- -0.25 * (1 + MDVI) / (1 - MDVI) * log((MDVI - DVI) / (MDVI - 0.09))
  dc <- pmin(0.0336 + 0.0477 / cos(SZA), 1)
  return(PPI / (0.5 / cos(SZA) * (1 - dc) + dc))
}
# toZ <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
VIspan = 0.2
fluxVI <- flux %>% 
  mutate(DVI_M = NBAR_NIR_M - NBAR_Red_M,
         DVI_raw = NIR_raw - Red_raw) %>% 
  mutate_at(vars(starts_with("DVI")), ~replace(., . < 0.01, NA)) %>% 
  mutate(NDVI_M = DVI_M / (NBAR_NIR_M + NBAR_Red_M),
         NDVI_raw = DVI_raw / (NIR_raw + Red_raw)) %>% 
  mutate(EVI2_M = 2.5 * DVI_M / (NBAR_NIR_M + 2.4 * NBAR_Red_M + 1),
         EVI2_raw = 2.5 * DVI_raw / (NIR_raw + 2.4 * Red_raw + 1)) %>% 
  group_by(ID) %>% 
  mutate(MDVI_M = max(quantile(DVI_M, 1, na.rm = T) + 0.005, 0.18),
         MDVI_raw = max(quantile(DVI_raw, 1, na.rm = T) + 0.005, 0.18)) %>% 
  mutate(PPI_M = calPPI(DVI_M, MDVI_M, SZAlocal),
         PPI_raw = calPPI(DVI_raw, MDVI_raw, SZAlocal)) %>% 
  
  group_modify(~{
    .x %>%
      mutate(NDVI_M_loeFiltered = loessFilter(., NDVI_M),
             NDVI_raw_loeFiltered = loessFilter(., NDVI_raw),
             EVI2_M_loeFiltered = loessFilter(., EVI2_M),
             EVI2_raw_loeFiltered = loessFilter(., EVI2_raw),
             PPI_M_loeFiltered = loessFilter(., PPI_M),
             PPI_raw_loeFiltered = loessFilter(., PPI_raw))
  }) %>% ungroup() 
  
save(fluxVI, file = "data/RData/fluxVI.RData")
write_csv(fluxVI %>% 
            select(ID, date, GPP, GPP_loess,
                   EVI2_M_loeFiltered, EVI2_raw_loeFiltered, EVI2_M, EVI2_raw, 
                   PPI_M_loeFiltered, PPI_raw_loeFiltered, PPI_M, PPI_raw, 
                   NDVI_M_loeFiltered, NDVI_raw_loeFiltered,NDVI_M, NDVI_raw),
          "data/2_timesat_inputs/Flux_sites_S2_VI_SCL45_100m_mean.csv")




