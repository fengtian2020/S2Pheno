
# libraries ----------------------------------------------------------------

library(sf)
library(tidyverse)
library(readr)
library(lubridate)


# data preparation --------------------------------------------------------

siteQA <- read_csv("data/0_ground_data_raw/phenocam/PhenoCam_sites.csv") %>% 
  filter(visual_check != "bad") %>% 
  mutate(ID = paste(site, veg_type, roi, sep = "_")) %>% 
  # filter(ID != "monteblanco_EB_1000", ID != "montenegro_SH_1000") %>%
  select(ID) %>% pull()

files <- list.files(path = "data/0_ground_data_raw/phenocam", recursive = T,
                    pattern = "*_3day.csv", 
                    full.names = T)
phenocam0 <- NULL
for (i in files) {
  ID <- str_remove(basename(i), "_3day.csv")
  skipn <- read_csv(i) %>% filter(str_detect(`#`, "^#")) %>% nrow
  tmp <- read_csv(i, skip = skipn + 1) %>% 
    select(date, GCC = gcc_75, #gcc_90, #midday_gcc, gcc_mean, 
           #smooth_gcc_75, # smooth_gcc_mean, smooth_gcc_50, smooth_gcc_90,
    ) %>% 
    mutate(ID = ID) %>% 
    filter(between(date, ymd("2017-4-1"), ymd("2020-3-31")))
  phenocam0 <- bind_rows(phenocam0, tmp)
}

phenocam0 <- phenocam0 %>% 
  drop_na() %>% 
  add_count(ID) %>% 
  filter(ID %in% siteQA,
         n > 120) %>% 
  select(-n)
sites <- phenocam0 %>% distinct(ID) %>% pull()
days <- seq.Date(ymd("2017-4-1"), ymd("2020-3-31"), 1)
phenocam <- tibble(ID = rep(sites, each = length(days)),
                   date = rep(days, length(sites))) %>% 
  left_join(phenocam0) %>% 
  group_by(ID) %>%
  mutate(index = row_number()) %>%
  group_modify(~{
    .x %>%
      mutate(GCC_loess = predict(loess(GCC ~ index, data = ., span = 0.2,
                                       family = "symmetric",
                                       na.action = na.exclude)))
  }) %>% ungroup() %>% 
  mutate(LC = str_split(ID, "_", simplify = T)[,2])


# Sentinel 2 at phenocam sites --------------------------------------------

level_key <- c("MX" = "Mixed forest",
               "SH" = "Mixed forest",
               "DN" = "Coniferous forest",
               "EN" = "Coniferous forest",
               "DB" = "Broad-leaved forest",
               "EB" = "Broad-leaved forest",
               "GR" = "Grassland",
               "WL" = "Wetland",
               "TN" = "Wetland",
               "AG" = "Agriculture")
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
phenocamS2 <- read_sf("data/1_gee_extracted_data/phenocam/phenocam_S2_NBAR_raw_201704_202003_10m_SCL45.shp") %>%
  mutate(date = ymd(date)) %>%
  as_tibble() %>%
  mutate(ID = paste(site, veg_type, str_pad(roi, 4, pad = "0"), sep = "_")) %>% 
  select(ID, date, SZAlocal, #SCL, 
         Red_raw = B4, NIR_raw = B8, 
         NBAR_Red_M = NBAR_Red, NBAR_NIR_M = NBAR_NIR,
         NBAR_Red_R = NBAR_Red_R, NBAR_NIR_R = NBAR_NIR_R) %>%
  mutate_at(vars(ends_with("raw")), function(x) x * 0.0001) %>% 
  mutate_at(vars(contains("Red"), contains("NIR")),
            ~replace(., . > 1, NA)) %>% 
  right_join(phenocam) %>% 
  mutate(LC = recode(LC, !!!level_key)) %>% 
  group_by(ID) %>% 
  mutate(index = row_number()) %>% 
  group_modify(~{
    .x %>%
      mutate(NBAR_NIR_M_loeFilter = loessFilter(., NBAR_NIR_M),
             NBAR_Red_M_loeFilter = loessFilter(., NBAR_Red_M),
             NIR_raw_loeFilter = loessFilter(., NIR_raw),
             Red_raw_loeFilter = loessFilter(., Red_raw)
      )
  }) %>% ungroup()

# calculate VIs -----------------------------------------------------------

calPPI <- function(DVI, MDVI, SZA) {
  PPI <- -0.25 * (1 + MDVI) / (1 - MDVI) * log((MDVI - DVI) / (MDVI - 0.09))
  dc <- pmin(0.0336 + 0.0477 / cos(SZA), 1)
  return(PPI / (0.5 / cos(SZA) * (1 - dc) + dc))
}
# toZ <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
VIspan = 0.2
phenocamVI <- phenocamS2 %>% 
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
  
save(phenocamVI, file = "data/RData/phenocamVI.RData")
write_csv(phenocamVI %>% 
            select(ID, date, GCC, GCC_loess,
                   EVI2_M_loeFiltered, EVI2_raw_loeFiltered, 
                   PPI_M_loeFiltered, PPI_raw_loeFiltered,
                   NDVI_M_loeFiltered, NDVI_raw_loeFiltered),
          "data/2_timesat_inputs/PhenoCam_sites_S2_VI_SCL45_10m.csv")





