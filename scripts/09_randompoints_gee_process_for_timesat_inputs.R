
library(sf)
library(tidyverse)
library(lubridate)

# read data and calculate VIs
data <- read_sf("data/1_gee_extracted_data/randompoints/randomPoints_S2_NBAR_raw_201704_202003_10m.shp") %>%
  mutate(date = ymd(date)) %>%
  as_tibble() %>%
  mutate(ID = as.character(geometry)) %>%
  select(ID, date, SZAlocal, SCL, 
         Red_raw = B4, NIR_raw = B8, 
         NBAR_Red_M = NBAR_Red, NBAR_NIR_M = NBAR_NIR) %>%
  filter(SCL %in% c(4, 5)) %>% select(-SCL) %>% 
  mutate_at(vars(ends_with("raw")), function(x) x * 0.0001) %>% 
  mutate_at(vars(contains("Red"), contains("NIR")),
            ~replace(., . > 1, NA))


IDlist = data %>% distinct(ID) %>% pull()
days = seq.Date(ymd("2017-4-1"), ymd("2020-3-31"), 1)
dataf <- tibble(ID = rep(IDlist, each = length(days)),
                date = rep(days, length(IDlist))) %>% 
  left_join(data) 

calPPI <- function(DVI, MDVI, SZA) {
  PPI <- -0.25 * (1 + MDVI) / (1 - MDVI) * log((MDVI - DVI) / (MDVI - 0.09))
  dc <- pmin(0.0336 + 0.0477 / cos(SZA), 1)
  return(PPI / (0.5 / cos(SZA) * (1 - dc) + dc))
}

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

randompointsVI <- dataf %>%
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
  mutate(index = row_number()) %>% 
  
  group_modify(~{
    .x %>%
      mutate(
        NDVI_M_loeFiltered = loessFilter(., NDVI_M),
        NDVI_raw_loeFiltered = loessFilter(., NDVI_raw),
        EVI2_M_loeFiltered = loessFilter(., EVI2_M),
        EVI2_raw_loeFiltered = loessFilter(., EVI2_raw),
        PPI_M_loeFiltered = loessFilter(., PPI_M),
        PPI_raw_loeFiltered = loessFilter(., PPI_raw),
      )
  }) %>% ungroup()

save(randompointsVI, file = "data/RData/randompointsVI.RData")

write_csv(randompointsVI %>% 
            select(ID, date,
                   EVI2_M_loeFiltered, EVI2_raw_loeFiltered,
                   PPI_M_loeFiltered, PPI_raw_loeFiltered,
                   NDVI_M_loeFiltered, NDVI_raw_loeFiltered),
          "data/2_timesat_inputs/Randomn_points_S2_VI_SCL45_10m.csv")




