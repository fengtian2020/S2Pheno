library(sf)
library(tidyverse)
library(lubridate)
library(ggpubr)

# prepare pep data --------------------------------------------------------

pep <- read_csv("data/0_ground_data_raw/pep725/pep725_smets0.csv") 
pep %>% filter(year %in% c(2018, 2019),
               phase_id %in% c(10, 11, 95, 100, 151, 182, 205, 223)) %>% distinct(genus) %>% arrange(genus) %>% 
  # add_count(phase_id) %>% 
  # distinct(phase_id, .keep_all = T) %>% 
  # select(phase_id, n) %>% 
  # arrange(phase_id) %>% 
  print(n = 100)

level_key <- c("Acer" = "Broad-leaved-tree",
               "Actinidia" = "Fruit-tree",
               "Aesculus" = "Broad-leaved-tree",
               "Alnus" = "Broad-leaved-tree",
               "Avena" = "Crop",
               "Beta" = "Crop",
               "Betula" = "Broad-leaved-tree",
               "Brassica" = "Crop",
               "Castanea" = "Fruit-tree",
               "Corylus" = "Fruit-tree",
               "Fagus" = "Broad-leaved-tree",
               "Fraxinus" = "Broad-leaved-tree",
               "Helianthus" = "Crop",
               "Hordeum" = "Crop",
               "Juglans" = "Fruit-tree",
               "Larix" = "Deciduous Coniferous-tree",
               "Malus" = "Fruit-tree",
               "meadow" = "Meadow",
               "Picea" = "Coniferous-tree",
               "Pinus" = "Coniferous-tree",
               "Populus" = "Broad-leaved-tree",
               "Prunus" = "Fruit-tree",
               "Punica" = "Fruit-tree",
               "Pyrus" = "Fruit-tree",
               "Quercus" = "Broad-leaved-tree",
               "Ribes" = "Fruit-tree",
               "Robinia" = "Broad-leaved-tree",
               "Rosmarinus" = "other",
               "Salix" = "Broad-leaved-tree",
               "Sambucus" = "Broad-leaved-tree",
               "Secale" = "Crop",
               "Solanum" = "Crop",
               "Sorbus" = "Broad-leaved-tree",
               "Syringa" = "Broad-leaved-tree",
               "Tilia" = "Broad-leaved-tree",
               "Triticum" = "Crop",
               "Vaccinium" = "Fruit-tree",
               "Vitis" = "Vineyards",
               "Zea" = "Crop")
pepi <- pep %>% 
  filter(year %in% c(2018, 2019),
         phase_id %in% c(10, 11, 95, 100, 151, 182, 205, 223)) %>% 
  select(s_id, provider_id, gss_id, 
         lon, lat, alt, alt_dem,
         genus, pheno = date, phase_id, cult_season) %>% 
  mutate(type = recode(genus, !!!level_key)) %>% 
  filter(type %in% c("Meadow", "Broad-leaved-tree", "Coniferous-tree")) 

pepi %>% add_count(genus) %>% distinct(genus, n) %>% arrange(genus)

write_csv(pepi, file = "data/0_ground_data_raw/pep725/pep725_in_use.csv")

# read data and calculate VIs ---------------------------------------------

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

##
calPPI <- function(DVI, MDVI, SZA) {
  PPI <- -0.25 * (1 + MDVI) / (1 - MDVI) * log((MDVI - DVI) / (MDVI - 0.09))
  dc <- pmin(0.0336 + 0.0477 / cos(SZA), 1)
  return(PPI / (0.5 / cos(SZA) * (1 - dc) + dc))
}
# toZ <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
VIspan = 0.2

# read data and calculate VIs
files <- list.files("data/1_gee_extracted_data/pep725", pattern = "*.shp", full.names = T)
formatpep <- function(tmp) {
  tmp %>% 
    mutate(date = ymd(S2date)) %>%
    as_tibble() %>%
    mutate(check = rowSums(across(c(B4, B8, NBAR_Red, NBAR_NIR)), na.rm = T)) %>% 
    filter(check != 0) %>% 
    select(s_id, date, LC, Red_raw = B4, NIR_raw = B8, 
           NBAR_Red_M = NBAR_Red, NBAR_NIR_M = NBAR_NIR, SZAlocal) %>%
    mutate_at(vars(ends_with("raw")), function(x) x * 0.0001) %>% 
    mutate_at(vars(contains("Red"), contains("NIR")),
              ~replace(., . > 1, NA))
}
pepPheno <- bind_rows(read_sf(files[1]) %>% mutate(LC = "Broad-leaved-tree") %>% formatpep()) %>%
  bind_rows(read_sf(files[2]) %>% mutate(LC = "Coniferous-tree") %>% formatpep()) %>%
  bind_rows(read_sf(files[3]) %>% mutate(LC = "Meadow") %>% formatpep())

pepPheno0 <- pepPheno %>% 
  mutate(ID = paste(LC, s_id, sep = "_")) %>% 
  select(ID, everything(), -c(s_id, LC))

sites <- pepPheno0 %>% distinct(ID) %>% pull()
days <- seq.Date(ymd("2017-4-1"), ymd("2020-3-31"), 1)
pepVI <- tibble(ID = rep(sites, each = length(days)),
                date = rep(days, length(sites))) %>% 
  left_join(pepPheno0) %>%  
  
  mutate(DVI_M = NBAR_NIR_M - NBAR_Red_M,
         DVI_raw = NIR_raw - Red_raw) %>% 
  mutate_at(vars(starts_with("DVI")), ~replace(., . < 0.01, NA)) %>% 
  mutate(NDVI_M = DVI_M / (NBAR_NIR_M + NBAR_Red_M),
         NDVI_raw = DVI_raw / (NIR_raw + Red_raw)) %>% 
  mutate(EVI2_M = 2.5 * DVI_M / (NBAR_NIR_M + 2.4 * NBAR_Red_M + 1),
         EVI2_raw = 2.5 * DVI_raw / (NIR_raw + 2.4 * Red_raw + 1)) %>% 
  group_by(ID) %>% 
  mutate(index = row_number()) %>% 
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

save(pepVI, file = "data/RData/pepVI4kmMean.RData")
write_csv(pepVI %>% 
            select(ID, date, 
                   EVI2_M_loeFiltered, EVI2_raw_loeFiltered,
                   PPI_M_loeFiltered, PPI_raw_loeFiltered,
                   NDVI_M_loeFiltered, NDVI_raw_loeFiltered),
          "data/2_timesat_inputs/PEP_S2_VI_SCL45_LC_4km_mean.csv")



