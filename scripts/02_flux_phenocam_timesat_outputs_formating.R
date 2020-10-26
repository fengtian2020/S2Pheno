library(sf)
library(tidyverse)
library(readr)
library(lubridate)


# FLUX --------------------------------------------------------------------

load("data/RData/fluxVI.RData")
# normalize to z-scores for visualization
forscale <- fluxVI %>%
  filter(year(date) == 2018) %>% 
  group_by(ID) %>%
  summarise(across(c(GPP_loess, contains("_raw_loeFiltered")),
                   list(mean = mean, sd = sd), na.rm = T, 
                   .names = "{col}_{fn}")) %>% 
  rename_with(~ paste0(str_extract(., "^[A-Z2]{3,4}"), 
                       str_extract(., "[a-z]{2,4}$")), -ID)


###  TIMESAT smoothing time series

fluxVItimesat <- read_csv("data/3_timesat_outputs/smooth_time_series/smoothed_Flux_sites_S2_VI_SCL45_100m_mean_newDL.csv") %>% 
  select(ID, date, 
         EVI2_M_loeFiltered_fit = EVI2_M_loeFiltered,
         EVI2_raw_loeFiltered_fit = EVI2_raw_loeFiltered,
         PPI_M_loeFiltered_fit = PPI_M_loeFiltered,
         PPI_raw_loeFiltered_fit = PPI_raw_loeFiltered,
         NDVI_M_loeFiltered_fit = NDVI_M_loeFiltered,
         NDVI_raw_loeFiltered_fit = NDVI_raw_loeFiltered,
  ) %>%
  right_join(fluxVI) %>%
  left_join(forscale) %>%
  mutate(LC = fct_relevel(factor(LC), "Broad-leaved forest",
                          "Coniferous forest", "Mixed forest",
                          "Agriculture", "Grassland", "Wetland")) %>%
  group_by(ID) %>% 
  mutate(across(c(GPP, GPP_loess), 
                ~ (.x - GPPmean)/GPPsd, .names = "{col}_Z")) %>% 
  mutate(across(c(contains("NDVI_raw"), contains("NDVI_M")),
                ~ (.x - NDVImean)/NDVIsd, .names = "{col}_Z")) %>% 
  mutate(across(c(contains("PPI_raw"), contains("PPI_M")),
                ~ (.x - PPImean)/PPIsd, .names = "{col}_Z")) %>%
  mutate(across(c(contains("EVI2_raw"), contains("EVI2_M")),
                ~ (.x - EVI2mean)/EVI2sd, .names = "{col}_Z")) %>% 
  ungroup() 


### TIMESAT phenometrics thresholds 

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
fluxLC <- read_sf("data/1_gee_extracted_data/flux/flux_drought2018_site_CORINE_LC_100m.shp") %>%
  as_tibble() %>% filter(year == 2018, !is.na(b1)) %>% 
  mutate(LC = recode(b1, !!!level_key)) %>% 
  select(ID, LC)

### amplitude thresholds
files <- Sys.glob("data/3_timesat_outputs/phenometrics/flux/*.csv") %>%
  str_subset("(GPP_loess|PPI_M_loeFiltered|EVI2_M_loeFiltered|NDVI_M_loeFiltered|PPI_raw_loeFiltered|EVI2_raw_loeFiltered|NDVI_raw_loeFiltered)")

fluxVIpheno <- NULL
for (i in files) {
  tmp <- str_extract(i, "\\d.\\d{3}_(\\w{3,4})_(\\w{1,3})")
  thres <- str_split(tmp, "_", simplify = TRUE)[1] 
  VIname <- str_split(tmp, "_", simplify = TRUE)[2]
  BRDF <- str_split(tmp, "_", simplify = TRUE)[3]
  tmp <- read_csv(i) %>% 
    transmute(ID, year,
              SOS = start_of_season,
              EOS = end_of_season,
              thres = thres, 
              BRDF = BRDF,
              VIname = VIname) %>% 
    left_join(fluxLC)
  if (tmp$BRDF[1] == "loe") {
    tmp <- mutate(tmp, BRDF = "raw") %>% bind_rows(mutate(tmp, BRDF = "M"))
  }
  fluxVIpheno <- bind_rows(fluxVIpheno, tmp)
}



# Phenocam ----------------------------------------------------------------

load("data/RData/phenocamVI.RData")

# normalize to z-scores for visualization
forscale <- phenocamVI %>%
  filter(year(date) == 2018) %>% 
  group_by(ID) %>%
  summarise(across(c(GCC_loess, contains("_raw_loeFiltered")),
                   list(mean = mean, sd = sd), na.rm = T, 
                   .names = "{col}_{fn}")) %>% 
  rename_with(~ paste0(str_extract(., "^[A-Z2]{3,4}"), 
                       str_extract(., "[a-z]{2,4}$")), -ID)


### TIMESAT fitting results 

phenocamVItimesat <- read_csv("data/3_timesat_outputs/smooth_time_series/smoothed_PhenoCam_sites_S2_VI_SCL45_10m_newDL.csv") %>% 
  select(ID, date, GCC_loess, 
         EVI2_M_loeFiltered_fit = EVI2_M_loeFiltered,
         EVI2_raw_loeFiltered_fit = EVI2_raw_loeFiltered,
         PPI_M_loeFiltered_fit = PPI_M_loeFiltered,
         PPI_raw_loeFiltered_fit = PPI_raw_loeFiltered,
         NDVI_M_loeFiltered_fit = NDVI_M_loeFiltered,
         NDVI_raw_loeFiltered_fit = NDVI_raw_loeFiltered) %>%
  right_join(phenocamVI %>% select(-GCC_loess)) %>% 
  left_join(forscale) %>% 
  mutate(LC = fct_relevel(factor(LC), "Broad-leaved forest",
                          "Coniferous forest", "Mixed forest",
                          "Agriculture", "Grassland", "Wetland")) %>%
  group_by(ID) %>% 
  mutate(across(c(GCC, GCC_loess), 
                ~ (.x - GCCmean)/GCCsd, .names = "{col}_Z")) %>% 
  mutate(across(c(contains("NDVI_raw"), contains("NDVI_M")),
                ~ (.x - NDVImean)/NDVIsd, .names = "{col}_Z")) %>% 
  mutate(across(c(contains("PPI_raw"), contains("PPI_M")),
                ~ (.x - PPImean)/PPIsd, .names = "{col}_Z")) %>%
  mutate(across(c(contains("EVI2_raw"), contains("EVI2_M")),
                ~ (.x - EVI2mean)/EVI2sd, .names = "{col}_Z")) %>% 
  ungroup() 


#### TIMESAT phenometrics thresholds

### amplitude thresholds
files <- Sys.glob("data/3_timesat_outputs/phenometrics/phenocam/*.csv") %>%
  str_subset("(GCC_loess|PPI_M_loeFiltered|EVI2_M_loeFiltered|NDVI_M_loeFiltered|PPI_raw_loeFiltered|EVI2_raw_loeFiltered|NDVI_raw_loeFiltered)")

phenocamVIpheno <- NULL
for (i in files) {
  tmp <- str_extract(i, "\\d.\\d{3}_(\\w{3,4})_(\\w{1,3})")
  thres <- str_split(tmp, "_", simplify = TRUE)[1]
  VIname <- str_split(tmp, "_", simplify = TRUE)[2]
  BRDF <- str_split(tmp, "_", simplify = TRUE)[3]
  tmp <- read_csv(i) %>% 
    transmute(ID, year,
              SOS = start_of_season,
              EOS = end_of_season,
              thres = thres, 
              BRDF = BRDF,
              VIname = VIname) %>% 
    mutate(LC = str_split(ID, "_", simplify = T)[,2])
  if (tmp$BRDF[1] == "loe") {
    tmp <- mutate(tmp, BRDF = "raw") %>% 
      bind_rows(mutate(tmp, BRDF = "M"))
  }
  phenocamVIpheno <- bind_rows(phenocamVIpheno, tmp)
}

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

phenocamVIpheno <- phenocamVIpheno %>% 
       mutate(LC = recode(LC, !!!level_key))


# FLUX and phenocam combination -------------------------------------------

## TIMESAT smoothed time series
fluxVItimesat1 <- fluxVItimesat %>% 
  mutate(validate = "GPP") %>% 
  rename_at(vars(contains("GPP")), ~str_replace(., "GPP", "Gvar"))
phenocamVItimesat1 <- phenocamVItimesat %>% 
  mutate(validate = "GCC") %>% 
  rename_at(vars(contains("GCC")), ~str_replace(., "GCC", "Gvar")) %>% 
  # completely mismatch GCC and VI, likely due to large uncertaities in GCC location
  filter(!(ID %in% c("monteblanco_SH_1000", "monteblanco_EB_1000", "montenegro_SH_1000",
                     "donanafuenteduque_WL_1000", "DE-Geb_AG_1000")))

VItimesat <- bind_rows(fluxVItimesat1, phenocamVItimesat1)
save(VItimesat, file = "data/RData/VItimesat_flux_phenocam.RData")

## TIMESAT phenometrics
fluxVIpheno1 <- fluxVIpheno %>% 
  mutate(validate = "GPP")
phenocamVIpheno1 <- phenocamVIpheno %>% 
  mutate(validate = "GCC") %>% 
  # completely mismatch GCC and VI, likely due to large uncertaities in GCC location
  filter(!(ID %in% c("monteblanco_SH_1000", "monteblanco_EB_1000", "montenegro_SH_1000",
                     "donanafuenteduque_WL_1000", "DE-Geb_AG_1000")))

VIpheno <- bind_rows(fluxVIpheno1, phenocamVIpheno1) 
save(VIpheno, file = "data/RData/VIpheno_flux_phenocam.RData")
