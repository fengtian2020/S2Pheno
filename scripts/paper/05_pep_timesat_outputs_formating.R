library(tidyverse)
library(lubridate)

files <- Sys.glob("data/3_timesat_outputs/phenometrics/pep725_new/*.csv") %>%
  str_subset("(PPI_M_loeFiltered|EVI2_M_loeFiltered|NDVI_M_loeFiltered|PPI_raw_loeFiltered|EVI2_raw_loeFiltered|NDVI_raw_loeFiltered)")

pepVIpheno <- NULL
for (i in files) {
  tmp <- str_extract(i, "\\d.\\d{2,3}_(\\w{3,4})_(\\w{1,3})")
  thres <- str_split(tmp, "_", simplify = TRUE)[1] %>% 
    str_extract("\\d{2}") %>% str_c("%")
  VIname <- str_split(tmp, "_", simplify = TRUE)[2]
  BRDF <- str_split(tmp, "_", simplify = TRUE)[3]
  tmp <- read_csv(i) %>% 
    separate(ID, c("type", "s_id"), sep = "_") %>% 
    transmute(s_id = as.numeric(s_id), type, year,
              SOS = start_of_season,
              EOS = end_of_season,
              thres = thres, 
              BRDF = BRDF,
              VIname = VIname) %>% 
    pivot_longer(cols = c(SOS, EOS), names_to = "metric",
                 values_to = "S2doy")
  pepVIpheno <- bind_rows(pepVIpheno, tmp)
}

raw <- read_csv("data/0_ground_data_raw/pep725/pep725_in_use.csv") %>% 
  mutate(year = year(pheno), PEPdoy = yday(pheno)) %>% 
  select(s_id, type, phase_id, year, PEPdoy) %>% 
  group_by(s_id, type, phase_id, year) %>% 
  summarise(PEPdoy = mean(PEPdoy)) %>% 
  ungroup()

pepS2 <- raw %>% 
  left_join(pepVIpheno) %>% 
  mutate(phasePair = paste(phase_id, metric, sep = "_")) %>% 
  select(-c(phase_id, metric)) %>% 
  filter(phasePair %in% c("10_SOS", "11_SOS", "223_SOS", "182_SOS",
                          "100_EOS", "151_EOS", "205_EOS", "95_EOS")) %>% 
  mutate(phasePair = factor(phasePair,
                            levels = c("10_SOS", "11_SOS", "223_SOS", "182_SOS",
                                       "100_EOS", "151_EOS", "205_EOS", "95_EOS"),
                            labels = c("Leaf emerged",
                                       "Leaf unfolded",
                                       "Leaf unfolded 50%",
                                       "Fresh green 25%",
                                       "Start of harvest",
                                       "Silage harvest",
                                       "Coloration  50%",
                                       "Leaf fallen 50%")),
         VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>% 
  
  # filter crop leaf emerged
  filter(!(type == "Crop")) #%>% (PEPdoy <= 150 | PEPdoy >= 260)
# for report, 70% sites
# filter(!(s_id %in% c(read_csv("Sen2Pheno/PEP/PEP_site_for_valiation.csv") %>% pull(s_id))))

save(pepS2, file = "data/RData/pepS2.RData")


