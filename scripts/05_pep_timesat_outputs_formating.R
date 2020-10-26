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















###############################################################
### TIMESAT phenometrics thresholds assessment
###############################################################
files <- Sys.glob("Sen2Pheno/TIMESAT/PEP/results_PEP_S2_VI_SCL45_LC_4km_mean_DL-SP_sos_eos/*.csv") %>%
  str_subset("(PPI_M_loeFiltered|EVI2_M_loeFiltered|NDVI_M_loeFiltered|PPI_raw_loeFiltered|EVI2_raw_loeFiltered|NDVI_raw_loeFiltered)")
#i <- files[1]
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

raw <- read_csv("Sen2Pheno/PEP/pep725_in_use.csv") %>% 
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
  filter(!(type == "Crop" & 
             phasePair == "Leaf emerged")) %>% #(PEPdoy <= 150 | PEPdoy >= 260)
  # for report, 70% sites
  filter(!(s_id %in% c(read_csv("Sen2Pheno/PEP/PEP_site_for_valiation.csv") %>% pull(s_id))))

save(pepS2, file = "Sen2Pheno/RData/pepS2report.RData")

# set.seed(2020)
# pepS2 %>% distinct(s_id, type) %>% 
#   filter(!(type %in% c("Fruit-tree", "Vineyards"))) %>% 
#   add_count(type) %>% distinct(s_id, type, n) %>% 
#   group_by(type) %>% 
#   sample_frac(size = .7) %>% 
#   ungroup() %>% 
#   distinct(s_id)

##################################
# statistics of pep records
pepS2 %>% filter(thres == "25%", BRDF == "raw", VIname == "NDVI", type != "Crop") %>% 
  ggplot(aes(PEPdoy, color = phasePair))+
  geom_freqpoly()+
  # scale_colour_discrete() +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(type)) +
  labs(x = "day of year") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        # legend.key.size = unit(1, "cm"),
        legend.position = "top",
  ) 
ggsave("Sen2Pheno/figures/PEP_stat_3_paper.pdf",
       width = 20, height = 10, units = "cm")



#################################
## comparison all

brdf = "M"
ind = "diff"

pepcmb <- function(brdf, ind) {
  threscomall <- pepS2 %>% 
    filter(BRDF == brdf) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, phasePair, thres, type) %>% 
    summarise(cor = cor(PEPdoy, S2doy, method = "spearman",
                        use = "pairwise.complete.obs"),
              diff = mean(PEPdoy - S2doy, na.rm = T),
              std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    # filter(cor > 0.2) %>% 
    pivot_longer(cols = c(cor, diff, std),
                 names_to = "indicator",
                 values_to = "value") %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) 
  
  if (ind == "cor") {
    threscomall %>% 
      filter(indicator == ind, value > 0.2) %>% 
      select(-indicator) %>% 
      
      ggplot(aes(x = thres, y = value, color = VIname)) +
      geom_line() +
      geom_point(size = 2) +
      geom_vline(xintercept = 25, linetype = 2, color = "grey60")+
      scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
      facet_grid(type ~ phasePair) +
      # facet_grid(phasePair ~ type) +
      theme_bw(base_size = 17)+
      labs(x = "Amplitude thresholds (%)", 
           y = "Spatial correlation") +
      theme(panel.grid = element_blank(),
            # panel.grid.major.x = element_blank(),
            legend.title = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 17)
            # legend.direction = "horizontal"
      )
  } else {
    threscomall %>% filter(indicator != "cor", type != "Crop") %>% 
      pivot_wider(names_from = indicator, values_from = value) %>% 
      ggplot(aes(x = thres, y = diff, color = VIname)) +
      # geom_ribbon(aes(ymin = diff - std, ymax = diff + std, 
      #                 fill = VIname), alpha = 0.1, outline.type = "none") +
      # geom_line() +
      # geom_point(size = 2) +
      geom_pointrange(aes(ymin = diff - std, ymax = diff + std),
                      fatten = 3,
                      position=position_dodge(width = 3))+
      geom_vline(xintercept = 25, linetype = 2, color = "grey60")+
      geom_hline(yintercept = 0, linetype = 2, color = "grey60")+
      scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
      # scale_fill_manual(values = c("tan3", "blue", "firebrick3")) +
      facet_grid(type ~ phasePair, scales = "free") +
      theme_bw(base_size = 20)+
      labs(x = "Amplitude thresholds (%)", 
           y = "Time difference (day)") +
      theme(panel.grid = element_blank(),
            # panel.grid.major.x = element_blank(),
            legend.title = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 20)
            # legend.direction = "horizontal"
      )
  }
}

# pepcmb("raw", "cor")
# ggsave("Sen2Pheno/figures/PEP_correlation_thres.pdf",
#        width = 45, height = 30, units = "cm")

pepcmb("raw", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_thres_3_report.pdf",
       width = 40, height = 25, units = "cm")
pepcmb("M", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_thres_3_report_NBAR.pdf",
       width = 40, height = 25, units = "cm")

pepcmb("M", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_thres_3_paper.pdf",
       width = 40, height = 25, units = "cm")
pepcmb("raw", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_thres_3_paper_raw.pdf",
       width = 40, height = 25, units = "cm")

##########
# bias 25% threshold

pepS2 %>% mutate(sel = paste0(thres, phasePair)) %>% #distinct(phasePair)
  filter(BRDF == "raw", VIname == "PPI", type != "Crop", 
         sel %in% c("25%Leaf unfolded", "25%Leaf emerged", "25%Leaf unfolded 50%",
                    "25%Fresh green 25%", "15%Leaf fallen 50%", "15%Coloration  50%")) %>% 
  select(type, phasePair, PEPdoy, S2doy) %>% 
  group_by(phasePair, type) %>% 
  summarise(diff = mean(PEPdoy - S2doy, na.rm = T),
            std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(phasePair = fct_rev(phasePair),
         upper = diff + std,
         lower = diff - std) %>% 
  
  ggplot(aes(x = diff, y = phasePair)) +
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin = lower, xmax = upper, color = type),
                  fatten = 4,
                  position=position_dodge(width = .3),
  )+
  scale_x_continuous(breaks = seq(-50, 70, 10))+
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PEP doy  -  Sentinel2 doy", y = "") +
  theme_bw(base_size = 15) +
  theme(legend.title = element_blank(),
        legend.position = c(.82, .8),
  ) 
# stat_cor(aes(color = type, label = ..r.label..), method = "spearman", #size = 4.2,
#          label.x.npc = "left", label.y.npc = "top")

ggsave("Sen2Pheno/figures/PEP_PPI_amplitude_25_15_diff_3_report.pdf",
       width = 20, height = 10, units = "cm")



#############################################################
# fixed thresholds

files <- Sys.glob("/projects/eko/fs5/CROSSDRO/Sen2Pheno/TIMESAT/PEP/fixed_thresholds/*.csv") %>% 
  str_subset("(PPI_M_loeFiltered|EVI2_M_loeFiltered|NDVI_M_loeFiltered|PPI_raw_loeFiltered|EVI2_raw_loeFiltered|NDVI_raw_loeFiltered)")
# files <- files[2:length(files)]
pepVIpheno <- NULL
for (i in files) {
  tmp <- str_extract(i, "\\d.\\d{3}_(\\w{3,4})_(\\w{1,3})")
  thres <- str_split(tmp, "_", simplify = TRUE)[1]
  # str_extract("\\d{2}") %>% str_c("%")
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

raw <- read_csv("Sen2Pheno/PEP/pep725_in_use.csv") %>% 
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
  filter(!(type == "Crop" & 
             phasePair == "Leaf emerged")) %>% #(PEPdoy <= 150 | PEPdoy >= 260)
  # for report, 70% sites
  filter(!(s_id %in% c(read_csv("Sen2Pheno/PEP/PEP_site_for_valiation.csv") %>% pull(s_id)))) %>% 
  mutate(thres = paste0(VIname, thres))

# tt <- pepS2 %>% 
#   filter(VIname == "PPI", thres %in% c("0.060", "0.120", "0.180")) %>% 
#   mutate(thres = recode(thres, 
#                         "0.060" = "0.061",
#                         "0.120" = "0.121",
#                         "0.180" = "0.181"))
# ff <- pepS2 %>% 
#   filter(!(VIname == "PPI" & thres %in% c("0.060", "0.120", "0.180"))) 
# 
# pepS2 <- bind_rows(tt, ff)

brdf = "M"
ind = "diff"

level_key <- c("NDVI0.020" = "1", "EVI20.015" = "1", "PPI0.060" = "1",
               "NDVI0.040" = "2", "EVI20.030" = "2", "PPI0.120" = "2",
               "NDVI0.060" = "3", "EVI20.045" = "3", "PPI0.180" = "3",
               "NDVI0.080" = "4", "EVI20.060" = "4", "PPI0.240" = "4",
               "NDVI0.100" = "5", "EVI20.075" = "5", "PPI0.300" = "5",
               "NDVI0.120" = "6", "EVI20.090" = "6", "PPI0.360" = "6",
               "NDVI0.140" = "7", "EVI20.105" = "7", "PPI0.420" = "7",
               "NDVI0.160" = "8", "EVI20.120" = "8", "PPI0.480" = "8",
               "NDVI0.180" = "9", "EVI20.135" = "9", "PPI0.540" = "9",
               "NDVI0.200" = "10", "EVI20.150" = "10", "PPI0.600" = "10")

pepcmb <- function(brdf, ind) {
  threscomall <- pepS2 %>% 
    mutate(thres = as.numeric(recode(thres, !!!level_key))) %>%
    filter(BRDF == brdf) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, phasePair, thres, type) %>% 
    summarise(cor = cor(PEPdoy, S2doy, method = "spearman",
                        use = "pairwise.complete.obs"),
              diff = mean(PEPdoy - S2doy, na.rm = T),
              std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    # filter(cor > 0.2) %>% 
    pivot_longer(cols = c(cor, diff, std),
                 names_to = "indicator",
                 values_to = "value")
  
  threscomall %>% filter(indicator != "cor", type != "Crop") %>% 
    pivot_wider(names_from = indicator, values_from = value) %>%
    ggplot(aes(x = thres, y = diff, color = VIname)) +
    geom_pointrange(aes(ymin = diff - std, ymax = diff + std),
                    fatten = 3,
                    position=position_dodge(width = 0.4)
    )+
    scale_x_continuous(breaks = seq(1,10,2))+
    geom_vline(xintercept = 5, linetype = 2, color = "grey60")+
    geom_hline(yintercept = 0, linetype = 2, color = "grey60")+
    scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
    facet_grid(type ~ phasePair, scales = "free") +
    theme_bw(base_size = 20)+
    labs(x = "VI fixed thresholds (level)", 
         y = "Time difference (day)") +
    theme(panel.grid = element_blank(),
          # panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 20)
          # legend.direction = "horizontal"
    )
}

# pepcmb("raw", "cor")
# ggsave("Sen2Pheno/figures/PEP_correlation_thres.pdf",
#        width = 45, height = 30, units = "cm")

pepcmb("raw", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_fixed_thres_3_report.pdf",
       width = 40, height = 25, units = "cm")
pepcmb("M", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_fixed_thres_3_report_NBAR.pdf",
       width = 40, height = 25, units = "cm")

pepcmb("M", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_fixed_thres_3_paper.pdf",
       width = 40, height = 25, units = "cm")
pepcmb("raw", "diff")
ggsave("Sen2Pheno/figures/PEP_bias_fixed_thres_3_paper_raw.pdf",
       width = 40, height = 25, units = "cm")

##########
# fixed threshold bias

pepS2 %>% mutate(sel = paste0(thres, phasePair)) %>% #filter(str_detect(thres, "PPI")) %>% distinct(thres)
  filter(BRDF == "raw", VIname == "PPI", type != "Crop", 
         sel %in% c("PPI0.300Leaf unfolded", "PPI0.300Leaf emerged", "PPI0.300Leaf unfolded 50%",
                    "PPI0.300Fresh green 25%", "PPI0.120Leaf fallen 50%", "PPI0.120Coloration  50%")) %>% 
  select(type, phasePair, PEPdoy, S2doy) %>% 
  group_by(phasePair, type) %>% 
  summarise(diff = mean(PEPdoy - S2doy, na.rm = T),
            std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(phasePair = fct_rev(phasePair),
         upper = diff + std,
         lower = diff - std) %>% 
  
  ggplot(aes(x = diff, y = phasePair)) +
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin = lower, xmax = upper, color = type),
                  fatten = 4,
                  position=position_dodge(width = .3),
  )+
  scale_x_continuous(breaks = seq(-50, 70, 10))+
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PEP doy  -  Sentinel2 doy", y = "") +
  theme_bw(base_size = 15) +
  theme(legend.title = element_blank(),
        legend.position = c(.82, .8),
  ) 
# stat_cor(aes(color = type, label = ..r.label..), method = "spearman", #size = 4.2,
#          label.x.npc = "left", label.y.npc = "top")

ggsave("Sen2Pheno/figures/PEP_PPI_fixedthres_25_15_diff_3_report.pdf",
       width = 20, height = 10, units = "cm")













