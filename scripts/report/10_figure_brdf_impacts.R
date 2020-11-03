library(tidyverse)
library(lubridate)
library(ggpubr)

# Impacts of BRDF correction on VI values ---------------------------------

load("data/RData/randompointsVI.RData")

forbrdfvalue <- randompointsVI %>% 
  filter(between(year(date), 2018, 2019)) %>% 
  drop_na() %>% 
  select(ID, date,
         NDVI_M, NDVI_raw,
         EVI2_M, EVI2_raw,
         PPI_M, PPI_raw) %>%
  pivot_longer(cols = -c(ID, date),
               names_to = c("VIname", "BRDF"),
               names_sep = "_",
               values_to = "VIval") %>%
  pivot_wider(names_from = "BRDF", values_from = "VIval", values_fn = list(VIval = mean)) %>%
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")))

forbrdfvalue  %>% 
  filter(!(VIname == "PPI" & (M > 3 | raw > 3 | M < 0 | raw < 0))) %>%
  filter(!(VIname == "EVI2" & (M > 1 | raw > 1))) %>%
  ggplot(aes(x = raw, y = M)) +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c(direction = -1) +
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  
  facet_wrap(vars(VIname), scales = "free") +
  labs(x = "TOC based VI", y = "NBAR based VI") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())+
  stat_regline_equation(aes(label =  paste(..eq.label..)), size = 4)+
  stat_cor(aes(label = ..rr.label..), method = "pearson",
           size = 4, color = "black",geom = "text",
           label.y.npc = 0.85,digits = 2)

ggsave("figures/report/figure_BRDF_impacts_value.pdf",
       width = 25, height = 8, units = "cm")



# BRDF impact on phenometrics ---------------------------------------------

files <- Sys.glob("data/3_timesat_outputs/phenometrics/randompoints/*.csv") %>%
  str_subset("0.25") %>% str_subset("M_loeFiltered|raw_loeFiltered")

data <- NULL
for (i in files) {
  VIname <- str_extract(i, "NDVI|EVI2|PPI")
  brdf <- str_sub(str_extract(i, "_raw|_M"), 2)
  tmp <- read_csv(i) %>%
    transmute(ID, year,
              SOS = start_of_season,
              EOS = end_of_season,
              VIname = VIname,
              BRDF = brdf) %>%
    filter(between(year, 2018, 2019)) %>%
    mutate(lat = as.numeric(str_extract(ID, " \\d{2}.\\d{1,}")))
  data <- bind_rows(data, tmp)
}


databrdf <- data %>% 
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric", 
               values_to = "metricValue") %>% 
  pivot_wider(names_from = BRDF, values_from = c(metricValue),
              values_fn = list(metricValue = mean)) %>% 
  mutate(metric = factor(metric, levels = c("SOS", "EOS"),
                         labels = c("SOS", "EOS")))

tbpsos <- databrdf %>% 
  filter(metric == "SOS") %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>% 
  
  ggplot(aes(raw, M)) +
  geom_bin2d(bins = 150) +
  scale_x_continuous(n.breaks = 4)+
  scale_y_continuous(n.breaks = 4)+
  scale_fill_viridis_c(direction = -1) +
  # geom_abline(slope = 1, intercept = 0, color = "grey40") +
  # facet_grid(metric ~ VIname, scales = "free") +
  facet_wrap(vars(VIname), scales = "fixed") +
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        # legend.position = "top",
        legend.text = element_text(size = 15)) +
  labs(y = "SOS from NBAR-based VI",
       x = "SOS from TOC-based VI")+
  stat_regline_equation(aes(label =  paste(..eq.label..)), size = 4)+
  stat_cor(aes(label = ..rr.label..), method = "pearson",
           size = 4, color = "black",geom = "text",
           label.y.npc = 0.85)

tbpeos <- databrdf %>% 
  filter(metric == "EOS") %>% 
  mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI"))) %>% 
  
  ggplot(aes(raw, M)) +
  geom_bin2d(bins = 150) +
  scale_x_continuous(n.breaks = 4)+
  scale_y_continuous(n.breaks = 4)+
  scale_fill_viridis_c(direction = -1) +
  # geom_abline(slope = 1, intercept = 0, color = "grey40") +
  # facet_grid(metric ~ VIname, scales = "free") +
  facet_wrap(vars(VIname), scales = "fixed") +
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        # legend.position = "top",
        legend.text = element_text(size = 15)) +
  labs(y = "EOS from NBAR-based VI",
       x = "EOS from TOC-based VI")+
  stat_regline_equation(aes(label =  paste(..eq.label..)), size = 4)+
  stat_cor(aes(label = ..rr.label..), method = "pearson",
           size = 4, color = "black",geom = "text",
           label.y.npc = 0.85, digits = 2)

ggarrange(tbpsos, tbpeos, labels = c("a", "b"),
          ncol = 1, common.legend = T, legend = "right",
          font.label = list(size = 17))

ggsave("figures/report/figure_BRDF_impacts_phenometrics.pdf",
       width = 25, height = 18, units = "cm")
