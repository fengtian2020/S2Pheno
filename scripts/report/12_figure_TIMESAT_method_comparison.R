library(sf)
library(tidyverse)
library(lubridate)
library(ggpubr)

# Sys.glob("data/3_timesat_outputs/smooth_time_series/*Flux*.csv") %>% 
filesGPP <- Sys.glob("data/3_timesat_outputs/method_comparison/*GPP*.csv")

GPPpheno <- NULL
for (i in filesGPP) {
  method <- str_extract(i, "DLSP|DL|SP")
  tmp <- read_csv(i) %>% 
    transmute(ID, year,
              SOS = start_of_season,
              EOS = end_of_season,
              method = method) %>% 
    filter(year == 2018) %>% 
    select(-year) %>% 
    pivot_longer(cols = c(SOS, EOS),
                 values_to = "GPP",
                 names_to = "metric")
  GPPpheno <- bind_rows(GPPpheno, tmp)
}

filesPPI <- Sys.glob("data/3_timesat_outputs/method_comparison/*PPI*.csv")

PPIpheno <- NULL
for (i in filesPPI) {
  method <- str_extract(i, "DLSP|DL|SP")
  tmp <- read_csv(i) %>% 
    transmute(ID, year,
              SOS = start_of_season,
              EOS = end_of_season,
              method = method) %>% 
    filter(year == 2018) %>% 
    select(-year) %>% 
    pivot_longer(cols = c(SOS, EOS),
                 values_to = "PPI",
                 names_to = "metric")
  PPIpheno <- bind_rows(PPIpheno, tmp)
}

GPPpheno %>% left_join(PPIpheno) %>% 
  mutate(closet = abs(GPP - PPI)) %>% 
  group_by(ID, method, metric) %>% 
  slice(which.min(closet)) %>% 
  ungroup() %>% 
  mutate(method = factor(method,
                         levels = c("SP", "DL", "DLSP")),
         metric = factor(metric, levels = c("SOS", "EOS"))) %>% 
  
  ggplot(aes(x = PPI, y = GPP, color = metric)) +
  geom_point() +
  scale_colour_brewer(palette = "Dark2") +
  geom_abline(slope = 1, intercept = 0, color = "grey40")+
  facet_wrap(vars(method)) +
  theme_bw(base_size = 15) +
  labs(x = "PPI DOY", y = "GPP DOY") +
  theme(legend.position = c(0.92, 0.2),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        title = element_text(color = "black"))+
  stat_cor(aes(label = ..r.label..), method = "pearson", size = 4.2,
           label.x.npc = "left", label.y.npc = "top")

ggsave("figures/report/figure_method_cmp_pheno_raw.pdf",
       width = 25, height = 9, units = "cm")

####################################################  

# LC fluxsite --------------------------------------------------------------

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

gpp <- read_csv("data/2_timesat_inputs/Flux_sites_S2_VI_SCL45_100m_mean.csv") %>% 
  select(ID, date, gpp = GPP_loess) %>% 
  filter(year(date) == 2018) %>% 
  left_join(fluxLC)

dlsp <- read_csv("data/3_timesat_outputs/smooth_time_series/smoothed_Flux_sites_S2_VI_SCL45_100m_mean_newDLSP.csv") %>% 
  select(ID, date, dlsp = PPI_raw_loeFiltered) %>% 
  left_join(fluxLC)

dl <- read_csv("data/3_timesat_outputs/smooth_time_series/smoothed_Flux_sites_S2_VI_SCL45_100m_mean_newDL.csv") %>% 
  select(ID, date, dl = PPI_raw_loeFiltered) %>% 
  left_join(fluxLC)

sp <- read_csv("data/3_timesat_outputs/smooth_time_series/smoothed_Flux_sites_S2_VI_SCL45_100m_mean_SP.csv") %>% 
  select(ID, date, sp = PPI_raw_loeFiltered) %>% 
  left_join(fluxLC)

cmp <- gpp %>% left_join(dlsp) %>% left_join(dl) %>% left_join(sp) %>% 
  pivot_longer(cols = c(dlsp, dl, sp),
               names_to = "method",
               values_to = "value") %>% 
  group_by(ID) %>% 
  mutate(across(where(is.numeric), 
                list(Z = ~ (.x - mean(.x, na.rm = T)) / sd(.x, na.rm = T)))) %>%
  ungroup() %>% 
  mutate(method = factor(method,
                         levels = c("sp", "dl", "dlsp"),
                         labels = c("SP", "DL", "DL-SP")))

# rmse <- cmp %>% mutate(diff = abs(gpp_Z - value_Z)) %>% 
#   group_by(method) %>% 
#   summarise(rmse = mean(diff))
# ggplot(cmp, aes(x = value_Z, y = gpp_Z)) +
#   geom_point(size = 0.1)+
#   geom_abline(slope = 1, intercept = 0, color = "grey40")+
#   facet_wrap(vars(method)) +
#   theme_bw(base_size = 15) +
#   labs(x = "PPI Z-score", y = "GPP Z-score") +
#   theme(panel.grid = element_blank()) +
#   stat_cor(aes(label = ..rr.label..), method = "pearson", #size = 4.2,
#            label.x.npc = "left", label.y.npc = "top")


data <- cmp %>%
  group_by(ID, LC, method) %>% 
  summarise(cor = cor(gpp_Z, value_Z,
                      method = "pearson",
                      use = "pairwise.complete.obs")^2) %>% 
  ungroup()

boxLC <- ggplot(data, aes(x = method, y = cor)) +
  geom_boxplot(aes(fill = method), width = 0.4, outlier.shape = NA) +
  geom_point(size = 1) +
  geom_line(aes(group = ID), color = "grey40", size = 0.3) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(vars(LC), scales = "free_y") +
  theme_bw(base_size = 15) +
  labs(x = "", y = expression("Temporal "~italic(R)^2~" per site")) +
  theme(legend.position = "top",
        legend.title = element_blank(), 
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        title = element_text(color = "black"))

#########
scaterLC <- ggplot(cmp, aes(x = value_Z, y = gpp_Z)) +
  geom_point(aes(color = LC), size = 0.2, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "grey40")+
  # annotate("text", x = 2.6, y = 3.0, angle = 52, label = "1:1") +
  geom_smooth(aes(color = LC), method = lm, se = F, size = 0.5) +
  xlim(-2, 3) +
  ylim(-2, 3) +
  scale_colour_brewer(palette = "Dark2") +
  # stat_poly_eq(formula = y~x, aes(label = ..rr.label..), size = 4.2, parse = TRUE) +
  facet_wrap(vars(method)) +
  theme_bw(base_size = 15) +
  labs(x = "PPI Z-score", y = "GPP Z-score") +
  theme(legend.position = "top",
        legend.title = element_blank(), 
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        title = element_text(color = "black"))+
  stat_cor(aes(label = ..rr.label..), method = "pearson", size = 4.2,
           label.x.npc = "left", label.y.npc = "top")


cmpbar <- cmp %>% 
  group_by(LC, method) %>% 
  summarise(cor = cor(gpp_Z, value_Z, 
                      method = "pearson",
                      use = "pairwise.complete.obs")^2)

barLC <- ggplot(cmpbar, aes(x = LC, y = cor, fill = method)) +
  geom_col(position = "dodge") +
  scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) +
  geom_text(aes(label = round(cor, 2)), 
            position = position_dodge(width = 1), vjust = -0.5) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  labs(x = "", y = expression("Spatiotemperal "~italic(R)^2~" per LC")) +  
  theme(legend.position = "top",
        legend.title = element_blank(), 
        legend.text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        title = element_text(color = "black"))

ggarrange(boxLC, 
          ggarrange(scaterLC, barLC, 
                    labels = c("b", "c"), 
                    ncol = 1,
                    font.label = list(size = 15)),
          ncol = 2, 
          labels = "a", 
          font.label = list(size = 15))

ggsave("figures/report/figure_method_cmp_trahectory_raw.pdf",
       width = 35, height = 20, units = "cm")

