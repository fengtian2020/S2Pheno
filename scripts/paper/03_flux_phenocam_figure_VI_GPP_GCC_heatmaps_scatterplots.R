library(tidyverse)
library(lubridate)
library(ggpubr)
library(ggpmisc)
library(ggh4x)

load("data/RData/VIpheno_flux_phenocam.RData")

VIpheno0 <- VIpheno %>% 
  filter(between(year, 2018, 2019), !(str_detect(ID, "^[:upper:]") & year == 2019))

brdfm = "M"
vali = "GCC"
phenoData <- function(brdfm, vali){
  GphenoVI <- 
    VIpheno0 %>% filter(BRDF == brdfm, validate == vali, VIname != vali) %>%  #, VIname == "PPI", ID == "borgocioffisouth_AG_1000", ID == "juncabalejo_WL_1000"
    select(-c(BRDF, validate)) %>% 
    rename(thres_VI = thres) %>% 
    pivot_longer(cols = c(SOS, EOS),
                 names_to = "metric",
                 values_to = "VI")
  Gpheno <-
    VIpheno0 %>% filter(BRDF == brdfm, VIname == vali) %>% #, ID =="borgocioffisouth_AG_1000", ID == "juncabalejo_WL_1000"
    select(-c(BRDF, validate, VIname)) %>%
    rename(thres_Gvar = thres) %>% 
    pivot_longer(cols = c(SOS, EOS),
                 names_to = "metric",
                 values_to = "Gvar") %>% 
    left_join(GphenoVI, by = c("ID", "LC", "metric", "year")) %>% 
    
    ### only reserve the closet match for each year, in case there are 
    ### double growing season in eithor GPP or VI
    mutate(closet = abs(Gvar - VI)) %>% #ifelse(abs(Gvar - VI) > 183, abs(Gvar - VI), abs(Gvar - VI)) %>% 
    group_by(ID, LC, thres_Gvar, thres_VI, metric, VIname, year) %>% 
    slice(which.min(closet)) %>% 
    filter(closet < 210) %>% # to remove peaking year differ between VI and GPP/GCC, located in December and January.
    ungroup() %>% 
    
    mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")),
           metric = factor(metric, levels = c("SOS", "EOS"))) %>% 
    mutate(diff = abs(Gvar-VI),
           thres_VI = as.numeric(str_extract(thres_VI, "\\d{2}")),
           thres_Gvar = as.numeric(str_extract(thres_Gvar, "\\d{2}")),
           # thres_VI_P = paste0(str_extract(thres_VI, "\\d{2}"), "%"),
           # thres_Gvar = paste0(str_extract(thres_Gvar, "\\d{2}"), "%"),
           validate = vali)
  return(Gpheno)
}

# example scatterplots ----------------------------------------------------

phenoData("M", "GPP") %>%
  # filter(thres_Gvar %in% c("15%", "25%", "45%"), thres_VI_V == 25, metric == "EOS", VIname == "PPI") %>%
  filter(thres_Gvar %in% c(15, 25, 45), thres_VI == 25, metric == "SOS", VIname == "PPI") %>%
  group_by(thres_Gvar) %>%
  summarise(`absolute bias` = mean(diff))

dat_bias <- data.frame(
  label = c("Mean absolute bias = 18.1 days",
            "Mean absolute bias = 13.9 days",
            "Mean absolute bias = 11.8 days"),
  thres_Gvar = c(15, 25, 45)
)

phenoData("M", "GCC") %>% #filter(thres_Gvar %in% c(45), thres_VI == 15, metric == "EOS", VIname == "PPI", ID == "juncabalejo_WL_1000") 
  filter(thres_Gvar %in% c(15, 25, 45), thres_VI == 15, metric == "EOS", VIname == "PPI") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_abline(slope = 1, intercept = 0, color = "grey70")+
  geom_point() +
  geom_text(aes(label = ID), hjust = 1)+
  # scale_x_continuous(limits = c(0, 270))+
  # scale_y_continuous(limits = c(0, 270))+
  facet_wrap(vars(thres_Gvar)) +
  theme_bw(base_size = 14)+
  labs(x = "PPI SOS at threshold 25% (DOY)",
       y = "GPP SOS (DOY)") +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))+
  # geom_text(data = dat_bias, mapping = aes(x = 0, y = 220, label = label),
  #           hjust = 0) +
  stat_cor(aes(label = ..r.label..), method = "pearson", #size = 4.2,
           label.x.npc = "left", label.y.npc = "top")
  
ggsave("figures/figure_VI_GPP_GCC_scatterplots.pdf",
       width = 25, height = 10, units = "cm")


# heatmap matrix NBAR ----------------------------------------------------------

# GPP
gppm <- phenoData("M", "GPP") %>% 
  group_by(thres_Gvar, thres_VI, VIname, metric) %>% 
  summarise(meanDiff = mean(diff, na.rm = T),
            cor = cor(Gvar, VI, method = "pearson"),
            pv = cor.test(Gvar, VI, method = "pearson", alternative = "two.sided")$p.value
            ) %>%
  ungroup()%>% 
  mutate(across(pv, ~ifelse(. < 0.05, 1, 0)))
gppmm <- gppm %>% 
  group_by(thres_VI, VIname, metric) %>%
  slice(which.min(meanDiff)) %>%
  ungroup()


gpp <- ggplot(data = gppm, aes(thres_VI, thres_Gvar)) +
  geom_tile(aes(fill = meanDiff)) + 
  geom_point(data = gppmm, aes(color = cor), size = 4) +
  geom_text(data = gppmm, aes(label = round(cor, 2)), vjust = 2) +
  geom_text(data = gppmm, aes(label = round(meanDiff)), vjust = -1) +
  scale_fill_distiller(name="Bias (Days)", direction = 1, limits = c(11, 91)) +
  scale_size(range = c(1, 7), 
             guide = guide_legend(
               reverse = T, title = expression(italic("R"))
             )) + 
  scale_color_distiller(palette = "YlOrRd", direction = 1, limits = c(0.6, 0.9), 
                       name = expression(italic("R"))) +
  # scale_colour_viridis_c(option = "plasma", name = expression(italic("R"))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "VI amplitude threshold (%)",
       y = "GPP amplitude threshold (%)") +
  facet_grid(metric ~ VIname) +
  # theme_classic(base_size = 15)
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        # panel.background = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_colourbar(order = 1),
         color = guide_colourbar(order = 2))
gpp
# GCC
gccm <- phenoData("M", "GCC") %>% 
  group_by(thres_Gvar, thres_VI, VIname, metric) %>% 
  summarise(meanDiff = abs(mean(diff, na.rm = T)),
            cor = cor(Gvar, VI, method = "pearson",
                      use = "pairwise.complete.obs"),
            pv = cor.test(Gvar, VI, method = "pearson", alternative = "two.sided", 
                          use = "pairwise.complete.obs")$p.value
  ) %>%
  ungroup() %>% 
  mutate(across(pv, ~ifelse(. < 0.05, 1, 0)))
gccmm <- gccm %>% 
  group_by(thres_VI, VIname, metric) %>%
  slice(which.min(meanDiff)) %>%
  ungroup()

gcc <- ggplot(data = gccm, aes(thres_VI, thres_Gvar)) +
  geom_tile(aes(fill = meanDiff)) + 
  geom_point(data = gccmm, aes(color = cor), size = 4) +
  geom_text(data = gccmm, aes(label = round(cor, 2)), vjust = 2) +
  geom_text(data = gccmm, aes(label = round(meanDiff)), vjust = -1) +
  scale_fill_distiller(name="Bias (Days)", direction = 1, limits = c(11, 91)) +
  scale_size(range = c(1, 7), 
             guide = guide_legend(
               reverse = T, title = expression(italic("R"))
             )) + 
  scale_color_distiller(palette = "YlOrRd", direction = 1, limits = c(0.58, 0.9),
                        name = expression(italic("R"))) +
  # scale_colour_viridis_c(option = "plasma", name = expression(italic("R"))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "VI amplitude threshold (%)",
       y = "GCC amplitude threshold (%)") +
  facet_grid(metric ~ VIname) +
  # theme_classic(base_size = 15)
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        # panel.background = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_colourbar(order = 1),
         color = guide_colourbar(order = 2))
gcc

ggarrange(gpp, gcc, 
          labels = c("a", "b"), 
          common.legend = T,
          legend = "right",
          font.label = list(size = 19),
          ncol = 1)
ggsave("figures/figure_VI_GPP_GCC_heatmatrix_NBAR.pdf",
       width = 36, height = 50, units = "cm")



# scatterplots NBAR -------------------------------------------------------
#GPP
gppmmb <- gppmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_min(both, with_ties = F) %>% 
  mutate(perf = "Best VI threshold")
gppmmw <- gppmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_max(both, with_ties = F) %>% 
  mutate(perf = "Worst VI threshold")
gppmmbw <- bind_rows(gppmmb, gppmmw) %>% 
  mutate(label = paste0("Bias = ", round(meanDiff), " d", "\n", "R = ", round(cor,2))) %>% 
  mutate(label1 = paste0("thres = ", thres_VI, "%"))
gppmm2 <- bind_rows(gppmmb, gppmmw) %>% 
  select(metric, VIname, thres_Gvar, thres_VI, perf) %>% 
  left_join(phenoData("M", "GPP") %>% select(metric, VIname, thres_Gvar, thres_VI, Gvar, VI))

gppsos <- gppmm2 %>% filter(metric == "SOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gppmmbw, metric == "SOS"), 
            aes(npcx = 0.03, npcy = 0.95, label = label), 
            hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gppmmbw, metric == "SOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI SOS (day of year)", y = "GPP SOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))
gppeos <- gppmm2 %>% filter(metric == "EOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gppmmbw, metric == "EOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gppmmbw, metric == "EOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI EOS (day of year)", y = "GPP EOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

# GCC
gccmmb <- gccmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_min(both, with_ties = F) %>% 
  mutate(perf = "Best VI threshold")
gccmmw <- gccmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_max(both, with_ties = F) %>% 
  mutate(perf = "Worst VI threshold")
gccmmbw <- bind_rows(gccmmb, gccmmw) %>% 
  mutate(label = paste0("Bias = ", round(meanDiff), " d", "\n", "R = ", round(cor,2))) %>% 
  mutate(label1 = paste0("thres = ", thres_VI, "%"))
gccmm2 <- bind_rows(gccmmb, gccmmw) %>% 
  select(metric, VIname, thres_Gvar, thres_VI, perf) %>% 
  left_join(phenoData("M", "GCC") %>% select(metric, VIname, thres_Gvar, thres_VI, Gvar, VI))

gccsos <- gccmm2 %>% filter(metric == "SOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gccmmbw, metric == "SOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gccmmbw, metric == "SOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI SOS (day of year)", y = "GCC SOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))
gcceos <- gccmm2 %>% filter(metric == "EOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  geom_text(aes(label = ))
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gccmmbw, metric == "EOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gccmmbw, metric == "EOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI EOS (day of year)", y = "GCC EOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggarrange(gppsos, gppeos, gccsos, gcceos,
          labels = c("a", "b", "c", "d"), 
          font.label = list(size = 19),
          ncol = 2, nrow = 2)

ggsave("figures/figure_VI_GPP_GCC_scatterplots_NBAR.pdf",
       width = 40, height = 27, units = "cm")


# heatmap matrix raw ----------------------------------------------------------

# GPP
gppm <- phenoData("raw", "GPP") %>% 
  group_by(thres_Gvar, thres_VI, VIname, metric) %>% 
  summarise(meanDiff = abs(mean(diff, na.rm = T)),
            cor = cor(Gvar, VI, method = "pearson",
                      use = "pairwise.complete.obs"),
            pv = cor.test(Gvar, VI, method = "pearson", alternative = "two.sided", 
                          use = "pairwise.complete.obs")$p.value
  ) %>%
  ungroup()%>% 
  mutate(across(pv, ~ifelse(. < 0.05, 1, 0)))
gppmm <- gppm %>% 
  group_by(thres_VI, VIname, metric) %>%
  slice(which.min(meanDiff)) %>%
  ungroup()


gpp <- ggplot(data = gppm, aes(thres_VI, thres_Gvar)) +
  geom_tile(aes(fill = meanDiff)) + 
  geom_point(data = gppmm, aes(color = cor), size = 4) +
  geom_text(data = gppmm, aes(label = round(cor, 2)), vjust = -1) +
  geom_text(data = gppmm, aes(label = round(meanDiff)), vjust = 2) +
  scale_fill_distiller(name="Bias (Days)", direction = 1, limits = c(11, 91)) +
  scale_size(range = c(1, 7), 
             guide = guide_legend(
               reverse = T, title = expression(italic("R"))
             )) + 
  scale_color_distiller(palette = "YlOrRd", direction = 1,limits = c(0.58, 0.9),
                        name = expression(italic("R"))) +
  # scale_colour_viridis_c(option = "plasma", name = expression(italic("R"))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "VI amplitude threshold (%)",
       y = "GPP amplitude threshold (%)") +
  facet_grid(metric ~ VIname) +
  # theme_classic(base_size = 15)
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        # panel.background = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_colourbar(order = 1),
         color = guide_colourbar(order = 2))
gpp

# GCC
gccm <- phenoData("raw", "GCC") %>% 
  group_by(thres_Gvar, thres_VI, VIname, metric) %>% 
  summarise(meanDiff = abs(mean(diff, na.rm = T)),
            cor = cor(Gvar, VI, method = "pearson",
                      use = "pairwise.complete.obs"),
            pv = cor.test(Gvar, VI, method = "pearson", alternative = "two.sided", 
                          use = "pairwise.complete.obs")$p.value
  ) %>%
  ungroup() %>% 
  mutate(across(pv, ~ifelse(. < 0.05, 1, 0)))
gccmm <- gccm %>% 
  group_by(thres_VI, VIname, metric) %>%
  slice(which.min(meanDiff)) %>%
  ungroup()

gcc <- ggplot(data = gccm, aes(thres_VI, thres_Gvar)) +
  geom_tile(aes(fill = meanDiff)) + 
  geom_point(data = gccmm, aes(color = cor), size = 4) +
  geom_text(data = gccmm, aes(label = round(cor, 2)), vjust = -1) +
  geom_text(data = gccmm, aes(label = round(meanDiff)), vjust = 2) +
  scale_fill_distiller(name="Bias (Days)", direction = 1, limits = c(11, 91)) +
  scale_size(range = c(1, 7), 
             guide = guide_legend(
               reverse = T, title = expression(italic("R"))
             )) + 
  scale_color_distiller(palette = "YlOrRd", direction = 1,limits = c(0.58, 0.9),
                        name = expression(italic("R"))) +
  # scale_colour_viridis_c(option = "plasma", name = expression(italic("R"))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "VI amplitude threshold (%)",
       y = "GCC amplitude threshold (%)") +
  facet_grid(metric ~ VIname) +
  # theme_classic(base_size = 15)
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        # panel.background = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_colourbar(order = 1),
         color = guide_colourbar(order = 2))
gcc

ggarrange(gpp, gcc, 
          labels = c("a", "b"), 
          common.legend = T,
          legend = "right",
          font.label = list(size = 19),
          ncol = 1)
ggsave("figures/figure_VI_GPP_GCC_heatmatrix_raw.pdf",
       width = 36, height = 50, units = "cm")


# scatterplots raw -------------------------------------------------------
#GPP
gppmmb <- gppmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_min(both, with_ties = F) %>% 
  mutate(perf = "Best VI threshold")
gppmmw <- gppmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_max(both, with_ties = F) %>% 
  mutate(perf = "Worst VI threshold")
gppmmbw <- bind_rows(gppmmb, gppmmw) %>% 
  mutate(label = paste0("Bias = ", round(meanDiff), " d", "\n", "R = ", round(cor,2))) %>% 
  mutate(label1 = paste0("thres = ", thres_VI, "%"))
gppmm2 <- bind_rows(gppmmb, gppmmw) %>% 
  select(metric, VIname, thres_Gvar, thres_VI, perf) %>% 
  left_join(phenoData("M", "GPP") %>% select(metric, VIname, thres_Gvar, thres_VI, Gvar, VI))

gppsos <- gppmm2 %>% filter(metric == "SOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gppmmbw, metric == "SOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gppmmbw, metric == "SOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI SOS (day of year)", y = "GPP SOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))
gppeos <- gppmm2 %>% filter(metric == "EOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gppmmbw, metric == "EOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gppmmbw, metric == "EOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI EOS (day of year)", y = "GPP EOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

# GCC
gccmmb <- gccmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_min(both, with_ties = F) %>% 
  mutate(perf = "Best VI threshold")
gccmmw <- gccmm %>% 
  mutate(both = meanDiff / cor) %>% 
  group_by(VIname, metric) %>% 
  slice_max(both, with_ties = F) %>% 
  mutate(perf = "Worst VI threshold")
gccmmbw <- bind_rows(gccmmb, gccmmw) %>% 
  mutate(label = paste0("Bias = ", round(meanDiff), " d", "\n", "R = ", round(cor,2))) %>% 
  mutate(label1 = paste0("thres = ", thres_VI, "%"))
gccmm2 <- bind_rows(gccmmb, gccmmw) %>% 
  select(metric, VIname, thres_Gvar, thres_VI, perf) %>% 
  left_join(phenoData("M", "GCC") %>% select(metric, VIname, thres_Gvar, thres_VI, Gvar, VI))

gccsos <- gccmm2 %>% filter(metric == "SOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gccmmbw, metric == "SOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gccmmbw, metric == "SOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI SOS (day of year)", y = "GCC SOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))
gcceos <- gccmm2 %>% filter(metric == "EOS") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_point(aes(color = VIname), size = 1.2, show.legend = F)+
  scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
  geom_smooth(method = "lm", size = 0.5, se = F) +
  geom_text_npc(data = filter(gccmmbw, metric == "EOS"), 
                aes(npcx = 0.03, npcy = 0.95, label = label), 
                hjust = 0, size = 4.5) +
  geom_text_npc(data = filter(gccmmbw, metric == "EOS"), 
                aes(npcx = 0.97, npcy = 0.03, label = label1), 
                hjust = 1, size = 4.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey40") +
  facet_grid(perf ~ VIname, scales = "fixed")+ 
  theme_bw(base_size = 15) +
  labs(x = "VI EOS (day of year)", y = "GCC EOS (day of year)") +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggarrange(gppsos, gppeos, gccsos, gcceos,
          labels = c("a", "b", "c", "d"), 
          font.label = list(size = 19),
          ncol = 2, nrow = 2)

ggsave("figures/figure_VI_GPP_GCC_scatterplots_raw.pdf",
       width = 40, height = 27, units = "cm")


