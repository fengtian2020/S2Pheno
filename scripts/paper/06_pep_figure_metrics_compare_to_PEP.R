
library(tidyverse)
library(ggpubr)

load("data/RData/pepS2.RData")


# plot scatterplots -------------------------------------------------------

brdf <- "M"

pepS2f3 <- pepS2 %>% 
  mutate(closet = abs(PEPdoy - S2doy)) %>% 
  group_by(s_id, type, year, thres, BRDF, VIname, phasePair) %>% 
  slice(which.min(closet)) %>% select(-closet) %>% 
  ungroup() %>% 
  mutate(pp = paste(phasePair, type, sep = " - "),
         pp = factor(pp, 
                     levels = rev(c("Leaf emerged - Broad-leaved-tree",
                                    "Leaf unfolded - Broad-leaved-tree",
                                    "Leaf unfolded 50% - Broad-leaved-tree",
                                    "Leaf emerged - Coniferous-tree",
                                    "Leaf unfolded - Coniferous-tree",
                                    "Leaf unfolded 50% - Coniferous-tree",
                                    "Fresh green 25% - Meadow",
                                    "Coloration 50% - Broad-leaved-tree",
                                    "Leaf fallen 50% - Broad-leaved-tree")),
                     labels = rev(c("Leaf emerged - Decidous",
                                    "Leaf unfolded - Decidous",
                                    "Leaf unfolded 50% - Decidous",
                                    "Leaf emerged - Evergreen",
                                    "Leaf unfolded - Evergreen",
                                    "Leaf unfolded 50% - Evergreen",
                                    "Fresh green 25% - Meadow",
                                    "       Coloration 50% - Decidous",
                                    "       Leaf fallen 50% - Decidous"))))

  # group_by(s_id, type, year, BRDF, VIname, phasePair) %>%
  # mutate(diffat = max(abs(PEPdoy - S2doy))) %>%
  # ungroup()
  # filter(s_id == 21, BRDF == "M") %>% arrange(phasePair, type, VIname) %>%  print(n = 10000)
  # filter(diffat < 120)

# pepS2f3 %>%
#   filter(BRDF == "M", thres == "50%", VIname == "PPI") %>% #, phasePair == "Leaf unfolded 50%"
#   select(-c(thres, BRDF, year)) %>%
#   group_by(VIname, type) %>%
#   ggplot(aes(x = S2doy, y = PEPdoy)) +
#   geom_point(alpha = 0.2, size = 1) +
#   geom_smooth(method = "lm", size = 0.8) +
#   facet_grid(type ~ phasePair, scales = "free")+
#   labs(x = "VI-derived SOS/EOS (day of year)",
#        y = "PEP725 phenophase (day of year)") +
#   theme_bw(base_size = 17)+
#   theme(panel.grid = element_blank(),
#         legend.title = element_blank(),mutate(pp = paste(phasePair, type, sep = " - "),
#         legend.position = "top",
#         legend.text = element_text(size = 17))
# ggsave("figures/figure_metrics_compare_to_PEP_scatterplots_outlierfiltered_NBAR.pdf",
#        width = 35, height = 20, units = "cm")


# parameters mean bias ------------------------------------------
brdf <- "M"
pepcmb <- function(brdf) {
  diffdfsos <- pepS2f3 %>% 
    filter(BRDF == brdf & !(pp %in% c("       Coloration 50% - Decidous",
                                      "       Leaf fallen 50% - Decidous"))) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, pp, thres) %>%add_count() %>% ungroup() %>%
    group_by(VIname, thres) %>% add_count() %>% ungroup() %>%
    mutate(weight = n / nn) %>%
    # filter(thres == "25%", VIname == "NDVI") %>% distinct(weight, phasePair, type) %>% summarise(w = sum(weight))
    group_by(VIname, thres, pp) %>% 
    summarise(diff = mean(PEPdoy - S2doy, na.rm = T),
              weight = first(weight),
              n = first(n)) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) 
  # number of paired PEP and VI statistics
  diffdfsos %>% group_by(pp) %>% summarise(n = mean(n))
  # 1 Fresh green 25% - Meadow       777. 
  # 2 Leaf unfolded 50% - Evergreen  194. 
  # 3 Leaf unfolded - Evergreen       12  
  # 4 Leaf emerged - Evergreen       766. 
  # 5 Leaf unfolded 50% - Decidous   264. 
  # 6 Leaf unfolded - Decidous      1100. 
  # 7 Leaf emerged - Decidous         99.0
  
  stddfsos <- pepS2f3 %>% 
    filter(BRDF == brdf & !(pp %in% c("       Coloration 50% - Decidous",
                                      "       Leaf fallen 50% - Decidous"))) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, pp, thres) %>% 
    summarise(std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}")))
  
  diffdfeos <- pepS2f3 %>% 
    filter(BRDF == brdf & (pp %in% c("       Coloration 50% - Decidous",
                                     "       Leaf fallen 50% - Decidous"))) %>% 
    select(-c(BRDF, year)) %>% 
    group_by(VIname, pp, thres) %>% add_count() %>% ungroup() %>%
    group_by(VIname, thres) %>% add_count() %>% ungroup() %>%
    mutate(weight = n / nn) %>%
    group_by(VIname, thres, pp) %>% 
    summarise(diff = mean(PEPdoy - S2doy, na.rm = T),
              weight = first(weight),
              n = first(n)) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) 
  # number of paired PEP and VI statistics
  diffdfeos %>% group_by(pp) %>% summarise(n = mean(n))
  # 1 "       Leaf fallen 50% - Decidous" 1141.
  # 2 "       Coloration 50% - Decidous"  1247.
  stddfeos <- pepS2f3 %>% 
    filter(BRDF == brdf & (pp %in% c("       Coloration 50% - Decidous",
                                     "       Leaf fallen 50% - Decidous"))) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, pp, thres) %>% 
    summarise(std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}")))
  
  # cordf <- pepS2f3 %>% 
  #   filter(BRDF == brdf) %>% 
  #   select(-c(BRDF, year)) %>%
  #   group_by(VIname, phasePair, thres, type) %>% 
  #   summarise(cor = cor(PEPdoy, S2doy, 
  #                       use = "pairwise.complete.obs")^2 * s - 150) %>% 
  #   ungroup() %>% 
  #   mutate(thres = as.numeric(str_extract(thres, "\\d{2}")))
  
  sosdata <- diffdfsos %>% left_join(stddfsos) 
  eosdata <- diffdfeos %>% left_join(stddfeos) 
  
  sosdatamean <- sosdata %>% 
    group_by(thres, VIname) %>% 
    summarise(across(c(diff, std), ~ sum(abs(.) * weight))) %>% ungroup() %>% 
    mutate(pp = "Weighted mean") %>%
    bind_rows(sosdata) %>% 
    mutate(pp = factor(pp, 
                       levels = rev(c("Leaf emerged - Decidous",
                                      "Leaf unfolded - Decidous",
                                      "Leaf unfolded 50% - Decidous",
                                      "Leaf emerged - Evergreen",
                                      "Leaf unfolded - Evergreen",
                                      "Leaf unfolded 50% - Evergreen",
                                      "Fresh green 25% - Meadow",
                                      "Weighted mean")),
                       labels = rev(c("Leaf emerged - Decidous\n(99)",
                                      "Leaf unfolded - Decidous\n(1100)",
                                      "Leaf unfolded 50% - Decidous\n(264)",
                                      "Leaf emerged - Evergreen\n(766)",
                                      "Leaf unfolded - Evergreen\n(12)",
                                      "Leaf unfolded 50% - Evergreen\n(194)",
                                      "Fresh green 25% - Meadow\n(777)",
                                      "Weighted mean"))))
  eosdatamean <- eosdata %>% 
    mutate(pp = factor(pp,
                       levels = rev(c("       Coloration 50% - Decidous",
                                      "       Leaf fallen 50% - Decidous")),
                       labels = rev(c("       Coloration 50% - Decidous\n(1247)",
                                      "       Leaf fallen 50% - Decidous\n(1141)"))))

  sosdatam <- sosdatamean %>% group_by(pp, VIname) %>% slice_min(abs(diff)) %>% ungroup()
  eosdatam <- eosdatamean %>% group_by(pp, VIname) %>% slice_min(abs(diff)) %>% ungroup()
  
  sos <- sosdatamean %>% 
    ggplot(aes(thres, pp)) +
    geom_tile(aes(fill = diff)) +
    geom_point(data = sosdatam, size = 3) +
    geom_text(aes(label = round(diff)), vjust = -1) +
    geom_text(aes(label = round(std)), color = "grey40", vjust = 2) +
    scale_fill_distiller(name="Bias (Days)", palette = "RdBu", limits = c(-55, 55)) +
    # coord_fixed() +
    scale_x_continuous(breaks = seq(5, 50, 5)) +
    labs(x = "VI amplitude threshold (%)",
         y = "") +
    facet_wrap(vars(VIname)) +
    # theme_classic(base_size = 15)
    theme(axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 15),
          panel.grid = element_blank(),
          # panel.background = element_blank(),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 13))
  
  eos <- eosdatamean %>% 
    ggplot(aes(thres, pp)) +
    geom_tile(aes(fill = diff)) +
    geom_point(data = eosdatam, size = 3) +
    geom_text(aes(label = round(diff)), vjust = -1) +
    geom_text(aes(label = round(std)), color = "grey40", vjust = 2) +
    scale_fill_distiller(name="Bias (Days)", palette = "RdBu", limits = c(-113, 113)) +
    # coord_fixed() +
    scale_x_continuous(breaks = seq(5, 50, 5)) +
    labs(x = "VI amplitude threshold (%)",
         y = "") +
    facet_wrap(vars(VIname)) +
    # theme_classic(base_size = 15)
    theme(axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 15),
          panel.grid = element_blank(),
          # panel.background = element_blank(),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 13))

  plot <- ggarrange(sos, eos, heights = c(7, 2.4),
            labels = c("a", "b"),
            legend = "right",
            font.label = list(size = 19),
            ncol = 1)
  return(plot)
}

pepcmb("M")
ggsave("figures/figure_metrics_compare_to_PEP_NBAR.pdf",
       width = 42, height = 25, units = "cm")


pepcmb("raw")
ggsave("figures/figure_metrics_compare_to_PEP_raw.pdf",
       width = 42, height = 25, units = "cm")

library(ggh4x)
# parameters of R2 ------------------------------------------
brdf <- "M"
pepcmb <- function(brdf, limt) {
  data <- pepS2f3 %>% 
    filter(BRDF == brdf, S2doy > -60) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, pp, thres) %>% 
    summarise(cor = cor(PEPdoy, S2doy),
              pv = cor.test(PEPdoy, S2doy, method = "pearson", alternative = "two.sided")$p.value) %>%
    ungroup() %>% 
    # mutate(across(cor, ~ifelse(pv < 0.05, cor, NA))) %>%
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) %>%
    mutate(pp = factor(pp, 
                       levels = rev(c("Leaf emerged - Decidous",
                                      "Leaf unfolded - Decidous",
                                      "Leaf unfolded 50% - Decidous",
                                      "Leaf emerged - Evergreen",
                                      "Leaf unfolded - Evergreen",
                                      "Leaf unfolded 50% - Evergreen",
                                      "Fresh green 25% - Meadow",
                                      "       Coloration 50% - Decidous",
                                      "       Leaf fallen 50% - Decidous")),
                       labels = rev(c("Leaf emerged - Decidous\n(99)",
                                      "Leaf unfolded - Decidous\n(1100)",
                                      "Leaf unfolded 50% - Decidous\n(264)",
                                      "Leaf emerged - Evergreen\n(766)",
                                      "Leaf unfolded - Evergreen\n(12)",
                                      "Leaf unfolded 50% - Evergreen\n(194)",
                                      "Fresh green 25% - Meadow\n(777)",
                                      "       Coloration 50% - Decidous\n(1247)",
                                      "       Leaf fallen 50% - Decidous\n(1141)"))))
  
  datam <- data %>% group_by(pp, VIname) %>% slice_max(cor) %>% ungroup()
  
  heatmap <- data %>% 
    ggplot(aes(thres, pp)) +
    geom_tile(aes(fill = cor)) +
    geom_point(data = datam, size = 3) +
    geom_text(data = datam, aes(label = round(cor, 2)), vjust = -1) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1*limt, limt), 
                         name = expression(italic("R"))) +
    scale_x_continuous(breaks = seq(5, 50, 5)) +
    labs(x = "VI amplitude threshold (%)",
         y = "") +
    facet_wrap(vars(VIname)) +
    theme(axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          strip.text = element_text(size = 15),
          panel.grid = element_blank(),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15))
  
  ##### scatter plots
  scatt <- pepS2f3 %>% 
    filter(BRDF == brdf) %>% 
    group_by(VIname, pp) %>% 
    # filter(thres == "50%", S2doy > -60, (pp != "Fresh green 25% - Meadow" & PEPdoy > 60)) %>%
    filter(thres == "50%", S2doy > -60,
           pp %in% c("Leaf unfolded 50% - Decidous", "Leaf unfolded 50% - Evergreen")) %>%
    mutate(pp = factor(pp,
                       levels = c("Leaf unfolded 50% - Decidous", "Leaf unfolded 50% - Evergreen"))) %>%
    ggplot(aes(x = S2doy, y = PEPdoy)) +
    geom_point(aes(color = VIname), show.legend = F, size = 1.2) +
    scale_color_manual(values = c("#FE6100", "blue", "#DC267F"))+
    # geom_smooth(method = "lm") +
    stat_cor(aes(label = ..r.label..), 
             method = "pearson", size = 4.5,
             label.x.npc = 0, label.y.npc = 0.98) +
    stat_cor(aes(label = ..p.label..), 
             method = "pearson", size = 4.7,
             label.x.npc = 0, label.y.npc = 0.82) +
    # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
    #          method = "pearson", size = 4.2,
    #          label.x.npc = 0, label.y.npc = 0.98) +
    scale_x_continuous(n.breaks = 4) +
    labs(x = "VI SOS (day of year)",
         y = "PEP725 phase (day of year)") +
    facet_nested(~ pp + VIname, scales = "free") +
    theme_bw(base_size = 17) +
    theme(panel.grid = element_blank(),
          # legend.text = element_text(size = 14),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))
  # scatt
  plot <- ggarrange(heatmap, scatt, 
                    heights = c(6, 3),
                    labels = c("a", "b"),
                    font.label = list(size = 19),
                    ncol = 1, nrow = 2)
    
  return(plot)
}

pepcmb("M", 0.65)
ggsave("figures/figure_metrics_compare_to_PEP_R2_NBAR1.pdf",
       width = 42, height = 28, units = "cm")


pepcmb("raw", 0.65)
ggsave("figures/figure_metrics_compare_to_PEP_R2_raw1.pdf",
       width = 42, height = 28, units = "cm")


  

    
  






