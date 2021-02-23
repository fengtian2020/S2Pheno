library(tidyverse)
library(lubridate)
library(ggpubr)

load("data/RData/VItimesat_flux_phenocam.RData")


VItimesat0 <- VItimesat %>% 
  filter(between(year(date), 2018, 2019), !(str_detect(ID, "^[:upper:]") & year(date) >= 2019))

VIcomp <- function(variable, brdfm) {
  if (brdfm == "raw") {
    data <- VItimesat0 %>% filter(validate == variable) %>%
      group_by(ID, LC) %>% 
      summarise(NDVI = cor(Gvar_Z, NDVI_raw_loeFiltered_Z,
                           method = "pearson",
                           use = "pairwise.complete.obs"),
                EVI2 = cor(Gvar_Z, EVI2_raw_loeFiltered_Z, 
                           method = "pearson",
                           use = "pairwise.complete.obs"),
                PPI = cor(Gvar_Z, PPI_raw_loeFiltered_Z, 
                          method = "pearson",
                          use = "pairwise.complete.obs")) %>% 
      pivot_longer(cols = c(NDVI, EVI2, PPI), 
                   names_to = "var", values_to = "cor") %>% 
      mutate(var = factor(var, levels = c("NDVI", "EVI2", "PPI")))
    
    VIselectZ <- VItimesat0 %>% filter(validate == variable) %>% 
      select(Gvar_Z, LC,
             NDVI_raw_loeFiltered_Z, 
             EVI2_raw_loeFiltered_Z, 
             PPI_raw_loeFiltered_Z) %>% 
      pivot_longer(cols = c(NDVI_raw_loeFiltered_Z, 
                            EVI2_raw_loeFiltered_Z, 
                            PPI_raw_loeFiltered_Z), 
                   names_to = "VInames", values_to = "VIvalueZ") %>% 
      mutate_at(vars(VInames), factor) %>% 
      mutate(VInames = fct_relevel(VInames, 
                                   "NDVI_raw_loeFiltered_Z", 
                                   "EVI2_raw_loeFiltered_Z", 
                                   "PPI_raw_loeFiltered_Z"))
  } else {
    data <- VItimesat0 %>% filter(validate == variable) %>%
      group_by(ID, LC) %>% 
      summarise(NDVI = cor(Gvar_Z, NDVI_M_loeFiltered_Z,
                           method = "pearson",
                           use = "pairwise.complete.obs"),
                EVI2 = cor(Gvar_Z, EVI2_M_loeFiltered_Z, 
                           method = "pearson",
                           use = "pairwise.complete.obs"),
                PPI = cor(Gvar_Z, PPI_M_loeFiltered_Z, 
                          method = "pearson",
                          use = "pairwise.complete.obs")) %>% 
      pivot_longer(cols = c(NDVI, EVI2, PPI), 
                   names_to = "var", values_to = "cor") %>% 
      mutate(var = factor(var, levels = c("NDVI", "EVI2", "PPI")))
    
    VIselectZ <- VItimesat0 %>% filter(validate == variable) %>% 
      select(Gvar_Z, LC,
             NDVI_M_loeFiltered_Z, 
             EVI2_M_loeFiltered_Z, 
             PPI_M_loeFiltered_Z) %>% 
      pivot_longer(cols = c(NDVI_M_loeFiltered_Z, 
                            EVI2_M_loeFiltered_Z, 
                            PPI_M_loeFiltered_Z), 
                   names_to = "VInames", values_to = "VIvalueZ") %>% 
      mutate_at(vars(VInames), factor) %>% 
      mutate(VInames = fct_relevel(VInames, 
                                   "NDVI_M_loeFiltered_Z", 
                                   "EVI2_M_loeFiltered_Z", 
                                   "PPI_M_loeFiltered_Z"))
  }
  
  boxLC <- ggplot(data, aes(x = var, y = cor)) +
    geom_boxplot(aes(fill = var), width = 0.4, outlier.shape = NA, size = 0.4) +
    geom_point(size = 1.2) +
    # geom_text(aes(label = ID))+
    # geom_line(aes(group = ID), color = "grey40", size = 0.4) +
    scale_fill_manual(values = c("#FE6100", "#648FFF", "#DC267F")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    # facet_wrap(vars(LC)) +
    theme_bw(base_size = 15) +
    labs(x = "", y = "Temporal correlation per site") +
    theme(legend.position = "top",
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          panel.grid = element_blank(),
          title = element_text(color = "black"),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))
  
  if (variable == "GPP") {
    level_key <- c("Mixed forest" = "Mixed forest (5)",
                   "Coniferous forest" = "Evergreen forest (15)",
                   "Broad-leaved forest" = "Decidous forest (6)",
                   "Grassland" = "Grassland (4)",
                   "Wetland" = "Wetland (5)",
                   "Agriculture" = "Agriculture (14)")
    withn <- function(x) recode(x, !!!level_key)
  } else {
    level_key <- c("Mixed forest" = "Mixed forest (1)",
                   "Coniferous forest" = "Evergreen forest (7)",
                   "Broad-leaved forest" = "Decidous forest (7)",
                   "Grassland" = "Grassland (8)",
                   "Wetland" = "Wetland (5)",
                   "Agriculture" = "Agriculture (5)")
    withn <- function(x) recode(x, !!!level_key)
  }
  boxLC <- boxLC +
    facet_wrap(vars(LC), labeller = as_labeller(withn))
  
  #########
  scaterLC <- ggplot(VIselectZ, aes(x = VIvalueZ, y = Gvar_Z)) +
    geom_point(aes(color = LC), size = 0.2, alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey40")+
    # annotate("text", x = 2.6, y = 3.0, angle = 52, label = "1:1") +
    geom_smooth(aes(color = LC), method = lm, se = F, size = 0.4) +
    xlim(-2, 3) +
    ylim(-2, 3) +
    scale_colour_brewer(palette = "Set1") +
    # stat_poly_eq(formula = y~x, aes(label = ..rr.label..), size = 4.2, parse = TRUE) +
    facet_wrap(vars(VInames), labeller = as_labeller(function(s) str_extract(s, "[A-Z2]+"))) +
    theme_bw(base_size = 15) +
    labs(x = "VI Z-score", y = paste(variable, "Z-score", sep = " ")) +
    theme(legend.position = "top",
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          panel.grid = element_blank(),
          title = element_text(color = "black"),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))+
    stat_cor(aes(label = ..r.label..), method = "pearson", size = 4.2,
             label.x.npc = "left", label.y.npc = "top")
  
  
  VIselectZbar <- VIselectZ %>% 
    group_by(LC, VInames) %>% 
    summarise(cor = cor(Gvar_Z, VIvalueZ, 
                        method = "pearson",
                        use = "pairwise.complete.obs"))
  
  barLC <- ggplot(VIselectZbar, aes(x = LC, y = cor, fill = VInames)) +
    geom_col(position = "dodge") +
    scale_y_continuous(limits = c(-0.2,1.05), expand = c(0,0)) +
    geom_text(aes(label = round(cor, 2)), 
              position = position_dodge(width = 1), vjust = -0.5) +
    scale_fill_manual(values = c("#FE6100", "#648FFF", "#DC267F"), 
                      labels = c("NDVI", "EVI2", "PPI")) +
    theme_bw(base_size = 15) +
    labs(x = "", y = "Spatiotemperal correlation") +
    theme(legend.position = "top",
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          title = element_text(color = "black"),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))
  if (variable == "GPP") {
    ggarrange(boxLC, 
              ggarrange(scaterLC, barLC, 
                        labels = c("b", "c"), 
                        ncol = 1,
                        font.label = list(size = 15)),
              ncol = 2, 
              labels = "a", 
              font.label = list(size = 15))
  } else {
    ggarrange(boxLC, 
              ggarrange(scaterLC, barLC, 
                        labels = c("e", "f"), 
                        ncol = 1,
                        font.label = list(size = 15)),
              ncol = 2, 
              labels = "d", 
              font.label = list(size = 15))
  }
  
}


ggarrange(VIcomp("GPP", "raw"), VIcomp("GCC", "raw"), ncol = 1)

ggsave("figures/figure_trajectory_consistency_raw.pdf",
       width = 35, height = 40, units = "cm")

ggarrange(VIcomp("GPP", "M"), VIcomp("GCC", "M"), ncol = 1)

ggsave("figures/figure_trajectory_consistency_NBAR.pdf",
       width = 35, height = 40, units = "cm")

