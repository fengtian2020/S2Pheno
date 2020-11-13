
library(tidyverse)

load("data/RData/pepS2report.RData")

# plot scatterplots -------------------------------------------------------

pepS2f3 <- pepS2 %>% 
  mutate(foroutlier = PEPdoy - S2doy) %>% 
  group_by(VIname, phasePair, thres, type, BRDF) %>% 
  mutate(across(foroutlier, ~ ifelse(abs(. - mean(., na.rm = T)) > sd(., na.rm = T) * 3, NA, 1))) %>% 
  filter(!is.na(foroutlier)) %>% 
  
  mutate(foroutlier = PEPdoy - S2doy) %>% 
  group_by(VIname, phasePair, thres, type, BRDF) %>% 
  mutate(across(foroutlier, ~ ifelse(abs(. - mean(., na.rm = T)) > sd(., na.rm = T) * 3, NA, 1))) %>% 
  filter(!is.na(foroutlier)) %>% 
  
  mutate(foroutlier = PEPdoy - S2doy) %>% 
  group_by(VIname, phasePair, thres, type, BRDF) %>% 
  mutate(across(foroutlier, ~ ifelse(abs(. - mean(., na.rm = T)) > sd(., na.rm = T) * 3, NA, 1))) %>% 
  filter(!is.na(foroutlier)) 

pepS2f3 %>%
  filter(BRDF == "raw", thres == "25%", VIname == "PPI") %>% #, phasePair == "Leaf unfolded 50%"
  select(-c(thres, BRDF, year)) %>%
  group_by(VIname, type) %>%
  ggplot(aes(x = S2doy, y = PEPdoy)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(method = "lm", size = 0.8) +
  facet_grid(type ~ phasePair, scales = "free")+
  labs(x = "VI-derived SOS/EOS (day of year)",
       y = "PEP725 phenophase (day of year)") +
  theme_bw(base_size = 17)+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 17))
ggsave("figures/report/figure_metrics_compare_to_PEP_scatterplots_outlierfiltered_raw.pdf",
       width = 35, height = 20, units = "cm")

# amplitude threshold -----------------------------------------------------


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
ggsave("figures/report/figure_metrics_compare_to_PEP_raw_amplitude.pdf",
       width = 40, height = 25, units = "cm")

pepcmb("M", "diff")
ggsave("figures/report/figure_metrics_compare_to_PEP_NBAR_amplitude.pdf",
       width = 35, height = 20, units = "cm")



# fixed threshold ---------------------------------------------------------

# load("data/RData/pepS2reportfix.RData")
# 
# pepcmb <- function(brdf) {
#   diffdf <- pepS2 %>% 
#     filter(BRDF == brdf) %>% 
#     select(-c(BRDF, year)) %>%
#     group_by(VIname, phasePair, thres, type) %>% 
#     summarise(diff = mean(PEPdoy - S2doy, na.rm = T)) %>% 
#     ungroup() %>% 
#     mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) 
#   
#   cordf <- pepS2 %>% 
#     filter(BRDF == brdf) %>% 
#     select(-c(BRDF, year)) %>%
#     group_by(VIname, phasePair, thres, type) %>% 
#     summarise(cor = cor(PEPdoy, S2doy, 
#                         use = "pairwise.complete.obs")^2 * 500 - 150) %>% 
#     ungroup() %>% 
#     mutate(thres = as.numeric(str_extract(thres, "\\d{2}")))
#   
#   stddf <- pepS2 %>% 
#     filter(BRDF == brdf) %>% 
#     select(-c(BRDF, year)) %>%
#     group_by(VIname, phasePair, thres, type) %>% 
#     summarise(std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
#     ungroup() %>% 
#     mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) %>% 
#     mutate(indicator = "diff")
#   
#   mainn <- left_join(cordf, diffdf) %>%  #, by = c("VIname", "phasePair", "thres", "type")
#     pivot_longer(cols = c(cor, diff),
#                  names_to = "indicator",
#                  values_to = "value") %>% 
#     left_join(stddf)
#   
#   plot <- mainn %>% 
#     ggplot(aes(x = thres, y = value, color = VIname)) +
#     
#     geom_vline(xintercept = 25, linetype = 2, color = "grey60", size = 0.3)+
#     geom_hline(yintercept = 0, linetype = 2, color = "grey60", size = 0.3)+
#     
#     geom_pointrange(aes(size = indicator, shape = indicator,
#                         ymin = value - std, ymax = value + std),
#                     fatten = 6,
#                     position=position_dodge(width = 6))+
#     scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
#     scale_size_manual(values = c(0.2, 0.15), guide = F) +
#     scale_shape_manual(values = c(1, 16),
#                        labels = c(expression(italic("R"^"2")),
#                                   "Time difference")) +
#     
#     scale_y_continuous("Time difference (day)",breaks = c(-100, 0, 100),
#                        sec.axis = 
#                          sec_axis(~(. + 150) / 500, breaks = c(0, 0.2, 0.4),
#                                   name = expression(italic("R"^"2"))))+
#     facet_grid(type ~ phasePair) +
#     labs(x = "VI amplitude threshold (%)") +
#     theme_bw(base_size = 17)+
#     theme(panel.grid = element_blank(),
#           legend.title = element_blank(),
#           legend.position = "top",
#           legend.text = element_text(size = 17))
#   return(plot)
# }
# 
# pepcmb("M")
# ggsave("figures/report/figure_metrics_compare_to_PEP_NBAR_fix.pdf",
#        width = 35, height = 20, units = "cm")
# 
# 
# pepcmb("raw")
# ggsave("figures/report/figure_metrics_compare_to_PEP_raw_fix.pdf",
#        width = 35, height = 20, units = "cm")

    
  






