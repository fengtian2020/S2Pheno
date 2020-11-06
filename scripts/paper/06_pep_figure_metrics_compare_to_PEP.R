
library(tidyverse)

load("data/RData/pepS2.RData")


# plot scatterplots -------------------------------------------------------

brdf <- "M"

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
  filter(BRDF == "M", thres == "25%", VIname == "PPI") %>% #, phasePair == "Leaf unfolded 50%"
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
ggsave("figures/figure_metrics_compare_to_PEP_scatterplots_outlierfiltered_NBAR.pdf",
       width = 35, height = 20, units = "cm")


# parameters of R2 and mean bias ------------------------------------------

pepcmb <- function(brdf) {
  diffdf <- pepS2f3 %>% 
    filter(BRDF == brdf) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, phasePair, thres, type) %>% 
    summarise(diff = mean(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) 
  
  s <- 300
  
  cordf <- pepS2f3 %>% 
    filter(BRDF == brdf) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, phasePair, thres, type) %>% 
    summarise(cor = cor(PEPdoy, S2doy, 
                        use = "pairwise.complete.obs")^2 * s - 150) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}")))
  
  stddf <- pepS2f3 %>% 
    filter(BRDF == brdf) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, phasePair, thres, type) %>% 
    summarise(std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) %>% 
    mutate(indicator = "diff")
  
  mainn <- left_join(cordf, diffdf) %>%  #, by = c("VIname", "phasePair", "thres", "type")
    pivot_longer(cols = c(cor, diff),
                 names_to = "indicator",
                 values_to = "value") %>% 
    left_join(stddf)
  
  plot <- mainn %>% 
    ggplot(aes(x = thres, y = value, color = VIname)) +
    
    geom_vline(xintercept = 25, linetype = 2, color = "grey60", size = 0.3)+
    geom_hline(yintercept = 0, linetype = 2, color = "grey60", size = 0.3)+
    
    geom_pointrange(aes(size = indicator, shape = indicator,
                        ymin = value - std, ymax = value + std),
                    fatten = 6,
                    position=position_dodge(width = 6))+
    scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
    scale_size_manual(values = c(0.2, 0.15), guide = F) +
    scale_shape_manual(values = c(1, 16),
                       labels = c(expression(italic("R"^"2")),
                                  "Time difference")) +
    
    scale_y_continuous("Time difference (day)",breaks = c(-100, 0, 100),
                       sec.axis = 
                         sec_axis(~(. + 150) / s, #breaks = c(0, 0.2, 0.4),
                                  name = expression(italic("R"^"2"))))+
    facet_grid(type ~ phasePair) +
    labs(x = "VI amplitude threshold (%)") +
    theme_bw(base_size = 17)+
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 17))
  return(plot)
}

pepcmb("M")
ggsave("figures/figure_metrics_compare_to_PEP_outlierfiltered_NBAR.pdf",
       width = 35, height = 20, units = "cm")


pepcmb("raw")
ggsave("figures/figure_metrics_compare_to_PEP_outlierfiltered_raw.pdf",
       width = 35, height = 20, units = "cm")

  
  

    
  






