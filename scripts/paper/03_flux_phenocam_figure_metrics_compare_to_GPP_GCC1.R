library(tidyverse)

load("data/RData/VIpheno_flux_phenocam.RData")


VIpheno0 <- VIpheno %>% 
  filter(between(year, 2018, 2019), !(str_detect(ID, "^[:upper:]") & year == 2019))

brdfm <- "M"
vali <- "GPP"
phenoData <- function(brdfm, vali){
  GphenoVI <- 
    VIpheno0 %>% filter(BRDF == brdfm, validate == vali, VIname != vali) %>% 
    select(-c(BRDF, validate)) %>% 
    rename(thres_VI = thres) %>% 
    pivot_longer(cols = c(SOS, EOS),
                 names_to = "metric",
                 values_to = "VI")
  Gpheno <-
    VIpheno0 %>% filter(BRDF == brdfm, VIname == vali) %>% 
    select(-c(BRDF, validate, VIname)) %>%
    rename(thres_Gvar = thres) %>% 
    pivot_longer(cols = c(SOS, EOS),
                 names_to = "metric",
                 values_to = "Gvar") %>% 
    left_join(GphenoVI, by = c("ID", "LC", "metric", "year")) %>% 
    
    ### only reserve the closet match for each year, in case there are 
    ### double growing season in eithor GPP or VI
    mutate(closet = abs(Gvar - VI)) %>% 
    group_by(ID, LC, thres_Gvar, thres_VI, metric, VIname, year) %>% 
    slice(which.min(closet)) %>% 
    filter(closet < 210) %>% # to remove peaking year differ between VI and GPP/GCC, located in December and January.
    ungroup() %>% 
    
    mutate(VIname = factor(VIname, levels = c("NDVI", "EVI2", "PPI")),
           metric = factor(metric, levels = c("SOS", "EOS"))) %>% 
    mutate(diff = abs(Gvar-VI),
           thres_VI = as.numeric(str_extract(thres_VI, "\\d{2}")),
           thres_Gvar = paste0(str_extract(thres_Gvar, "\\d{2}"), "%"),
           validate = vali)
  return(Gpheno)
}
phenoThres <- function(brdfm) {
  data <- bind_rows(phenoData(brdfm, "GPP"), phenoData(brdfm, "GCC")) %>% 
    mutate(validate = factor(validate, levels = c("GPP", "GCC"))) %>%
    
    group_by(thres_Gvar, metric, thres_VI, VIname, validate) %>% 
    summarise(meanDiff = mean(diff, na.rm = T),
              cor = cor(Gvar, VI, method = "pearson",
                        use = "pairwise.complete.obs")) %>%
    ungroup() %>% 
    group_by(VIname, metric, thres_VI, validate) %>% 
    slice(which.min(abs(meanDiff))) %>% 
    ungroup() %>% 
    pivot_longer(cols = c(meanDiff, cor),
                 names_to = "measure",
                 values_to = "val") %>% 
    mutate(thres_Gvar = ifelse(measure == "cor", NA, thres_Gvar))
    
  plotR <- data %>% filter(measure == "cor") %>% 
    ggplot(aes(x = thres_VI, color = VIname)) +
    geom_point(aes(y = val), size = 1.8) +
    geom_line(aes(y = val), size = .6) +
    scale_x_continuous(breaks = seq(10,50,10))+
    scale_color_manual("VI", values = c("#FE6100", "#648FFF", "#DC267F")) +
    scale_y_continuous(expression(paste(italic("R")))) +
    facet_nested(. ~ validate + metric) +
    theme_bw(base_size = 17)+
    labs(x = "VI amplitude threshold (%)", title = "Spatial correlation") +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 17),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.2),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))
  
  plotB <- data %>% filter(measure == "meanDiff") %>% 
    ggplot(aes(x = thres_VI, color = VIname)) +
    geom_point(aes(y = val), size = 1.8) +
    geom_line(aes(y = val), size = .6) +
    scale_x_continuous(breaks = seq(10,50,10))+
    scale_color_manual("VI", values = c("#FE6100", "#648FFF", "#DC267F")) +
    scale_y_continuous("Day") +
    facet_nested(. ~ validate + metric) +
    theme_bw(base_size = 17)+
    labs(x = "VI amplitude threshold (%)", title = "Mean absolute bias") +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 17),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.2),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))
  
  ggarrange(plotB, plotR,
            labels = c("a", "b"),
            common.legend = T,
            ncol = 1, nrow = 2,
            font.label = list(size = 19))
}

phenoThres("M")
ggsave("figures/figure_phenometrics_compare_to_GPP_GCC_NBAR1.pdf",
       width = 25, height = 25, units = "cm")

phenoThres("raw")
ggsave("figures/figure_phenometrics_compare_to_GPP_GCC_raw1.pdf",
       width = 25, height = 25, units = "cm")
