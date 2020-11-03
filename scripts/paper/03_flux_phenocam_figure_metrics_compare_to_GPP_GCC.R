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
  bind_rows(phenoData(brdfm, "GPP"), phenoData(brdfm, "GCC")) %>% 
    mutate(validate = factor(validate, levels = c("GPP", "GCC"))) %>%
    
    group_by(thres_Gvar, metric, thres_VI, VIname, validate) %>% 
    summarise(meanDiff = mean(diff, na.rm = T),
              sdDiff = sd(diff, na.rm = T),
              cor = cor(Gvar, VI, method = "pearson",
                        use = "pairwise.complete.obs")*100) %>%
    ungroup() %>% 
    group_by(VIname, metric, thres_VI, validate) %>% 
    slice(which.min(abs(meanDiff))) %>% 
    ungroup() %>% 
    pivot_longer(cols = c(meanDiff, cor),
                 names_to = "measure",
                 values_to = "val") %>% 
    mutate(thres_Gvar = ifelse(measure == "cor", NA, thres_Gvar)) %>% 
    
    ggplot(aes(x = thres_VI, color = VIname)) +
    geom_vline(xintercept = 25, color = "grey60", linetype = 2)+
    geom_point(aes(y = val, shape = measure), size = 2.5) +
    geom_line(aes(y = val, linetype = measure), size = .8) +
    geom_text(aes(y = val, label = thres_Gvar),
              position = position_nudge(y = -2))+
    scale_x_continuous(breaks = seq(5,50,5))+
    # geom_point(aes(y = cor*100)) +
    # geom_line(aes(y = cor*100), linetype = 2) +
    scale_color_manual("VI", values = c("tan3", "blue", "firebrick3")) +
    scale_shape_manual(values = c(1, 16),
                       labels = c("Spatial correlation", "Mean absolute bias")) +
    scale_linetype_manual(values = c(2, 1),
                          labels = c("Spatial correlation", "Mean absolute bias")) +
    scale_y_continuous("Mean absolute bias (day)", breaks = seq(0, 80, 20),
                       sec.axis = 
                         sec_axis(~./100,
                                  breaks = seq(0, .8, .2),
                                  name = expression(paste("Spatial correlation (",italic("R"),")", sep = ""))))+
    # facet_wrap(vars(metric)) +
    facet_grid(rows = vars(metric), cols = vars(validate), scales = "free")+
    theme_bw(base_size = 17)+
    labs(x = "VI amplitude threshold (%)") +
    theme(panel.grid = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 17))
}

phenoThres("M")
ggsave("figures/figure_phenometrics_compare_to_GPP_GCC_NBAR.pdf",
       width = 25, height = 25, units = "cm")

phenoThres("raw")
ggsave("figures/figure_phenometrics_compare_to_GPP_GCC_raw.pdf",
       width = 25, height = 25, units = "cm")
