library(tidyverse)
library(lubridate)
library(ggpubr)

load("data/RData/VIpheno_flux_phenocam.RData")

VIpheno0 <- VIpheno %>% 
  filter(between(year, 2018, 2019), !(str_detect(ID, "^[:upper:]") & year == 2019))

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

phenoData("raw", "GPP") %>% 
  filter(thres_Gvar %in% c("15%", "25%", "45%"), thres_VI == 25, metric == "SOS", VIname == "PPI") %>% 
  group_by(thres_Gvar) %>% 
  summarise(`absolute bias` = mean(diff))

dat_bias <- data.frame(
  label = c("Mean absolute bias = 20.0 days", 
            "Mean absolute bias = 13.5 days", 
            "Mean absolute bias = 11.3 days"),
  thres_Gvar = c("15%", "25%", "45%")
)

phenoData("raw", "GPP") %>% 
  filter(thres_Gvar %in% c("15%", "25%", "45%"), thres_VI == 25, metric == "SOS", VIname == "PPI") %>% 
  ggplot(aes(x = VI, y = Gvar)) +
  geom_abline(slope = 1, intercept = 0, color = "grey70")+
  geom_point() +
  scale_x_continuous(limits = c(0, 270))+
  scale_y_continuous(limits = c(0, 270))+
  facet_wrap(vars(thres_Gvar)) +
  theme_bw(base_size = 14)+
  labs(x = "PPI SOS at threshold 25% (DOY)",
       y = "GPP SOS (DOY)") +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))+
  geom_text(data = dat_bias, mapping = aes(x = 0, y = 220, label = label),
            hjust = 0) +
  stat_cor(aes(label = ..r.label..), method = "pearson", #size = 4.2,
           label.x.npc = "left", label.y.npc = "top")
  
ggsave("figures/report/figure_VI_GPP_GCC_scatterplots_raw.pdf",
       width = 25, height = 10, units = "cm")
