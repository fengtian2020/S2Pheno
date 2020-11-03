
library(tidyverse)

# amplitude threshold -----------------------------------------------------

load("data/RData/pepS2report.RData")

pepS2 %>% mutate(sel = paste0(thres, phasePair)) %>% #distinct(phasePair)
  filter(BRDF == "raw", VIname == "PPI", type != "Crop", 
         sel %in% c("25%Leaf unfolded", "25%Leaf emerged", "25%Leaf unfolded 50%",
                    "25%Fresh green 25%", "15%Leaf fallen 50%", "15%Coloration  50%")) %>% 
  select(type, phasePair, PEPdoy, S2doy) %>% 
  group_by(phasePair, type) %>% 
  summarise(diff = mean(PEPdoy - S2doy, na.rm = T),
            std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(phasePair = fct_rev(phasePair),
         upper = diff + std,
         lower = diff - std) %>% 
  
  ggplot(aes(x = diff, y = phasePair)) +
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin = lower, xmax = upper, color = type),
                  fatten = 4,
                  position=position_dodge(width = .3),
  )+
  scale_x_continuous(breaks = seq(-50, 70, 10))+
  scale_color_brewer(palette = "Dark2") +
  labs(x = "PEP doy  -  Sentinel2 doy", y = "") +
  theme_bw(base_size = 15) +
  theme(legend.title = element_blank(),
        legend.position = c(.82, .8),
  ) 
# stat_cor(aes(color = type, label = ..r.label..), method = "spearman", #size = 4.2,
#          label.x.npc = "left", label.y.npc = "top")

ggsave("Sen2Pheno/figures/PEP_PPI_amplitude_25_15_diff_3_report.pdf",
       width = 20, height = 10, units = "cm")