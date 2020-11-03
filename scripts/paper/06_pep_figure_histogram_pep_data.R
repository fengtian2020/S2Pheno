
library(tidyverse)
load("data/RData/pepS2.RData")

# statistics of pep records
pepS2 %>% filter(thres == "25%", BRDF == "raw", VIname == "NDVI", type != "Crop") %>% 
  ggplot(aes(PEPdoy, color = phasePair))+
  geom_freqpoly()+
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(vars(type)) +
  labs(x = "Day of year", y = "Count") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        # legend.key.size = unit(1, "cm"),
        legend.position = "top",
  ) 
ggsave("figures/figure_histogram_pep_data.pdf",
       width = 20, height = 10, units = "cm")


