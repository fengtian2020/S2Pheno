library(tidyverse)
library(lubridate)

###############################
# example time series of one pep site
load("data/RData/pepVI4kmMean.RData")

write_csv(pepVI %>% filter(str_detect(ID, "8212")) %>% 
            select(ID, date,PPI_raw_loeFiltered),
          "data/2_timesat_inputs/PEP_S2_VI_SCL45_LC_4km_mean_site_8212.csv")

timesat <- read_csv("data/3_timesat_outputs/smooth_time_series/smoothed_PPI_raw_loeFiltered_PEP_S2_VI_SCL45_LC_4km_mean_site_8212_DL-SP.csv") %>% 
  rename(PPI_fit = PPI_raw_loeFiltered)


pepEX <- pepVI %>% filter(str_detect(ID, "8212")) %>% 
  left_join(timesat) %>% 
  separate(ID, c("type", "site"), sep = "_") %>% 
  select(type, date, PPI = PPI_raw_loeFiltered, PPI_fit) %>% 
  filter(between(date, ymd("2019-01-01"), ymd("2019-12-31")))

pepEX %>% 
  ggplot(aes(x = date, color = type)) +
  geom_point(aes(y = PPI)) +
  geom_line(aes(y = PPI_fit)) +
  scale_color_brewer(palette = "Dark2")+
  theme_bw(base_size = 20) +
  labs(x = "")+
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.grid = element_blank()) 

ggsave("figures/figure_example_PEP_site_curves.pdf",
       width = 18, height = 18, units = "cm")

