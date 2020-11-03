
setwd("D:/S2Pheno")

library(tidyverse)
library(ggpubr)

files <- list.files("data/timesat/transects", full.names = T, recursive = T)

data <- NULL
for (i in files) {
  trans <- str_extract(i, "Scandic|Spain")
  brdf <- str_extract(i, "(?<=_)raw|M(?=_)")
  thres <- str_extract(i, "(?<=_)\\d\\.\\d{2,3}") %>% as.numeric()
  type <- str_extract(i, "fixed|amplitude")
  tmp <- read_csv(i) %>%
    transmute(ID, year,
              SOS = start_of_season,
              EOS = end_of_season,
              trans = trans,
              BRDF = brdf,
              thres = thres,
              type = type) %>%
    filter(between(year, 2019, 2019)) %>%
    mutate(lat = as.numeric(str_extract(ID, " \\d{2}.\\d{1,}")))
  data <- bind_rows(data, tmp)
}

demfiles <- list.files("data/S2sites", pattern = "*DEM.shp", full.names = T)
library(sf)
transDEM <- NULL
for (i in demfiles) {
  trans <- str_extract(i, "Scandic|Spain")
  tmp <- read_sf(i) %>%
    as_tibble() %>%
    mutate(ID = as.character(geometry),
           trans = trans) %>%
    select(ID, DEM = AVE_DSM, trans)
  transDEM <- bind_rows(transDEM, tmp)
}



datar <- data %>% 
  pivot_longer(cols = c(SOS, EOS),
               names_to = "metric", 
               values_to = "pheno") %>% 
  filter((type == "amplitude" & metric == "SOS" & thres == 0.25) |
           (type == "amplitude" & metric == "EOS" & thres == 0.15) |
           (type == "fixed" & metric == "SOS" & thres == 0.3) |
           (type == "fixed" & metric == "EOS" & thres == 0.12)) %>% 
  left_join(transDEM) %>% 
  mutate(metric = factor(metric, levels = c("SOS", "EOS"),
                         labels = c("SOS", "EOS")),
         trans = factor(trans, levels = c("Scandic", "Spain"),
                        labels = c("Northern Europe", "Southern Europe")),
         type = factor(type, labels = c("Amplitude threshold",
                                        "Fixed threshold")))


axisscale1 <- 8
north <- datar %>% filter(BRDF == "M", trans == "Northern Europe") %>% 
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) + #
  # scale_shape_manual("Land cover type", values = seq(0, 6)) +
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.15, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale1, color = "DEM"), size = 0.5, alpha = 0.1, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale1, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(0, 400),
                     sec.axis = sec_axis(trans = ~.*axisscale1, 
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(type), scales = "free") +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

axisscale2 <- 20
south <- datar %>% filter(BRDF == "M", trans == "Southern Europe") %>% 
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) +
  scale_shape_manual("Land cover type", values = seq(0, 6)) +
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale2, color = "DEM"), size = 0.5, alpha = 0.15, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale2, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(-100, 500), breaks = seq(0, 500, 100),
                     sec.axis = sec_axis(trans = ~.*axisscale2, 
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(type), scales = "free") +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

ggarrange(north, south, common.legend = T,
          labels = c("a", "b"), ncol = 1)

ggsave("figures/report/figure_transect_PPI_NBAR.pdf",
       width = 25, height = 25, units = "cm")





axisscale1 <- 8
north <- datar %>% filter(BRDF == "raw", trans == "Northern Europe") %>% 
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) + #
  # scale_shape_manual("Land cover type", values = seq(0, 6)) +
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.15, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale1, color = "DEM"), size = 0.5, alpha = 0.1, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale1, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(0, 400),
                     sec.axis = sec_axis(trans = ~.*axisscale1, 
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(type), scales = "free") +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

axisscale2 <- 20
south <- datar %>% filter(BRDF == "raw", trans == "Southern Europe") %>% 
  ggplot(aes(x = lat))+
  geom_point(aes(y = pheno, color = metric), size = 0.5, alpha = 0.3, shape = 1) +
  scale_shape_manual("Land cover type", values = seq(0, 6)) +
  geom_smooth(aes(y = pheno, color = metric), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  geom_point(aes(y = DEM/axisscale2, color = "DEM"), size = 0.5, alpha = 0.15, shape = 1) +
  geom_smooth(aes(y = DEM/axisscale2, color = "DEM"), method="loess", span = 0.1, size = 0.7, se = F,
              method.args = list(family = "symmetric")) +
  scale_color_manual("", breaks = c("SOS", "EOS", "DEM"),
                     values = c(RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2], "grey40")) +
  scale_y_continuous(name = "Day of year", limits = c(-100, 500), breaks = seq(0, 500, 100),
                     sec.axis = sec_axis(trans = ~.*axisscale2, 
                                         name = "Elevation (km)",
                                         breaks = seq(0, 4000, 1000),
                                         labels = seq(0, 4, 1))) +
  facet_wrap(vars(type), scales = "free") +
  theme_bw(base_size = 15) +
  labs(x = "Latitude")

ggarrange(north, south, common.legend = T,
          labels = c("a", "b"), ncol = 1)

ggsave("figures/report/figure_transect_PPI_raw.pdf",
       width = 25, height = 25, units = "cm")

