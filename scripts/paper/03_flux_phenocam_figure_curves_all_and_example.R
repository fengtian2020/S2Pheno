library(tidyverse)
library(lubridate)

load("data/RData/VItimesat_flux_phenocam.RData")

#time series comparison between GPP/GCC and VIs
VItimesat1 <- VItimesat %>% 
  filter(between(date, ymd("2017-10-01"), ymd("2019-12-31")), 
         !(str_detect(ID, "^[:upper:]") & year(date) >= 2019))

VItimesatPointZ <- VItimesat1 %>% 
  select(ID, LC ,date, validate, 
         Gvar_nbar = Gvar_raw_Z,
         Gvar_toc = Gvar_raw_Z,
         NDVI_nbar = NDVI_M_loeFiltered_Z,
         EVI2_nbar = EVI2_M_loeFiltered_Z, 
         PPI_nbar = PPI_M_loeFiltered_Z,
         NDVI_toc = NDVI_raw_loeFiltered_Z, 
         EVI2_toc = EVI2_raw_loeFiltered_Z, 
         PPI_toc = PPI_raw_loeFiltered_Z) %>% 
  pivot_longer(cols = -c(ID, LC, date, validate), 
               names_to = c("name", "BRDF"), 
               names_sep = "_",
               values_to = "pointZ") %>%
  mutate_if(is.character, factor)

VItimesatZ <- VItimesat1 %>% 
  select(ID, LC, date, validate,
         Gvar_nbar = Gvar_Z,
         Gvar_toc = Gvar_Z,
         NDVI_nbar = NDVI_M_loeFiltered_fit_Z,
         EVI2_nbar = EVI2_M_loeFiltered_fit_Z, 
         PPI_nbar = PPI_M_loeFiltered_fit_Z,
         NDVI_toc = NDVI_raw_loeFiltered_fit_Z,
         EVI2_toc = EVI2_raw_loeFiltered_fit_Z, 
         PPI_toc = PPI_raw_loeFiltered_fit_Z) %>% 
  pivot_longer(cols = -c(ID, LC, date, validate), 
               names_to = c("name", "BRDF"), 
               names_sep = "_",
               values_to = "lineZ") %>%
  mutate_if(is.character, factor) %>% 
  full_join(VItimesatPointZ) %>% 
  mutate(BRDF = factor(BRDF, levels = c("nbar", "toc"), 
                       labels = c("NBAR", "TOC")))

########
tsplot <- function(df, lc, co. = co) {
  ggplot(data = df %>% filter(LC == lc, BRDF == "NBAR"),
         aes(x = date, col = name)) +
    geom_point(aes(y = pointZ), size = 0.3, show.legend = T) +
    geom_line(aes(y = lineZ), size = 0.5, show.legend = T) +
    scale_colour_manual(values = c("grey20", "#FE6100", "blue", "#DC267F")) +
    scale_x_date(date_labels = "%b") +
    facet_wrap(vars(ID), ncol = co., scales = "free_y") +
    labs(y = "Z-score", x = "", title = lc) +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = rel(1.1)),
          legend.position = "top", 
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(size = 0.4),
          strip.background = element_rect(size = 0.4))
}

VItimesatZGCC <- VItimesatZ %>% filter(validate == "GCC") %>% 
  mutate(name = factor(name, levels = c("Gvar", "NDVI", "EVI2", "PPI"),
                       labels = c("GCC", "NDVI", "EVI2", "PPI"))) %>%
  arrange(name)
co <- 3
LCbf <- tsplot(VItimesatZGCC, "Broad-leaved forest", co)
LCcf <- tsplot(VItimesatZGCC, "Coniferous forest", co)
LCmf <- tsplot(VItimesatZGCC, "Mixed forest", co)
LCa <- tsplot(VItimesatZGCC, "Agriculture", co)
LCg <- tsplot(VItimesatZGCC, "Grassland", co)
LCw <- tsplot(VItimesatZGCC, "Wetland", co)

library(ggpubr)
ggarrange(LCbf, LCcf, LCmf, LCa, LCg, LCw, 
          labels = c("a", "b", "c", "d", "e", "f"), 
          font.label = list(size = 17),
          common.legend = F, heights = c(1.5, 1.5, 0.8, 1.1, 1.5, 1.1),
          ncol = 1)
ggsave("figures/figure_GCC_curves_NBAR5.pdf",
       width = 22, height = 80, units = "cm")


VItimesatZGPP <- VItimesatZ %>% filter(validate == "GPP") %>% 
  mutate(name = factor(name, levels = c("Gvar", "NDVI", "EVI2", "PPI"),
                       labels = c("GPP", "NDVI", "EVI2", "PPI"))) %>%
  arrange(name)

co <- 4
LCbf <- tsplot(VItimesatZGPP, "Broad-leaved forest", co)
LCcf <- tsplot(VItimesatZGPP, "Coniferous forest", co)
LCmf <- tsplot(VItimesatZGPP, "Mixed forest", co)
LCa <- tsplot(VItimesatZGPP, "Agriculture", co)
LCg <- tsplot(VItimesatZGPP, "Grassland", co)
LCw <- tsplot(VItimesatZGPP, "Wetland", co)

ggarrange(LCbf, LCcf, LCmf, LCa, LCg, LCw, 
          labels = c("a", "b", "c", "d", "e", "f"), 
          font.label = list(size = 17),
          common.legend = F, heights = c(1.25, 2, 1.25, 2, 0.85, 1.25),
          ncol = 1)
ggsave("figures/figure_GPP_curves_NBAR5.pdf",
       width = 22, height = 80, units = "cm")




# example curves at DE-HoH site -------------------------------------------

VItimesatZ %>%
  filter(str_detect(ID, "DE-HoH"), BRDF == "NBAR") %>% 
  mutate(across(c(validate, name), as.character)) %>%
  mutate(name = ifelse(name == "Gvar", "GPP / GCC", name),
         name = factor(name, levels = c("GPP / GCC", "NDVI", "EVI2", "PPI"))) %>%
  mutate(validate = factor(validate, levels = c("GPP", "GCC"))) %>%
  
  ggplot(aes(x = date, col = name)) +
  geom_point(aes(y = pointZ), size = 1, show.legend = T) +
  geom_line(aes(y = lineZ), size = 0.5, show.legend = T) +
  # geom_text(aes(y = pointZ, label = date)) +
  scale_colour_manual(values = c("grey20", "#FE6100", "blue", "#DC267F")) +
  scale_size_manual(values = c(0.8, 0.8, 0.8, 0.8))+
  scale_x_date(date_labels = "%b %Y") +
  facet_wrap(vars(validate), ncol = 1, scales = "free_y") +
  labs(y = "Z-score", x = "") +
  theme_bw(base_size = 17) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        legend.position = "top", 
        legend.text = element_text(size = 17),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggsave("figures/figure_VI_GPP_GCC_example_curves_NBAR5.pdf",
       width = 20, height = 22, units = "cm")




# GCC filtered pixels -----------------------------------------------------


VItimesatZ %>% filter(validate == "GCC") %>% 
  filter(ID %in% c("monteblanco_SH_1000", "monteblanco_EB_1000", "montenegro_SH_1000",
                   "donanafuenteduque_WL_1000", "DE-Geb_AG_1000", "borgocioffinorth_AG_1000")) %>% 
  mutate(name = factor(name, levels = c("Gvar", "NDVI", "EVI2", "PPI"),
                       labels = c("GCC", "NDVI", "EVI2", "PPI"))) %>%
  arrange(name) %>% 
  filter(BRDF == "NBAR") %>% 
  ggplot(aes(x = date, col = name)) +
  geom_point(aes(y = pointZ), size = 0.3, show.legend = T) +
  geom_line(aes(y = lineZ), size = 0.5, show.legend = T) +
  scale_colour_manual(values = c("grey20", "#FE6100", "blue", "#DC267F")) +
  scale_x_date(date_labels = "%b") +
  facet_wrap(vars(ID), ncol = 3, scales = "free_y") +
  labs(y = "Z-score", x = "") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        legend.position = "top", 
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggsave("figures/figure_VI_GCC_excluded.pdf",
       width = 25, height = 15, units = "cm")



# SE-Svb -----------------------------------------------------


VItimesatZ %>% filter(validate == "GPP") %>% 
  filter(ID %in% c("SE-Svb")) %>% 
  mutate(name = factor(name, levels = c("Gvar", "NDVI", "EVI2", "PPI"),
                       labels = c("GPP", "NDVI", "EVI2", "PPI"))) %>%
  arrange(name) %>% 
  filter(BRDF == "NBAR") %>% 
  ggplot(aes(x = date, col = name)) +
  geom_point(aes(y = pointZ)) +
  geom_line(aes(y = lineZ), size = 0.5, show.legend = T) +
  scale_colour_manual(values = c("grey20", "#FE6100", "blue", "#DC267F")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %d") +
  labs(y = "Z-score", x = "Year 2017−2018") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        legend.position = "top", 
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggsave("figures/figure_VI_SE_Svb.pdf",
       width = 25, height = 15, units = "cm")

VItimesatPoint <- VItimesat1 %>% 
  select(ID, LC ,date, validate, 
         RED_nbar = NBAR_Red_M_loeFilter,
         NIR_nbar = NBAR_NIR_M_loeFilter,
         Gvar_nbar = Gvar_raw,
         NDVI_nbar = NDVI_M_loeFiltered,
         EVI2_nbar = EVI2_M_loeFiltered, 
         PPI_nbar = PPI_M_loeFiltered) %>% 
  pivot_longer(cols = -c(ID, LC, date, validate), 
               names_to = c("name", "BRDF"), 
               names_sep = "_",
               values_to = "point") %>%
  mutate_if(is.character, factor)



# VItimesatPoint %>% drop_na() %>% 
#   mutate(ndvi = (NIR_nbar - RED_nbar) / (NIR_nbar + RED_nbar)) %>% 
#   ggplot(aes(x = NDVI_nbar, y = ndvi)) +
#   geom_point()


VItimesat <- VItimesat1 %>% 
  select(ID, LC, date, validate,
         Gvar_nbar = Gvar,
         RED_nbar = NBAR_Red_M_loeFilter,
         NIR_nbar = NBAR_NIR_M_loeFilter,
         NDVI_nbar = NDVI_M_loeFiltered_fit,
         EVI2_nbar = EVI2_M_loeFiltered_fit, 
         PPI_nbar = PPI_M_loeFiltered_fit) %>% 
  pivot_longer(cols = -c(ID, LC, date, validate), 
               names_to = c("name", "BRDF"), 
               names_sep = "_",
               values_to = "line") %>%
  mutate_if(is.character, factor) %>% 
  full_join(VItimesatPoint) %>% 
  mutate(BRDF = factor(BRDF, levels = c("nbar", "toc"), 
                       labels = c("NBAR", "TOC")))



VItimesat %>% filter(validate == "GPP") %>% 
  # filter(between(date, ymd("2018-07-31"), ymd("2018-08-03"))) %>% 
  filter(ID %in% c("SE-Svb"), name != "Gvar") %>% 
  mutate(name = factor(name, levels = c("NDVI", "EVI2", "PPI", "RED", "NIR"))) %>%
  arrange(name) %>% 
  filter(BRDF == "NBAR") %>% 
  ggplot(aes(x = date, col = name)) +
  geom_point(aes(y = point), show.legend = T) +
  # geom_text(aes(y = point, label = date)) +
  # geom_line(aes(y = lineZ), size = 0.5, show.legend = T) +
  scale_colour_manual(values = c("#FE6100", "blue", "#DC267F", "grey40", "black")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %d") +
  labs(y = "VI & reflectance", x = "Year 2017-2018") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        legend.position = "top", 
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggsave("figures/figure_VI_SE_Svb1.pdf",
       width = 25, height = 15, units = "cm")

VItimesat %>% filter(validate == "GPP") %>% 
  filter(between(date, ymd("2018-08-10"), ymd("2018-9-10"))) %>%
  filter(ID %in% c("SE-Svb"), name != "Gvar") %>% 
  mutate(name = factor(name, levels = c("NDVI", "EVI2", "PPI", "RED", "NIR"))) %>%
  arrange(name) %>% 
  filter(BRDF == "NBAR") %>% drop_na() %>% 
  ggplot(aes(x = date, y = point, col = name)) +
  geom_point(aes(size = name)) +
  # geom_text(aes(y = point, label = date)) +
  # geom_line(aes(y = lineZ), size = 0.5, show.legend = T) +
  geom_smooth(method = "lm", se = F, size= 0.5) +
  scale_colour_manual(values = c("#FE6100", "blue", "#DC267F", "black", "grey50")) +
  scale_size_manual(values = c(2,2,2,1.2,1.2)) +
  # scale_x_date(date_breaks = "2 month", date_labels = "%b %d") +
  labs(y = "VI & reflectance", x = "") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        legend.position = "top", 
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.4),
        strip.background = element_rect(size = 0.4))

ggsave("figures/figure_VI_SE_Svb2.pdf",
       width = 13, height = 15, units = "cm")
