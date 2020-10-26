library(sf)
library(tidyverse)
library(rnaturalearth)

load("data/RData/VItimesat_flux_phenocam.RData")
load("data/RData/VIpheno_flux_phenocam.RData")


# European-centric ETRS89 Lambert Azimuthal Equal-Area projection
euEqualArea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

europe <- ne_countries(continent = "europe", returnclass = "sf", scale = 50) %>%
  st_crop(xmin = -10, ymin = 34, xmax = 41, ymax = 71) %>%
  filter(admin != "Russia")

fluxsites <- VIpheno %>% filter(validate == "GPP") %>% distinct(ID) %>% pull()

flux <- read_sf("data/1_gee_raw_reflectance/flux/flux_drought2018_site_S2_NBAR_raw_201704_202003_100m_SCL45_mean.shp") %>% 
  st_centroid() %>% 
  as_tibble() %>%
  distinct(ID, .keep_all = T) %>%
  filter(ID %in% fluxsites) %>% 
  mutate(type = "FLUX",
         lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2]) %>% 
  select(ID, lon, lat, type)

phenocamsites <- VIpheno %>% filter(validate == "GCC") %>% distinct(ID) %>% pull()

phenocam <- read_sf("data/1_gee_raw_reflectance/phenocam/phenocam_S2_NBAR_raw_201704_202003_10m_SCL45.shp") %>%
  as_tibble() %>%
  mutate(ID = paste(site, veg_type, str_pad(roi, 4, pad = "0"), sep = "_")) %>% 
  distinct(ID, .keep_all = T) %>% 
  filter(ID %in% phenocamsites) %>% 
  mutate(type = "PhenoCam",
         lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2]) %>% 
  select(ID, lon, lat, type) 

pep <- read_csv("data/0_ground_data_raw/pep725/pep725_in_use.csv") %>%
  mutate(ID = paste0("PEP_", s_id)) %>%
  distinct(ID, lon, lat) %>%
  mutate(type = "PEP725")

# for the project
# set.seed(2020)
# pep <- read_csv("data/pep725/pep725_in_use.csv") %>%
#   mutate(ID = paste0("PEP_", s_id)) %>%
#   distinct(ID, lon, lat) %>%
#   mutate(type = "PEP725") %>%
#   sample_frac(0.7)
# 
# pep70 <- pep %>% pull(ID)
# pep30 <- read_csv("Sen2Pheno/PEP/pep725_in_use.csv") %>% 
#   mutate(ID = paste0("PEP_", s_id)) %>% 
#   distinct(ID, lon, lat) %>% 
#   filter(!(ID %in% c(pep70))) %>% 
#   separate(ID, sep = "_", c("p", "s_id")) %>% 
#   select(s_id, lon, lat, -p)
# write_csv(pep30, "data/pep725/PEP_site_for_valiation.csv")

transect0 <- read_sf("data/1_gee_extracted_data/transects/TransPoints_Scandic_S2_DEM.shp") %>%
  as_tibble() %>%
  mutate(type = "Transect",
         lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2],
         ID = paste0("tran_", lat)) %>% 
  distinct(lat, .keep_all = T) %>% 
  select(ID, lon, lat, type)

transect <- read_sf("data/1_gee_extracted_data/transects/TransPoints_Spain_S2_DEM.shp") %>%
  as_tibble() %>%
  mutate(type = "Transect",
         lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2],
         ID = paste0("tran_", lat)) %>% 
  distinct(lat, .keep_all = T) %>% 
  select(ID, lon, lat, type) %>% 
  bind_rows(transect0)

randomn <- read_sf("data/1_gee_extracted_data/randompoints/randomn_points_location.shp") %>% 
  as_tibble() %>% 
  mutate(type = "Random points",
         lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2])

ground <- flux %>%
  bind_rows(phenocam) %>%
  bind_rows(pep) %>%
  bind_rows(transect) %>%
  bind_rows(randomn) %>% 
  mutate(type = factor(type,
                       levels = c("Random points", "PEP725", "FLUX", "PhenoCam", "Transect"),
                       labels = c("Random points (10 000 pixels)",
                                  "PEP725 (1321 sites)",
                                  "Eddy covariance (49 sites)",
                                  "PhenoCam (35 sites)",
                                  "Transects (4 500 pixels)"))) %>%
  st_as_sf(coords = c("lon", "lat"))
st_crs(ground) <- 4326

ggplot(data = ground) +
  geom_sf(data = europe,
          fill = "grey85",
          col = "grey85",
          size = 0.3) +
  geom_sf(aes(color = type, size = type, shape = type))+
  scale_shape_manual(values = c(20, 20, 20, 3, 20))+
  scale_color_manual(values = c("grey70", "grey20", "red", "blue", "orange")) +
  scale_size_manual(values = c(0.01, 0.5, 2, 3, 0.1)) +
  theme_bw() +
  coord_sf(crs = euEqualArea) +
  theme(panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.201, 0.87),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, 'cm'))

ggsave("figures/figure_study_area_and_data.pdf",
       width = 20, height = 18, units = "cm")


ggsave("figures/report_figure_study_area_and_data.pdf",
       width = 20, height = 18, units = "cm")