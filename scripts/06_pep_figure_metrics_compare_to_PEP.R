
library(tidyverse)

load("data/RData/pepS2.RData")

brdf = "M"
ind = "diff"

pepcmb <- function(brdf, ind) {
  threscomall <- pepS2 %>% 
    filter(BRDF == brdf) %>% 
    select(-c(BRDF, year)) %>%
    group_by(VIname, phasePair, thres, type) %>% 
    summarise(cor = cor(PEPdoy, S2doy, method = "pearson",
                        use = "pairwise.complete.obs")^2,
              diff = mean(PEPdoy - S2doy, na.rm = T),
              std = sd(PEPdoy - S2doy, na.rm = T)) %>% 
    ungroup() %>% 
    pivot_longer(cols = c(cor, diff, std),
                 names_to = "indicator",
                 values_to = "value") %>% 
    mutate(thres = as.numeric(str_extract(thres, "\\d{2}"))) 
  
  if (ind == "cor") {
    threscomall %>% 
      filter(indicator == ind) %>% 
      select(-indicator) %>% 
      
      ggplot(aes(x = thres, y = value, color = VIname)) +
      geom_line() +
      # geom_vline(xintercept = 25, linetype = 2, color = "grey60")+
      # geom_hline(yintercept = c(0, 0.2), linetype = 2, color = "grey60")+
      geom_point(size = 2) +
      scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
      facet_grid(type ~ phasePair) +
      # facet_grid(phasePair ~ type) +
      theme_bw(base_size = 17)+
      labs(x = "Amplitude thresholds (%)", 
           y = expression(paste("Spatial ",italic("R"^"2"), sep = "")))+
      theme(panel.grid = element_blank(),
            # panel.grid.major.x = element_blank(),
            legend.title = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 17)
            # legend.direction = "horizontal"
      )
  } else {
    threscomall %>% filter(indicator != "cor", type != "Crop") %>% 
      pivot_wider(names_from = indicator, values_from = value) %>% 
      ggplot(aes(x = thres, y = diff, color = VIname)) +
      # geom_ribbon(aes(ymin = diff - std, ymax = diff + std, 
      #                 fill = VIname), alpha = 0.1, outline.type = "none") +
      # geom_line() +
      # geom_point(size = 2) +
      geom_vline(xintercept = 25, linetype = 2, color = "grey60")+
      geom_hline(yintercept = 0, linetype = 2, color = "grey60")+
      geom_pointrange(aes(ymin = diff - std, ymax = diff + std),
                      fatten = 3,
                      position=position_dodge(width = 3))+
      scale_color_manual(values = c("tan3", "blue", "firebrick3")) +
      # scale_fill_manual(values = c("tan3", "blue", "firebrick3")) +
      facet_grid(type ~ phasePair, scales = "free") +
      theme_bw(base_size = 20)+
      labs(x = "Amplitude thresholds (%)", 
           y = "Time difference (day)") +
      theme(panel.grid = element_blank(),
            # panel.grid.major.x = element_blank(),
            legend.title = element_blank(),
            legend.position = "top",
            legend.text = element_text(size = 20)
            # legend.direction = "horizontal"
      )
  }
}


pepcmb("M", "cor")
ggsave("figures/figure_metrics_compare_to_PEP_correlation_NBAR.pdf",
       width = 45, height = 30, units = "cm")

pepcmb("M", "diff")
ggsave("figures/figure_metrics_compare_to_PEP_NBAR.pdf",
       width = 40, height = 25, units = "cm")


pepcmb("raw", "cor")
ggsave("figures/figure_metrics_compare_to_PEP_correlation_raw.pdf",
       width = 45, height = 30, units = "cm")

pepcmb("raw", "diff")
ggsave("figures/figure_metrics_compare_to_PEP_raw.pdf",
       width = 40, height = 25, units = "cm")
