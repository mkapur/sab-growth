require(dplyr)
require(ggplot2)
require(gridExtra)
scenarios <- read.csv('inputs/scenarios.csv',na.strings = 'NA') ## manual file

## average ## of age X fish per dataset

a6df <- list.files("genData", full.names = T) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  filter(Age == 6)
dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds


## panel plot of one gendata for each scenario
scens <-  unique(scenarios$DESC)
scens.title <- scens
levels(scens.title) <-  c( "(b) Break at 25 deg.", "(c) Break at 49 deg.",
                     "(d) Low Contrast at 25 deg.", "(e) Overlap 20-25 deg.","(a) No Breaks")
plist <- list()
for(q in 1:length(scens)){
  # num <- floor(runif(1,1,100))
  scen <- scens[q]
  tempdf <- read.csv(paste0("genData/",scen,"_",num,".csv")) %>% filter(Age == 6); cat(q,nrow(tempdf),"\n")
  col2 <- ifelse(unique(tempdf$REG)[2] == 'R3', 'seagreen','slategrey') ## UNIQUE COLOR DEPENDING ON REGION
  
  levels(tempdf$REG) <-  c('Regime 1', ifelse(levels(tempdf$REG)[2] == 'R2', 'Regime 2','Regime 3'))  
  plist[[q]] <- ggplot(tempdf, aes(x = Latitude_dd, y = Length_cm, color = REG)) +
    theme_classic() +
    theme(legend.position = 'none') +
    scale_y_continuous(limits = c(100,400)) +
    scale_color_manual(values = c('goldenrod',col2), breaks = c('Regime 1','Regime 2','Regime 3')) +
    labs(x = 'Latitude',y = 'Length of Age-6 Fish (cm)', title = scens.title[q]) +
    annotate("text",x = 10, y = 390, label = paste0("n = ",nrow(tempdf))) +
    geom_point()
  
}
lay <- rbind(c(1,1,1,1),
             c(2,2,3,3),
             c(4,4,5,5))
grid.arrange(grobs = plist, layout_matrix = lay) %>% ggsave(plot = .,  file = paste0(getwd(),"/plots/scenarios.png"), width = 6, height = 8, units = 'in', dpi = 480)

