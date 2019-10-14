library(DescTools)
## SAB STUFF ## ----
## average # of age X SAB per sex in fulldat
load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
# full_data %>% filter(Age %in% c(4,6,10,30)) %>%
#   group_by(Age, Sex) %>% summarise(n=n()) %>% write.csv(.,paste0("./output_data/n_sab_sex_age.csv"),row.names=F)
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

## Figure 8 SAB map for a6 and a10 breaks, likely redone by Elliot ----
# usa <- map_data("world")
#   dat <- all_data %>% filter(Age == 30 & Sex == 'M')
#     ggplot() + 
#     geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
#     coord_quickmap() +
#     scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
#     scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
#     theme_minimal() +
#     theme(panel.grid.major = element_blank(),
#           axis.title =element_blank(), 
#           legend.position = c(0.1,0.2),
#           legend.key.size = unit(0.75, "cm"),
#           legend.text = element_text(size = 12)) +
#     geom_hline(yintercept = c(36,50),lwd = 1.1, linetype = 'dashed', col = 'red') +
#     geom_vline(xintercept =  -130, lwd = 1.1, linetype = 'dashed', col = 'red') +
#       ## ecosystem break
#     geom_vline(xintercept =  -145, lwd = 1.1, linetype = 'dashed', col = 'blue') +
#     geom_rect(fill = 'white', aes(xmin = -180, xmax = -130.1, ymin = 30, ymax = 49.9)) +
#     geom_point(data = dat, aes(x = Longitude_dd, y = Latitude_dd, size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
#     scale_fill_viridis_c(guide = "legend") +
#     labs(fill = paste0("Length (cm)"),
#          size = paste0("Length (cm)"),
#          title = "GAM- and Ecosystem-based Regions with Age 30 Fish") +
#       geom_label(aes(x = c(-155,-138,-120,-120,-120), 
#                     y = c(65,65,65,40,32),
#                     label = c(paste('Region ',5:1)) ),size = 6, col = 'black',fill = 'white',show.legend   =FALSE)
#     ggsave(plot = last_plot(),  file = paste0("./figures/sab_zones.png"), width = 8, height = 6, units = 'in', dpi = 480)
#     
# df1<-all_data %>% filter(Age %in% c(4,6,30)) %>% 
#   group_by(REG, Age) %>% summarise(n = n()) 
# 
# df1 %>%
#   group_by(Age) %>% summarise(mn = mean(n))

## predicts and parest for sab ----
## parest: see if actually different. PHASE 1 = pre-merge, PHASE 2- temporal merges
usedate <- "2019-10-04" #Sys.Date() ## USE TODAY

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

for(phase in c("phase1","phase2")){
  parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_",usedate,"_",phase,'.csv')) %>%
    filter(variable ==  "Linf") %>%
    mutate(source = 'Estimated') %>%
    mutate(REG2 = gsub("_.*", "\\1", REG),
           REG3 = as.numeric(substr(REG,2,2)),
           REG4 = sub('_([^_]*)$', '',REG),
           Sex = gsub(".*_", "\\1", REG),
           lwr = value - 1.96 * sd,
           upr = value + 1.96 * sd,
           matchcol = 'black')
  # parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
  parest <- parest[order(parest$REG2,parest$Sex),]

  # parest <- rbind(parest, read.csv("./input_data/true_sab_vals.csv"))
  # levels(parest$REG) <- c("ALL","R1","AK","R2","BC","R1","WC")
  # parest$REG <- factor(parest$REG ,levels=c("ALL","R1","WC","R2","BC","R3","AK",'R4'))

  ## check if CIs overlap
  if(phase == 'phase1'){
    for(i in seq(1,nrow(parest),2)) { ## will go by sex
      parest$match[i:(i+1)] <- ifelse(c(parest$lwr[i], parest$upr[i]) %overlaps%  c(parest$lwr[i+1], parest$upr[i+1]),
                                      paste0('OVERLAP ', parest$REG2[i+1]), 'NO OVERLAP')
    }

  } else if(phase == 'phase2') {
    for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
      if(parest$REG3 < 5){
        tmp <- subset(parest, Sex == parest$Sex[i] & REG3 == parest$REG3[i]+1)
        parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
      } ## end R < 5
    } ## end parest rows
  } ## end else phase2

  parest$matchcol <- parest$match != 'NO OVERLAP'
  parest$Sex <- factor(parest$Sex)
  levels(parest$Sex) <- c('Females','Males')
  ggplot(parest, aes(x = REG4, y = value, col = matchcol))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'right',
          legend.background = element_blank(),
          axis.text = element_text(size = 10,angle = 45),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,100)) +
    scale_color_manual(values = c('black','red')) +
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatiotemporal x Sex Stratum', y = "",
         col = ifelse(phase == 'phase1', 'CI overlap within Region + Sex','CI Overlap Adjacent Region x Sex'),
         title = paste0(phase," Linf Estimates")) +
    facet_wrap(~Sex )

ggsave(plot = last_plot(),  file = paste0("./figures/sab_parest_",Sys.Date(),"_",phase,".png"), width = 10, height = 8, units = 'in', dpi = 480)
write.csv(parest, file = paste0("./GAM_output/overlap_",Sys.Date(),"_",phase,".csv"),row.names=F)


## Figure 9 SAB Fits ----
ypreds <- read.csv(paste0("./GAM_output/SAB_predicts_",usedate,"_",phase,".csv"))
ypreds$gamREG <- paste0('GAM-Detected ',ypreds$gamREG)
levels(ypreds$REG) <- c('Alaska','British Columbia','US West Coast')
levels(ypreds$Sex) <- c('Females','Males')
levels(ypreds$Period) <- c('pre-2010','2010-2018','All Years')
for(i in 1:nrow(ypreds)){
  ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1,
                              'All Years', paste(ypreds$Period[i]))
}

ggplot(ypreds, aes(x = Age, y = Predicted, col = REG, linetype = Period )) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size=10))+
  
  scale_alpha(guide = 'none') +
  guides(shape = FALSE) +
  
  scale_y_continuous(limits = c(0,110)) +
  scale_x_continuous(limits = c(0,65)) +
  scale_colour_manual(values = c('black','grey22','grey44','grey77'),
                      guide = guide_legend(override.aes = list(
                        shape = c(15:17),
                        linetype = c("solid", "dashed", "dotted")))) +
  
  geom_point(alpha = 0.2, aes(y = Length_cm, shape = REG, col = REG)) +
  geom_line(lwd = 1.1, col = 'black')+
 
  
  labs(y = 'Length (cm)', x= 'Age (years)', col = "Actual Data Source") +
  # scale_linetype_manual(values=c("solid", "dashed", "dotted"))+
  facet_wrap(~gamREG + Sex , ncol = 4)

ggsave(plot = last_plot(),
       file = paste0("./figures/sab_fits_",Sys.Date(),"_",phase,".png"),
       width = 10, height = 12, units = 'in', dpi = 520)


fd_summary_gamREG <- ypreds %>%
  # filter(Age < 75) %>%
  group_by(Age, Sex, gamREG,Period) %>%
  dplyr::summarise(meanL = mean(Length_cm), sdmeanL = sd(Length_cm), meanPred = mean(Predicted))

ggplot(fd_summary_gamREG, aes(x = Age, col = gamREG, group = gamREG)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size=10))+
  scale_alpha(guide = 'none') +
  scale_y_continuous(limits = c(0,110)) +
  scale_x_continuous(limits = c(0,65)) +

  scale_colour_manual(values = c('grey66','grey44',cbbPalette[1:5],'grey22'),
                      guide = guide_legend(override.aes = list(
                        linetype = c(NA,NA,rep(1,5),NA),
                        shape = c(15,16,rep(NA,5),17)))) +
  geom_point(data = ypreds, alpha = 0.5, 
             aes(x = Age, y = Length_cm, pch = REG, col = REG))+
  guides(pch = FALSE) +

  
  geom_line(aes(y = meanPred, col = gamREG), lwd = 1.1)+
  labs(y = 'Length (cm)', x= 'Age (years)', col = "") +
  facet_wrap(~Sex +Period, ncol = 4)

ggsave(plot = last_plot(),
       file = paste0("./figures/STARstyle_sab_fits_",Sys.Date(),"_",phase,".png"),
       width = 10, height = 7, units = 'in', dpi = 520)

cat(phase," done \n")
} ## end phase



# all_data %>% filter(Age %in% c(4,6,30)) %>% group_by(Age,Sex,REG) %>% summarise(n = n())


## Make T3 to avoid copypaste problems
## Right now the L1/L2 is reported using age 4, but consider re-doing for 0.5

LR <- function(Linf, k, t0, AREF){
  lreg <- Linf*(1-exp(-k*(AREF - t0)))
  return(lreg)
}
## Get numbers by category for final phase.
ypreds <- read.csv(paste0("./GAM_output/SAB_predicts_",usedate,"_phase2.csv"))
samplesize <- ypreds %>% group_by(cREG) %>% summarise(N = n()) %>% select(N) %>% as.vector()


sabpe <- read.csv("./GAM_output/SAB_parEst_gam_2019-10-04_phase2.csv")
sabpe <- sabpe %>%
  mutate(REG = as.character(REG)) %>%
  select(-sd) %>%
  filter(variable != 'Sigma') %>%
  mutate(
         Region = substr(REG,1,2),
         Period =  sapply(strsplit(REG, "_"), function(x) x[2]),
         Sex = sub('(^[^_]+_[^_]+)_(.*
                   )$', '\\2', REG) ) %>%
  tidyr::spread(key = variable, value = value) %>%
  mutate(   L1_0.5 = round(LR(Linf, k, t0, AREF = 0.5),2)) %>%
  mutate(L1 = ifelse(L1_0.5 <0, L1, L1_0.5)) %>%
  select(Region, Sex, Period, Linf, k, t0, L1, L2) 

cbind(sabpe, ypreds %>% group_by(cREG) %>% summarise(N = n()) %>% select(N)) %>% 
  select(Region, Sex, Period, N, Linf, k, t0, L1, L2) %>%

write.csv(., file = paste0("./figures/table3_",Sys.Date(),".csv"),row.names = F)

