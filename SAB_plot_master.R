
## SAB STUFF ## ----

## average # of age 6 fish per dataset----
# a6df <- list.files("./IBM_output/datasets", full.names = T) %>%
#   lapply(read.csv) %>%
#   bind_rows() %>%
#   filter(Age == 6)
# dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds

## average # of age X SAB per sex in fulldat
# load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
# full_data %>% filter(Age %in% c(4,6,10,30)) %>% 
#   group_by(Age, Sex) %>% summarise(n=n()) %>% write.csv(.,paste0("./output_data/n_sab_sex_age.csv"),row.names=F)


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
# for(phase in c("phase1","phase2")){
#   parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_2019-04-15_",phase,'.csv')) %>% 
#     filter(variable ==  "Linf") %>% 
#     mutate(source = 'Estimated') %>%
#     mutate(REG2 = gsub("_.*", "\\1", REG),
#            REG3 = as.numeric(substr(REG,2,2)),
#            REG4 = sub('_([^_]*)$', '',REG),
#            Sex = gsub(".*_", "\\1", REG),
#            lwr = value - 1.96 * sd,
#            upr = value + 1.96 * sd,
#            matchcol = 'black')
#   # parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
#   parest <- parest[order(parest$REG2,parest$Sex),]
#   
#   # parest <- rbind(parest, read.csv("./input_data/true_sab_vals.csv")) 
#   # levels(parest$REG) <- c("ALL","R1","AK","R2","BC","R1","WC")
#   # parest$REG <- factor(parest$REG ,levels=c("ALL","R1","WC","R2","BC","R3","AK",'R4'))
#   
#   ## check if CIs overlap
#   if(phase == 'phase1'){
#     for(i in seq(1,nrow(parest),2)) { ## will go by sex
#       parest$match[i:(i+1)] <- ifelse(c(parest$lwr[i], parest$upr[i]) %overlaps%  c(parest$lwr[i+1], parest$upr[i+1]),
#                                       paste0('OVERLAP ', parest$REG2[i+1]), 'NO OVERLAP')
#     }
#     
#   } else if(phase == 'phase2') {
#     for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
#       if(parest$REG3 < 5){
#         tmp <- subset(parest, Sex == parest$Sex[i] & REG3 == parest$REG3[i]+1)
#         parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
#       } ## end R < 5
#     } ## end parest rows
#   } ## end else phase2
#   
#   parest$matchcol <- parest$match != 'NO OVERLAP'
#   parest$Sex <- factor(parest$Sex)
#   levels(parest$Sex) <- c('Females','Males')
#   ggplot(parest, aes(x = REG4, y = value, col = matchcol))+
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           legend.position = 'right',
#           legend.background = element_blank(),
#           axis.text = element_text(size = 10,angle = 45),
#           axis.title = element_text(size = 10),
#           legend.text = element_text(size = 10),
#           strip.text = element_text(size=14))+
#     scale_y_continuous(limits = c(0,100)) +
#     scale_color_manual(values = c('black','red')) +
#     geom_point() +
#     geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
#     labs(x = 'Spatiotemporal x Sex Stratum', y = "", 
#          col = ifelse(phase == 'phase1', 'CI overlap within Region + Sex','CI Overlap Adjacent Region x Sex'),
#          title = paste0(phase," Linf Estimates")) +
#     facet_wrap(~Sex )
# 
# ggsave(plot = last_plot(),  file = paste0("./figures/sab_parest_",Sys.Date()-2,"_",phase,".png"), width = 10, height = 8, units = 'in', dpi = 480)
# write.csv(parest, file = paste0("./GAM_output/overlap_",Sys.Date()-2,"_",phase,".csv"),row.names=F)


## Figure 9 SAB Fits ----
# ypreds <- read.csv(paste0("./GAM_output/SAB_predicts_2019-04-15_",phase,".csv"))
# 
# ypreds$gamREG <- paste0('Region ',ypreds$gamREG)
# levels(ypreds$Sex) <- c('Females','Males')
# levels(ypreds$Period) <- c('pre-2010','2010-Present','All Years')
# for(i in 1:nrow(ypreds)){
#   ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1, 
#                               'All Years', paste(ypreds$Period[i]))
# }
# 
# ggplot(ypreds, aes(x = Age, y = Predicted, col = REG, linetype = Period )) +
#   theme_classic() +
#   theme(panel.grid = element_blank(),
#         legend.position = 'right',
#         legend.background = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         strip.text = element_text(size=10))+
#   scale_alpha(guide = 'none') +
#   scale_y_continuous(limits = c(0,110)) +
#   scale_color_brewer(palette =  'Accent')+
#   geom_point(alpha = 0.5, aes(y = Length_cm)) +
#   geom_line(lwd = 1.1, col = 'black')+
#   labs(y = 'Length (cm)', x= 'Age (years)', col = "Actual Data Source") +
#   scale_linetype_manual(values=c("solid", "dashed", "dotted"))+
#   facet_wrap(~gamREG + Sex , ncol = 4)
# 
# ggsave(plot = last_plot(),  
#        file = paste0("./figures/sab_fits_",Sys.Date()-2,"_",phase,".png"), 
#        width = 10, height = 12, units = 'in', dpi = 520)
# cat(phase," done \n")
# } ## end phase



# all_data %>% filter(Age %in% c(4,6,30)) %>% group_by(Age,Sex,REG) %>% summarise(n = n())
