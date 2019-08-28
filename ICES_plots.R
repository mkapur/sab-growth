## rebooting some manuscript figures in theme_black
require(ggplot2)
require(kaputils)
require(dplyr)
source("C:/Users/mkapur/Dropbox/kaputils/R/theme_black.R")

## blank black map -- three regions
usa <- map_data("world") 

xlims <- list(c(-180,-140),c(-133,-125),c(-140,-115))
ylims <- list(c(50,70),c(45,55),c(30,50))
plist <- list();idx <- 1
for(i in 1:3){
# plist[[idx]] <-
  ggplot() +
  geom_polygon(data = usa, 
               aes(x = long, y = lat, group = group),
               fill = 'grey66') +
  coord_quickmap() +
  scale_x_continuous(expand = c(0,0), limits = xlims[[i]], 
                     breaks = seq(xlims[[i]][1],xlims[[i]][2],10),
                     labels = paste(seq(xlims[[i]][1],xlims[[i]][2],10), "°W")) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(ylims[[i]][1],ylims[[i]][2]), 
                     breaks =  seq(ylims[[i]][1],ylims[[i]][2],10), 
                     labels =  paste(seq(ylims[[i]][1],ylims[[i]][2],10), "°N"))  +
  theme_black(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        axis.title =element_blank(),
        legend.position = c(0.1,0.2),
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size = 12))
# idx <- idx+1
ggsave(plot = last_plot(),
       file = paste0("./figures/BLACK_zone",i,".png"),
       width = 10, height = 10, units = 'in', dpi = 480)

}
lay <- cbind(1,2,3)
grid.arrange(grobs = plist, layout_matrix = lay)  %>%


## black map w breaks (not currents)
usa <- map_data("world") 
  dat <- all_data %>% filter(Age == 30 & Sex == 'M')
    ggplot() +
    geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = 'grey66') +
    coord_quickmap() +
    scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
    scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
    theme_black(base_size = 20) +
    theme(panel.grid.major = element_blank(),
          axis.title =element_blank(),
          legend.position = c(0.1,0.2),
          legend.key.size = unit(0.75, "cm"),
          legend.text = element_text(size = 12)) +
    geom_hline(yintercept = c(36,50),lwd = 1.1, linetype = 'dashed', col = 'red') +
    geom_vline(xintercept =  -130, lwd = 1.1, linetype = 'dashed', col = 'red') +
      ## ecosystem break
    geom_vline(xintercept =  -145, lwd = 1.1, linetype = 'dashed', col = 'blue') +
    geom_rect(fill = 'black', aes(xmin = -180, xmax = -130.1, ymin = 30, ymax = 49.9)) +
    # geom_point(data = dat, aes(x = Longitude_dd, y = Latitude_dd, size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
    # scale_fill_viridis_c(guide = "legend") +
    labs(fill = paste0("Length (cm)"),
         size = paste0("Length (cm)"),
         title = "GAM- and Ecosystem-based Regions") +
      geom_label(aes(x = c(-155,-138,-120,-120,-120),
                    y = c(65,65,65,40,32),
                    label = c(paste('Region ',5:1)) ),size = 8, col = 'black',fill = 'white',show.legend   =FALSE)
    ggsave(plot = last_plot(),  file = paste0("./figures/BLACK_sab_zones.png"), width = 18, height = 12, units = 'in', dpi = 480)


## black neatplot from sab with faceted curves
ypreds <- read.csv(paste0("C:/Users/mkapur/Dropbox/UW/sab-growth/GAM_output/SAB_predicts_2019-04-15_phase2.csv"))

ypreds$gamREG <- paste0('Region ',ypreds$gamREG)
levels(ypreds$Sex) <- c('Females','Males')
levels(ypreds$Period) <- c('pre-2010','2010-Present','All Years')
for(i in 1:nrow(ypreds)){
  ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1,
                              'All Years', paste(ypreds$Period[i]))
}

fd_summary_gamREG <- ypreds %>%
  filter(Age < 31) %>%
  group_by(Age, Sex, gamREG,Period) %>%
  dplyr::summarise(meanL = mean(Length_cm), sdmeanL = sd(Length_cm), meanPred = mean(Predicted))

ggplot(fd_summary_gamREG, aes(x = Age, col = gamREG, group = gamREG)) +
  theme_black(base_size = 20) +
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  scale_alpha(guide = 'none') +
  scale_y_continuous(limits = c(0,110)) +
  scale_color_brewer(palette =  'Accent')+
  geom_point( aes(y = meanL)) +
  geom_line(aes(y = meanPred, col = gamREG), lwd = 1.1)+
  labs(y = 'Length (cm)', x= 'Age (years)', col = "GAM-Detected Region") +
  # scale_linetype_manual(values=c("solid", "dashed", "dotted"))+
  facet_wrap(~Sex +Period, ncol = 4)
# 
ggsave(plot  = last_plot(),
       file = "C:/Users/mkapur/Dropbox/UW/sab-growth/figures/VBGF_meanL_periodXregion_black.png",
       height = 8, width = 16, unit = 'in', dpi = 520)


phase <- "phase2"
parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_2019-04-15_",phase,'.csv')) %>%
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
  theme_black(base_size = 20) +
  # theme(panel.grid = element_blank(),
  #       legend.position = 'right',
  #       legend.background = element_blank(),
  #       axis.text = element_text(size = 10,angle = 45),
  #       axis.title = element_text(size = 10),
  #       legend.text = element_text(size = 10),
  #       strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c('white','red')) +
  geom_point() +
  geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
  labs(x = 'Spatiotemporal x Sex Stratum', y = "",
       col = ifelse(phase == 'phase1', 'CI overlap within Region + Sex','CI Overlap Adjacent Region x Sex'),
       title = paste0(phase," Linf Estimates")) +
  facet_wrap(~Sex )

ggsave(plot = last_plot(),  file = paste0("./presentation/sab_parest_",Sys.Date(),".png"), width = 10, height = 8, units = 'in', dpi = 480)
# write.csv(parest, file = paste0("./GAM_output/overlap_",Sys.Date()-2,"_",phase,".csv"),row.names=F)

## smoother plots in black ----
## Fit GAMS to ID breakpoints ----
load(paste0("./input_data/gam_data_sab_0315.rda")) ## all_data -- made using gam_dataprep and 15k subsample
source(paste0(getwd(),"/functions/Deriv.R"))
require(mgcv)
for(AA in c(4,6,30)){
  for(SS in c("F","M")){
    dat <- all_data %>% filter( Age == AA  & Sex == SS)
    # dat <- all_data %>% filter( Age == 6  | Age == 10  & Sex == SS)
    
    mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd),data = dat)

    
    ## get & eval derivatives ----
    pdat <- dat
    pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
    p2 <- predict(mod, newdata = pdat) ## raw predicts
    pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2],se2_lon = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,1])
    
    df.res <- df.residual(mod)
    crit.t <- qt(0.025, df.res, lower.tail = FALSE)
    ## variances are additive just FYI
    pdat <- transform(pdat,
                      upper = predLen + (crit.t * (se2_lat+se2_lon+se2_yr)),
                      lower = predLen - (crit.t * (se2_lat+se2_lon+se2_yr)))
    
    breaksdf <- list()
    Terms <- c("Year","Latitude_dd","Longitude_dd")
    for(t in 1:length(Terms)){
      Term <- Terms[t]
      
      
      newD <- data.frame(
        Year = seq(1981, 2017, length = 100),
        
        Latitude_dd = seq(0, 64, length = 100),
        Longitude_dd = seq(-186, -116, length = 100)
      )## whatever isn't modeled will just get ignored
      
      xlims <- data.frame(
        Year = c(1990, 2017),
        
        Latitude_dd = c(30, 60),
        Longitude_dd = c(-150, -116      ))
      
      m2.d <- Deriv(mod, newdata = newD)
      m2.dci <- confint(m2.d, term = Term)
      
      crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
      crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails
      
      
      ## identify where mean crosses zero or falls out of bounds (NAs where it does)
      m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
                                d = m2.d[[Term]]$deriv,
                                upper = m2.dci[[Term]]$upper,
                                lower = m2.dci[[Term]]$lower,
                                crit.eval = crit.eval,
                                eval = 0)
      
      pix <- !is.na(m2.dsig.zeros)
      
      vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
      breaksdf[[t]] <- sort(c(unique(round(vals)))) ## get rounded unique
      ## fill NAs in bdf for binding
      for(i in 1:length(breaksdf)){
        if (length(breaksdf[[i]]) == 0){ ## fill NA for empty
          breaksdf[[i]] <- NA
        }
      }## end breaksdf
      
      
    } ## end terms
    
    plist <- list()
    idx <- 1
    ## extract plotting stuff
    pd <- plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c(-10,ifelse(t != 3,10,500)))
    for(t in 1:length(Terms)){
      
      ## plot smooth
      temp0 <- pd[t] ## get this smoother
      temp <- data.frame(cbind(temp0[[1]]$x,temp0[[1]]$fit, temp0[[1]]$se)); names(temp) = c('x','fit','se')
      plist[[idx]] <-  ggplot(temp, aes(x = x, y = fit)) +
        theme_black(base_size = 18) +
        theme(panel.grid = element_blank())+
        geom_vline(xintercept = ifelse(is.na(breaksdf[[t]]),-999,breaksdf[[t]]), lwd = 1.1, linetype = 'dashed', col = 'red') +
        # scale_x_continuous(expand = c(0,0), limits = c(ifelse(t == 2, 30, min(temp$x)),max(temp$x))) +
        scale_x_continuous(expand = c(0,0), limits = c(xlims[1,t],xlims[2,t])) +
        scale_y_continuous(limits = c(ifelse(t == 3, -15, min(temp$fit)*3), ifelse(t == 3, 10, max(temp$fit)*3)))  +
        
        geom_line(lwd = 1.1, col = 'white') +
        geom_line(aes(y= fit-se), linetype = 'dashed', col = 'white') +
        geom_line(aes(y= fit+se), linetype = 'dashed', col = 'white') +
        geom_rug(sides = 'b') +
        labs(x = Terms[t], y = "smoother", title = paste0(letters[idx],") ",'Smoother for ', temp0[[1]]$xlab))
      
      ## plot deriv
      idx <- idx+1
      CI <- confint(m2.d, term = Terms[t])
      m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv, CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
      plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
        theme_black(base_size = 18) +
        theme(panel.grid = element_blank())+
        geom_line(lwd = 1.1, col = 'white') +
        geom_hline(yintercept = 0, col = 'grey22') +
        geom_line(aes(y= upper), linetype = 'dashed', col = 'white') +
        geom_line(aes(y= lower), linetype = 'dashed', col = 'white') +
        geom_vline(xintercept = ifelse(is.na(breaksdf[[t]]),-999,breaksdf[[t]]), lwd = 1.1, linetype = 'dashed', col = 'red') +
        # scale_x_continuous(expand = c(0,0), limits = c(ifelse(t == 2, 30, min(m2.dtemp$x)),max(m2.dtemp$x))) +
        scale_x_continuous(expand = c(0,0), limits = c(xlims[1,t],xlims[2,t])) +
        
        labs(x = Terms[t], y = "f'(x)", title =  paste0(letters[idx],") ",'First Derivative for ', temp0[[1]]$xlab))
      idx <- idx+1
    } ## end terms plotting loop
  
    
    
    lay <- rbind(c(1,2),
                 c(3,4),
                 c(5,6))
    grid.arrange(grobs = plist, layout_matrix = lay)  %>%
      ggsave(plot = .,  file = paste0("./figures/BLACK_sab_smooth_map_",AA,SS,".png"), width = 12, height = 15, units = 'in', dpi = 480)
  } ## end SS
} ## end AA

## fake age/length/weight plot for illustration
Linf <- 64
t0 <- -2
k <- 0.4
len <- len2  <- wt <- wt2<- NULL
wa <-1.3528e-2
wb <- 3.42971
for(a in 1:30){
  len[a] <- Linf*(1-exp(-k*(a-t0)))
  len2[a] <- Linf*0.97*(1-exp(-k*0.75*(a-t0*1.2)))
  
  wt[a] <- wa*len[a]*wb
  wt2[a] <- wa*len2[a]*wb
  
}

df<- data.frame(cbind(len,len2,wt,wt2, 'age' = 1:30)) %>%
  melt(id = c('age'))
df2 <- data.frame(cbind(len,len2,wt,wt2, 'age' = 1:30)) %>%
  melt(id = c('len'))
plist <- list()
plist[[1]] <- ggplot(subset(df, variable == 'len'  | variable == 'len2'), 
                     aes(x = age, y = value, col = variable, group= variable)) +
  theme_black(base_size = 20) +
  geom_line(lwd = 1.1) +
  labs(y = 'length') +
  scale_color_brewer(palette =  'Accent')+
  scale_y_continuous(limits = c(50,65))+
  theme(legend.position = 'none')
plist[[2]] <- ggplot(subset(df2, variable == 'wt'| variable == 'wt2'),
                     aes(x = len, y = value*10, col = variable, group= variable)) +
  theme_black(base_size = 20) +
  geom_line(lwd = 1.1) +
  labs(y = 'weight', x = 'length') +
  scale_color_brewer(palette =  'Accent')+
  # scale_x_continuous(limits = c(10,75))+
  scale_y_continuous(limits = c(15,30))+
  
  theme(legend.position = 'none')

lay <- rbind(c(1,2))
grid.arrange(grobs = plist, layout_matrix = lay)  %>%
  ggsave(plot = .,  file = paste0("./figures/BLACK_lenweightage.png"),
  width = 10, height = 7, units = 'in', dpi = 480)
