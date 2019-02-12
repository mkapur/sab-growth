# gs <- lapply(1:6, function(ii) 
#   grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)),
#            textGrob(ii)))
# grid.arrange(grobs=gs, ncol=4, 
#              top="top label", bottom="bottom\nlabel", 
#              left="left label", right="right label")
# grid.rect(gp=gpar(fill=NA))
lay <- rbind(c(1,2),
             c(1,3),
             c(1,4),
             c(1,5),
             c(1,6))
grid.arrange(grobs = gs, layout_matrix = lay)

plot_grid(breakhist,plist[[3]],plist[[3]], ncol = 2, align = 'h', axis = 'l')
