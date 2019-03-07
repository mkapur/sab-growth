makeMod <- function(scenario,dat){
  mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd), data = dat) ## WTF is this

  ## plotting model
  # png(paste0(outdir,"/boot_",b,"_mod.png"))
  # layout(matrix(1:3, ncol = 1))
  # plot(mod, scale = 0)
  # pacf(resid(mod), lag.max = 36, main = "pACF")
  # graphics.off()
  
  ## now try with some AR structures and check AIC
  # mod1 <- gam(Length_cm ~   s(Latitude_dd),
  #             correlation = corAR1(form = ~ 1|Year, p = 1),
  #             data = dat)
  # mod2 <- gam(Length_cm ~ s(Year, bs = "cc") +
  #               s(Latitude_dd),
  #             correlation = corAR1(form = ~ 1|Year, p = 2),
  #             data = dat)
  # mod3 <- gam(Length_cm ~  s(Year, bs = "cc") +
  #               s(Latitude_dd),
  #             correlation = corAR1(form = ~ 1|Year, p = 3),
  #            data = dat)
  
  cat("built model \n")
  
  # ms <- MuMIn::model.sel(mod,mod1,mod2,mod3) ## no improvement
  # if(ms[1]$correlation != ''){ ## if model without correlation picked
  #   stop("model without correlation picked, hand inspect \n")
# }
  png(
    file = paste0(outdir,"/boot_",b,  "_gamcheck.png"
    ),
    height = 6,
    width = 8,
    units = 'in',
    res = 500
  )
  layout(matrix(1:4, ncol = 2))
  gam.check(mod)
  graphics.off()
  cat("plotted gam check \n")
  
  return(mod)
}
