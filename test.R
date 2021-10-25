#### using data from promact ####

NumIter = 10000
luxdata = seq(1, 28, 1)
time = seq(0,801,1)
PT_time_points = promact$time


#### light output ####

PT_heights = promact$X100mM
lux100 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X40mM
lux40 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X20mM
lux20 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X10mM
lux10 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X5mM
lux5 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X2mM
lux2 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X1mM
lux1 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X0.5mM
lux0.5 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X0.1mM
lux0.1 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X0.05mM
lux0.05 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)
PT_heights = promact$X0.01mM
lux0.01 = lux_extended_sub_final(time,PT_time_points,PT_heights,luxparams,starts)

# data inaccuracy error FIXED?
# add on empty 0s
# lux0.5fix = c(lux0.5$Lux, numeric(133))

luxoutput = data.frame(t(rbind(lux100$Lux, lux40$Lux, lux20$Lux, lux10$Lux, lux5$Lux, lux2$Lux, lux1$Lux, lux0.5$Lux, lux0.1$Lux, lux0.05$Lux, lux0.01$Lux)))
luxoutput = luxoutput/max(luxoutput)

ggplot(luxoutput, aes(x = time)) +
  geom_line(aes(y = X1, color = "100mM")) +
  geom_line(aes(y = X2, color = "40mM")) +
  geom_line(aes(y = X3, color = "20mM")) +
  geom_line(aes(y = X4, color = "10mM")) +
  geom_line(aes(y = X5, color = "5mM")) +
  geom_line(aes(y = X6, color = "2mM")) +
  geom_line(aes(y = X7, color = "1mM")) +
  geom_line(aes(y = X8, color = "0.5mM")) +
  geom_line(aes(y = X9, color = "0.1mM")) +
  geom_line(aes(y = X10, color = "0.05mM")) +
  geom_line(aes(y = X11, color = "0.01mM")) +
  ggtitle("Light Output from Promoter Activity") +
  ylab("Light Output (normalized)") +
  xlab("Time")


#### MCMC ####

# temp = MCMC_Prom_Activity(PT_time_points,PT_heights,NumIter=NumIter,frac=0.35,NORMFAC = 10000,Beta=1000,Sigma=0.01,Const_h=1)
# 
# luxdata = rbind(luxdata, temp[[2]][NumIter,])
# 
# luxdata = data.frame(t(luxdata))
# 
# #PLOT 5B
# 
# ggplot(luxdata, aes(x = luxdata)) +
#   geom_line(aes(y = X, color = "100mM")) +
#   geom_line(aes(y = X.1, color = "40mM")) +
#   geom_line(aes(y = X.2, color = "20mM")) +
#   geom_line(aes(y = X.3, color = "10mM")) +
#   geom_line(aes(y = X.4, color = "5mM")) +
#   geom_line(aes(y = X.5, color = "2mM")) +
#   geom_line(aes(y = X.6, color = "1mM")) +
#   geom_line(aes(y = X.7, color = "0.5mM")) +
#   geom_line(aes(y = X.8, color = "0.1mM")) +
#   geom_line(aes(y = X.9, color = "0.05mM")) +
#   geom_line(aes(y = X.10, color = "0.01mM")) +
#   ggtitle("Luciferase Activty") +
#   ylab("Inferred Promoter Activity") +
#   xlab("Time")







