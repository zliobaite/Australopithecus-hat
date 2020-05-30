# 2019 11 26 I.Zliobaite

# this is to accompany my commentary to Du et al. (2020)
# the script is not optimized for efficieny and contains a lot of "hard-coded" values, sorry about that
# the random seed should allow to reproduce the runs exactly as reported in my commentary

tm = c(4.15,4.145,4.095,4.09,4,3.97,3.875,3.85,3.82,3.795,3.745,3.68,3.665,3.645,3.63,3.59,3.58,3.45,3.44,3.425,3.4,3.375,3.36,3.33,3.27,3.2685,3.22,3.16,3.05,3.035)
tm_sd <- sd(tm)
tm_dif <- diff(-tm)

t_min <- min(tm)
t_max <- max(tm)

param_ag <- 3
param_hf <- (param_ag - 1)/2


mov_agerage <- function(data,lg) {
  dd <- length(data)
  data_av <- c()
  for (sk in 1:(dd - lg + 1)){
    data_av <- c(data_av,mean(data[sk:(sk+lg-1)]))
  }
  return(data_av)
}
#test <- mov_agerage(c(1,2,3,4,5),3)

sstp <- (t_max - t_min)/29
real_quantiles <- rev(seq(t_min,t_max,sstp))
rr_dif_index <- mov_agerage(real_quantiles,2)


tm_mov <- mov_agerage(tm_dif,param_ag)
tm_dif_index <- mov_agerage(tm,2)

ind_range <- (1+param_hf):(29-param_hf)
  
pdf('fig_observed.pdf',width = 7, heigh = 4)
plot(tm_dif_index[ind_range],tm_mov,type = 'l',lwd = 3,ylim = c(0,0.09),xlim=rev(range(tm_dif_index[ind_range])),xlab = 'Age, Ma', ylab = 'Mean age gap between observations, MA')
dev.off()


set.seed(1981)

data_rn <- c()
data_rn_raw <- c()

for (sk in 1:1000){
  rn <- runif(28, min = t_min, max = t_max)
  rn <- sort(rn,decreasing = TRUE)
  rn <- c(t_max,rn,t_min)
  rn_dif <- diff(-rn)
  
  rn_mov <- mov_agerage(rn_dif,param_ag)
  rn_dif_index <- tm_dif_index
  data_rn <- rbind(data_rn,rn_mov)
  data_rn_raw <- rbind(data_rn_raw,rn)
}

rn_mn <- apply(data_rn,2,mean)
rn_sd <- apply(data_rn,2,sd)

rn_raw <- apply(data_rn_raw,2,mean)



pdf('fig_both.pdf',width = 7, heigh = 4)
plot(NA,NA,xlim=rev(range(tm_dif_index[ind_range])),ylim = c(0,0.09),xlab = 'Age, Ma', ylab = 'Mean age gap between observations, Ma',main = '(c) Expected age gaps with constant abundance')
lines(rn_dif_index[ind_range],rn_mn,type = 'l',lwd = 2,lty=5,col = 'grey50')
lines(rn_dif_index[ind_range],rn_mn+rn_sd,type = 'l',lty=5,col = 'grey50')
lines(rn_dif_index[ind_range],rn_mn-rn_sd,type = 'l',lty=5,col = 'grey50')
lines(tm_dif_index[ind_range],tm_mov,type = 'l',lwd = 4,col = 'black')
dev.off()

stp <- 1/29
uniform_quantiles <- seq(0,1,stp)

pdf('fig_uniformity_flat.pdf',width = 4.5, height = 4.5)
plot(uniform_quantiles,rn_raw,ylim = rev(range(rn_raw)),pch=16,xlab = 'Uniform quantiles',ylab = 'Age, Ma',main = '(e) Expectation with constant abundance')
lines(c(0,1),c(t_max,t_min))
dev.off()

pdf('fig_uniformity_observed.pdf',width = 4.5, height = 4.5)
plot(uniform_quantiles,tm,ylim = rev(range(tm)),pch=16,xlab = 'Uniform quantiles',ylab = '',main = '(f) Data from Du et al. (2020)')
lines(c(0,1),c(t_max,t_min))
dev.off()


data_rr <- c()
data_rr_raw <- c()
slot_length <- (t_max - t_min)/29

for (sk in 1:1000){
rr1 <- runif(4, min = t_min, max = t_min+8*slot_length)
rr2 <- runif(4, min = t_min+8*slot_length, max = t_min+12*slot_length)
rr3 <- runif(4, min = t_min+12*slot_length, max = t_min+14*slot_length)
rr4 <- runif(4, min = t_min+14*slot_length, max = t_min+15*slot_length)
rr5 <- runif(4, min = t_min+15*slot_length, max = t_min+17*slot_length)
rr6 <- runif(4, min = t_min+17*slot_length, max = t_min+21*slot_length)
rr7 <- runif(4, min = t_min+21*slot_length, max = t_min+28*slot_length)
rr <- c(rr1,rr2,rr3,rr4,rr5,rr6,rr7)

rr <- sort(rr,decreasing = TRUE)
rr <- c(t_max,rr,t_min)
rr_dif <- diff(-rr)

rr_mov <- mov_agerage(rr_dif,param_ag)
data_rr <- rbind(data_rr,rr_mov)
data_rr_raw <- rbind(data_rr_raw,rr)
}

rr_mn <- apply(data_rr,2,mean)
rr_sd <- apply(data_rr,2,sd)

rr_raw <- apply(data_rr_raw,2,mean)


pdf('fig_peak.pdf',width = 7, heigh = 4)
plot(NA,NA,ylim = c(0,0.09),xlim=rev(range(tm_dif_index[ind_range])),xlab = 'Age, Ma', ylab = "", main = '(d) Expected age gaps with waxing and waning')
lines(rr_dif_index[ind_range],rr_mn,type = 'l',lwd = 2,lty=5,col = 'grey50')
lines(rr_dif_index[ind_range],rr_mn+rr_sd,type = 'l',lty=5,col = 'grey50')
lines(rr_dif_index[ind_range],rr_mn-rr_sd,type = 'l',lty=5,col = 'grey50')
lines(tm_dif_index[ind_range],tm_mov,type = 'l',lwd = 4,col='black')
dev.off()

pdf('fig_uniformity_peak.pdf',width = 4.5, height = 4.5)
plot(uniform_quantiles,rr_raw,ylim = rev(range(rr_raw)),pch=16,xlab = 'Uniform quantiles',ylab = '',main = '(g) Expectation with waxing and waning')
#plot(uniform_quantiles,rr_raw,ylim = rev(range(rr_raw)),pch=16,xlab = 'Uniform quantiles',ylab = "Mean age gap between observations, Ma")
lines(c(0,1),c(t_max,t_min))
dev.off()

pdf('fig_abund_peak.pdf',width = 7, heigh = 3)
plot(NA,NA,ylim = c(0,4.1),xlim=rev(range(t_min,t_max)),xlab = 'Age, Ma', ylab = '',yaxt='n',main = '(b) Waxing and waning abundance')
lines(c(t_min,t_min+8*slot_length),c(4/8,4/8))
polygon(c(t_min,t_min+8*slot_length,t_min+8*slot_length,t_min),c(4/8,4/8,0,0),col="black",border = FALSE)
lines(c(t_min+8*slot_length,t_min+12*slot_length),c(4/4,4/4))
polygon(c(t_min+8*slot_length,t_min+12*slot_length,t_min+12*slot_length,t_min+8*slot_length),c(4/4,4/4,0,0),col="black",border = FALSE)
lines(c(t_min+12*slot_length,t_min+14*slot_length),c(4/2,4/2))
polygon(c(t_min+12*slot_length,t_min+14*slot_length,t_min+14*slot_length,t_min+12*slot_length),c(4/2,4/2,0,0),col="black",border = FALSE)
lines(c(t_min+14*slot_length,t_min+15*slot_length),c(4/1,4/1))
polygon(c(t_min+14*slot_length,t_min+15*slot_length,t_min+15*slot_length,t_min+14*slot_length),c(4/1,4/1,0,0),col="black",border = FALSE)
lines(c(t_min+15*slot_length,t_min+17*slot_length),c(4/2,4/2))
polygon(c(t_min+15*slot_length,t_min+17*slot_length,t_min+17*slot_length,t_min+15*slot_length),c(4/2,4/2,0,0),col="black",border = FALSE)
lines(c(t_min+17*slot_length,t_min+21*slot_length),c(4/4,4/4))
polygon(c(t_min+17*slot_length,t_min+21*slot_length,t_min+21*slot_length,t_min+17*slot_length),c(4/4,4/4,0,0),col="black",border = FALSE)
lines(c(t_min+21*slot_length,t_min+28*slot_length),c(4/8,4/8))
polygon(c(t_min+21*slot_length,t_min+28*slot_length,t_min+28*slot_length,t_min+21*slot_length),c(4/8,4/8,0,0),col="black",border = FALSE)
dev.off()

pdf('fig_abund_flat.pdf',width = 7, heigh = 3)
plot(NA,NA,ylim = c(0,4.1),xlim=rev(range(t_min,t_max)),xlab = 'Age, Ma', ylab = 'Abundance',yaxt='n', main = '(a) Constant abundance')
lines(c(t_min,t_max),c(28/28,28/28))
polygon(c(t_min,t_max,t_max,t_min),c(28/28,28/28,0,0),col="black",border = FALSE)
dev.off()


# this is noise experiment, not reported

data_rv <- c()
data_rv_raw <- c()
slot_length <- (t_max - t_min)/29

for (sk in 1:1000){
  rv1 <- runif(4, min = t_min, max = t_min+8*slot_length) + rnorm(4,0,tm_sd)
  rv2 <- runif(4, min = t_min+8*slot_length, max = t_min+12*slot_length) + rnorm(4,0,tm_sd)
  rv3 <- runif(4, min = t_min+12*slot_length, max = t_min+14*slot_length) + rnorm(4,0,tm_sd)
  rv4 <- runif(4, min = t_min+14*slot_length, max = t_min+15*slot_length) + rnorm(4,0,tm_sd)
  rv5 <- runif(4, min = t_min+15*slot_length, max = t_min+17*slot_length) + rnorm(4,0,tm_sd)
  rv6 <- runif(4, min = t_min+17*slot_length, max = t_min+21*slot_length) + rnorm(4,0,tm_sd)
  rv7 <- runif(4, min = t_min+21*slot_length, max = t_min+28*slot_length) + rnorm(4,0,tm_sd)
  rv <- c(rv1,rv2,rv3,rv4,rv5,rv6,rv7)
  
  rv <- sort(rv,decreasing = TRUE)
  rv <- c(t_max + rnorm(1,0,tm_sd),rv,t_min + rnorm(1,0,tm_sd))
  rv_dif <- diff(-rv)
  
  rv_mov <- mov_agerage(rv_dif,param_ag)
  data_rv <- rbind(data_rv,rv_mov)
  data_rv_raw <- rbind(data_rv_raw,rv)
}

rv_mn <- apply(data_rv,2,mean)
rv_sd <- apply(data_rv,2,sd)

rv_raw <- apply(data_rv_raw,2,mean)

pdf('fig_noise.pdf',width = 7, heigh = 4)
plot(NA,NA,ylim = c(0,0.09),xlim=rev(range(tm_dif_index[ind_range])),xlab = 'Age, Ma', ylab = 'Mean age gap between observations, Ma')
lines(rr_dif_index[ind_range],rv_mn,type = 'l',lwd = 2,lty=5,col = 'grey50')
lines(rr_dif_index[ind_range],rv_mn+rv_sd,type = 'l',lty=5,col = 'grey50')
lines(rr_dif_index[ind_range],rv_mn-rv_sd,type = 'l',lty=5,col = 'grey50')
lines(tm_dif_index[ind_range],tm_mov,type = 'l',lwd = 3,col='black')
dev.off()

pdf('fig_uniformity_noise.pdf',width = 4.5, height = 4.5)
plot(uniform_quantiles,rv_raw,ylim = rev(range(rv_raw)),pch=16,xlab = 'Uniform quantiles',ylab = 'Age, Ma')
lines(c(0,1),c(t_max,t_min))
dev.off()