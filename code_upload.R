library(dlnm)
library(splines) 
library(lubridate)
library(haven)
library(dplyr)
library(purrr)  
library(rgl)
library(mgcv)
library(ggplot2)

rsvmet  <- readRDS("data/rsvmet_exampledata.rds")

rsvmet <- rsvmet %>%
  mutate(year = as.factor(format(date, "%Y"))) %>%
  mutate(month = as.factor(format(date, "%m"))) 


#0 Define an equation for calculating QAIC (Quasi-Akaike Information Criterion)
fqaic <- function(model) {
  loglik <- sum(dpois(model$y, model$fitted.values, log = TRUE))
  phi <- summary(model)$dispersion
  k <- length(coef(model))
  qaic <- (-2 * loglik / phi) + 2 * k  
  return(qaic)
}


## 1.1 Cross-basis function of temperature----------------------
cb.tem = crossbasis(rsvmet$avg_tem, lag=14,   
                    argvar = list(fun = "ns", df =  4),
                    arglag = list(fun = "ns", df =  3))

modeltem = glm(rsvpos ~ cb.tem + ns(rhu, df = 3)  + ns(ssd, df = 3)+ ns(win, df = 4) + year + month + DOW + holiday,
                family = quasipoisson(), data = rsvmet) 

fqaic(modeltem)
summary(modeltem)

pred.tem <- crosspred(cb.tem, modeltem, by =0.1, cumul = TRUE, cen = median(rsvmet$avg_tem, na.rm = TRUE))

#Cumulative Effect Plot
cralltem <- crossreduce(cb.tem,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),type="overall")
plot(cralltem,xlab="tem",ylab="RR",col=2,lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1),ylim=c(0, 8))

#3d plot
plot(pred.tem,ticktype='detailed',border='#3366FF',xlab="tem",ylab="lag(days)",zlab="RR",col='#99FFCC',
     shade = 0.1,cex.lab=1.1,cex.axis=1.1,lwd=1,theta = 20, phi = 25,ltheta = -35)


### 1.1.1 Find the lag and temperature corresponding to the maximum RR value
# Extract the numeric part from column names and convert to numeric
# Extract the RR value matrix
matRRfit <- pred.tem$matRRfit    
matRRlow <- pred.tem$matRRlow    
matRRhigh <- pred.tem$matRRhigh 
numeric_colnames <- as.numeric(gsub("lag", "", colnames(matRRfit)))

# Find the maximum value and its position
max_rr <- max(matRRfit, na.rm = TRUE)  # Maximum RR value
max_pos <- which(matRRfit == max_rr, arr.ind = TRUE)  # Position of the maximum value

# Extract the corresponding lag day and temperature
max_temp <- as.numeric(rownames(matRRfit)[max_pos[1]])  # Corresponding temperature
max_lag <- numeric_colnames[max_pos[2]]  # Corresponding lag day

# Extract the corresponding lower and upper bounds
max_rr_low <- matRRlow[max_pos[1], max_pos[2]]  # Lower bound
max_rr_high <- matRRhigh[max_pos[1], max_pos[2]]  # Upper bound

# Create a result dataframe
result_df <- data.frame(
  Max_RR = max_rr,
  Lower_Bound = max_rr_low,
  Upper_Bound = max_rr_high,
  Temperature = max_temp,
  Lag_Day = max_lag
)

print(result_df)


### 1.1.2 View RR values at different lags for a specific temperature
# Determine the index of the target temperature
target_temp <- 29.6
temp_index <- which.min(abs(as.numeric(rownames(pred.tem1$matRRfit)) - target_temp))

# Extract RR matrix and confidence interval matrices
matRRfit <- pred.tem$matRRfit    # Instantaneous effect matrix
matRRlow <- pred.tem$matRRlow    # Lower bound matrix
matRRhigh <- pred.tem$matRRhigh  # Upper bound matrix

# Retrieve values for different lag days corresponding to the target temperature
rr_values <- matRRfit[temp_index, ]      # RR values
rr_low <- matRRlow[temp_index, ]         # Lower bound
rr_high <- matRRhigh[temp_index, ]       # Upper bound

# Create a dataframe to display results
result_df <- data.frame(
  Lag_Day = as.numeric(gsub("lag", "", colnames(matRRfit))),  # Lag days
  RR = rr_values,
  Lower_Bound = rr_low,
  Upper_Bound = rr_high
)

# Display the results
print(result_df)

## 2.1 Cross-basis function of relative humidity----------------------
cb.rhu = crossbasis(rsvmet$rhu, lag=14, 
                    argvar = list(fun = "ns", df = 3),
                    arglag = list(fun = "ns", df =3))

modelrhu = glm(rsvpos ~ cb.rhu  + ns(avg_tem, df = 3) + ns(win, df = 3)+ year + month + DOW + holiday,
                family = quasipoisson(), data = rsvmet)

fqaic(modelrhu)

summary(modelrhu)         

pred.rhu <- crosspred(cb.rhu, modelrhu, by = 0.5 , cumul = TRUE, cen = median(rsvmet$rhu, na.rm = TRUE))  

#Cumulative Effect Plot
crallrhu <- crossreduce(cb.rhu,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),type="overall")
plot(crallrhu,xlab="rhu",ylab="RR",col=2,lwd=2,cex.lab=1.2,cex.axis=1.2,mar=c(1,2,0,1),ylim=c(0,4))

#3d plot
plot(pred.rhu,ticktype='detailed',border='#3366FF',xlab="Relative humidity",ylab="lag(days)",zlab="RR",col='#99FFCC',
     shade = 0.1,cex.lab=1.0,cex.axis=0.9,lwd=1,theta = 20, phi = 25,ltheta = -35) 

### 2.1.1 Find the lag and value corresponding to the maximum RR for rhu
# Extract the RR value matrix
matRRfit <- pred.rhu1$matRRfit    # Central RR value matrix
matRRlow <- pred.rhu1$matRRlow    # RR lower bound matrix
matRRhigh <- pred.rhu1$matRRhigh  # RR upper bound matrix

# Extract lag day information from column names
numeric_colnames <- as.numeric(gsub("lag", "", colnames(matRRfit)))

# Find the maximum RR value and its position
max_rr <- max(matRRfit, na.rm = TRUE)  # Maximum RR value
max_pos <- which(matRRfit == max_rr, arr.ind = TRUE)  # Position of the maximum value

# Extract the corresponding humidity value and lag day
max_rhu <- as.numeric(rownames(matRRfit)[max_pos[1]])  # Corresponding humidity value
max_lag <- numeric_colnames[max_pos[2]]  # Corresponding lag day

# Extract the lower and upper bounds corresponding to the maximum value
max_rr_low <- matRRlow[max_pos[1], max_pos[2]]  # Lower bound
max_rr_high <- matRRhigh[max_pos[1], max_pos[2]]  # Upper bound

# Create a result dataframe
result_rhu <- data.frame(
  Max_RR = max_rr,
  Lower_Bound = max_rr_low,
  Upper_Bound = max_rr_high,
  Relative_Humidity = max_rhu,
  Lag_Day = max_lag
)

# Display the result
print(result_rhu)

### 2.1.2 Find the matRR values for a specific rhu at different lags
# Assume the target relative humidity (rhu) is 50%
target_rhu <- 50

# Find the row index corresponding to the target rhu
rhu_index <- which.min(abs(as.numeric(rownames(pred.rhu1$matRRfit)) - target_rhu))

# Extract matRRfit (instantaneous effect matrix) and its confidence interval matrices
matRRfit <- pred.rhu1$matRRfit    # Instantaneous effect matrix
matRRlow <- pred.rhu1$matRRlow    # Lower bound matrix
matRRhigh <- pred.rhu1$matRRhigh  # Upper bound matrix

# Extract RR values and confidence intervals corresponding to the target rhu
rr_values <- matRRfit[rhu_index, ]      # Instantaneous RR values
rr_low <- matRRlow[rhu_index, ]         # RR lower bound
rr_high <- matRRhigh[rhu_index, ]       # RR upper bound

# Create a result dataframe
result_df <- data.frame(
  Lag_Day = as.numeric(gsub("lag", "", colnames(matRRfit))),  # Lag days
  RR = rr_values,
  Lower_Bound = rr_low,
  Upper_Bound = rr_high
)

# Display the result
print(result_df)


## 3 Calculation of AF and AN----------------------
# AF-cold
attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),
       range=c(-15,median(rsvmet$avg_tem, na.rm = TRUE)),dir="forw")*100
quantile(attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,sim=T,nsim=5000,
                cen=median(rsvmet$avg_tem, na.rm = TRUE),range=c(-15,median(rsvmet$avg_tem, na.rm = TRUE))),c(0.025,0.975))

# AF-extreme cold
attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),
       range=c(-15,quantile(rsvmet$avg_tem, probs = 0.01)),dir="forw")*100
quantile(attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,sim=T,nsim=5000,
                cen=median(rsvmet$avg_tem, na.rm = TRUE), range=c(-15,quantile(rsvmet$avg_tem, probs = 0.01))),c(0.025,0.975))*100

# AF-moderate cold
attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),
       range=c(quantile(rsvmet$avg_tem, probs = 0.01),median(rsvmet$avg_tem, na.rm = TRUE)),dir="forw")*100
quantile(attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,sim=T,nsim=5000,
                cen=median(rsvmet$avg_tem, na.rm = TRUE)),c(0.025,0.975))*100

# AN-cold
attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),
       range=c(-15,median(rsvmet$avg_tem, na.rm = TRUE)),dir="forw",type="an")
quantile(attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,sim=T,nsim=5000,type="an",
                cen=median(rsvmet$avg_tem, na.rm = TRUE),range=c(-15,median(rsvmet$avg_tem, na.rm = TRUE))),c(0.025,0.975))

# AN-extreme cold
attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),
       range=c(-15,quantile(rsvmet$avg_tem, probs = 0.01)),dir="forw",type="an")
quantile(attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,sim=T,nsim=5000,type="an",
                cen=median(rsvmet$avg_tem, na.rm = TRUE), range=c(-15,quantile(rsvmet$avg_tem, probs = 0.01))),c(0.025,0.975))

# AN-moderate cold
attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,cen=median(rsvmet$avg_tem, na.rm = TRUE),
       range=c(quantile(rsvmet$avg_tem, probs = 0.01),median(rsvmet$avg_tem, na.rm = TRUE)),dir="forw",type="an")
quantile(attrdl(rsvmet$avg_tem,cb.tem,rsvmet$rsvpos,modeltem,sim=T,nsim=5000,type="an",
                cen=median(rsvmet$avg_tem, na.rm = TRUE)),c(0.025,0.975))

# AF-wetness
attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),
       range=c(cen=median(rsvmet$rhu, na.rm = TRUE),100),dir="forw")*100
quantile(attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,sim=T,nsim=5000,
                cen=median(rsvmet$rhu, na.rm = TRUE),range=c(median(rsvmet$rhu, na.rm = TRUE)),100),c(0.025,0.975))*100

# AF-extreme wetness
attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),
       range=c(quantile(rsvmet$rhu, probs = 0.99),100),dir="forw")*100
quantile(attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,sim=T,nsim=5000,
                cen=median(rsvmet$rhu, na.rm = TRUE), range=c(quantile(rsvmet$rhu, probs = 0.99),100)),c(0.025,0.975))*100

# AF-moderate wetness
attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),
       range=c(median(rsvmet$rhu, na.rm = TRUE),quantile(rsvmet$rhu, probs = 0.99)),dir="forw")*100
quantile(attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,sim=T,nsim=5000,
                cen=median(rsvmet$rhu, na.rm = TRUE), range=c(median(rsvmet$rhu, na.rm = TRUE)),quantile(rsvmet$rhu, probs = 0.99)),c(0.025,0.975))*100

# AN-wetness
attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),type="an",
       range=c(cen=median(rsvmet$rhu, na.rm = TRUE),100),dir="forw")
quantile(attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,sim=T,nsim=5000,type="an",
                cen=median(rsvmet$rhu, na.rm = TRUE),range=c(median(rsvmet$rhu, na.rm = TRUE)),100),c(0.025,0.975))

# AN-extreme wetness
attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),type="an",
       range=c(quantile(rsvmet$rhu, probs = 0.99),100),dir="forw")
quantile(attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,sim=T,nsim=5000,type="an",
                cen=median(rsvmet$rhu, na.rm = TRUE), range=c(quantile(rsvmet$rhu, probs = 0.99),100)),c(0.025,0.975))

# AN-moderate wetness
attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,cen=median(rsvmet$rhu, na.rm = TRUE),type="an",
       range=c(median(rsvmet$rhu, na.rm = TRUE),quantile(rsvmet$rhu, probs = 0.99)),dir="forw")
quantile(attrdl(rsvmet$rhu,cb.rhu,rsvmet$rsvpos,modelrhu,sim=T,nsim=5000,type="an",
                cen=median(rsvmet$rhu, na.rm = TRUE), range=c(median(rsvmet$rhu, na.rm = TRUE)),quantile(rsvmet$rhu, probs = 0.99)),c(0.025,0.975))

## 4 Plotting -------
### 4.1 Cumulative ----------------------
# Start plotting to a TIFF file
tiff("plot/cumR1.tiff", width = 10, height = 5, units = "in", res = 330)
# Set three plots to be displayed in one row
par(mfrow=c(1, 2), mar = c(5, 5, 3, 1.5))  
plot(cralltem1, xlab="Mean temperature (°C)", ylab="RR", col="#1f90de", ci.arg = list(col= "#D3EAFA"),
     lwd=1.8, cex.lab=1.35, cex.axis=1.2, mar=c(1,2,0,1), ylim=c(0, 7))
# Add character "(A)" in the top-left corner
text(x = -17, y = 8.11, labels = "(A)", pos = 4, cex = 1.3, col = "black", xpd = TRUE) 
plot(crallrhu1, xlab="Relative humidity (%)", ylab="RR", col="#1f90de", ci.arg = list(col="#D3EAFA"),
     lwd=2, cex.lab=1.3, cex.axis=1.2, mar=c(1,2,0,1), ylim=c(0, 5))
# Add character "(B)" in the top-left corner
text(x = 10, y = 5.8, labels = "(B)", pos = 4, cex = 1.3, col = "black", xpd = TRUE) 
# Close the TIFF file device
dev.off()

### 4.2 3D ----------------------
# Start plotting to a TIFF file
tiff("plot/3DR1.tiff", width = 9.5, height = 5, units = "in", res = 320)
par(mfrow=c(1, 2), mar=c(1.5,4.5,1.5,2))  # 1 row, 2 columns
plot(pred.tem1, ticktype='detailed', border='#3366FF', xlab="Mean temperature (°C)", ylab="Lag days", zlab="RR", col='#99FFCC',
     shade = 0.1, cex.lab=1.0, cex.axis=1, lwd=1, theta = 20, phi = 25, ltheta = -35) 
# Add character "(A)" in the top-left corner
text(x = -0.8, y = 0.6, labels = "(A)", adj = c(0.5, 0.5), col = 'black', cex=1.3, xpd = TRUE) 
plot(pred.rhu1, ticktype='detailed', border='#3366FF', xlab="Relative humidity (%)", ylab="Lag days", zlab="RR", col='#99FFCC',
     shade = 0.1, cex.lab=1.0, cex.axis=1, lwd=1, theta = 20, phi = 25, ltheta = -35) 
text(x = -0.8, y = 0.6, labels = "(B)", adj = c(0.5, 0.5), col = 'black', cex=1.3, xpd = TRUE)  
# Close the TIFF file device
dev.off()

### 4.3 Variable-specific ----------------------
temperature_1th <- quantile(rsvmet$avg_tem, probs = 0.01) # 1st percentile temperature
temperature_10th <- quantile(rsvmet$avg_tem, probs = 0.1) # 10th percentile temperature
temperature_99th <- quantile(rsvmet$avg_tem, probs = 0.99) # 99th percentile temperature
temperature_90th <- quantile(rsvmet$avg_tem, probs = 0.9) # 90th percentile temperature
rhu_1th <- quantile(rsvmet$rhu, probs = 0.01) # 1st percentile relative humidity
rhu_10th <- quantile(rsvmet$rhu, probs = 0.1) # 10th percentile relative humidity
rhu_99th <- quantile(rsvmet$rhu, probs = 0.99) # 99th percentile relative humidity
rhu_90th <- quantile(rsvmet$rhu, probs = 0.9) # 90th percentile relative humidity
pred.tem10 <- crosspred(cb.tem, modeltem, by = 0.5, cumul = TRUE, cen = median(rsvmet$avg_tem, na.rm = TRUE))
pred.tem1 <- crosspred(cb.tem, modeltem, by = 0.12, cumul = TRUE, cen = median(rsvmet$avg_tem, na.rm = TRUE))
pred.tem90 <- crosspred(cb.tem, modeltem, by = 0.2, cumul = TRUE, cen = median(rsvmet$avg_tem, na.rm = TRUE))
pred.tem99 <- crosspred(cb.tem, modeltem, by = 0.5, cumul = TRUE, cen = median(rsvmet$avg_tem, na.rm = TRUE))
pred.rhuspe <- crosspred(cb.rhu, modelrhu, by = 0.1, cumul = TRUE, cen = median(rsvmet$rhu, na.rm = TRUE))
#########################################

plot(pred.rhuspe, "slices", type = "p", pch = 16, col ="#DB726B", cex = 1.4, var = 99,
     ci = "bars", ci.arg = list(col = "#ffc2b9", lwd = 2.5, lty = 1), 
     ylab = " ", xlab = "",  # Disable default xlab
     main = "", cex.main = 1)

# Start plotting
tiff("plot/varspecific1R1.tiff", width =12, height = 6, units = "in", res = 320)
par(mgp = c(3, 0.8, 0), mfrow=c(2, 4), mar=c(3.8,3.4,3,1))  # 2 rows, 4 columns

# Extremely low temperature
plot(pred.tem1, "slices", type = "p", pch = 16, col = '#519bd7', cex = 1.4, var = 0.54,
     ci.arg = list(col = "#add1ef", lwd = 2.5, lty=1),
     ci = "bars", ylab = " ", xlab = "", main = " ", cex.main = 1)
mtext("Lag (days)", side = 1, line = 2.2, cex = 0.7)
mtext("RR", side = 2, line = 2.15, cex = 0.7)
text(x = -3.7, y = 1.311, labels = "(A)", pos = 4, cex = 1.4, col = "black", xpd = TRUE)
title("Extremely low temperature 
(0.6°C)", line = 0, cex.main = 1.2, font.main = 2)

# Low temperature
plot(pred.tem10, "slices", type = "p", pch = 16, col = '#519bd7', cex = 1.4, var = 5.5,
     ci.arg = list(col = "#add1ef", lwd = 2.5, lty=1),
     ci = "bars", ylab = " ", xlab = "", main = " ", cex.main = 1)
mtext("Lag (days)", side = 1, line = 2.2, cex = 0.7)
mtext("RR", side = 2, line = 2.2, cex = 0.7)
text(x = -3.7, y = 1.378, labels = "(B)", pos = 4, cex = 1.4, col = "black", xpd = TRUE)
title("Low temperature 
(5.5°C)", line = 0, cex.main = 1.15, font.main = 2)

# High temperature
plot(pred.tem90, "slices", type = "p", pch = 16, col = "#DB726B", cex = 1.4, var = 29.6, 
     ci.arg = list(col = "#ffc2b9", lwd = 2.5, lty=1),
     xlab = "", ylim = c(0.75, 1.4), ci = "bars", ylab = " ", main = " ", cex.main = 1)
mtext("Lag (days)", side = 1, line = 2.2, cex = 0.7)
mtext("RR", side = 2, line = 2.15, cex = 0.7)
text(x = -3.7, y = 1.496, labels = "(C)", pos = 4, cex = 1.4, col = "black", xpd = TRUE)
title("High temperature 
(29.6°C)", line = 0, cex.main = 1.2, font.main = 2)

# Extremely high temperature
plot(pred.tem99, "slices", type = "p", pch = 16, col = "#DB726B", cex = 1.4, var = 33.5, 
     ci = "bars", ci.arg = list(col = "#ffc2b9", lwd = 2.5, lty = 1), 
     ylab = " ", xlab = "", main = "", cex.main = 1)
mtext("Lag (days)", side = 1, line = 2.2, cex = 0.7)
mtext("RR", side = 2, line = 2.2, cex = 0.7)
text(x = -3.7, y = 1.807, labels = "(D)", pos = 4, cex = 1.4, col = "black", xpd = TRUE)
title("Extremely high temperature
(33.5°C)", line = 0, cex.main = 1.15, font.main = 2)

# Extremely low relative humidity
plot(pred.rhuspe, "slices", type = "p", pch = 16, col = '#519bd7', cex = 1.4, var = 43,
     ci = "bars", ci.arg = list(col = "#add1ef", lwd = 2.5, lty = 1), 
     ylab = " ", xlab = "", main = "", cex.main = 1)
mtext("Lag (days)", side = 1, line = 2.2, cex = 0.7)
mtext("RR", side = 2, line = 2.2, cex = 0.7)
text(x = -3.7, y = 1.118, labels = "(E)", pos = 4, cex = 1.4, col = "black", xpd = TRUE)
title("Extremely low relative humidity 
(42.8%)", line = 0, cex.main = 1.15, font.main = 2)     

# Low relative humidity
plot(pred.rhuspe, "slices", type = "p", pch = 16, col = '#519bd7', cex = 1.4, var = 58,
     ci = "bars", ci.arg = list(col = "#add1ef", lwd = 2.5, lty = 1), 
     ylab = " ", xlab = "", main = "", cex.main = 1)
mtext("Lag (days)", side = 1, line = 2.2, cex = 0.7)
mtext("RR", side = 2, line = 2.2, cex = 0.7)
text(x = -3.7, y = 1.057, labels = "(F)", pos = 4, cex = 1.4, col = "black", xpd = TRUE)
title("Low relative humidity 
(58.0%)", line = 0, cex.main = 1.15, font.main = 2)

# Close TIFF file device
dev.off()
