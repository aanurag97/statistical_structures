#Code for Assignment 1 of Statistical Structures in Data course of
#PGDBA 2022-24 by Anurag Shukla (22BM6JP08)

#Importing required Libraries
library(DescTools)
library(modeest)
library(moments)
#setwd("C:/Users/Anurag Shukla/Documents/pgdba/study/ISI/SSD/assignment1")
#Generating Data and corresponding Histograms
set.seed(08)
X_poisson <- rpois(1000, lambda=3)
hist(X_poisson, main="Histogram of Poisson distribution", col = 'skyblue3')
set.seed(08)
X_geo <- rgeom(1000, prob = 0.2)
hist(X_geo, main="Histogram of Geometric distribution", col = 'skyblue3')
set.seed(08)
Y_exp <- rexp(1000, rate = 1)
hist(Y_exp, main="Histogram of Exponential distribution", col = 'skyblue3')
set.seed(08)
Y_gamma <- rgamma(1000, shape = 5)
hist(Y_gamma, main="Histogram of Gamma distribution", col = 'skyblue3')


#Question 1 - X as Poisson and Y as Gamma Distribution
#Data generated for all 4 distributions
#Frequency Distribution Table for data exported as CSV file
write.csv(table(X_geo), file="X_geo_table.csv")
write.csv(table(X_poisson), file="X_poisson_table.csv")
range(Y_exp)
breaks <- seq(0, 7, by=1)
fd_exp <- cut(Y_exp, breaks, right=F) #defining intervals
write.csv(table(fd_exp), file="Y_exp_table.csv")
range(Y_gamma)
breaks <- seq(0, 18, by=2)
fd_gamma <- cut(Y_gamma, breaks, right=F) #defining intervals
write.csv(table(fd_gamma), file="Y_gamma_table.csv")

#Computing Data Measures for generated Data
#Mean
X_geo_mean <- mean(X_geo)
X_poisson_mean <- mean(X_poisson)
Y_exp_mean <- mean(Y_exp)
Y_gamma_mean <- mean(Y_gamma)
#Median
X_geo_median <- median(X_geo)
X_poisson_median <- median(X_poisson)
Y_exp_median <- median(Y_exp)
Y_gamma_median <- median(Y_gamma)
#Mode
X_geo_mode <- Mode(X_geo)
X_poisson_mode <- Mode(X_poisson)
Y_exp_mode <- mlv(Y_exp, method = "meanshift")
Y_gamma_mode <- mlv(Y_gamma, method = "meanshift")
#Range
range_geo <- range(X_geo)
range_poisson <- range(X_poisson)
range_gamma <- range(Y_gamma)
range_exp <- range(Y_exp)
#Interquartile Range
iqr_geo <- IQR(X_geo)
iqr_poisson <- IQR(X_poisson)
iqr_gamma <- IQR(Y_gamma)
iqr_exp <- IQR(Y_exp)
#Standard Deviation
sd_geo <- sd(X_geo)
sd_poisson <- sd(X_poisson)
sd_exp <- sd(Y_exp)
sd_gamma <- sd(Y_gamma)
#Skewness
skew_geo <- skewness(X_geo)
skew_poisson <- skewness(X_poisson)
skew_exp <- skewness(Y_exp)
skew_gamma <- skewness(Y_gamma)
#Kurtosis
kurt_geo <- kurtosis(X_geo)
kurt_poisson <- kurtosis(X_poisson)
kurt_exp <- kurtosis(Y_exp)
Kurt_gamma <- kurtosis(Y_gamma)

#Box and Whisker Plots
boxplot(X_poisson, Y_gamma,
        main = 'Box-and-whisker plot',
        names = c('Poisson','Gamma'), col = c('green','orange'), ylab = 'Variable values')
boxplot(X_geo, X_poisson, Y_gamma, Y_exp, 
        main = 'Box-and-whisker Plot of all variables',
        names = c('Geomertric','Poisson','Gamma','Exponential'),
        col = c('red','yellow','orange','green'))


#Question 2 - Fitting Normal Distribution on generated Data
#Used Poisson and Gamma Distributions
#Poisson Distribution
set.seed(08)
poi_norm <- rnorm(1000, mean = X_poisson_mean, sd = sd_poisson)
approx_poi_norm <- round(poi_norm, digits = 0) #discrete normal distribution
fd_poi_norm <- table(approx_poi_norm)
#range(poi_norm)
#breaks <- seq(-3.5,10.5, by=1)
#fd_poi_norm <- cut(poi_norm, breaks, right = F)
write.csv(fd_poi_norm, file="poi_norm_table.csv")
plot(ecdf(X_poisson), col = "blue", main = "ECDF comparisson for Poisson")
plot(ecdf(approx_poi_norm), col = "red", add = T)
plot(ecdf(poi_norm), col = "green", add = T)
qqnorm(X_poisson)
qqline(X_poisson)
#Gamma Distribution
set.seed(08)
gamma_norm <- rnorm(1000, mean = Y_gamma_mean, sd = sd_gamma)
range(gamma_norm)
breaks <- seq(-4,18, by=2)
fd_gamma_norm <- cut(gamma_norm, breaks, right = F)
write.csv(table(fd_gamma_norm), file="gamma_norm_table.csv")
plot(ecdf(Y_gamma), col = "blue", main = "ECDF comparisson for Gamma")
plot(ecdf(gamma_norm), col = "red", add = T)
qqnorm(Y_gamma)
qqline(Y_gamma)

#Kolmogorov-Smirnov test for goodness-of-fit with density plots
poi_fit <- ks.test(X_poisson,approx_poi_norm)
gamma_fit <- ks.test(Y_gamma,gamma_norm)
seq_poi <- seq(min(X_poisson), max(X_poisson), length(1000))
set.seed(08)
d_poi <- dnorm(seq_poi, mean = X_poisson_mean, sd = sd_poisson)
hist(X_poisson, prob = T, main = "Poisson with fitted Normal curve")
lines(seq_poi, d_poi, col = "red", lwd = 2)
seq_gamma = seq(min(Y_gamma), max(Y_gamma), length(1000))
set.seed(08)
d_gamma <- dnorm(seq_gamma, mean = Y_gamma_mean, sd = sd_gamma)
hist(Y_gamma, prob = T, ylim = c(0,0.20), main = "Gamma with fitted Normal curve")
lines(seq_gamma, d_gamma, col = "red", lwd = 2)


#Question 3 - Regression on Gamma(Y1) and Exponential(Y2) Data
#Computing regressiong line
regress <- lm(Y_gamma ~ Y_exp)
plot(Y_exp, Y_gamma, pch = "+", main = 'Scatter Plot with Regression Line')
abline(regress, lwd = 3, col = 'red')

#Sum of Squares
sse <- sum((fitted(regress)-Y_gamma)^2)
ssr <- sum((fitted(regress)-mean(Y_gamma))^2)
tss <- sse + ssr
#F-Test and Anova Table
MSR <- ssr/(summary(regress)$df[1] - 1)
MSE <- sse/summary(regress)$df[2]
fstat <- MSR/MSE
anova(regress)
write.csv(anova(regress), file="anova.csv")

#Coefficient of Determination (R-squared value)
rsquared <- ssr/tss
par(mfrow = c(2,2))
plot(regress)
par(mfrow = c(1,2))
plot(predict(regress), residuals(regress), pch = 20, col = 'red',
     xlab = 'Predicted Values', ylab = 'Residuals')
abline(regress, lwd = 3, col = 'blue')
plot(predict(regress), rstudent(regress), pch = 20, col = 'green',
     xlab = 'Predicted Values', ylab = 'Studentized Residuals')

#Influential points w.r.t Cooke's Distance
cooks_dist <- cooks.distance(regress)
influential_points <- cooks_dist[(cooks_dist > (3 * mean(cooks_dist, na.rm = TRUE)))]
write.csv(influential_points, file="influential_points.csv")
