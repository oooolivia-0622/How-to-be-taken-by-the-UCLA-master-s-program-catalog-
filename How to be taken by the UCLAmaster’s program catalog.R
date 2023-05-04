# Add the packages we need 
library(carData)
library(car)
library(olsrr)
library(moments)
library(MASS)
library(boot)
library(relaimpo)
library(ISLR)
library(leaps)
library(broom)
library(robustbase)
library(qpcR)

#Set working path
#setwd('C:/Users/(userName)/Desktop')

#Read the data 
data <- read.csv("Admission.csv")
data <- data.frame(data)

#Data processing

##Preview the correlation coefficient between each regressors and response
corr_variables <- cor(data,use = "pairwise.complete.obs")
corr_variables <- data.frame(corr_variables)
View(corr_variables)
##Generate correlation matrix diagram (function comes from car package)
spm(data,smooth=list(lty.smooth=2, spread = F),main="Matrix of correlations between variables")

#Change the variable name to facilitate subsequent operations

##x1 ~ GRE_Score,x2 ~ TOEFL_Score,x3 ~ University_Rating,x4 ~ SOP,x5 ~ LOR,x6 ~ cGPA,z1 ~ Research,y ~ Chance of dmit
names(data) <-c ("x1","x2","x3","x4","x5","x6","z1","y")
attach(data)
##Because z1 is a categorical variable, it is set as a factor
data$z1 <- as.factor(z1) # set z1 to be a factor
str(z1)
summary(z1)

#Define some functions we need 

##Visualize normality histogram
residplot <- function(fit, nbreaks=10) {  
  z <- rstudent(fit) 
  hist(z, breaks=nbreaks, freq=FALSE, 
       xlab="Studentized Residual", 
       main="Distribution of Errors") 
  rug(jitter(z), col="brown") 
  curve(dnorm(x, mean=mean(z), sd=sd(z)), 
        add=TRUE, col="blue", lwd=2) 
  lines(density(z)$x, density(z)$y, 
        col="red", lwd=2, lty=2) 
  legend("topright", 
         legend = c( "Normal Curve", "Kernel Density Curve"), 
         lty=1:2, col=c("blue","red"), cex=.7) 
} 
##Plot high leverage points
hat.plot <- function(fit){
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  plot(hatvalues(fit),main = "Index Plot of Hat Values")
  abline(h=c(2,3)*p/n,col="red",lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}

#Build a preliminary model
model1 <- lm(y~x1+x2+x3+x4+x5+x6+z1,data = data)
summary(model1)
par(mfrow = c(2,2))
plot(model1)

#Model diagnosis

##Test linearity
crPlots(model1) 
##Test Normality
par(mfrow = c(1,2))
qqPlot(model1,labels=row.names(states),id.method="identity",simulate=TRUE,main="Q-Q Plot") 
residplot(model1) 
##Test Independence of residuals
durbinWatsonTest(model1)
##Test Homoscedasticity
par(mfrow = c(1,1))
ncvTest(model1)  
spreadLevelPlot(model1)
##Test multicollinearity
##Judgment criteria: if the vif value is less than 10, there is no multicollinearity
vif(model1)

#Test abnormal observations on preliminary model
##Outlier
plot(x=fitted(model1),y=rstudent(model1))
abline(h=3,col="red",lty=2)
abline(h=-3,col="red",lty=2)
which(abs(rstudent(model1))>3)
outlierTest(model1)
##High leverage point
hat.plot(model1)
##Hige influence point
cutoff <- 4/(nrow(500-length(model1$coefficients)-2)) 
plot(model1,which=4,cook.levels = cutoff)
abline(h=cutoff,lty=2,col="red")

#Variable selection
##Obtain the results of two variable selections and build models for comparison

##Stepwise Regression
stepwise_model <- step(model1,direction = 'both')
summary(stepwise_model)
par(mfrow = c(2,2))
plot(stepwise_model)
##All-Subsets Regression
leaps <- regsubsets(y ~.,data=data, nvmax =10) 
summary.leaps <- summary(leaps)
par(mfrow=c(2,2))
plot(summary.leaps$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
plot(summary.leaps$adjr2 ,xlab="Number of Variables ",ylab="Adjusted RSq",type="l")
points (6,summary.leaps$adjr2[6], col="red",cex=2,pch=20)
plot(summary.leaps$cp ,xlab="Number of Variables ",ylab="Cp",type='l')
points (6,summary.leaps$cp [6],col="red",cex=2,pch=20)
plot(summary.leaps$bic ,xlab="Number of Variables ",ylab="BIC",type='l')
points(6,summary.leaps$bic [6],col="red",cex=2,pch =20)
coef(leaps ,6)
ols_step_all_possible(model1)
##The model after all-subsets regression
all_subsets_model<-lm(y~x1+x2+x3+x5+x6+z1,data = data)
summary(all_subset_model)
par(mfrow = c(2,2))
plot(all_subset_model)

#Do Transformation

##If we use powerTransform(), we need to ensure all data >0.
min(y)
min(x1)
min(x2)
min(x3)
min(x4)
min(x5)
min(x6)
min(z1+1)
##Transform regressors using function in R
summary(powerTransform(x1))
summary(powerTransform(x2))
summary(powerTransform(x3))
summary(powerTransform(x5))
summary(powerTransform(x6))
summary(powerTransform(z1+1))
##Transform response by box-cox transformation 
bc_model <- boxcox(y~x1+x2+x3+x4+x5+x6+z1,data = data)
lambda <- bc_model$x
lik <- bc_model$y
bc <- cbind(lambda,lik)
bc[order(-lik),]
¦Ë=2 # calculate ¦Ë=2, the best value for y in transformation 

#Remodel after variable selection and transformation
model2<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data)
summary(model2)
par(mfrow = c(2,2))
plot(model2)

#Test abnormal observations on the new model and processing

##Outliers
par(mfrow = c(1,1))
plot(x=fitted(model2),y=rstudent(model2))
abline(h=3,col="red",lty=2)
abline(h=-3,col="red",lty=2)
which(abs(rstudent(model2))>3)
outlierTest(model2)
##High leverage point
hat.plot(model2)
##High influence point
cutoff <- 4/(nrow(500-length(model2$coefficients)-2))
plot(model2,which=4,cook.levels = cutoff)
abline(h=cutoff,lty=2,col="red")
##High influence point Test
##Prove that the influence point in this model have no influence on the model
olsrr::ols_plot_dffits(model2)
olsrr::ols_plot_dfbetas(model2)
olsrr::ols_plot_diagnostics(model2)
##Dealing with outliers and high leverage points
###Remove and re-detect outliers and high leverage values for the first time
data_out_0<-data[-c(10,11,53,65,66,67,118,437),]
model2_0<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_0)
which(abs(rstudent(model2_0))>3)
outlierTest(model2_0)
hat.plot(model2_0)
###The second time
data_out_1<-data_out_0[-c(39,57,62,63),]
model2_1<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_1)
which(abs(rstudent(model2_1))>3)
outlierTest(model2_1)
hat.plot(model2_1)
###The third time
data_out_2<-data_out_1[-c(39,40,106),]
model2_2<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_2)
which(abs(rstudent(model2_2))>3)
outlierTest(model2_2)
hat.plot(model2_2)
###The forth time
data_out_3<-data_out_2[-c(39,40,106),]
model2_3<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_3)
which(abs(rstudent(model2_3))>3)
outlierTest(model2_3)
hat.plot(model2_3)
###The fifth time
data_out_4<-data_out_3[-c(55,79),]
model2_4<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_4)
which(abs(rstudent(model2_4))>3)
outlierTest(model2_4)
hat.plot(model2_5)
###The sixth time
data_out_5<-data_out_4[-c(309),]
model2_5<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_5)
which(abs(rstudent(model2_5))>3)
outlierTest(model2_5)
hat.plot(model2_5)

#Final model
##Use the last deleted outlier data
model3<-lm((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_5)
summary(model3)
par(mfrow = c(2,2))
plot(model3)

#Final model diagnosis

##Test linearity
crPlots(model3) 
##Test Normality
par(mfrow = c(1,2))
qqPlot(model3) 
residplot(model3) 
##Test Error Independence
durbinWatsonTest(model3)
##Test Homoscedasticity
par(mfrow = c(1,1))
ncvTest(model3)  
spreadLevelPlot(model3)
##Test VIF
vif(model3)

#Conclusion

##Values of various statistics
glance(model3)
leaps_final <-regsubsets((y^¦Ë-1)/¦Ë~x1+log(x2)+x3+x5+x6+z1,data = data_out_5, nvmax =10) 
summary_leaps_final<-summary(leaps_final)
View(summary_leaps_final)
PRESS(model3,verbose = TRUE)
##Relative weight between independent variables
###Which variable is obtained from the graph has the greatest impact on the result
crlm <- calc.relimp(model3, type =  "car", rela = TRUE )
crlm
par(mfrow = c(1,1))
plot(crlm)