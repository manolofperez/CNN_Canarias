library(abc)
setwd("/Users/manolo/Desktop/dartr_Mario/enviar/CNN/")

######################################################
#Euphorbia                                           #
######################################################

###################Model comparison###################

sust<-read.table("Euphorbia/3rdRound/Results/Pred_testSet_CalibratedModelPredictions.txt")
models<-rep(1:3,each=5000)
emp<-read.table("Euphorbia/3rdRound/Results/Pred_Emp_CalModelPredictions.txt")
emp<-apply(emp, 2, FUN = median)

#cv.modsel <- cv4postpr(models, sust, nval=10, tol =c(.05,.01,.005,.002,.001), method="rejection")
cv.modsel <- cv4postpr(models, sust, nval=100, tol =.05, method="rejection")
Rej.05<-postpr(emp, models, sust, tol = .05, method = "rejection")
summary(cv.modsel)
summary(Rej.05)

#################Parameter estimation#################

parameters<-read.table("Euphorbia/3rdRound/Results/parSSHBackCol.txt")
colnames(parameters)=c("Theta","Tm","Tb","T1","T2","T3","T4","T5","T6","T7","Ne","NGC","NTor","NToc","NLP","NLG","NEH","NC","m","FoundedSizeRatio")
parameters=within(parameters,rm(Theta,T7))
parameters$NGC=parameters$NGC*parameters$Ne
parameters$NTor=parameters$NTor*parameters$Ne
parameters$NToc=parameters$NToc*parameters$Ne
parameters$NLP=parameters$NLP*parameters$Ne
parameters$NLG=parameters$NLG*parameters$Ne
parameters$NEH=parameters$NEH*parameters$Ne
parameters$NC=parameters$NC*parameters$Ne

sust<-read.table("Euphorbia/3rdRound/Results/testSet_ParameterPredictions.txt")
emp<-read.table("Euphorbia/3rdRound/Results/Emp_ParametersPredictions.txt")
emp<-apply(emp, 2, FUN = median)

#cv.parest <- cv4abc(parameters, sust, nval=10, tol =c(.05,.01,.005,.002,.001), method = "rejection")
cv.parest <- cv4abc(parameters, sust, nval=100, tol =c(.005), method = "rejection")

REJ.parest.005 <- abc(emp, parameters, sust, tol = .005,method = "rejection")

summary(cv.parest)

summary(REJ.parest.005)

######################################################
#Kleinia                                             #
######################################################

###################Model comparison###################

sust<-read.table("Kleinia/3rdRound/Results/Pred_testSet_CalibratedModelPredictions.txt")
models<-rep(1:3,each=5000)
emp<-read.table("Kleinia/3rdRound/Results/Pred_Emp_CalModelPredictions.txt")
emp<-apply(emp, 2, FUN = median)

#cv.modsel <- cv4postpr(models, sust, nval=10, tol =c(.05,.01,.005,.002,.001), method="rejection")
cv.modsel <- cv4postpr(models, sust, nval=100, tol =.05, method="rejection")
Rej.05<-postpr(emp, models, sust, tol = .05, method = "rejection")
summary(cv.modsel)
summary(Rej.05)

#################Parameter estimation#################

parameters<-read.table("Kleinia/3rdRound/Results/parEastWest.txt")
colnames(parameters)=c("Theta","Tm","Tb","T1","T2","T3","T4","T5","T6","T7","Ne","NGC","NTor","NToc","NLP","NLG","NEH","NC","m","FoundedSizeRatio")
parameters=within(parameters,rm(Theta,Tm,Tb,NC,m))
parameters$NGC=parameters$NGC*parameters$Ne
parameters$NTor=parameters$NTor*parameters$Ne
parameters$NToc=parameters$NToc*parameters$Ne
parameters$NLP=parameters$NLP*parameters$Ne
parameters$NLG=parameters$NLG*parameters$Ne
parameters$NEH=parameters$NEH*parameters$Ne

sust<-read.table("Kleinia/3rdRound/Results/testSet_ParameterPredictions.txt")
emp<-read.table("Kleinia/3rdRound/Results/Emp_ParametersPredictions.txt")
emp<-apply(emp, 2, FUN = median)

#cv.parest <- cv4abc(parameters, sust, nval=10, tol =c(.05,.01,.005,.002,.001), method = "rejection")
cv.parest <- cv4abc(parameters, sust, nval=100, tol =c(.005), method = "rejection")



REJ.parest.005 <- abc(emp, parameters, sust, tol = .005,method = "rejection")

summary(cv.parest)

summary(REJ.parest.005)
