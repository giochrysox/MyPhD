rm(list=ls())

#loading different libraries
library(ROCR)
library(randomForest)
library(fingerprint)
library(party)
library(rcdk)
library(kernlab)
library(e1071)
source('Classification Summary Function.R')


#Adding Outcome column 
dataset<- read.csv("C:/users/xpyso/Desktop/TS22/AR.csv")
colnames(dataset)[1447]<-"Ids"
colnames(dataset)[1446]<-"Activity"

## Make a random train and test set 80% and 20%
ind<-sample(2,nrow(dataset),replace=TRUE,prob=c(0.8,0.2))
trainset<-dataset[ind==1,]
testset<-dataset[ind==2,]

## Modeling with two algorithms Naive Bayes, Random Forest  
rf_model<-randomForest(as.factor(Activity)~.,data = trainset,ntree=500,proximity=TRUE)
nb_model<-naiveBayes(as.factor(Activity)~.,data = trainset)


#For Naive Bayes model
predict_nb<-predict(nb_model,newdata=testset,type="raw")
pred.nb<-prediction(predict_nb[,2],testset$Activity)
nb.auc<-performance(pred.nb,'tpr','fpr')

#For Random Forest Model
#table(predict(rf), trainset$Outcome)
#print(rf)
#predict_rf<-predict(rf_model,newdata=testset)
#table(predict_rf, testset$Outcome)
predict_rf<-predict(rf_model,newdata=testset,type="prob")
pred.rf<-prediction(predict_rf[,2],testset$Activity)
rf.auc<-performance(pred.rf,'tpr','fpr')

#Plotting the ROC Curve of three algorithm

rfauc<-performance(pred.rf,"auc")@y.values[[1]]
nbauc<-performance(pred.nb,"auc")@y.values[[1]]

Methods<- c('Random Forest','Naive Bayes')
AUCScore<-c(rfauc,nbauc)
data.frame(Methods,AUCScore)
#PLot and adding legend 
plot(rf.auc,col='red',lty=1,main='ROC Curve Comparison of Random Forest V/s Naive Bayes')
plot(nb.auc,col='green',add=TRUE,lty=2)

L<-list(bquote("Random Forest"== .(rfauc)), bquote("Naive Bayes"== .(nbauc)))
legend("bottomright",legend=sapply(L, as.expression),col=c('red','green'),lwd=2,bg="gray",pch=14,text.font=2,cex=0.6)

#Calculate Classification Statistics for each algorithms and put it in a dataframe
rf.label<-predict(rf_model,testset)
rfdata<-classificationsummary(rf.label,testset$Activity,'Active','Inactive')

nb.label<-predict(nb_model,testset)
nbdata<-classificationsummary(nb.label,testset$Activity,'Active','Inactive')


Methods<- c('Random Forest','Naive Bayes')
AUC<-c(rfauc,nbauc)
Accuracy<-c(rfdata$accuracy,nbdata$accuracy)
Sensitivity<-c(rfdata$sensitivity,nbdata$sensitivity)
Specificity<-c(rfdata$specificity,nbdata$specificity)
F1Score <- c(rfdata$F1score,nbdata$F1Score)
data.frame(Methods,Accuracy,Sensitivity,Specificity,AUC,F1Score)
varImpPlot(rf_model,cex=0.5)

