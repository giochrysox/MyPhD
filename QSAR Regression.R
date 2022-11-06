# Remove all the history variables 
rm(list=ls()) 
#----------------------------------------------------
library(QSARdata)
library(caret)
library(gplots)
library(latticeExtra)
library(lattice)
library(ggplot2)

#----------------------------------------------------
# Function to compute RMSE and R2 
r2se<- function (obs,pred){
  rmse<-(mean((obs -pred)^2))^0.5 
  ssr<-sum((obs - pred)^2)
  sst<-sum((obs - mean(obs))^2)
  R2<-1-(ssr/sst)
  output<-list(RMSE=rmse,RSquared=R2)
  return(output)
}

#----------------------------------------------------
# Function to plot qsar results with ggplot2 
plotsar <-function(results){
  
  myplot<-ggplot(results,aes(x=observed,y=predicted,color=factor(set),shape=factor(set)))+geom_point()+scale_colour_manual(values=c("blue", "yellow"))+
    theme(axis.ticks.y=element_line(size=1),axis.text.y=element_text(size=16),axis.ticks.length=unit(0.25,"cm"), 
          panel.background = element_rect(fill = "black", colour = NA))+
    theme(axis.title=element_text(size=12,face="bold",colour = 'white'),plot.background = element_rect(colour = 'black', fill = 'black'))+
    theme( legend.title = element_text(size = 12, face = "bold", hjust = 0, colour = 'white'),
           legend.background=element_rect(color="black",fill="black"),legend.text=element_text(size=12,color="white",face="bold"),
           axis.title.x = element_text(size = 12, colour = 'white', vjust = 1) ,axis.title.y = element_text(size = 12, colour = 'white'))+
    ggtitle("Predicted v/s test set QSAR results")+
    theme(plot.title=element_text(lineheight=.8,face="bold",color="white",size=14))+stat_smooth(span = 0.9)
  
  
  return(myplot)
}


#----------------------------------------------------
# Loading data from QSARdata package.
# If you find error in loading the package download the tar file and 
# install via local install 


q.data<-read.csv("C:/Users/xpyso/Desktop/TS23/Karthikeyan-MP.csv")
q.data<-na.omit(q.data) # Omit rows for which data is not present

#----------------------------------------------------
# Data cleaning
descs <- q.data[, !apply(q.data, 2, function(x) any(is.na(x)) )]
#Near constant columns
descs <- descs[, !apply( descs, 2, function(x) length(unique(x)) == 1 )]

r2 <- which(cor(descs[2:185])^2 > .29, arr.ind=TRUE)

r2 <- r2[ r2[,1] > r2[,2] , ]
d <- descs[, -unique(r2[,2])]


#----------------------------------------------------
# Normalizing the data 
# Power transform is a useful data transformation technique used to stabilize variance, make the data more normal distribution-like, 
# improve the validity of measures # of association such as the Pearson correlation between variables and for other data stabilization procedures

T <- preProcess(d[,2:dim(d)[2]],method = "BoxCox")
data <-predict(T,d[,2:dim(d)[2]])

#----------------------------------------------------
# Create training and test set
ind<-sample(2,nrow(data),replace=TRUE,prob=c(0.8,0.2))
s<-dim(data)[2]
trainset<-data[ind==1,2:s]
testset<-data[ind==2,2:s]

#----------------------------------------------------
# Use OLS to model the data
y.test <- testset$MeltingPoint
q.ols <- lm(MeltingPoint ~ . , data=trainset)

#----------------------------------------------------
# Summarizing results
summary(q.ols)
pred.ols.train<-predict(q.ols,trainset) #predict train set
pred.ols.test<-predict(q.ols,testset) #predict test set

rlmValues <- data.frame(obs = y.test,pred = pred.ols.test)
r2se(rlmValues$obs,rlmValues$pred)

#----------------------------------------------------
# Plot the results 
y<-xyplot(trainset$MeltingPoint ~ pred.ols.train, 
          type = c('g', 'p'),  xlab = "Predicted", ylab = "Observed",
          panel = function(x,y, ...){ 
            panel.xyplot(x,y,...)
            panel.lines(x, predict(q.ols), col = 'black', lwd = 2) 
          } 
) 


y+as.layer(xyplot(y.test ~ pred.ols.test,pch=17,col="black"))
