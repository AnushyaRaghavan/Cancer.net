setwd("C:\\Users\\DELL-PC\\Documents\\R\\breast-cancer-wisconsin-R-main")
breastAll <- read.csv("BreastCancer.csv")

breast <- breastAll[,-1]

library(caret)

set.seed(1992)#Ensuring the  same sample at each progeam run
intrain<-createDataPartition(y=breast$Class,p=0.7,list=FALSE)# Partitioning the daata

training   <- breast[ intrain , ] # Training 
validation <- breast[-intrain , ] # validaation

library(rpart)

dtree <- rpart(Class ~ ., data=training, method="class") # Model Buillding
#uses ginis impurity index and recursive approach

library(rpart.plot)

# Graphical view of tree model
rpart.plot(dtree , type=4, extra=1)

# Predicting on validation data
dtree.pred <- predict(dtree, validation, type="class")
dtree.pred[1]  
dtree.pred[7]
dtree.pred[19]
dtree.pred

# Model Evaluation
tbl <- table(dtree.pred, validation$Class,dnn=c("Predicted", "Actual"))#dimension names,gives predicted and actual values of all values in prediction set
tbl
confusionMatrix(tbl)          #this will print actual and predicted values of all.
# cross (diagonal) values are predicted and wrong, matching values are predicted and right


##### Conditional Inference Tree ########

library(party)
fit.ctree <- ctree(Class~., data=training) # model building

plot(fit.ctree, main="Conditional Inference Tree",
     type="simple")

plot(fit.ctree, main="Conditional Inference Tree",
     type="extended")

ctree.pred <- predict(fit.ctree, validation, type="response")
ctree.perf <- table( ctree.pred,validation$Class ,
                    dnn=c("Predicted", "Actual"))

confusionMatrix(ctree.perf)

library(randomForest)
rfmodel <- randomForest(Class ~ ., data = training, importance = TRUE)
rfmodel
rfmodelpred <- predict(rfmodel,validation,type="class")
rfmodelpred

rfpredtbl <- table(rfmodelpred,validation$Class,dnn=c("Predicted", "Actual"))
rfpredtbl
confusionMatrix(rfpredtbl)


