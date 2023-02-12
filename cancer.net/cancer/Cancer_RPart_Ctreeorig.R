breastAll <- read.csv("BreastCancer.csv")

breast <- breastAll[,-1]

library(caret)
set.seed(1992)
intrain<-createDataPartition(y=breast$Class,p=0.7,list=FALSE)

training   <- breast[ intrain , ]
validation <- breast[-intrain , ]

library(rpart)
dtree <- rpart(Class ~ ., data=training,  method="class")

dtree$cptable
plotcp(dtree)

optim_CP <- dtree$cptable[which.min(dtree$cptable[,"xerror"]),"CP"]
optim_CP

dtree.pruned <- prune(dtree, cp= 0.01775148 ) 


library(rpart.plot)

rpart.plot(dtree , type=4, extra=1,main="Before Pruning")

rpart.plot(dtree.pruned, type = 4, extra = 1,main="After Pruning")

library(rattle)
fancyRpartPlot(dtree.pruned)


## Before Pruning ###
dtree.pred <- predict(dtree, validation, type="class")
tbl <- table(dtree.pred, validation$Class, dnn=c("Predicted", "Actual"))
confusionMatrix(tbl)


### After Pruning ###
dtree.pred1 <- predict(dtree.pruned, validation, type="class")
tbl_pruned <- table(dtree.pred1, validation$Class, dnn=c("Predicted", "Actual"))
confusionMatrix(tbl_pruned)


##### Conditional Inference Tree ########

library(party)
fit.ctree <- ctree(Class~., data=training)

plot(fit.ctree, main="Conditional Inference Tree", type="simple")

plot(fit.ctree, main="Conditional Inference Tree", type="extended")

ctree.pred <- predict(fit.ctree, validation, type="response")
ctree.perf <- table( ctree.pred,validation$Class , dnn=c("Predicted", "Actual"))

confusionMatrix(ctree.perf)


