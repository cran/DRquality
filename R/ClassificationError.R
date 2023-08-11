ClassificationError=function(OutputDistances,Cls,k=5){
#acc=ClassificationError(OutputDistances,Cls)$Error
#Methods are evaluated by k-nearest neighbor classification accuracy (with k = 5), that is, each sample in the
#visualization is classified by majority vote of its k nearest neighbors in the visualization, and the
#classification is compared to the ground truth label. [Venna 2010]
#INPUT
# OutputDistances[1:n,1:n]
# Cls[1:n]
# Optional
# k                   number of k nearest neighbors, in Venna 2010 set to 5 (here default)
# OUTPUT
# Error               1-Accuracy[1]
# Accuracy[1]
# NewCls
#author: MT 07/2016
n=nrow(OutputDistances)
NewCls=Cls*NaN
#requireRpackage('FastKNN')
requireNamespace('FastKNN')
for(j in 1:n){
  votes=FastKNN::k.nearest.neighbors(j,OutputDistances,k)
  uniqueClasses <- sort(na.last=T,unique(Cls[votes]))
  numberOfClasses <- length(uniqueClasses)
  countPerClass <- rep(0, numberOfClasses) #
  for(i in 1: numberOfClasses) {
    inClassI <- sum(Cls[votes] == uniqueClasses[i]) # counts all occurances of uniqueClass[i] in cls
    countPerClass[i] = inClassI # updates countPerClass[i] to the number of occurances of uniqueClasses[i] in cls.
  }
  ind=which.max(countPerClass)
  
  NewCls[j]=uniqueClasses[ind]
}

accuracy=FCPS::ClusterAccuracy(Cls, NewCls)

return(list(Error=1-accuracy,Accuracy=accuracy,KNNCls=NewCls))
}