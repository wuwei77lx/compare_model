# predict peak-gene results and truth peak-gene results:A data frame with n rows and 3 variables.
# column1:{gene};column2:{peak};column3:{score}

#AUPRC
#' @param predict Model for predicting peak-gene
#' @param truth truth peak-gene
library(PRROC)# load the required libraries
AUP<-function(predict,truth){
  data<-merge(truth,predict,by=c("gene","peak"),all.x=T)
  data$score.y[is.na(data$score.y)]=0
  pr_curve_data <- pr.curve(scores.class0 = abs(as.numeric(data$score.y)), weights.class0 = data$score.x, curve = TRUE)
  return(pr_curve_data)
  
}


#Early Precision
#' @param predict Model for predicting peak-gene
#' @param truth truth peak-gene
EPR<-function(predict,truth){
  k<-length(truth$score[truth$score==1])
  data<-merge(truth,predict,by=c("gene","peak"),all.x=T)
  data$score.y<-as.numeric(data$score.y)
  data$score.y[is.na(data$score.y)]=0
  num<-sum(data$score.y!=0.0)
  if(num>=k){
    values<-abs(data$score.y)[order(abs(data$score.y), decreasing = TRUE)[k]]
    data1<-data[abs(data$score.y)>=values,]
    data1$score.y<-1
    early_precision=sum(data1$score.x==data1$score.y)/dim(data1)[1]
  }else{
    data1<-data[abs(data$score.y)>0,]
    data1$score.y<-1
    early_precision=sum(data1$score.x==data1$score.y)/dim(data1)[1]
    print(num)
  }
  print(early_precision)
  return(early_precision)
}