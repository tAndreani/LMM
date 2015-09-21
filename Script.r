#Load data
line_specific <- read.table("merged.isoforms.fpkm.tracking.final.no.30.and.70.correct.txt",header=T,sep="\t")
reference <- read.table("CufflinksQuantification_parsed_column_renamed_ordered.csv ",header=T,sep="\t")

#Extract common transcript
transcript_reference <- reference$Transcript_ID
trancript_line_specific <- line_specific$Transcript_ID
reference_TrId_in_LineSpecific <- reference[which(reference$Transcript_ID%in%trancript_line_specific),]
line_specific_TrId_in_Reference <- line_specific[which(line_specific$Transcript_ID%in%transcript_reference),]

#Sort
library("plyr")
sort_line_specific <- arrange(reference_TrId_in_LineSpecific,Transcript_ID)
sort_reference <- arrange(line_specific_TrId_in_Reference,Transcript_ID)

#Perform Correlation
A <- sort_line_specific[,c(3:80)]
B <- sort_reference[,c(3:80)]
#explore your data
plot(A[,1],B[,1],xlim=c(1,10000),ylim=c(1,10000),main="Methods Comparison r=0.92",xlab="Sample1 Method1",ylab="Sample1 Method 2 ")
#correlate
correlation <- cor(A, B)
correlation_diagonal <- diag(cor(A,B))
boxplot(correlation_diagonal)


#Consider different source of variablity i.e. the samples and the methods
.libPaths("/nfs/users/rg/tandreani/R/x86_64-redhat-linux-gnu-library/3.2")
library("lme4")
library("reshape2")

matrix1=sort_reference
matrix2=sort_line_specific
for (i in 1:nrow(matrix1)){
  m1 <- melt(matrix1[i,2:80])
  m2 <- melt(matrix2[i,2:80])
  m1=cbind(m1[,2:3],rep(1,78))
  m2=cbind(m2[,2:3],rep(2,78))
  names=c("sample","express","method")
  colnames(m1)=names
  colnames(m2)=names
  df=rbind(m1,m2)
  df=df[,c(2,3,1)]
  tx=as.character(matrix1[i,2])
  print(head(df))
  # fit the lmm
  res <- lmer(express~(1|method)+(1|sample),data=df)
  # variance components
  var_random_effect<-as.numeric(VarCorr(res))
  var_residual<-attr(VarCorr(res),"sc")^2
  # save v
  total_variance <- rbind(c(var_random_effect,var_residual),deparse.level=1)
  col_names_lmm=c("method","sample","residual")
  colnames(total_variance)=col_names_lmm
  row_names_lmm=tx
  rownames(total_variance)=row_names_lmm
  if(i==1){
    output=total_variance
  }else{
    output=rbind(output,total_variance)
  }
}
