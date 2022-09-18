
#Loading necessary libraries. You will need to install these into your R environment using either install.packages or bioconductor. You can google each one.
library(caret)
library(gplots)
library(RColorBrewer)
library(doParallel)
library(sva)
library(GEOquery)
library(limma)
library(affy)
library(MLeval)
library(WebGestaltR)
library(glmnet)
library(plyr)
library(mda)
library(readr)
library(minfi)
library(EpiDISH)


#Small helper functions

#Column vars
colVars <-function(X)
{
  return(rowVars(t(X)))
}

#Pearson correlation with vector v for every row in x returned as sorted vector
getCors <-function(x,v)
{
  res=apply(x,1,function(c) cor(x=c,y=v))
  return(res)
}

#Convert between beta and M values (make sure values are reasonable)
getM <-function(betas2)
{
  betas=betas2
  rmin=0.00001
  rmax=0.99999
  betas[betas<rmin]<-rmin
  betas[betas>rmax]<-rmax
  return(log(betas/(1-betas),2))
}

getMeanGaps <- function(vals)
{
  meangaps=apply(vals,1, function(x){
    sorted=x[order(x)]
    shifted=c(sorted[-1],sorted[length(sorted)])
    return(mean(shifted-sorted))
  })
  return(meangaps)
}

adjustSVA <- function(mat,age)
{
  n.sv = num.sv(mat,age,method="leek")
  svobj = sva(mat,age,model.matrix(~1,data=mat),n.sv=n.sv)
  return(svobj)
}

#Adjust for cell type composition by adjusting the betas. Returns adjusted betas. Takes a bit of time.
#Batch correction happens prior to this. From Jones, Meaghan J. (2015). [Methods in Molecular Biology]  || Adjusting for Cell Type Composition in DNA Methylation Data Using a Regression-Based Approach. , (Chapter 262), –.         doi:10.1007/7651_2015_262     
#adjustCellType <- function(betas,composition)
#{
#  beta.lm<-apply(betas, 1, function(x)
#    {
#    composition[colnames(betas),]->specificComp
##    lm(x~IC+Epi,specificComp)}
#    )
#  cat("Modeled each beta as a linear model of the cell compositions\n")
#  residuals<-t(sapply(beta.lm,function(x)residuals
#                      (summary(x))))
#  cat("Calculated residuals\n")
##  colnames(residuals)<-colnames(betas)
#  adj.betas<-residuals+matrix(apply(betas, 1, mean),
#                              nrow=nrow(residuals), ncol=ncol(residuals))
#  cat("adjusted betas calculated as mean+residuals\n")
#  return(adj.betas)
#}

#Vars
#setwd("Documents/Tally/iscience/revisions/BuccalComparison/") #Set this to whatever your working dir is. 
#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)  #may be required for some big GEO files
## SETUP
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
#stopCluster(cl)
#I'm commenting this out because it is a pain to read in this giant csv file. Instead I recommend just loading in the Rds object using readRDS below
#gse137688<-read_csv("GSE137688_AverageBeta_MavanMoms.csv",col_names = TRUE)
meta137688=read.delim("metaGSE137688.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
#saveRDS(gse137688,"gse137688.Rds")

gse137688<-readRDS("gse137688.Rds")

#Convert to data.frame, which is easier to work with (may not be necessary)
df137688=as.data.frame(gse137688)

#You may or may not need to set this. If you are having trouble loading in GEO datasets you may need to.
#readr::local_edition(1)

#Load in the smoking dataset
gse<-getGEO("GSE94876",GSEMatrix = TRUE)
ages <- as.numeric(sub("age: ","",gse$GSE94876_series_matrix.txt.gz$characteristics_ch1.2))
race <- sub("race: ","",gse$GSE94876_series_matrix.txt.gz$characteristics_ch1.1)
smoker<- sub("class: ","",gse$GSE94876_series_matrix.txt.gz$characteristics_ch1)
df94876=as.data.frame(gse$GSE94876_series_matrix.txt.gz)

#Load in the saliva GSE111223
GSE111223=getGEO("GSE111223")# doesn't work to get data only meta
ages2 <- as.numeric(sub("age: ","",GSE111223$GSE111223_series_matrix.txt.gz$characteristics_ch1.1))
PD <- rep("Healthy",length(ages2))
PD[grep("Parkinson",GSE111223$GSE111223_series_matrix.txt.gz$`disease state:ch1`)]<-"PD"
#instead we will process from idats.
targets=data.frame("Basename"=substr(list.files(path="GSE111223/",pattern="Red.idat"),start=1,stop=30))
#library(R.utils)
#for(f in list.files(path="GSE111223/",pattern=".gz")){
#  gunzip(paste0("GSE111223/",f))  
#}
RGset=read.metharray.exp(base = "GSE111223/idat",targets=targets)
noob=preprocessNoob(RGset)
gsnoob=mapToGenome(noob)
#Annotate the sites
anno=getAnnotation(gsnoob)
#Get betas
df111223=t(getBeta(noob))

#Load in GSE40279 (blood)
#=getGEO("GSE40279")
#df40279=as.data.frame(GSE40279$GSE40279_series_matrix.txt.gz)
#agesblood=as.numeric(GSE40279$GSE40279_series_matrix.txt.gz$`age (y):ch1`)
#bloodstatus=rep("Healthy",length(agesblood))

#Load in GSE50586
GSE50586<-getGEO("GSE50586")
df50586=as.data.frame(GSE50586$GSE50586_series_matrix.txt.gz)
#colnames(df50586)<-substr(colnames(df50586),31,999)
#colnames(df94876)<-substr(colnames(df94876),31,999)
df50586=df50586[,order(colnames(df50586))]
df94876=df94876[,order(colnames(df94876))]
df137688=df137688[,order(colnames(df137688))]
df111223=df111223[,order(colnames(df111223))]

#Find the common probes between all 4 datasets loaded (they should overlap almost completely)
commonProbes=intersect(colnames(df50586),intersect(colnames(df94876),intersect(colnames(df137688),colnames(df111223))))
#Sanity checks
a=df50586[,colnames(df50586)%in%commonProbes]
b=df94876[,colnames(df94876)%in%commonProbes]
c=df137688[,colnames(df137688)%in%commonProbes]
d=df111223[,colnames(df111223)%in%commonProbes]
#e=df40279[,colnames(df40279)%in%commonProbes]
all(colnames(a)==colnames(b))
all(colnames(a)==colnames(c))
all(colnames(a)==colnames(d))
all(colnames(a)==colnames(e))
dim(a)
dim(b)
dim(c)
dim(d)
#dim(e)

#Combine by common probes:
df_comb=rbind(df50586[,colnames(df50586)%in%commonProbes],df94876[,colnames(df94876)%in%commonProbes],df137688[,colnames(df137688)%in%commonProbes],df111223[,colnames(df111223)%in%commonProbes])
#Combine the ages and status effects as well to correspond to each of the samples
ages_comb=c(as.numeric(GSE50586$GSE50586_series_matrix.txt.gz$`age:ch1`),as.numeric(ages),as.numeric(meta137688$Age),as.numeric(ages2))
status_comb=c(rep("DS",10),rep("Healthy",10),smoker,rep("Healthy",length(meta137688$Age)),PD)
status_comb[status_comb=="Non-Tobacco User"]<-"Healthy"

c2=getCors(ms[,status_comb=="Healthy"],ages_comb[status_comb=="Healthy"])
#hist(c2,breaks=30)
#Get most variable sites
v=rowVars(ms[,status_comb=="Healthy"])
names(v)<-rownames(ms) 
#hist(v,breaks=30)

#Get sites that don't have gaps in their values (they are relatively evenly distributed)
g=getMeanGaps(ms[,status_comb=="Healthy"])

#filter corelated CpGs that are changing and are well distributed
c=c2[v>0.1&g<0.02]

#now let's sort by correlation magnitude
c=c[order(abs(c),decreasing = TRUE)]


#Identify cell type
df_comb[is.na(df_comb)]<-0
df_mat=t(as.matrix(df_comb))
betas=df_mat[rowVars(df_mat)>0.001,]
#One unbiased method for identifying confounding factors:
#out = Tsisal(betas, K =NULL, knowRef = NULL, possibleCellNumber = 2:7)
#colnames(out$estProp)<-c("A","B","C")
#A biased method specific for identifying 3 cell types:
#outEpi <-epidish(beta.m=df_mat,ref.m = centEpiFibIC.m,method = "RPC")$estF
out2 <- hepidish(beta.m = df_mat, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')

out2[out2<0]<-0
out2<-cbind(out2,rowSums(out2[,c(2,3,4,5,6,7,9)]))
celltypecol=rgb(out2[,c(1,8,10)])



#Setting batches for batch correction with COMBAT
#batches=c(rep(1,20),rep(2,120),rep(3,250),rep(4,length(ages2)))

#run batch correction
#combated=ComBat(t(df_comb),batch = batches,mod = ages_comb)
#msc=getM(combated)
#msc=msc[!is.na(rowSums(msc)),]
#run cell type correction
#corrected=adjustCellType(combated,as.data.frame(out2))

#now let's just check to see how variable the betas are between health conditions after all of the corrections. Can we find any diff CpGs?
#dmpSmokersVsHealthy=dmpFinder(corrected[,-c(1:10)],pheno=status_comb[-c(1:10)],type="categorical")
#dmpHealthyVsOther=dmpFinder(corrected,pheno=status_comb,type="categorical")

#convert to M values for ML
#ms=getM(corrected)
ms=getM(t(df_comb))
ms=ms[!is.na(rowSums(ms)),]

set.seed(42)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000000000, 19) 
seeds[[11]]<-sample(1000000000,1)
####Train control (settings for ML, check out ?trainControl)
control<- trainControl(method="repeatedcv",number = 10, repeats=1, verboseIter = TRUE,selectionFunction = "best",predictionBounds = c(0,100), allowParallel = TRUE,savePredictions = "final")
control100<- trainControl(method="repeatedcv",number = 10, repeats=100, verboseIter = TRUE,selectionFunction = "best",predictionBounds = c(0,100), allowParallel = TRUE,savePredictions = "final")

eGrid <- expand.grid(.alpha = (1:19) * 0.05, 
                     .lambda = (1:19) * 0.05)

common=read.table("CommonCpGs130.txt")
dfmat130=df_mat[rownames(df_mat)%in% common[,1],]
#dfmat130=combated[rownames(combated)%in% common[,1],]
ms130=getM(dfmat130)
dfall=as.data.frame(t(ms130))
dfall$Epi=out2[,1]
dfall$Fib=out2[,2]
dfall$B=out2[,3]
dfall$NK=out2[,4]
dfall$CD4T=out2[,5]
dfall$CD8T=out2[,6]
dfall$Mono=out2[,7]
dfall$Neutro=out2[,8]
#ms130=ms[row.names(ms) %in% common[,1],]
all_r2_130_200 = vector()
all_mae_130_200 = vector()
all_rmse_130_200 = vector()
all_p_ds_130_200 = vector()
all_p_ms_130_200 = vector()
for(i in 1:100){
df130=as.data.frame(t(ms130[,status_comb=="Healthy"]))
df130$Age=as.numeric(ages_comb[status_comb=="Healthy"])
df130$Epi=out2[status_comb=="Healthy",1]
df130$Fib=out2[status_comb=="Healthy",2]
df130$B=out2[status_comb=="Healthy",3]
df130$NK=out2[status_comb=="Healthy",4]
df130$CD4T=out2[status_comb=="Healthy",5]
df130$CD8T=out2[status_comb=="Healthy",6]
df130$Mono=out2[status_comb=="Healthy",7]
df130$Neutro=out2[status_comb=="Healthy",8]
#impute missing values
impute=preProcess(df130, method = c("medianImpute"))
df130=predict(impute,df130)
ss<-seq(1,431,1)
#random 300
s200<-sample(ss,200)
#same random 300 
#df130=df130[s300,]
#same random 200
df130=df130[s200,]
#df130=df130[s100,]
#set.seed(13)
model130 <- train(Age~.,data=df130,method="lm",trControl=control, preProcess=c("center","scale"))
all_r2_130_200 <- append(all_r2_130_200, model130$results$Rsquared)
all_mae_130_200 <- append(all_mae_130_200, model130$results$MAE)
all_rmse_130_200 <- append(all_rmse_130_200, model130$results$RMSE)
pred=predict(model130,dfall)
DAHealthy=model130$pred$pred-model130$pred$obs
DADS=pred[status_comb=="DS"]-ages_comb[status_comb=="DS"]
DASNUFF=pred[status_comb=="Moist Snuff User"]-ages_comb[status_comb=="Moist Snuff User"]    
all_p_ms_130_200 <- append(all_p_ms_130_200,wilcox.test(DAHealthy,DASNUFF)$p.value)
all_p_ds_130_200 <- append(all_p_ds_130_200, wilcox.test(DAHealthy,DADS)$p.value)
}
all_r2_130_300 = vector()
all_mae_130_300 = vector()
all_rmse_130_300 = vector()
all_p_ds_130_300 = vector()
all_p_ms_130_300 = vector()
for(i in 1:100){
  df130=as.data.frame(t(ms130[,status_comb=="Healthy"]))
  df130$Age=as.numeric(ages_comb[status_comb=="Healthy"])
  df130$Epi=out2[status_comb=="Healthy",1]
  df130$Fib=out2[status_comb=="Healthy",2]
  df130$B=out2[status_comb=="Healthy",3]
  df130$NK=out2[status_comb=="Healthy",4]
  df130$CD4T=out2[status_comb=="Healthy",5]
  df130$CD8T=out2[status_comb=="Healthy",6]
  df130$Mono=out2[status_comb=="Healthy",7]
  df130$Neutro=out2[status_comb=="Healthy",8]
  #impute missing values
  impute=preProcess(df130, method = c("medianImpute"))
  df130=predict(impute,df130)
  ss<-seq(1,431,1)
  #random 300
  s300<-sample(ss,300)
  #same random 300 
  #df130=df130[s300,]
  #same random 200
  df130=df130[s300,]
  #df130=df130[s100,]
  #set.seed(13)
  model130 <- train(Age~.,data=df130,method="lm",trControl=control, preProcess=c("center","scale"))
  all_r2_130_300 <- append(all_r2_130_300, model130$results$Rsquared)
  all_mae_130_300 <- append(all_mae_130_300, model130$results$MAE)
  all_rmse_130_300 <- append(all_rmse_130_300, model130$results$RMSE)
  pred=predict(model130,dfall)
  DAHealthy=model130$pred$pred-model130$pred$obs
  DADS=pred[status_comb=="DS"]-ages_comb[status_comb=="DS"]
  DASNUFF=pred[status_comb=="Moist Snuff User"]-ages_comb[status_comb=="Moist Snuff User"]    
  all_p_ms_130_300 <- append(all_p_ms_130_300,wilcox.test(DAHealthy,DASNUFF)$p.value)
  all_p_ds_130_300 <- append(all_p_ds_130_300, wilcox.test(DAHealthy,DADS)$p.value)
}



dfall=as.data.frame(t(ms))
dfall$Epi=out2[,1]
dfall$Fib=out2[,2]
dfall$B=out2[,3]
dfall$NK=out2[,4]
dfall$CD4T=out2[,5]
dfall$CD8T=out2[,6]
dfall$Mono=out2[,7]
dfall$Neutro=out2[,8]


all_mae_200 = vector()
all_r2_200 = vector()
all_rmse_200 = vector()
all_p_ds_200 = vector()
all_p_ms_200 = vector()
for(i in 1:100){
  df=as.data.frame(t(ms[rownames(ms)%in% names(c)[1:10000],status_comb=="Healthy"]))
  df$Age=as.numeric(ages_comb[status_comb=="Healthy"])
  df$Epi=out2[status_comb=="Healthy",1]
  df$Fib=out2[status_comb=="Healthy",2]
  df$B=out2[status_comb=="Healthy",3]
  df$NK=out2[status_comb=="Healthy",4]
  df$CD4T=out2[status_comb=="Healthy",5]
  df$CD8T=out2[status_comb=="Healthy",6]
  df$Mono=out2[status_comb=="Healthy",7]
  df$Neutro=out2[status_comb=="Healthy",8]
  ss<-seq(1,431,1)
  s200<-sample(ss,200)
  df=df[s200,]
  model <- train(Age~.,data=df,method="glmnet",trControl=control, importance=TRUE,tuneGrid=eGrid, preProcess=c("center","scale"))
  all_r2_200<-append(all_r2_200,model$results$Rsquared[bestIndex])
  all_rmse_200<-append(all_rmse_200,model$results$RMSE[bestIndex])
  all_mae_200<-append(all_mae_200,model$results$MAE[bestIndex])
  pred=predict(model,dfall)
  DAHealthy=model$pred$pred-model$pred$obs
  DADS=pred[status_comb=="DS"]-ages_comb[status_comb=="DS"]
  DASNUFF=pred[status_comb=="Moist Snuff User"]-ages_comb[status_comb=="Moist Snuff User"]
  all_p_ds_200<-append(all_p_ds_200,wilcox.test(DAHealthy,DADS)$p.value)
  all_p_ms_200<-append(all_p_ms_200,wilcox.test(DAHealthy,DASNUFF)$p.value)
  sink("unbiased_200.txt")
  print(all_mae_200)
  print(all_r2_200)
  print(all_rmse_200)
  print(all_mae_200)
  print(all_p_ds_200)
  print(all_p_ms_200)
  sink()
  print(i)
}

all_mae_300 = vector()
all_r2_300 = vector()
all_rmse_300 = vector()
all_p_ds_300 = vector()
all_p_ms_300 = vector()
for(i in 1:100){
  df=as.data.frame(t(ms[rownames(ms)%in% names(c)[1:10000],status_comb=="Healthy"]))
  df$Age=as.numeric(ages_comb[status_comb=="Healthy"])
  df$Epi=out2[status_comb=="Healthy",1]
  df$Fib=out2[status_comb=="Healthy",2]
  df$B=out2[status_comb=="Healthy",3]
  df$NK=out2[status_comb=="Healthy",4]
  df$CD4T=out2[status_comb=="Healthy",5]
  df$CD8T=out2[status_comb=="Healthy",6]
  df$Mono=out2[status_comb=="Healthy",7]
  df$Neutro=out2[status_comb=="Healthy",8]
  ss<-seq(1,431,1)
  s300<-sample(ss,300)
  df=df[s300,]
  model <- train(Age~.,data=df,method="glmnet",trControl=control, importance=TRUE,tuneGrid=eGrid, preProcess=c("center","scale"))
  all_r2_300<-append(all_r2_300,model$results$Rsquared[bestIndex])
  all_rmse_300<-append(all_rmse_300,model$results$RMSE[bestIndex])
  all_mae_300<-append(all_mae_300,model$results$MAE[bestIndex])
  pred=predict(model,dfall)
  DAHealthy=model$pred$pred-model$pred$obs
  DADS=pred[status_comb=="DS"]-ages_comb[status_comb=="DS"]
  DASNUFF=pred[status_comb=="Moist Snuff User"]-ages_comb[status_comb=="Moist Snuff User"]
  all_p_ds_300<-append(all_p_ds_300,wilcox.test(DAHealthy,DADS)$p.value)
  all_p_ms_300<-append(all_p_ms_300,wilcox.test(DAHealthy,DASNUFF)$p.value)
  sink("unbiased_300_2.txt")
  print(all_mae_300)
  print(all_r2_300)
  print(all_rmse_300)
  print(all_mae_300)
  print(all_p_ds_300)
  print(all_p_ms_300)
  sink()
  print(i)
}

#all_p_ms_200<-scan('all_p_ms_200.txt', what=numeric(), sep=" ", quiet=TRUE)
#all_p_ds_200<-scan('all_p_ds_200.txt', what=numeric(), sep=" ", quiet=TRUE)
#all_r_200<-scan('all_r_200.txt', what=numeric(), sep=" ", quiet=TRUE)
#all_mae_200<-scan('all_mae_200.txt', what=numeric(), sep=" ", quiet=TRUE)

all_r_200 = vector()
all_r_300 = vector()
for (i in all_r2_200){all_r_200<-append(all_r_200,sqrt(i))}
for (i in all_r2_300){all_r_300<-append(all_r_300,sqrt(i))}
all_p_ms_200_log= vector()
all_p_ds_200_log= vector()
all_p_ms_300_log= vector()
all_p_ds_300_log= vector()
for (i in all_p_ms_200){all_p_ms_200_log<-append(all_p_ms_200_log,-log(i))}
for (i in all_p_ms_300){all_p_ms_300_log<-append(all_p_ms_300_log,-log(i))}
for (i in all_p_ds_200){all_p_ds_200_log<-append(all_p_ds_200_log,-log(i))}
for (i in all_p_ds_300){all_p_ds_300_log<-append(all_p_ds_300_log,-log(i))}

DS_plog_boxes=list("DS 200"=all_p_ds_200_log,"DS 300"=all_p_ds_300_log)
DS_p_boxes=list("DS 200"=all_p_ds_200,"DS 300"=all_p_ds_300)

pdf("DA_DS_pvals_random_sampling_log")
boxplot(DS_plog_boxes,main="Delta Age DS p-values",xlab = "Condition and number of samples", ylab = "-log(p-values)", col=c("Blue","orange"))
dev.off()

pdf("DA_DS_pvals_random_sampling")
boxplot(DS_p_boxes,main="Delta Age DS p-values",xlab = "Condition and number of samples", ylab = "p-values", col=c("Blue","orange"))
dev.off()

MS_plog_boxes=list("MS 200"=all_p_ms_200_log,"MS 300"=all_p_ms_300_log)
MS_p_boxes=list("MS 200"=all_p_ms_200,"MS 300"=all_p_ms_300)

pdf("DA_MS_pvals_random_sampling_log")
boxplot(MS_plog_boxes,main="Delta Age MS p-values",xlab = "Condition and number of samples", ylab = "-log(p-values)", col=c("Blue","orange"))
dev.off()

pdf("DA_MS_pvals_random_sampling")
boxplot(MS_p_boxes,main="Delta Age MS p-values",xlab = "Condition and number of samples", ylab = "p-values", col=c("Blue","orange"))
dev.off()


r_boxes = list("200"=all_r_200,"300"=all_r_300)
pdf("10K_R_random_sampling")
boxplot(r_boxes,main="10K Rvalues",xlab = "number of samples", ylab = "R values", col=c("Blue","orange"))
dev.off()


mae_boxes = list("200"=all_mae_200,"300"=all_mae_300)
pdf("10K_MAE_random_sampling")
boxplot(mae_boxes,main="10K MAE values",xlab = "number of samples", ylab = "MAE values", col=c("Blue","orange"))
dev.off()

all_r_130_200 = vector()
all_r_130_300 = vector()
for (i in all_r2_130_200){all_r_130_200<-append(all_r_130_200,sqrt(i))}
for (i in all_r2_130_300){all_r_130_300<-append(all_r_130_300,sqrt(i))}
all_p_ms_130_200_log= vector()
all_p_ds_130_200_log= vector()
all_p_ms_130_300_log= vector()
all_p_ds_130_300_log= vector()
for (i in all_p_ms_130_200){all_p_ms_130_200_log<-append(all_p_ms_130_200_log,-log(i))}
for (i in all_p_ms_130_300){all_p_ms_130_300_log<-append(all_p_ms_130_300_log,-log(i))}
for (i in all_p_ds_130_200){all_p_ds_130_200_log<-append(all_p_ds_130_200_log,-log(i))}
for (i in all_p_ds_130_300){all_p_ds_130_300_log<-append(all_p_ds_130_300_log,-log(i))}

DS_plog_130_boxes=list("DS 200"=all_p_ds_130_200_log,"DS 300"=all_p_ds_130_300_log)
DS_p_130_boxes=list("DS 200"=all_p_ds_130_200,"DS 300"=all_p_ds_130_300)

pdf("DA_DS_pvals_random_sampling_log_130")
boxplot(DS_plog_130_boxes,main="Delta Age DS p-values",xlab = "Condition and number of samples", ylab = "-log(p-values)", col=c("Blue","orange"))
dev.off()

pdf("DA_DS_pvals_random_sampling_130")
boxplot(DS_p_130_boxes,main="Delta Age DS p-values",xlab = "Condition and number of samples", ylab = "p-values", col=c("Blue","orange"))
dev.off()

MS_plog_130_boxes=list("MS 200"=all_p_ms_130_200_log,"MS 300"=all_p_ms_130_300_log)
MS_p_130_boxes=list("MS 200"=all_p_ms_130_200,"MS 300"=all_p_ms_130_300)

pdf("DA_MS_pvals_random_sampling_log_130")
boxplot(MS_plog_130_boxes,main="Delta Age MS p-values",xlab = "Condition and number of samples", ylab = "-log(p-values)", col=c("Blue","orange"))
dev.off()

pdf("DA_MS_pvals_random_sampling_130")
boxplot(MS_p_130_boxes,main="Delta Age MS p-values",xlab = "Condition and number of samples", ylab = "p-values", col=c("Blue","orange"))
dev.off()


r_130_boxes = list("200"=all_r_130_200,"300"=all_r_130_300)
pdf("130_R_random_sampling_130")
boxplot(r_130_boxes,main="Common130 Rvalues",xlab = "number of samples", ylab = "R values", col=c("Blue","orange"))
dev.off()


mae_boxes = list("200"=all_mae_130_200,"300"=all_mae_130_300)
pdf("130_MAE_random_sampling_130")
boxplot(mae_boxes,main="Common130 MAE values",xlab = "number of samples", ylab = "MAE values", col=c("Blue","orange"))
dev.off()


save.image("df_comb_na.RData")
