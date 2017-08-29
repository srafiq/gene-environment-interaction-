library(gdata)
library(car)


ped <-read.table(file="biobank_chr1_part1.ped")


genos<-ped[, -c(1:5)]

genos<-data.frame( genos[1], mapply( paste0, genos[-1][c(T,F)], genos[-1][c(F,T)])) 

genos<-genos[, -c(1)]


genos_modified<-lapply(genos, FUN=function(x) recode(x, "22=0;12=1;11=2"))
replace_00<-function(x){
x[x == "00"] <- NA
x
}


genos_modified[] <- lapply(genos_modified, replace_00)


write.table(genos_modified, file="genos_modified")


genos_mod<-as.matrix(read.table("genos_modified", header=TRUE, sep=" "))

genos_mod<-as.matrix(genos_mod)


snps<-read.table(file="biobank_chr22_part1.map")


snps_ids<-snps$V2

colnames(snps_ids)<-NULL

snp_id<-gsub("\\,.*", "", snps_ids)


colnames(genos_mod)<-snp_id



#####Lets reattach the Ids to the genotypes


ids<-ped[, c(1)]


genos_ready<-cbind(ids, genos_mod)


###merge with phenotype data

pheno<-read.table(file="pheno_final.txt", header=TRUE)



###merge pheno and genotype 


file1<-merge(pheno, genos_ready, by="ids")


write.csv(file1,file="file1.csv")



##########gxe scan 


library(lmtest)

data<-read.csv(file="file1.csv", header=T, sep=",")

data[is.na(data)]=-1



p_values=c()
for( i in 5:ncol(data))
{
  model1<-glm(data$mi~data[,colnames(data)[i]]+data$dash, data=data,na.action = na.exclude )
  temp=data
  data["IxDASH"]=data$dash*data[,colnames(data)[i]]
  model2=glm(data$mi~ data[,colnames(data)[i]]+data$dash+data$IxDASH,na.action = na.exclude)
  data=temp
  print(lrtest(model2,model1))
  model=lrtest(model2,model1)
  p_values=c(p_values,model$`Pr(>Chisq)`[2])
  
  
  
}
treshold=0.05
idx=which(p_values<treshold)
#idx  are    the ,,i''s  where    p_values  are smaller  then  treshold
print(idx)
file=as.data.frame(cbind(colnames(data[5:ncol(data)])[idx],as.numeric(p_values[idx])))
 
colnames(file)=c("Predictor","Corresponding p-value")
write.csv(file,'chr1_1.csv',row.names = F)

q()
n



  
