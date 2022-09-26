library(strawr)
library(tidyverse)
library(dbscan)
library(InteractionSet)
library(raster)
#library(devtools)
library(HiCcompare)
library(hictoolsr)
#library(terra)

#wt_loops <- read.file("WT_5kbLoops.txt")
#fs_loops <- read.file("FS_5kbLoops.txt")
setwd("/Users/phanstiel2/MK/Work/")

## Merge loops and convert to GInteractions

loops <- mergeBedpe(bedpeFiles = c("WT_5kbLoops.txt", "FS_5kbLoops.txt"), res = 10e3) |> as_ginteractions()
write.csv(loops,"loops")
setwd("/Users/phanstiel2/MK/Work")
head(loops)

#hicFiles = commandArgs(trailingOnly=TRUE)
hicFiles[1]
hicFiles <- c("GSM4259900_HEK_HiC_NUP_IDR_FS_A9_1_1_inter_30.hic")

mh_index <- function(buffer, loop, inner){
  m=(buffer*2)+1
  center <- buffer+1
  M <- matrix(data=NA,nrow=m,ncol=m)
  for(j in 1:m){
    l=buffer+j
    for(i in 1:m){
      k=(m+1)-j
      if((i <= (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-k-i
      }
      if((i <= (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-j-i
      }
      if((i > (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-i-j
      }
      if((i > (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-i-k
      }
    }
  }
  inner_val <- which(M == inner)
  new <- loop[inner_val]
  return(new)
}

head(loops)

loops_filtered <- filterBedpe(bedpe = loops, res = res, buffer = buffer)


normalized=data.frame()
for(j in 1:length(loops_filtered)){
  l <- calcApa(bedpe = loops_filtered [j], hic = hicFiles[1], norm = "NONE", res = 10000, buffer = 10, 
               filter = FALSE) 
  for(i in 0:10){
    normalized[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}

head(normalized)
colnames(normalized) <- c("0","1","2","3","4","5","6","7","8","9","10")
write.csv(normalized,"MH_obs_all.csv",quote=FALSE)

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero=data.frame()
for(i in 1:20){
  for(j in 1:11){
    a <- normalized[i,1]-normalized[i,11]
    b <- normalized[i,j]-normalized[i,11]
    zero[i,j] <- round((b/a),2)
  }
}
head(zero)
colnames(zero) <- c("0","1","2","3","4","5","6","7","8","9","10")


Group <- factor(0:10)
normalized1 <- na.omit(zero)
Mean<- as.vector(colMeans(normalized1))

###Plotting ########

df <- data.frame(Group,Mean)
df
library(ggplot2)
p <- ggplot(df,aes(Group,Mean,group=1))+geom_smooth(color=c("black"),se=F)+xlab("Distance from the center") + ylab("Signal relative to the center")+theme_classic() +scale_y_continuous(limits = c(0, 1))
p + theme(
  axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=16),
  axis.title.x = element_text(color="#993333", size=24, face="bold"),
  axis.title.y = element_text(color="#993333", size=24, face="bold")
)