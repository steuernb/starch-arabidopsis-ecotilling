
file <- read.table("output.txt", sep="\t", header=FALSE)

svg("output.svg")
par(mfrow=c(3,2))

numSNPs1 <-  file$V3/max(c(file$V3, file$V4))
numSNPs2 <-  file$V4/max(c(file$V3, file$V4))

for(chr in levels(file$V1)){
  
  x<- file$V2[file$V1==chr]
  y1<- file$V5[file$V1==chr] + file$V6[file$V1==chr]
  y2<- file$V5[file$V1==chr]
  
  y<-file$V7[file$V1==chr]/50000
  
  plot(x, (y2/y1), type = "h", xlab=chr, xlim = c(0, 31000000), ylim=c(0,1), ylab="freq")
  if(chr=="Chr1"){
    abline(v=11919790, col ="red") 
  }
  
  if(chr=="Chr5"){
    abline(v=15932714, col="orange")
  }
  
 }
dev.off()
