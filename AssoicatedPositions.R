library(dplyr)

path_source <- commandArgs(trailingOnly = TRUE)[1]
path_files<-'results/weights.txt'
path<-paste0(path_source,path_files)
lines <- readLines(path)

# Split each line into separate elements and store them as a list
weights<-as.list(as.numeric(lines))

plot_manhattan<-function(weights){
  betas<-as.numeric(weights)
  dataf<-data.frame(betas)
  p<-length(betas)
  d2 <- mahalanobis( x = as.matrix(betas), center = mean(betas), cov = var(betas) )
  f.stat <- (p - 1) * d2 / ( p-1 )
  neglog10pval <- -pf(f.stat, 1, p-1, lower.tail = FALSE, log.p = TRUE)
  dataf$P= 10^-(neglog10pval)
  dataf$P=p.adjust(dataf$P,method = "BH")
  dataf$position=seq(1,length(neglog10pval),1)
  dataf$BP<-(dataf$position)
  dataf$SNP=dataf$position
  
  
  # Assign 'CHR' values based on 'pos' ranges
  dataf$CHR[dataf$pos >= 1 & dataf$pos <= 2] <- 1
  dataf$CHR[dataf$pos >= 3 & dataf$pos <= 108] <- 2
  dataf$CHR[dataf$pos >= 109 & dataf$pos <= 233] <- 3
  dataf$CHR[dataf$pos >= 234 & dataf$pos <= 326] <- 4
  dataf$CHR[dataf$pos >= 327 & dataf$pos <= 420] <- 5
  dataf$CHR[dataf$pos >= 421 & dataf$pos <= 605] <- 6
  dataf$CHR[dataf$pos >= 606 & dataf$pos <= 711] <- 7
  dataf$CHR[dataf$pos >= 712 & dataf$pos <= 788] <- 8
  dataf$CHR[dataf$pos >= 789 & dataf$pos <= 938] <- 9
  dataf$CHR[dataf$pos >= 939 & dataf$pos <= 1104] <- 10
  dataf$CHR[dataf$pos >= 1105 & dataf$pos <= 1195] <- 11
  dataf$CHR[dataf$pos >= 1196 & dataf$pos <= 1323] <- 12
  return(dataf)
}
df<-plot_manhattan(weights)
df_<-df[df$P<0.05,]
print('Positions: ')
print(df_$position)
print('--------------------------------------------------------------------------')
print('Exons: ')
print(unique(df_$CHR))
