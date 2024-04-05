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
  return(dataf)
}
df<-plot_manhattan(weights)
write.csv(df[,c(2,3)],paste0(path_source,'/results/PostionalPvals.csv'))
df_<-df[df$P<0.05,]
print('Positions: ')
print(df_$position)
