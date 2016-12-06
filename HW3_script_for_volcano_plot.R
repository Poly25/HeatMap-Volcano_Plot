library("GEOquery")
gse <- getGEO ("GSE3541")[[1]]
gse # overview

mat <- exprs(gse) # assay data matrix
pData(gse)["title"] # phano data matix 
colnames(mat) <- c("control_sample1" , "control_sample2" ,"control_sample3","strain_sample1","strain_sample2","strain_sample3")
mat
# heat Map 
heatmap(mat , scale = "row", Rowv = NA , Colv = NA ,labRow = NA , main = "Epithelial cells - Heat Map" )
#par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)  #Sets "nice" graphical defaults

#for t_test -> p-values 
T1 <- log2(mat[,1:3]) # all features(all rows) of control samples - fisrt three columns in mat matrix
T2 <- log2(mat[,4:6]) # all features(all rows) of strain  samples - last three columns in mat matrix

p_values <- c(1:nrow(mat))

for (i in 1:nrow(mat)){
  # p_value of each feature 
  # T1[i,1:3],T2[i,1:3]) distribution of each feature (control and strain sample respectively)
  p_values[i] <- t.test(T1[i,1:3],T2[i,1:3])$p.value 
}

tail(p_values) #six last values of p_value

# mean of each row 
mT1 <- rowMeans(T1)
mT2 <- rowMeans(T2)
FC  <- mT2 - mT1 #Fold Change

head(mT1)
head(mT2)
head(FC)

#volcano plot
plot (FC , -log10(p_values) , pch =19 ,
      xlim = c(-4,4) , ylim = c(0,4) , # limits at x and y axis respectively
      col = ifelse(FC>1.5 & p_values <= 0.05 ,"red",ifelse(FC< -1.5 & p_values<=0.05,"green" ,"black")),
      xlab="log2 fold change", ylab="-log10 p-value") #Set axis labels

     
     