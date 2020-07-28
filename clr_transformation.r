install.packages('compositions')

library('compositions')

#Open the coverage file of the station of interest
genes <- read.csv("bowtie/S_bowtie.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
#Show the top of the dataframe 
head(genes)

#Set the X column as index 
rownames(genes) <- genes$X
genes$X <- NULL

#Perform centered log ratio transformation on the coverage file 
transformed <- clr(genes)

#Replace the coverage by their transformed values 
#S
genes$SO1 <- transformed[, 1]
genes$SO2 <- transformed[, 2]
genes$SR1 <- transformed[, 3]
genes$SR2 <- transformed[, 4]

#P
#genes$PO1 <- transformed[, 1]
#genes$PO2 <- transformed[, 2]
#genes$PO3 <- transformed[, 3]
#genes$PR1 <- transformed[, 4]
#genes$PR2 <- transformed[, 5]
#genes$PR3 <- transformed[, 6]

#CD
#genes$CDO1 <- transformed[, 1]
#genes$CDO2 <- transformed[, 2]
#genes$CDO3 <- transformed[, 3]
#genes$CDI1 <- transformed[, 4]
#genes$CDI2 <- transformed[, 5]
#genes$CDI3 <- transformed[, 6]
#genes$CDR2 <- transformed[, 7]
#genes$CDR3 <- transformed[, 8]

#Show the top of the created dataframe 
head(genes)

#Save the dataframe as csv 
write.csv(genes,"bowtie/S_bowtie_clr.csv")

