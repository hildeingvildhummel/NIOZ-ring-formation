#install.packages("devtools")
#devtools::install_github("ggloor/ALDEx_bioc")

library('ALDEx2')
library(dplyr)

#Read the csv file
df_1 <- read.csv("bowtie/C_bowtie.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
#Set the first column as index
rownames(df_1) <- df_1$X
df_1$X <- NULL
head(df_1)

#If CD samples are given, select partial dataframe
df <- df_1[, c(4, 5, 6, 7, 8)]
#Otherwise, select the complete dataframe
#df <- df_1


head(df)
#Convert the dataframe to the correct format
dt <- as.matrix(df)
#Select the conditions corresponding to the columns of the dataframe
#conds <- c('O', 'O', 'R', 'R')
#conds <- c('O', 'O', 'O', 'R', 'R', 'R')
conds <- c('O', 'O', 'O', 'R', 'R')

#Perform the ALDEx2 differenital abundance test
res <- aldex(dt, conds, denom = "iqlr")
#Order by the P-value of the Welch's t-test
res <- res[order(res$we.ep),]

#Select the results with a P-value below 0.05
top <- res[res$we.ep < 0.05, ]
top

#Save the dataframe with the significant values to a csv
write.csv(top,"Top_values_C_O_R.csv", row.names = TRUE)

#If CD samples is given continue with the part below

#Select the In and Ring samples
df <- df_1[, c(1,2,3,7,8)]

#Convert dataframe to the correct format 
dt <- as.matrix(df)
#Select the corresponding conditions
conds <- c('I', 'I', 'I', 'R', 'R')

#Perform Aldex2 differential abundance test 
res <- aldex(dt, conds, denom = "iqlr")
#Order results by the P-values of the Welch's t-test
res <- res[order(res$we.ep),]

#Select the values with a P-value below 0.05
top <- res[res$we.ep < 0.05, ]
top

#Save the significant results to a csv
write.csv(top,"Top_values_C_I_R.csv", row.names = TRUE)

#Select the final combination, In and Out samples 
df <- df_1[, c(1,2,3,4,5,6)]

#Convert to correct format
dt <- as.matrix(df)
#Slecht the corresponding conditions
conds <- c('I', 'I', 'I', 'O', 'O', 'O')

#Perform Aldex2 differential abundance test
res <- aldex(dt, conds, denom = "iqlr")
#Order the results by the P-value of the Welch's t-test
res <- res[order(res$we.ep),]

#Select results with a P-value below 0.05
top <- res[res$we.ep < 0.05, ]
top

#Save the significant results to a csv
write.csv(top,"Top_values_C_I_O.csv", row.names = TRUE)
