#install.packages('vegan')
#install.packages(c("cluster.datasets"), dependencies = TRUE)
#install.packages("FactoMineR")
#install.packages("factoextra")
#install.packages("labdsv")
#install.packages("hilldiv")

library('vegan')
library(cluster.datasets)
library(readr)
library(ggplot2)
library(gplots)
library("FactoMineR")
library("factoextra")
library(dplyr)
library(labdsv)
library(hilldiv)


#Read the 16s/18s csv file of the Archaea
df_A <- read.csv("16s_18s/files/OTU_table_Archaea.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
rownames(df_A) <- df_A$X
df_A$X <- NULL
colnames(df_A) <- c('CDI1', 'CDI2', 'CDI3', 'CDO1', 'CDO2', 'CDO3', 'CDR2', 'CDR3','PR1', 'PR2', 'PR3', 'PO1', 'PO2', 'PO3', 'SR1', 'SR2', 'SO1', 'SO2')
head(df_A)

#Read the 16s/18s csv file of the Bacteria
df_B <- read.csv("16s_18s/files/OTU_table_Bacteria.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
rownames(df_B) <- df_B$X
df_B$X <- NULL
colnames(df_B) <- c('CDI1', 'CDI2', 'CDI3', 'CDO1', 'CDO2', 'CDO3', 'CDR2', 'CDR3','PR1', 'PR2', 'PR3', 'PO1', 'PO2', 'PO3', 'SR1', 'SR2', 'SO1', 'SO2')
head(df_B)

#Read the 16s/18s csv file of the Eukaryota
df_E <- read.csv("16s_18s/files/OTU_table_Eukaryota.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
rownames(df_E) <- df_E$X
df_E$X <- NULL
colnames(df_E) <- c('CDI1', 'CDI2', 'CDI3', 'CDO1', 'CDO2', 'CDO3', 'CDR2', 'CDR3','PR1', 'PR2', 'PR3', 'PO1', 'PO2', 'PO3', 'SR1', 'SR2', 'SO1', 'SO2')
head(df_E)

df <- rbind(df_A, df_B, df_E)


#Normalize counts by TSS
#df_tss <- tss(df_E)
#df_tss <- tss(df_A)
#df_tss <- tss(df_B)
df_tss <- tss(df)


#Convert the dataframe to the right format
dt <- as.table(as.matrix(df_tss))

#CA plot for CD samples
res.ca <- CA(dt, ncp = 2, graph = TRUE)
pl <- ordiellipse(dt, c('C','C', 'C', 'C', 'C', 'C', 'C', 'C', 'P', 'P', 'P', 'P', 'P', 'P', 'S', 'S', 'S', 'S'), kind="se", conf=0.95, lwd=2, draw = "polygon", 
                  col="skyblue", border = "blue")
#CA information plots
get_eigenvalue(res.ca)
fviz_eig(res.ca)
fviz_ca_row(res.ca)
fviz_ca_col(res.ca)
fviz_ca_biplot(res.ca)

#CA plot for CD samples
res.ca <- CA(dt[,grepl("C", colnames(dt))], ncp = 2, graph = TRUE)
#CA information plots
get_eigenvalue(res.ca)
fviz_eig(res.ca)
fviz_ca_row(res.ca)
fviz_ca_col(res.ca)
fviz_ca_biplot(res.ca)

#CA plot for P samples
res.ca <- CA(dt[,grepl("P", colnames(dt))], ncp = 2, graph = TRUE)
#CA information plots for P samples
get_eigenvalue(res.ca)
fviz_eig(res.ca)
fviz_ca_row(res.ca)
fviz_ca_col(res.ca)
fviz_ca_biplot(res.ca)

#CA plot for S samples
res.ca <- CA(dt[,grepl("S", colnames(dt))], ncp = 2, graph = TRUE)
#CA information plots for S samples
get_eigenvalue(res.ca)
fviz_eig(res.ca)
fviz_ca_row(res.ca)
fviz_ca_col(res.ca)
fviz_ca_biplot(res.ca)

#Environmental info
env <- data.frame(site=c('Charles', 'Charles', 'Charles', 'Charles', 'Charles', 'Charles', 'Charles', 'Charles', 'Pascal', 'Pascal', 'Pascal', 'Pascal', 'Pascal','Pascal', 'Silvain', 'Silvain', 'Silvain', 'Silvain'), place=c('I', 'I', 'I', 'O', 'O', 'O', 'R', 'R', 'O', 'O', 'O', 'R', 'R', 'R', 'O', 'O', 'R', 'R'), salt=c(18, 18, 18, 18, 18, 18, 18, 18, 13, 13, 13, 13, 13, 13, 13, 13, 13,13))
env_table <- as.table(as.matrix(env))
rownames(env) <- sample_order <- c('CDI1', 'CDI2', 'CDI3', 'CDO1', 'CDO2', 'CDO3', 'CDR2', 'CDR3','PR1', 'PR2', 'PR3', 'PO1', 'PO2', 'PO3', 'SR1', 'SR2', 'SO1', 'SO2')

#Check influence of environmental arguments
#CD samples
envfit(t(dt[, 1:8]), env[1:8, ])

#P samples
envfit(t(dt[, 9:14]), env[9:14, ])
#S samples
envfit(t(dt[, 15:18]), env[15:18, ])


#Shannon
H <- diversity(t(df_tss), index = 'shannon')
H
#Richness
S <- specnumber(t(df_tss))
#Evenness
J <- H/log(S)
J

#Richness estimation
r <- estimateR(t(df_E), plot = TRUE)

sample_order <- c('CDI1', 'CDI2', 'CDI3', 'CDR2', 'CDR3', 'CDO1', 'CDO2', 'CDO3', 'PR1', 'PR2', 'PR3', 'PO1', 'PO2', 'PO3', 'SR1', 'SR2', 'SO1', 'SO2')

#Create a dataframe only containing the Chao1 index
chao_df <- data.frame(colnames(r), r[2,])
colnames(chao_df) <- c('sample', 'chao1')
plotting_chao <- chao_df[sample_order, ]


# lock in factor level order
plotting_chao$sample <- factor(plotting_chao$sample, levels = plotting_chao$sample)

#Plot the choa
p <-ggplot(data=plotting_chao[1:8,], aes(x=sample, y=chao1)) +
  geom_bar(stat="identity")
p

q <-ggplot(data=plotting_chao[9:18,], aes(x=sample, y=chao1)) +
  geom_bar(stat="identity")
q

k <-ggplot(data=plotting_chao, aes(x=sample, y=chao1)) +
  geom_bar(stat="identity")
k

#Plot ACE
ace_df <- data.frame(colnames(r), r[4,])
colnames(ace_df) <- c('sample', 'ACE')
plotting_ace = ace_df[sample_order, ]

# lock in factor level order
plotting_ace$sample <- factor(plotting_ace$sample, levels = plotting_ace$sample)

#Pllot the choa
p <-ggplot(data=plotting_ace[1:8,], aes(x=sample, y=ACE)) +
  geom_bar(stat="identity")
p

q <-ggplot(data=plotting_ace[9:18,], aes(x=sample, y=ACE)) +
  geom_bar(stat="identity")
q

k <-ggplot(data=plotting_ace, aes(x=sample, y=ACE)) +
  geom_bar(stat="identity")
k

#Plot the evenness
evenness_df <- data.frame(J)
colnames(evenness_df) = c('evenness')
evenness_df$sample = sample_order
evenness_df$sample = factor(evenness_df$sample, levels = evenness_df$sample)

p <-ggplot(data=evenness_df[1:8, ], aes(x=sample, y=evenness)) +
  geom_bar(stat="identity")
p

q <-ggplot(data=evenness_df[9:18, ], aes(x=sample, y=evenness)) +
  geom_bar(stat="identity")
q

k <-ggplot(data=evenness_df, aes(x=sample, y=evenness)) +
  geom_bar(stat="identity")
k

#Plot the shannon
shannon_df = data.frame(H)
colnames(shannon_df) = c('Shannon')
shannon_df$sample = sample_order
shannon_df$sample = factor(shannon_df$sample, levels = shannon_df$sample)

p <-ggplot(data=shannon_df[1:8, ], aes(x=sample, y=Shannon)) +
  geom_bar(stat="identity")
p

q <-ggplot(data=shannon_df[9:18, ], aes(x=sample, y=Shannon)) +
  geom_bar(stat="identity")
q

k <-ggplot(data=shannon_df, aes(x=sample, y=Shannon)) +
  geom_bar(stat="identity")
k
