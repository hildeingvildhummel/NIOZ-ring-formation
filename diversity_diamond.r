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

#Read the csv file created by DIAMOND annotation
df_A <- read.table("DIAMOND/CD_counts_genus.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
#df_A <- read.table("lineage_P.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
#df_A <- read.table("lineage_S.csv",header = TRUE, fileEncoding = "UTF-8-BOM", sep=",")
head(df_A)
#Sum the values if they are annotated the same
df <- aggregate(. ~ column, transform(df_A, column = column), sum)
#Set the annotation as index
rownames(df) <- df$column
df$column <- NULL
head(df)

#Normalize counts by TSS
df_tss <- tss(df)



#Convert dataframe to the right format
dt <- as.table(as.matrix(df_tss))

#CA plot
res.ca <- CA(dt, ncp = 2, graph = TRUE)
#CA information plots
get_eigenvalue(res.ca)
fviz_eig(res.ca)
fviz_ca_row(res.ca)
fviz_ca_col(res.ca)
fviz_ca_biplot(res.ca)


#Environmental info abount the samples 
#CD
#env <- data.frame(place = c('I', 'I', 'I', 'O', 'O', 'O', 'R', 'R'))
#P
env <- data.frame(place = c('O', 'O', 'O', 'R', 'R', 'R'))
#S
#env <- data.frame(place = c('O', 'O', 'R', 'R'))
env_table <- as.table(as.matrix(env))
rownames(env) <- sample_order <- colnames(df)

#Check influence of environmental arguments
envfit(t(dt), env)


#Shannon
H <- diversity(t(df_tss), index = 'shannon')
H
#Richness
S <- specnumber(t(df))
#Evenness
J <- H/log(S)


#Richness estimation
r <- estimateR(t(df), plot = TRUE)
r

#Convert Chao1 index results to dataframe
sample_order <- colnames(df)

chao_df <- data.frame(colnames(r), r[2,])
colnames(chao_df) <- c('sample', 'chao1')
plotting_chao <- chao_df[sample_order, ]


# lock in factor level order
plotting_chao$sample <- factor(plotting_chao$sample, levels = plotting_chao$sample)

#Plot the choa
p <-ggplot(data=plotting_chao, aes(x=sample, y=chao1)) +
  geom_bar(stat="identity")
p


#Plot the shannon
shannon_df = data.frame(H)
colnames(shannon_df) = c('Shannon')
shannon_df$sample = sample_order
shannon_df$sample = factor(shannon_df$sample, levels = shannon_df$sample)

p <-ggplot(data=shannon_df, aes(x=sample, y=Shannon)) +
  geom_bar(stat="identity")
p

