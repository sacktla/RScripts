
##Read file
cluster1 <- read.table(file.choose(), sep = "\t")
cluster2 <- read.table(file.choose(), sep = "\t")
cluster3 <- read.table(file.choose(), sep = "\t")
cluster4 <- read.table(file.choose(), sep = "\t")
cluster5 <- read.table(file.choose(), sep = "\t")
cluster6 <- read.table(file.choose(), sep = "\t")
cluster7 <- read.table(file.choose(), sep = "\t")
cluster8 <- read.table(file.choose(), sep = "\t")
cluster9 <- read.table(file.choose(), sep = "\t")

##Transform to matrix
processedCluster1 <- as.data.frame(cluster1[4:nrow(cluster1), 6:ncol(cluster1)], stringAsFactors = F)
head(processedCluster1)
charCluster1 <- lapply(processedCluster1, FUN = as.character)
head(charCluster1)
numericCluster1 <- lapply(charCluster1, FUN = as.numeric)
head(numericCluster1)
transformedCluster1 <- matrix(unlist(numericCluster1), ncol = 35, byrow = F)
head(transformedCluster1)
dim(transformedCluster1)

#Cluster2
ncol(cluster2)-5
processedCluster2 <- as.data.frame(cluster2[4:nrow(cluster2), 6:ncol(cluster2)], stringAsFactors = F)
head(processedCluster2)
charCluster2 <- lapply(processedCluster2, FUN = as.character)
head(charCluster2)
numericCluster2 <- lapply(charCluster2, FUN = as.numeric)
head(numericCluster2)
transformedCluster2 <- matrix(unlist(numericCluster2), ncol = 16, byrow = F)
head(transformedCluster2)
dim(transformedCluster2)

#Cluster3
colNum <- ncol(cluster3)-5
processedCluster3 <- as.data.frame(cluster3[4:nrow(cluster3), 6:ncol(cluster3)], stringAsFactors = F)
charCluster3 <- lapply(processedCluster3, FUN = as.character)
numericCluster3 <- lapply(charCluster3, FUN = as.numeric)
transformedCluster3 <- matrix(unlist(numericCluster3), ncol = colNum, byrow = F)
dim(transformedCluster3)

#Cluster4
colNum <- ncol(cluster4)-5
processedCluster4 <- as.data.frame(cluster4[4:nrow(cluster4), 6:ncol(cluster4)], stringAsFactors = F)
charCluster4 <- lapply(processedCluster4, FUN = as.character)
numericCluster4 <- lapply(charCluster4, FUN = as.numeric)
transformedCluster4 <- matrix(unlist(numericCluster4), ncol = colNum, byrow = F)
dim(transformedCluster4)

#Cluster5
colNum <- ncol(cluster5)-5
processedCluster5 <- as.data.frame(cluster5[4:nrow(cluster5), 6:ncol(cluster5)], stringAsFactors = F)
charCluster5 <- lapply(processedCluster5, FUN = as.character)
numericCluster5 <- lapply(charCluster5, FUN = as.numeric)
transformedCluster5 <- matrix(unlist(numericCluster5), ncol = colNum, byrow = F)
dim(transformedCluster5)

#Cluster6
colNum <- ncol(cluster6)-5
processedCluster6 <- as.data.frame(cluster6[4:nrow(cluster6), 6:ncol(cluster6)], stringAsFactors = F)
charCluster6 <- lapply(processedCluster6, FUN = as.character)
numericCluster6 <- lapply(charCluster6, FUN = as.numeric)
transformedCluster6 <- matrix(unlist(numericCluster6), ncol = colNum, byrow = F)
dim(transformedCluster6)

#Cluster7
colNum <- ncol(cluster7)-5
processedCluster7 <- as.data.frame(cluster7[4:nrow(cluster7), 6:ncol(cluster7)], stringAsFactors = F)
charCluster7 <- lapply(processedCluster7, FUN = as.character)
numericCluster7 <- lapply(charCluster7, FUN = as.numeric)
transformedCluster7 <- matrix(unlist(numericCluster7), ncol = colNum, byrow = F)
dim(transformedCluster7)

#Cluster8
colNum <- ncol(cluster8)-5
processedCluster8 <- as.data.frame(cluster8[4:nrow(cluster8), 6:ncol(cluster8)], stringAsFactors = F)
charCluster8 <- lapply(processedCluster8, FUN = as.character)
numericCluster8 <- lapply(charCluster8, FUN = as.numeric)
transformedCluster8 <- matrix(unlist(numericCluster8), ncol = colNum, byrow = F)
dim(transformedCluster8)

#Cluster9
colNum <- ncol(cluster9)-5
processedCluster9 <- as.data.frame(cluster9[4:nrow(cluster9), 6:ncol(cluster9)], stringAsFactors = F)
charCluster9 <- lapply(processedCluster9, FUN = as.character)
numericCluster9 <- lapply(charCluster9, FUN = as.numeric)
transformedCluster9 <- matrix(unlist(numericCluster9), ncol = colNum, byrow = F)
dim(transformedCluster9)

##Create the entire file with all the points
fullDataPoints <- cbind(transformedCluster1, transformedCluster2,
		 transformedCluster3,transformedCluster4, transformedCluster5,
			transformedCluster6,transformedCluster7,transformedCluster8,
			transformedCluster9)

dataFrameTotalCluster <- data.frame(geneNames, fullDataPoints)
head(dataFrameTotalCluster)



##Calculate SD for one gene for cluster 8
sd(transformedCluster8[1,])
##Calculate mean for one gene for cluster 8
mean(transformedCluster8[1,])

nrow(transformedCluster8)
ncol(transformedCluster8)
##Calculate SD for each gene in one cluster and the sd of that vector of SD
sdVector <- c()
for(i in 1:nrow(transformedCluster8)){
	geneVector <- c()	
	for(j in 1:ncol(transformedCluster8)){
		geneVector[j] <- transformedCluster8[i,j]
	} 
	sdGene <- sd(geneVector)
	sdVector[i] <- sdGene
}


#print out Standard Deviation vector
sdVector

#Calculate standard Deviation and mean of Standard Deviation vector

sd(sdVector)
mean(sdVector)

##Calculate the mean for all genes for each one of the clusters

calculateGeneMean <- function(matrixFile){

	orderedMeans <- c()
	for (i in 1:nrow(matrixFile)){
		orderedMeans[i] <- mean(matrixFile[i,])	
	}
	
	return(orderedMeans)
}

cluster1Means <- calculateGeneMean(transformedCluster1)
cluster2Means <- calculateGeneMean(transformedCluster2)
cluster3Means <- calculateGeneMean(transformedCluster3)
cluster4Means <- calculateGeneMean(transformedCluster4)
cluster5Means <- calculateGeneMean(transformedCluster5)
cluster6Means <- calculateGeneMean(transformedCluster6)
cluster7Means <- calculateGeneMean(transformedCluster7)
cluster8Means <- calculateGeneMean(transformedCluster8)
cluster9Means <- calculateGeneMean(transformedCluster9)

##Put all the means vectors in a matrix

matrixOfMeans <- function(...){
	meanMatrix <- cbind(...)
	
	return(meanMatrix)
	
}

meansMatrix <- matrixOfMeans(cluster1Means,cluster2Means, cluster3Means, 
			cluster4Means,cluster5Means, cluster6Means, cluster7Means,
			 cluster8Means,cluster9Means)

 
##Obtain the gene names, and create a data frame with row names as well as 
##The means matrix values

geneNames = cluster1[4:nrow(cluster1),3]

meansDataFrame <- data.frame(geneNames, meansMatrix)
##Confirm you have a squared Data frame
sapply(meansDataFrame, FUN = length)

##Calculate variance among these

#Obtain the names of the meansDataFrame
names(meansDataFrame)


variance = c()
for(i in 1:nrow(meansDataFrame)){
	geneValues <- c()
	for(j in 2:ncol(meansDataFrame)){
		geneValues[j-1] <- meansDataFrame[i,j]
	}
	variance[i] <- var(geneValues)	
}

#Confirm length
length(variance)

#See the variance vector
variance

#Combine the variance vector with the names

varianceDataFrame <- data.frame(geneNames, variance)

#Sort from highest to lowest variance

sortedValues <- varianceDataFrame[order(-varianceDataFrame$variance),]
tail(sortedValues)
#see sorted values
head(sortedValues)


##Means data frame

head(meansDataFrame)
sapply(meansDataFrame, class)


##I need to write a function that takes the means Data frame. Compares
#externally the means( mean1 - mean2)^2... Then create a matrix with these
#results. Normalize the results. and add them up row wise

calcMeanDiffSquare <- function(meanMatrix,i,j){
	meanDiffSquareVector <- c()
	vectorA <- c(meanMatrix[,i])
	vectorB <- c(meanMatrix[,j])
	for (num in 1:length(vectorA)){
		meanDiffSquareVector[num] <- (vectorA[num]-vectorB[num])^2
	}

	return(meanDiffSquareVector)
}
diffA <- calcMeanDiffSquare(meansDataFrame,2,3)
diffB <- calcMeanDiffSquare(meansDataFrame,3,4)
diffC <- calcMeanDiffSquare(meansDataFrame,4,5)
diffD <- calcMeanDiffSquare(meansDataFrame,5,6)
diffE <- calcMeanDiffSquare(meansDataFrame, 6,7)
diffF <- calcMeanDiffSquare(meansDataFrame,7,8)
diffG <- calcMeanDiffSquare(meansDataFrame,8,9)
diffH <- calcMeanDiffSquare(meansDataFrame,9,10)

##Obtaining the sum of the vectors, this will be used to normalize the data
sumA <- sum(diffA)
sumB <- sum(diffB)
sumC <- sum(diffC)
sumD <- sum(diffD)
sumE <- sum(diffE)
sumF <- sum(diffF)
sumG <- sum(diffG)
sumH <- sum(diffH)

##Create normalized vectors
normalizeDiffSquareVector <- function(vector, sum){
	normalizedVector <- c()
	for (i in 1:length(vector)){
		normalizedVector[i] <- vector[i]/max(vector)
	}

	return(normalizedVector)
}

normDiffA <- normalizeDiffSquareVector(diffA, sumA)
normDiffB <- normalizeDiffSquareVector(diffB, sumB)
normDiffC <- normalizeDiffSquareVector(diffC, sumC)
normDiffD <- normalizeDiffSquareVector(diffD, sumD)
normDiffE <- normalizeDiffSquareVector(diffE, sumE)
normDiffF <- normalizeDiffSquareVector(diffF, sumF)
normDiffG <- normalizeDiffSquareVector(diffG, sumG)
normDiffH <- normalizeDiffSquareVector(diffH, sumH)

##Create data frame with the values

normalizedDiffSquareMatrix <- cbind(normDiffA,normDiffB,
						normDiffC, normDiffD, normDiffE, normDiffF,
						normDiffG, normDiffH)

diffSquareDataFrame <- data.frame(geneNames, normalizedDiffSquareMatrix)

##create a vector with the sum of the normalized data. and rank the genes
sumDiffSquareVector <- c()
for (i in 1:nrow(diffSquareDataFrame)){
	sumDiffSquareVector[i] <- sum(diffSquareDataFrame[i,2:ncol(diffSquareDataFrame)])

}
##Create  data frame with the sum of Diff square vector and the gene names
length(sumDiffSquareVector)

sumDiffSquareDataFrame <- data.frame(geneNames, sumDiffSquareVector)
head(sumDiffSquareDataFrame)

#Reorder from high to low

orderedSumDiffSquare <- sumDiffSquareDataFrame[order(-sumDiffSquareDataFrame$sumDiffSquareVector),]

##Get headers and create data frame
clusterHeaders1 <- as.data.frame(cluster1[1,6:ncol(cluster1)], stringAsFactors = F)
clusterHeaders1A <- as.vector(t(clusterHeaders1))

clusterHeaders2 <- as.data.frame(cluster2[1,6:ncol(cluster2)], stringAsFactors = F)
clusterHeaders2A <- as.vector(t(clusterHeaders2))

clusterHeaders3 <- as.data.frame(cluster3[1,6:ncol(cluster3)], stringAsFactors = F)
clusterHeaders3A <- as.vector(t(clusterHeaders3))

clusterHeaders4 <- as.data.frame(cluster4[1,6:ncol(cluster4)], stringAsFactors = F)
clusterHeaders4A <- as.vector(t(clusterHeaders4))

clusterHeaders5 <- as.data.frame(cluster5[1,6:ncol(cluster5)], stringAsFactors = F)
clusterHeaders5A <- as.vector(t(clusterHeaders5))

clusterHeaders6 <- as.data.frame(cluster6[1,6:ncol(cluster6)], stringAsFactors = F)
clusterHeaders6A <- as.vector(t(clusterHeaders6))

clusterHeaders7 <- as.data.frame(cluster7[1,6:ncol(cluster7)], stringAsFactors = F)
clusterHeaders7A <- as.vector(t(clusterHeaders7))

clusterHeaders8 <- as.data.frame(cluster8[1,6:ncol(cluster8)], stringAsFactors = F)
clusterHeaders8A <- as.vector(t(clusterHeaders8))

clusterHeaders9 <- as.data.frame(cluster9[1,6:ncol(cluster9)], stringAsFactors = F)
clusterHeaders9A <- as.vector(t(clusterHeaders9))


totalFileHeaders <- c("Gene_Name",clusterHeaders1A, clusterHeaders2A, clusterHeaders3A,
			clusterHeaders4A, clusterHeaders5A, clusterHeaders6A,
			clusterHeaders7A, clusterHeaders8A, clusterHeaders9A)
totalFileHeaders
colnames(dataFrameTotalCluster) <- totalFileHeaders
head(dataFrameTotalCluster)


##Plot values
##This is taking into consideration the variance among means
##install ggplot2 package

install.packages("ggplot2")
library(ggplot2)

##try to do the qplot
head(dataFrameTotalCluster)
###MATRIX TOTAL CLUSTER (GENES IN X AND PATIENTS IN Y)####
tDataFrameTotalClusters <- t(dataFrameTotalCluster[,2:ncol(dataFrameTotalCluster)])
head(tDataFrameTotalClusters)
dataFrameTotalCluster[1,]
colnames(tDataFrameTotalClusters) <- dataFrameTotalCluster[,1]

###TABLED TOTAL CLUSTER DATA FRAME###
dataFrameClusters <- as.data.frame(as.table(tDataFrameTotalClusters))

###DATA FRAME OF TOTAL CLUSTER(GENES IN X AND PATIENTS IN Y)###
totalClusterDataFrameTransposed<- as.data.frame(tDataFrameTotalClusters)
head(totalClusterDataFrameTransposed)

##ADD AN EXTRA COLUMN INDICATING CLUSTER NUMBER###
cluster1Name <- rep("Cluster1", 35)
cluster2Name <- rep("Cluster2", 16)
cluster3Name <- rep("Cluster3", 51)
cluster4Name <- rep("Cluster4", 12)
cluster5Name <- rep("Cluster5", 28)
cluster6Name <- rep("Cluster6", 9)
cluster7Name <- rep("Cluster7", 78)
cluster8Name <- rep("Cluster8", 69)
cluster9Name <- rep("Cluster9", 11)

clusterIndex <- c(cluster1Name, cluster2Name, cluster3Name, cluster4Name,
			cluster5Name, cluster6Name, cluster7Name, cluster8Name,
			cluster9Name)
length(clusterIndex)

totalClusterDataFrameTransposed[,"Cluster Type"] <- clusterIndex
names(totalClusterDataFrameTransposed)[237] <- "ClusterType"
head(totalClusterDataFrameTransposed)
tail(totalClusterDataFrameTransposed)

##PLOTTING WITH GGPLOT2 USING THE VARIANCE ACROSS MEANS ANALYSIS
##genes are FOS, MMP3, EGR1, PLG, COL1A1, MMP9

#FOS
qplot(seq_along(totalClusterDataFrameTransposed$FOS), 
	totalClusterDataFrameTransposed$FOS, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "FOS EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#MMP3
qplot(seq_along(totalClusterDataFrameTransposed$MMP3), 
	totalClusterDataFrameTransposed$MMP3, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "MMP3 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#EGR1
qplot(seq_along(totalClusterDataFrameTransposed$EGR1), 
	totalClusterDataFrameTransposed$EGR1, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "EGR1 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#PLG
qplot(seq_along(totalClusterDataFrameTransposed$PLG), 
	totalClusterDataFrameTransposed$PLG, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "PLG EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#COL1A1
qplot(seq_along(totalClusterDataFrameTransposed$COL1A1), 
	totalClusterDataFrameTransposed$COL1A1, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "COL1A1 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()	

#MMP9
qplot(seq_along(totalClusterDataFrameTransposed$MMP9), 
	totalClusterDataFrameTransposed$MMP9, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "MMP9 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()


##Analysis of data using all clusters 

##Function that will compare all clusters
##Recycle outter cluster comparison and just run it as many times as you need it

calcMeanDiffSquare <- function(meanMatrix,i,j){
	meanDiffSquareVector <- c()
	vectorA <- c(meanMatrix[,i])
	vectorB <- c(meanMatrix[,j])
	for (num in 1:length(vectorA)){
		meanDiffSquareVector[num] <- (vectorA[num]-vectorB[num])^2
	}

	return(meanDiffSquareVector)
}

meanDiff12 <- calcMeanDiffSquare(meansDataFrame,1,2)
meanDiff13 <- calcMeanDiffSquare(meansDataFrame,1,3)
meanDiff14 <- calcMeanDiffSquare(meansDataFrame,1,4)
meanDiff15 <- calcMeanDiffSquare(meansDataFrame,1,5)
meanDiff16 <- calcMeanDiffSquare(meansDataFrame,1,6)
meanDiff17 <- calcMeanDiffSquare(meansDataFrame,1,7)
meanDiff18 <- calcMeanDiffSquare(meansDataFrame,1,8)
meanDiff19 <- calcMeanDiffSquare(meansDataFrame,1,9)
meanDiff23 <- calcMeanDiffSquare(meansDataFrame,2,3)
meanDiff24 <- calcMeanDiffSquare(meansDataFrame,2,4)
meanDiff25 <- calcMeanDiffSquare(meansDataFrame,2,5)
meanDiff26 <- calcMeanDiffSquare(meansDataFrame,2,6)
meanDiff27 <- calcMeanDiffSquare(meansDataFrame,2,7)
meanDiff28 <- calcMeanDiffSquare(meansDataFrame,2,8)
meanDiff29 <- calcMeanDiffSquare(meansDataFrame,2,9)
meanDiff34 <- calcMeanDiffSquare(meansDataFrame,3,4)
meanDiff35 <- calcMeanDiffSquare(meansDataFrame,3,5)
meanDiff36 <- calcMeanDiffSquare(meansDataFrame,3,6)
meanDiff37 <- calcMeanDiffSquare(meansDataFrame,3,7)
meanDiff38 <- calcMeanDiffSquare(meansDataFrame,3,8)
meanDiff39 <- calcMeanDiffSquare(meansDataFrame,3,9)
meanDiff45 <- calcMeanDiffSquare(meansDataFrame,4,5)
meanDiff46 <- calcMeanDiffSquare(meansDataFrame,4,6)
meanDiff47 <- calcMeanDiffSquare(meansDataFrame,4,7)
meanDiff48 <- calcMeanDiffSquare(meansDataFrame,4,8)
meanDiff49 <- calcMeanDiffSquare(meansDataFrame,4,9)
meanDiff56 <- calcMeanDiffSquare(meansDataFrame,5,6)
meanDiff57 <- calcMeanDiffSquare(meansDataFrame,5,7)
meanDiff58 <- calcMeanDiffSquare(meansDataFrame,5,8)
meanDiff59 <- calcMeanDiffSquare(meansDataFrame,5,9)
meanDiff67 <- calcMeanDiffSquare(meansDataFrame,6,7)
meanDiff68 <- calcMeanDiffSquare(meansDataFrame,6,8)
meanDiff69 <- calcMeanDiffSquare(meansDataFrame,6,9)
meanDiff78 <- calcMeanDiffSquare(meansDataFrame,7,8)
meanDiff79 <- calcMeanDiffSquare(meansDataFrame,7,9)
meanDiff89 <- calcMeanDiffSquare(meansDataFrame,8,9)


##put all the counts into a matrix

innerDiffSquareMatrix <- cbind(meanDiff12 ,meanDiff13 ,meanDiff14 ,meanDiff15 ,meanDiff16 ,meanDiff17 ,meanDiff18, 
					meanDiff19 ,meanDiff23 ,meanDiff24 ,meanDiff25 ,meanDiff26 ,meanDiff27 ,meanDiff28, 
					meanDiff29 ,meanDiff34 ,meanDiff35 ,meanDiff36 ,meanDiff37 ,meanDiff38 ,meanDiff39, 
					meanDiff45 ,meanDiff46 ,meanDiff47 ,meanDiff48 ,meanDiff49 ,meanDiff56 ,meanDiff57, 
					meanDiff58 ,meanDiff59 ,meanDiff67 ,meanDiff68 ,meanDiff69 ,meanDiff78 ,meanDiff79, 
					meanDiff89) 
dim(innerDiffSquareMatrix)

##transformedCluster1..9
##Calculate 1- cv(geneAclusterX)
##define a function that calculates the

##Function to calculate Coefficient of Variation

CV <- function(x){
	mean <- mean(x)
	std <- sd(x)
	coeffOfVar <- (std/mean)*100

	return(coeffOfVar)

}

##function to calculate from a gene cluster

calcCVcluster <- function(clusterName){
	cvVector <- c()
	for (i in 1:nrow(clusterName)){
		cvVector[i] <- 1 - CV(clusterName[i,])
	}

	return(cvVector)

}
head(transformedCluster1)
cluster1CoeffOfVarVector <- calcCVcluster(transformedCluster1)
cluster2CoeffOfVarVector <- calcCVcluster(transformedCluster2)
cluster3CoeffOfVarVector <- calcCVcluster(transformedCluster3)
cluster4CoeffOfVarVector <- calcCVcluster(transformedCluster4)
cluster5CoeffOfVarVector <- calcCVcluster(transformedCluster5)
cluster6CoeffOfVarVector <- calcCVcluster(transformedCluster6)
cluster7CoeffOfVarVector <- calcCVcluster(transformedCluster7)
cluster8CoeffOfVarVector <- calcCVcluster(transformedCluster8)
cluster9CoeffOfVarVector <- calcCVcluster(transformedCluster9)

rowBasedClusterCV <- rbind(cluster1CoeffOfVarVector, cluster2CoeffOfVarVector,
				cluster3CoeffOfVarVector, cluster4CoeffOfVarVector,
				cluster5CoeffOfVarVector, cluster6CoeffOfVarVector,
				cluster7CoeffOfVarVector, cluster8CoeffOfVarVector,
				cluster9CoeffOfVarVector)
head(rowBasedClusterCV)

##Matrix Operations
#######################################################################################


##Anova test
##Even though this needs to be done in a normally distibuted data set. I know that 
##some of the genes/cluster will not be complying but this should tell us whether
##across clusters they are from the same mean or different
##H0- These belong in the same group
##Reject at <.05
##Find the smallest pvalue. These will be your most different genes across clusters

anovaTest <- function(dataFrame, geneValue, clusterName){
	summary <- summary(aov(dataFrame$geneValue~dataFrame$clusterName))
	pValue <- summary[[1]]$'Pr(>F)'[1]

	return(pValue)
	
}


totalAnovaTest <- function(dataFrame, clusterNameHeader){
	pValues <- c()
	for (i in 1:(ncol(dataFrame)-1)){
		geneName <- colnames(dataFrame)[i]
		pValues[i] <- anovaTest(dataFrame,geneName,clusterNameHeader)
	}
	
	returnPvalues
}


##Needs to convert it with the gene name(use gene name vector) and convert it to a 
##data frame.

##totalClusterDataFrameTransposed

pValuesAnova <- totalAnovaTest(totalClusterDataFrameTransposed, ClusterType)

names("<- paste("Q", 1:ncol(totalClusterDataFrameTransposed)+1)
geneNames
pValues <- c()
for(i in 1:(ncol(totalClusterDataFrameTransposed)-1)){
	x <- c(unlist(totalClusterDataFrameTransposed[i]))
	tempSummary <- summary(aov(x ~ totalClusterDataFrameTransposed$ClusterType))
	pValues[i] <- tempSummary[[1]]$'Pr(>F)'[1]
}
pValues

##Create data frame with pvalues and gene names

pValueTable <- data.frame(geneNames,pValues)
pValueTable

##reorder from low to high

orderedPvalueTable <- pValueTable[order(pValueTable$pValues),]
head(orderedPvalueTable)

##Plot JUNB, TGFBR2, PIM1, IGF1, PTGS2, FAS (Highly variant)


##JUNB
qplot(seq_along(totalClusterDataFrameTransposed$JUNB), 
	totalClusterDataFrameTransposed$JUNB, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "JUNB EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

##TGFBR2
qplot(seq_along(totalClusterDataFrameTransposed$TGFBR2), 
	totalClusterDataFrameTransposed$TGFBR2, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "TGFBR2 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#PIM1
qplot(seq_along(totalClusterDataFrameTransposed$PIM1), 
	totalClusterDataFrameTransposed$PIM1, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "PIM1 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

##IGF1
qplot(seq_along(totalClusterDataFrameTransposed$IGF1), 
	totalClusterDataFrameTransposed$IGF1, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "IGF1 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#PTGS2
qplot(seq_along(totalClusterDataFrameTransposed$PTGS2), 
	totalClusterDataFrameTransposed$PTGS2, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "PTGS2 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#FAS
qplot(seq_along(totalClusterDataFrameTransposed$FAS), 
	totalClusterDataFrameTransposed$FAS, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "FAS EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

##Plot MST1R, TFRC, CDH1, TUBB, PGK1, DEK (low variant)

#DEK
qplot(seq_along(totalClusterDataFrameTransposed$DEK), 
	totalClusterDataFrameTransposed$DEK, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "DEK EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#PGK1
qplot(seq_along(totalClusterDataFrameTransposed$PGK1), 
	totalClusterDataFrameTransposed$PGK1, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "PGK1 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#TUBB
qplot(seq_along(totalClusterDataFrameTransposed$TUBB), 
	totalClusterDataFrameTransposed$TUBB, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "TUBB EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#CDH1
qplot(seq_along(totalClusterDataFrameTransposed$CDH1), 
	totalClusterDataFrameTransposed$CDH1, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "CDH1 EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#TFRC
qplot(seq_along(totalClusterDataFrameTransposed$TFRC), 
	totalClusterDataFrameTransposed$TFRC, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "TFRC EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()

#MST1R
qplot(seq_along(totalClusterDataFrameTransposed$MST1R), 
	totalClusterDataFrameTransposed$MST1R, 
	col = totalClusterDataFrameTransposed$ClusterType,
	xlab = "Samples", ylab = "Expression", main = "MST1R EXPRESSION")+ 
	stat_summary(fun.y = mean, geom = "point") +
	geom_smooth()


############################################################################################################
##PCA ANALYSIS

testPCA <- prcomp(totalClusterDataFrameTransposed[,1:(ncol(totalClusterDataFrameTransposed)-1)], center = TRUE, scale. = TRUE)
dim(testPCA)
testPCA
install.packages("devtools")
