## code modified from  https://github.com/JasonMathewsUCB/ElbowPoint
## A quick way of finding the elbow is to draw a line from the first to the last point of the curve and then find the data point that is farthest away from that line.

## Functions to Calculate ElbowPoint
library("pracma")

## functions
### Function that generates  Scalar Product at dimension level
dotProduct <- function(x,y){
    a <- matrix()
    x <- x
    y <- y
    for(i in 1:nrow(x)){
        a[i] <- ((x[i,1] * y[i,1]) + (x[i,2] * y[i,2]))

    }
    return (as.matrix(a))
}

ScalarMatrix <- function(x,y){
    a <- matrix(,length(x),2)
    x <- as.numeric(x)
    for(i in 1:length(x)){
        a[i,1] <- x[i] * y[1,1]
        a[i,2] <- x[i] * y[1,2]
    }
    return (a)
}


## Generate Test Data
test_data_allCoord <- function() {
    Trials <- seq(20,80,5)
    curve1 <- c(0.1, 0.08,	0.155555556,	0.214285714,	0.23,	        0.237037037,
                0.251428571,	0.254545455,	0.262962963,	0.267692308,	0.272727273,
                0.273333333,	0.273076923,	0.27394958,	    0.277037037,	0.278947368,
                0.281176471,	0.283597884,	0.285167464,	0.286956522,	0.287301587,
                0.288,	    0.288963211,	0.289506173,	0.289714286,	0.290185676,
                0.291358025,	0.291705069,	0.292672414,	0.293737374,	0.294497154,
                0.295714286,	0.296969697,	0.297933227,	0.299548872,	0.301424501,
                0.303783784,	0.306033376,	0.308669109,	0.311395349,	0.313747228,
                0.316190476,	0.318907988,	0.321663443,	0.324444444,	0.327417924,
                0.330382979,	0.333333333,	0.336420722,	0.339169811,	0.34204793,	0.344755245,
                0.347439353,	0.350097466,	0.352727273,	0.355326877,	0.357894737,	0.360429621,
                0.362930563,	0.365502646,	0.368032787,	0.370421836,	0.372775373,	0.375,	0.377375566,
                0.37962231,	0.382004264,	0.384341342,	0.386634461,	0.388806262,	0.390940236,	0.392962963,
                0.395025234,	0.397051597,	0.399042735,	0.401132578,	0.403311688,	0.405571383,
                0.407903674)
    m <- length(Trials)
    curve <- curve1[1:m]
    data.frame(Trials, curve)
}

# elbow function
elbow <- function(allCoord=NULL, test=FALSE) { # setting test to TRUE means get me a demo!
    if (test) allCoord <- test_data_allCoord()

    nPoints <- nrow(allCoord)
    Trials <- allCoord[,1]
    firstPoint = allCoord[1,]

    # get vector between first and last point - this is the line
    lineVec = tail(allCoord,n=1) - firstPoint
    allCoord <- as.matrix(allCoord)
    firstPoint = repmat(as.matrix(firstPoint), nPoints,1)

    # normalize the line vector
    Denominator = as.numeric(sqrt(lineVec[1]^2  + lineVec[2]^2))
    lineVecN = lineVec/Denominator

    # find the distance from each point to the line:
    # vector between all points and first point
    vecFromFirst = bsxfun("-", allCoord, firstPoint);
    scalarProduct = dotProduct(vecFromFirst, repmat(as.matrix(lineVecN),nPoints,1))
    vecFromFirstParallel = as.numeric(scalarProduct) * lineVecN

    vecFromFirstParallel <- ScalarMatrix(scalarProduct, lineVecN )
    NormVectorLine <- function(x){
        a <- matrix()
        for(i in 1:nrow(x)){
            a[i] <- sqrt(x[i,1]^2  + x[i,2]^2)
        }
        return(as.matrix(a))
    }
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = NormVectorLine(vecToLine)
    BestPointIndex <- which(distToLine == max(distToLine), arr.ind = TRUE)
    BestPointIndexRow <- BestPointIndex[1,1]
    idxOfBestPoint = c(allCoord[BestPointIndexRow,1],
                       allCoord[BestPointIndexRow,2])

    ### Plot to generate ElbowPoint
    plot(allCoord,type = "o")
         # col=ifelse(c(Trials,curve) == idxOfBestPoint , "red", "blue")
         # , cex = ifelse(c(Trials,curve) == idxOfBestPoint ,1.5, 1)
         # ,pch = ifelse(c(Trials,curve) == idxOfBestPoint ,19, 1))
    abline(v=idxOfBestPoint[1])

    idxOfBestPoint
}



