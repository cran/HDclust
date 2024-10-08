## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.show='hold'----------------------------------------------------------
library(HDclust)

## -----------------------------------------------------------------------------
# variable structure with one variable block of 3 variables and two mixture components, which corresponds to GMM model
Vb<- vb(1, dim=3, numst=2)
show(Vb)

## -----------------------------------------------------------------------------
# variable structure with two blocks. Dimensionality of data is 40. First block contains 10 variable with 3 mixture components. Second block has 30 variables with 5 mixture components. Variable order is natural.
Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
show(Vb)

## ----results='hide'-----------------------------------------------------------
# by default number of states in each block is varied from 1 to 9
data("faithful")
Vb <- vb(1, dim=2, numst=1)
set.seed(12345)
modelBIC <- hmmvbBIC(faithful, VbStructure=Vb)

## -----------------------------------------------------------------------------
show(modelBIC)

## ----fig1, fig.height = 5, fig.width = 5--------------------------------------
plot(modelBIC, xlab='# mixture components per block', ylab='BIC')

## ----results='hide'-----------------------------------------------------------
# user-provided configurations for the number of states in each block
data("sim3")
Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(1,1), varorder=list(c(1:10),c(11:40)))
set.seed(12345)
configs = list(c(1,2), c(3,5))
modelBIC <- hmmvbBIC(sim3[,1:40], VbStructure=Vb, configList=configs)

## -----------------------------------------------------------------------------
show(modelBIC)

## -----------------------------------------------------------------------------
# we just illustrate the loglikelihood of the first 10 data points beacause the data set is huge
getLoglikehd(modelBIC)[1:10]

## -----------------------------------------------------------------------------
data("sim3")
set.seed(12345)
Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
hmmvb <- hmmvbTrain(sim3[1:40], VbStructure=Vb)
show(hmmvb)

## -----------------------------------------------------------------------------
Vb <- getVb(hmmvb) # variable block structure
hmmChain <- getHmmChain(hmmvb) # list with HMM models
diagCov <- getDiagCov(hmmvb) # indicator whether covariance matrices are diagonal 
bic <- getBIC(hmmvb) # BIC value

# below we show HMM-VB parameters for the first variable block
numst <- getNumst(hmmChain[[1]]) # number of mixture components in variable block
prenumst <- getPrenumst(hmmChain[[1]]) # number of mixture components in the previous variable block
hmmParam <- getHmmParam(hmmChain[[1]]) # list with priors, transition probabilities, means, covariance matrices and other parameters of all states of HMM


## -----------------------------------------------------------------------------
data("sim2")
set.seed(12345)
hmmvb <- hmmvbTrain(sim2[1:5], searchControl=vbSearchControl(nperm=5), nthread=1)
show(hmmvb)

## ----results='hide'-----------------------------------------------------------
modelBIC <- hmmvbBIC(sim2[1:5], VbStructure=getVb(hmmvb), numst=1:15, nthread=1)

## -----------------------------------------------------------------------------
show(modelBIC)

## ----fig.show='hold', results='hide'------------------------------------------
data("faithful")
VbStructure <- vb(nb=1, dim=2, numst=1)
set.seed(12345)
modelBIC <- hmmvbBIC(faithful, VbStructure)

## -----------------------------------------------------------------------------
clust <- hmmvbClust(faithful, bicObj=modelBIC)
show(clust)

## ----fig2, fig.height = 5, fig.width = 5--------------------------------------
plot(faithful[,1], faithful[,2], xlab='eruptions', ylab='waiting', col=getClsid(clust))

## -----------------------------------------------------------------------------
clustParam <- getClustParam(clust)
mode <- clustParam$mode # a matrix with cluster modes
vseq <- clustParam$vseq # A list with integer vectors representing distinct Viterbi sequences for the dataset

## ----fig.show='hold', results='hide'------------------------------------------
# If variable block structure is unknown
data("sim2")
set.seed(12345)

# find variable block structure
hmmvb <- hmmvbTrain(sim2[,1:5], searchControl=vbSearchControl(nperm=5), nthread=1)

# refine number of states in variable block structure by model selection
modelBIC <- hmmvbBIC(sim2[,1:5], VbStructure=getVb(hmmvb), numst=1:15, nthread=1)

## -----------------------------------------------------------------------------
clust <- hmmvbClust(data=sim2[,1:5], bicObj=modelBIC)
show(clust)

## ----fig3, fig.height = 5, fig.width = 5--------------------------------------
palette(c(palette(), "purple", "brown")) # extend palette to show all 10 clusters
plot(clust)

## ----fig4, fig.height = 5, fig.width = 5--------------------------------------
plot(clust, method='PCA')

## ----fig.show='hold'----------------------------------------------------------
data("sim3")

# split data set in two halves
X1 = sim3[1:500,]
X2 = sim3[501:1000,]

set.seed(12345)
Vb <- vb(2, 40, c(10,30), c(3,5), list(c(1:10),c(11:40)))

# train HMM-VB on dataset X1 and cluster data
hmmvb <- hmmvbTrain(X1[1:40], VbStructure=Vb)
clust1 <- hmmvbClust(X1[1:40], model=hmmvb)
show(clust1)

# cluster data set X2 using HMM-VB and cluster parameters for dataset X1
clust2 <- hmmvbClust(X2[1:40], model=hmmvb, rfsClust=getClustParam(clust1))

show(clust2)

## ----fig.show='hold'----------------------------------------------------------
clust2 <- hmmvbClust(X2[1:40], model=hmmvb, rfsClust=getClustParam(clust1), control=clustControl(getlikelh = TRUE))
# we just illustrate the loglikelihood of the first 10 data points beacause the data set is huge
getLoglikehd(clust2)[1:10]

## ----fig.show='hold'----------------------------------------------------------
data("sim3")

set.seed(12345)
Vb <- vb(2, 40, c(10,30), c(3,5), list(c(1:10),c(11:40)))

# train HMM-VB
hmmvb <- hmmvbTrain(sim3[1:40], VbStructure=Vb)

# find all density modes
modes <- hmmvbFindModes(sim3[1:40], model=hmmvb)
show(modes)

# cluster density modes
mergedModes <- clustModes(modes, cutree.args = list(h=1.0))

show(mergedModes)

## ----fig5, fig.height = 5, fig.width = 5--------------------------------------
plot(mergedModes)

