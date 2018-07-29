#' @rdname getNb
#' @exportMethod getNb
#' @examples
#' # accessing nb in instance of class VB
#' Vb <- vb(2, dim=10, bdim=c(4,6), numst=c(3,11), varorder=list(c(1:4),c(5:10)))
#' getNb(Vb)
#' 
setMethod("getNb","VB", function(object){return(object@nb)})

#' @rdname getNb   
#' @examples
#' # accessing nb in instance of class HMMVB 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getNb(hmmvb)  
#'     
#' @exportMethod getNb
setMethod("getNb","HMMVB", function(object){return(object@VbStructure@nb)})

#' @rdname getDim
#' @exportMethod getDim
#' @examples
#' # accessing dim in instance of class VB
#' Vb <- vb(nb=2, dim=10, bdim=c(4,6), numst=c(3,11), varorder=list(c(1:4),c(5:10)))
#' getDim(Vb)
#' 
#' @exportMethod getDim
setMethod("getDim","VB", function(object){return(object@dim)})

#' @rdname getDim   
#' @examples
#' # accessing dim in instance of class HMM 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getDim(getHmmChain(hmmvb)[[1]])  
#'     
#' @exportMethod getDim
setMethod("getDim","HMM", function(object){return(object@dim)})

#' @rdname getDim  
#' @examples
#' # accessing dim in instance of class HMMVB 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getDim(hmmvb)  
#'     
#' @exportMethod getDim
setMethod("getDim","HMMVB", function(object){return(object@VbStructure@dim)})

#' @rdname getBdim
#' @exportMethod getBdim
#' @examples
#' # accessing bdim in instance of class VB
#' Vb <- vb(2, dim=10, bdim=c(4,6), numst=c(3,11), varorder=list(c(1:4),c(5:10)))
#' getBdim(Vb)
#' 
#' @exportMethod getBdim
setMethod("getBdim","VB", function(object){return(object@bdim)})

#' @rdname getBdim   
#' @examples
#' # accessing bdim in instance of class HMMVB 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getBdim(hmmvb)  
#'     
#' @exportMethod getBdim
setMethod("getBdim","HMMVB", function(object){return(object@VbStructure@bdim)})

#' @rdname getNumst
#' @exportMethod getNumst
#' @examples
#' # accessing numst in instance of class VB
#' Vb <- vb(2, dim=10, bdim=c(4,6), numst=c(3,11), varorder=list(c(1:4),c(5:10)))
#' getNumst(Vb)
#' 
#' @exportMethod getNumst
setMethod("getNumst","VB", function(object){return(object@numst)})

#' @rdname getNumst   
#' @examples
#' # accessing getNumst in instance of class HMM 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getNumst(getHmmChain(hmmvb)[[1]])  
#'     
#' @exportMethod getNumst
setMethod("getNumst","HMM", function(object){return(object@numst)})

#' @rdname getNumst  
#' @examples
#' # accessing numst in instance of class HMMVB 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getNumst(hmmvb)  
#'     
#' @exportMethod getNumst
setMethod("getNumst","HMMVB", function(object){return(object@VbStructure@numst)})

#' @rdname getNumst  
#' @examples
#' \donttest{
#' # accessing numst used in model selection in instance of class HMMVBBIC 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' modelBIC <- hmmvbBIC(sim3[,1:40], VbStructure=Vb)
#' getNumst(modelBIC)}
#'     
#' @exportMethod getNumst
setMethod("getNumst","HMMVBBIC", function(object){return(object@numst)})

#' @rdname getVarorder
#' @exportMethod getVarorder
#' @examples
#' # accessing varorder in instance of class VB
#' Vb <- vb(2, dim=10, bdim=c(4,6), numst=c(3,11), varorder=list(c(1:4),c(5:10)))
#' getVarorder(Vb)
#' 
#' @exportMethod getVarorder
setMethod("getVarorder","VB", function(object){return(object@varorder)})

#' @rdname getVarorder  
#' @examples
#' # accessing varorder in instance of class HMMVB 
#' data("sim3")
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' set.seed(12345)
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' getVarorder(hmmvb)  
#'     
#' @exportMethod getVarorder
setMethod("getVarorder","HMMVB", function(object){return(object@VbStructure@varorder)})

#' @exportMethod getPrenumst
setMethod("getPrenumst","HMM", function(object){return(object@prenumst)})

#' @exportMethod getHmmParam
setMethod("getHmmParam","HMM", function(object){
  return(list(a=object@a, a00=object@a00, mean=object@mean, sigma=object@sigma,
              sigmaInv=object@sigmaInv, sigmaDetLog=object@sigmaDetLog))
         })

#' @exportMethod getVb
setMethod("getVb","HMMVB", function(object){return(object@VbStructure)})

#' @exportMethod getHmmChain
setMethod("getHmmChain","HMMVB", function(object){return(object@HmmChain)})

#' @rdname getBIC  
#' @exportMethod getBIC
setMethod("getBIC","HMMVB", function(object){return(object@BIC)})

#' @rdname getBIC
#' @exportMethod getBIC
setMethod("getBIC","HMMVBBIC", function(object){return(object@BIC)})

#' @exportMethod getConfigList
setMethod("getConfigList","HMMVBBIC", function(object){return(object@configList)})

#' @exportMethod getOptHMMVB
setMethod("getOptHMMVB","HMMVBBIC", function(object){return(object@optHMMVB)})

#' @exportMethod getDiagCov
setMethod("getDiagCov","HMMVB", function(object){return(object@diagCov)})

#' @exportMethod getClustParam
setMethod("getClustParam","HMMVBclust", function(object){return(object@clustParam)})

#' @exportMethod getClsid
setMethod("getClsid","HMMVBclust", function(object){return(object@clsid)})

#' @exportMethod getSize
setMethod("getSize","HMMVBclust", function(object){return(object@size)})

#' @exportMethod show
setMethod("show","VB",
          function(object){
            cat("------------------------\n")
            cat("Variable block structure\n")
            cat("------------------------\n\n")
            
            cat("Data dimensionality =",object@dim,"\n\n")
            cat("Number of variable blocks =",object@nb,"\n\n")
            cat("Dimensionality of variable blocks:",paste(object@bdim, sep=" "),"\n\n")
            cat("Number of states in variable blocks:",paste(object@numst, sep=" "),"\n\n")
            cat("Variable order in variable blocks:\n")
            
            for (i in 1:object@nb){
              cat("Block",i,":",paste(object@varorder[[i]], sep=" "),"\n")
            }
            
            cat("------------------------\n")
          }
)

#' @exportMethod show
setMethod("show","HMM",
          function(object){
            cat("---\n")
            cat("HMM\n")
            cat("---\n\n")
            
            cat("Variable block dimensionality =",object@dim,"\n\n")
            cat("Number of states =",object@numst,"\n\n")
            cat("State probabilities:",paste(object@a00, sep=" "),"\n\n")
            cat("Transition probability matrix:\n"); print(object@a)
            
            cat("\nMeans:\n"); print(object@mean)
            
            cat("\nCovariance matrices:\n");
            for (i in 1:object@numst){
              cat("\n[,,", i, "]\n", sep = "")
              print(object@sigma[[i]])
            }
            cat("---\n")
          }
)

#' @exportMethod show
setMethod("show","HMMVB",
          function(object){
            cat("--------------------------------------\n")
            cat("Hidden Markov Model on Variable Blocks\n")
            cat("--------------------------------------\n\n")
            print(object@VbStructure)
            
            cat("\nBIC =",object@BIC,"\n\n")
            cat("Covariance matrices are diagonal =",object@diagCov,"\n\n")
            
            cat("To show parameters of HMMs, access elements in HmmChain list.\n")
          }
)

#' @exportMethod  
setMethod("show","HMMVBclust",
          function(object){
            cat("------------------------------------------------------\n")
            cat("Clustering with Hidden Markov Model on Variable Blocks\n")
            cat("------------------------------------------------------\n\n")
            cat("Number of clusters =",object@clustParam$ncls,"\n\n")
            cat("Cluster sizes:",paste(object@size, sep=" "),"\n\n")
          }
)

#' @exportMethod  
setMethod("show","HMMVBBIC",
          function(object){
            if (length(object@numst) > 0){
              cat("------------------------------------------------------\n")
              cat("Optimal number of states:",object@numst[which.min(object@BIC)],"with BIC:",min(object@BIC),"\n")
              cat("------------------------------------------------------\n\n")
            }
            else{
              cat("------------------------------------------------------\n")
              cat("Optimal configuration:",paste(object@configList[[which.min(object@BIC)]], sep=" "),"with BIC:",min(object@BIC),"\n")
              cat("------------------------------------------------------\n\n")
              
            }
              
          }
)

#' @exportMethod  
setMethod("plot", signature(x="HMMVBclust",y="missing"),
          function(x, y, ...){
            if(dim(x@data)[2]>2)
              Y = as.data.frame(prcomp(x@data)$x)
            else
              Y = data.frame(PC1=x@data[,1],PC2=x@data[,2])
            
            plot(Y[,1:2], col = x@clsid, ...)
          }
)