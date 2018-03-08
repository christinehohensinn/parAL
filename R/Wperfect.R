#'generate data matrices according to LLTM
#'
#'generate data matrices according to the original weight matrix
#'
#'The original data matrix and the original weight matrix are taken and an LLTM is estimated. With the resulting estimated parameters and the original data matrix new random data matrices are generated according to the original weight matrix. For each generated data matrix, the Rasch model and an LLTM is estimated with the original weight matrix. The correlation between the item parameters of the Rasch model and the item parameters reconstructed from the parameters estimated by the LLTM is calculated. With this procedure some kind of an upper benchmark of the goodness-of-fit of the weight matrix is generated.
#'
#'@param datmat data file
#'@param weightmat original weight matrix
#'@param repl number of weight matrices that are generated
#'@param keepmat if TRUE, the generated random data matrices are kept in an three-dimensional array
#'
#'@return 
#'\item{convLLTM}{convergence of the estimated LLTMs of the generated data matrices}
#'\item{convRM}{convergence of the estimated Rasch models of the generated data matrices}
#'\item{lltmpar}{item parameters reconstructed from the estimated parameters of the LLTM of the generated weight matrices}
#'\item{raschpar}{estimated Rasch model item parameters for the original data set}
#'\item{corlr}{correlation of the re-constructed item parameters by the LLTM and the estimated item parameters of the original Rasch model}
#'\item{sumcorlr}{short summary descriptive statistics of the correlations}
#'\item{randommat}{array containing the generated data matrices (if argument keepmat is set TRUE)}
#'
#' @importFrom pcIRT simDRM
#'
#' @author Christine Hohensinn
#' 
#' @keywords random data matrix, linear logistic test model
#' @references 
#' Baghaei, P. & Hohensinn, C. (submitted). A Method of Q-Matrix Validation for the Linear Logistic Test Model.
#' 
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in educational research. Acta Psychologica, 37, 359-374.
#'
#' @examples
#'
#' #generate ten data sets according to the original weight matrix
#'  
#' data(exampledata)
#' data(exampleweight)
#' 
#' Wperfect(exdat, orig.weight, repl=10)
#'
#' @export Wperfect
#' @rdname wperfect


Wperfect <- function(datmat, weightmat, repl=50, keepmat=FALSE){
  datmat <- as.matrix(datmat)
  weightmat <- as.matrix(weightmat)
  if(!all(datmat %in% c(0,1))){stop("data matrix has to contain only 0 and 1s")}
  if(!all(weightmat %in% c(0,1))){stop("weight matrix has to contain only 0 and 1s")}
  if(repl > 500){warning("The defined number of replications is high. The process may take some time!")}

#LLTM with original data
  lltmo <- DRM(datmat, desmat=weightmat)

  #objects for results
  convRM <- rep(NA, length=repl)
  convLLTM <- rep(NA, length=repl)
  lltmpar <- matrix(NA, nrow=nrow(weightmat), ncol=repl)
  rmpar <- matrix(NA, nrow=nrow(weightmat), ncol=repl)
  corlr <- rep(NA, length=repl)
  
  if(keepmat){
    randommat <- array(NA, dim=c(nrow(datmat), ncol(datmat),repl))
  } else {
    randommat <- NULL
  }

for (m in 1:repl){
    dat <- simDRM(lltmo$itempar, persons=nrow(datmat))
    rm.data <- DRM(dat$datmat)
    lltz <- DRM(dat$datmat, desmat=weightmat)
    lltzbeta.norm <- lltz$itempar-mean(lltz$itempar)
    
    corr <- cor(rm.data$itempar, lltzbeta.norm)
    
    #save results
    convLLTM[m] <- lltz$convergence
    convRM[m] <- rm.data$convergence
    lltmpar[,m] <- lltzbeta.norm
    rmpar[,m] <- rm.data$itempar
    corlr[m] <- corr

    if(keepmat){
      randommat[,,m] <- dat$datmat
    } else {
      randommat <- NULL
    }
    
    rm(lltz, lltzbeta.norm, rm.data, corr)
    
  }
  
  sumcorlr <- c(min(corlr), median(corlr), mean(corlr), quantile(corlr, probs=c(0.05,0.95)), max(corlr))
  names(sumcorlr) <- c("Min.", "5%", "Median", "Mean", "95%", "Max.")
  
res <- list(convLLTM=convLLTM, convRM=convRM, lltmpar=lltmpar, raschpar=rmpar, corlr=corlr, sumcorlr=sumcorlr, randommat=randommat)
  

class(res) <- "wperfect"

res

}
  
