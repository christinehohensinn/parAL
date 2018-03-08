#'generate random weight matrices
#'
#'generate random weight matrices according to the dimensions of the original weight matrix
#'
#'According to an original data file with a defined weight matrix to estimate a linear logistic test model (LLTM), random weight matrices with the same dimensions of the original one are generated. For each generated weight matrix, an LLTM is estimated using the function \code{DRM} of the R package pcIRT. A Rasch model is estimated for the original data set using the function \code{DRM} from pcIRT. The item parameters estimated by the Rasch model are correlated with the item parameters reconstructed from the parameters estimated by the LLTM.
#'
#'
#'@param datmat data file
#'@param weightmat original weight matrix
#'@param repl number of weight matrices that are generated
#'@param prop01 vector providing the proportions of Zeros and Ones of the weight matrix to be generated
#'@param keepmat if TRUE, the generated random weight matrices are kept in an three-dimensional array
#'
#'@return 
#'\item{conv}{convergence of the estimated LLTMs of the generated weight matrices}
#'\item{raschpar}{estimated Rasch model item parameters for the original data set}
#'\item{lltmpar}{item parameters reconstructed from the estimated parameters of the LLTM of the generated weight matrices}
#'\item{corlr}{correlation of the re-constructed item parameters by the LLTM and the estimated item parameters of the original Rasch model}
#'\item{sumcorlr}{short summary descriptive statistics of the correlations }
#'\item{randommat}{array containing the generated weight matrices (if argument keepmat is set TRUE)}
#'
#' @importFrom pcIRT DRM
#'
#' @author Christine Hohensinn
#' 
#' @keywords random weight matrix, linear logistic test model
#' @references 
#' Baghaei, P. & Hohensinn, C. (submitted). A Method of Q-Matrix Validation for the Linear Logistic Test Model.
#' 
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in educational research. Acta Psychologica, 37, 359-374.
#'
#' @examples
#'
#' #generate three random weight matrices for example data set
#' data(exampledata)
#' data(exampleweight)
#' 
#' Wrandom(exdat, orig.weight, repl=3)
#'
#' @export Wrandom
#' @rdname wrandom



Wrandom <- function(datmat, weightmat, repl=100, prop01=c(0.5,0.5), keepmat=FALSE){
  
  datmat <- as.matrix(datmat)
  weightmat <- as.matrix(weightmat)
  if(!all(datmat %in% c(0,1))){stop("data matrix has to contain only 0 and 1s")}
  if(!all(weightmat %in% c(0,1))){stop("weight matrix has to contain only 0 and 1s")}
  if(repl > 1000){warning("The defined number of replications is high. The process may take some time!")}
  
#objects for results
conv <- rep(NA, length=repl)
lltmpar <- matrix(NA, nrow=nrow(weightmat), ncol=repl)
corlr <- rep(NA, length=repl)

if(keepmat){
  randommat <- array(NA, dim=c(nrow(weightmat), ncol(weightmat),repl))
} else {
  randommat <- NULL
}

#Rasch model with original data
rasch <- DRM(datmat)

#generate random weight matrices
for(i in 1:repl){
    qz <- matrix(sample(c(0,1), size=nrow(weightmat)*ncol(weightmat), replace=TRUE, prob = prop01), ncol=ncol(weightmat))
    
    lltz <- DRM(datmat, desmat=qz)
    lltzbeta.norm <- lltz$itempar-mean(lltz$itempar)
    
    #correlation
    korr <- cor(rasch$itempar, lltzbeta.norm)
    
    #save results
    conv[i] <- lltz$convergence
    lltmpar[,i] <- lltzbeta.norm
    corlr[i] <- korr
    
    if(keepmat){
      randommat[,,i] <- qz
    } else {
      randommat <- NULL
    }

    rm(lltz, lltzbeta.norm, qz)
}

sumcorlr <- c(min(corlr), median(corlr), mean(corlr), quantile(corlr, probs=c(0.95)), max(corlr))
names(sumcorlr) <- c("Min.", "Median", "Mean", "95%", "Max.")

res <- list(conv=conv, lltmpar=lltmpar, raschpar=rasch$itempar, corlr=corlr, sumcorlr=sumcorlr, randommat=randommat)


class(res) <- "wrandom"

res
}

