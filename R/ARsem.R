ARsem <-
function (model,data,order=0){

#require(lavaan,quietly=TRUE)

upper.triangle <- function (x)
{
    y <- x
    y[row(x) > col(y)] <- 0
    return(y)
}

    data <- as.matrix(data)
    nvar <- ncol(data)
    data.lagged <- embed(data, dimension=(order+1))

    # new column names and indices
	v.names <- rep(colnames(data), (order+1))
	prefix  <- rep(0:order, each=nvar)
	colnames(data.lagged) <- paste(v.names, prefix, sep="_")
    
    # create AR matrix
	model <- t(matrix(model,nvar))
	#if(any(rowSums(model)==0)) stop('Some variables in data are not part of the model')

	m <- upper.triangle(matrix(1,order+1,order+1))
	M <- m%x%model
	M <- M + m%x%diag(nvar)-diag((order+1)*nvar)
	colnames(M) <- colnames (data.lagged)
	rownames(M) <- colnames (M)

	regr <- logical(0)
	  for (i in 1:nrow(M)){
	    
	    for (j in 1:ncol(M)){
	      
	      ifelse(M[i,j] == 1, r <- paste(rownames(M)[i], '~', colnames(M)[j], '\n', sep= ' '), r <- logical(0))
	  
	         regr<-c(regr,r)
	            
		    
				}
			      }
	    regr <- paste(regr,collapse=" ")#regr	 

            # merge for lavaan 0.3-1
           # regr <- lavaan3_merge(regr)
             
	   # regr
	    fit <- sem(regr,data=as.data.frame(data.lagged)) # Warning blijft verschijnen
	    
	    fit	  
                  }

