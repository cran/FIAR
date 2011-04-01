spm_eigen <-
function (x) 
{
z <-  .Call("La_rs", x, TRUE, PACKAGE = "base")$values
z
#return(list(values = z$values[order(z$values,decreasing=TRUE)]))
}

