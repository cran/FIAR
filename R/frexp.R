frexp <- function(x) {
  if (x == 0) 0
  else floor(log10(abs(x)))
}
