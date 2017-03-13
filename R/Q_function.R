require(gsl)

Q_approx <- function(x, eps1 = 10^{-30},eps2 = 10^{-7}) {

  if(x <= 1) {
    return(expint_E1(x)*exp(x))
  } else {
  f_p <- eps1
  C_p <- eps1
  D_p <- 0
  Delta <- 2 + eps2
  j <- 1
  while( abs(Delta-1)>= eps2) {
    j <- j+1
    D_c <- x + 2*j -1 -(j-1)^{2}*D_p
    C_c <- x + 2*j -1 -(j-1)^{2}/C_p
    D_c <- 1/D_c
    Delta <- C_c*D_c
    f_c <- f_p*Delta
    f_p <- f_c
    C_p <- C_c
    D_p <- D_c
  }
  return(1/(x+1 + f_c))
  }
}
