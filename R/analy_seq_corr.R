### Analytical solution for sequence correlation functions. ###
### Created on Jun 30, 2015 ###
### Mingzhi Lin ###

## Generic variables:
## n: population size
## u: mutation rate per site per generation
## r: recombination rate per site per generation
## f: fragment size (averaged)
## q: selection rate
## a: size of alphabet set
## x: distance between two sites

CalcD <- function(n, u, r, f, q, a = 4) {
  M = 2*u*matrix(c(-1, 1/(a-1), 1, -1/(a-1)),
             nrow = 2, ncol = 2, byrow = T)
  I = diag(2)
  b = (q+2*r*f/n)
  P0 = matrix(c(1, 0), nrow = 2, ncol = 1)
  x = (b * solve(b * I - M)) %*% P0
  x[2]
}

CalcP2 <- function(n, u, r, f, q2, q3, q4, x, a = 4) {
  alpha = a/(a - 1);
  d = CalcD(n, u, r, f, q2, a)
  h1 = r*f*(1 - exp(-x/f))
  h2 = r*f - h1

  # S matrix elements.
  a1 = -(2 * h2 + 4 * h1 * (n - 1))/n - q2;
  a2 = 4 * h1 * (n - 2)/n;
  a3 = 0;
  b1 = 2 * (h1 + h2)/n + 1/3*q3;
  b2 = -(6 * h2 / n + 2 * h1) - q3;
  b3 = 2 * h1 * (n - 3)/n;
  c1 = 0;
  c2 = 8 * (h1 + h2)/n + 2/3*q4;
  c3 = -12 * (h1 + h2)/n - q4;

  # Construct S matrix
  S = matrix(
    c(a1, a2, a3, b1, b2, b3, c1, c2, c3),
    nrow = 3,
    ncol = 3,
    byrow = TRUE)

  B0 = matrix(
    c(2*h2/n + q2, 0, 0),
    ncol = 1,
    nrow = 3)
  B1 = matrix(
    c(4*h1/n, 2/3*q3+4*(h1+h2)/n, 1/3*q4+4*(h1+h2)/n),
    ncol = 1,
    nrow = 3)

  P0 = matrix(c(1, 0, 0), ncol = 1, nrow = 3)
  P1 = matrix(c((1-d), d, 0), ncol = 1, nrow = 3)

  I = Matrix::Diagonal(3)
  Mu = 4*u*matrix(
    c(-1, alpha/(2*a), 0, 1, -alpha/2, alpha/a, 0, 1/2, -alpha/a),
    ncol = 3, nrow = 3, byrow = T)
  K1 = kronecker(S, I)
  K2 = kronecker(I, Mu)
  K = K1 + K2

  Sol = -solve(K) %*% (kronecker(B1, P1) + kronecker(B0, P0))
  Sol[3]
}

CalcCt <- function(n, u, r, f, q2, q3, q4, x, a = 4) {
  d = CalcD(n, u, r, f, q2, a)
  p2 = CalcP2(n, u, r, f, q2, q3, q4, x, a)
  p2 - d^2
}
