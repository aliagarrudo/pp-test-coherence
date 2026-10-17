# Example 3: coherence calibration in the normal--normal model
# X_i | theta ~ N(theta, sigma2), theta ~ N(0,1)
# A: theta = 0, B: |theta| <= theta0

# Prior probability of B
pi_B <- function(theta0) {
  2 * pnorm(theta0) - 1
}

# Posterior probability of B
g_n <- function(x, n, sigma2, theta0) {
  sigma <- sqrt(sigma2)
  
  pnorm(
    ((n + sigma2) * theta0 - n * x) /
      (sigma * sqrt(n + sigma2))
  ) -
    pnorm(
      (-(n + sigma2) * theta0 - n * x) /
        (sigma * sqrt(n + sigma2))
    )
}

# Bayes factors
BF_A <- function(x, n, sigma2) {
  sigma <- sqrt(sigma2)
  
  sqrt(n + sigma2) / sigma *
    exp(
      -n^2 * x^2 /
        (2 * sigma2 * (n + sigma2))
    )
}

BF_B <- function(x, n, sigma2, theta0) {
  piB <- pi_B(theta0)
  g <- g_n(x, n, sigma2, theta0)
  
  (1 - piB) / piB * g / (1 - g)
}

# ------------------------------------------------------------
# Fixed K_B: minimum K_A ensuring coherence
# ------------------------------------------------------------

x_B <- function(K_B, n, sigma2, theta0) {
  
  stopifnot(
    K_B > 0,
    K_B < BF_B(0, n, sigma2, theta0)
  )
  
  f <- function(x) {
    BF_B(x, n, sigma2, theta0) - K_B
  }
  
  upper <- 1
  while (f(upper) > 0) {
    upper <- 2 * upper
  }
  
  uniroot(f, c(0, upper))$root
}

K_A_star <- function(K_B, n, sigma2, theta0) {
  xb <- x_B(K_B, n, sigma2, theta0)
  BF_A(xb, n, sigma2)
}

# ------------------------------------------------------------
# Fixed K_A: boundary value for K_B
# ------------------------------------------------------------

x_A <- function(K_A, n, sigma2) {
  sigma <- sqrt(sigma2)
  
  stopifnot(
    K_A > 0,
    K_A < BF_A(0, n, sigma2)
  )
  
  sqrt(
    2 * sigma2 * (n + sigma2) / n^2 *
      log(
        sqrt(n + sigma2) /
          (sigma * K_A)
      )
  )
}

K_B_boundary <- function(K_A, n, sigma2, theta0) {
  xa <- x_A(K_A, n, sigma2)
  BF_B(xa, n, sigma2, theta0)
}

# ------------------------------------------------------------
# Numerical specification used in Example 3
# ------------------------------------------------------------

n <- 1
sigma2 <- 1
theta0 <- 1/2

piB <- pi_B(theta0)

# Direct calibration
K_B <- (1 - piB) / piB
xb <- x_B(K_B, n, sigma2, theta0)
KA_star <- K_A_star(K_B, n, sigma2, theta0)

# Reverse calibration
K_A <- 1/2
xa <- x_A(K_A, n, sigma2)
KB_boundary <- K_B_boundary(K_A, n, sigma2, theta0)

c(
  pi_B = piB,
  BF_B_0 = BF_B(0, n, sigma2, theta0),
  x_B = xb,
  K_A_star = KA_star,
  x_A = xa,
  K_B_boundary = KB_boundary
)

# ------------------------------------------------------------
# Figure 1: Bayes factors and incoherent region
# ------------------------------------------------------------

K_common <- 1/2

# Boundaries under the common threshold K_A = K_B = 1/2
x_B_common <- x_B(
  K_common,
  n,
  sigma2,
  theta0
)

x_A_common <- x_A(
  K_common,
  n,
  sigma2
)

x_grid <- seq(
  -2.5,
  2.5,
  length.out = 1500
)

pdf(
  "normal_normal.pdf",
  width = 10,
  height = 4.5
)

par(
  mar = c(4, 4.5, 1, 13)
)

plot(
  x_grid,
  BF_A(x_grid, n, sigma2),
  type = "n",
  xlim = c(-2.5, 2.5),
  ylim = c(0.35, 1.82),
  xlab = "x",
  ylab = "Bayes Factor",
  cex.lab = 1.1
)

# Regions in which B is rejected but A is not
rect(
  -x_A_common,
  0.35,
  -x_B_common,
  1.82,
  col = adjustcolor(
    "gray50",
    alpha.f = 0.18
  ),
  border = NA
)

rect(
  x_B_common,
  0.35,
  x_A_common,
  1.82,
  col = adjustcolor(
    "gray50",
    alpha.f = 0.18
  ),
  border = NA
)

# Bayes factors
lines(
  x_grid,
  BF_A(x_grid, n, sigma2),
  col = "red",
  lwd = 2
)

lines(
  x_grid,
  BF_B(x_grid, n, sigma2, theta0),
  col = "blue",
  lwd = 2
)

# Adaptive thresholds
abline(
  h = KA_star,
  col = "red",
  lty = 2,
  lwd = 2
)

abline(
  h = K_B,
  col = "blue",
  lty = 2,
  lwd = 2
)

# Boundaries associated with adaptive K_B
abline(
  v = c(-xb, xb),
  col = "gray40",
  lty = 3,
  lwd = 1
)

# Common threshold K_A = K_B = 1/2
abline(
  h = K_common,
  col = "black",
  lty = 4,
  lwd = 1.5
)

# Legend outside the plotting region
legend(
  "topright",
  inset = c(-0.38, 0),
  xpd = NA,
  legend = c(
    expression(BF[A](x)),
    expression(BF[B](x)),
    expression(K[A]^"*"),
    expression(K[B]),
    expression(
      paste(
        K[A],
        " = ",
        K[B],
        " = ",
        frac(1, 2)
      )
    ),
    "Incoherent region for K = 1/2"
  ),
  col = c(
    "red",
    "blue",
    "red",
    "blue",
    "black",
    adjustcolor(
      "gray50",
      alpha.f = 0.40
    )
  ),
  lty = c(
    1,
    1,
    2,
    2,
    4,
    NA
  ),
  lwd = c(
    2,
    2,
    2,
    2,
    1.5,
    NA
  ),
  pch = c(
    NA,
    NA,
    NA,
    NA,
    NA,
    15
  ),
  pt.cex = 1.2,
  cex = 0.9,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------
# Figure 2: K_A^*(K_B) and coherent/incoherent regions
# ------------------------------------------------------------

# Maximum values of the Bayes factors
K_B_max <- BF_B(
  0,
  n,
  sigma2,
  theta0
)

K_A_max <- BF_A(
  0,
  n,
  sigma2
)

# Minimum K_A ensuring coherence for a fixed K_B
K_A_star_fun <- function(KB) {
  
  if (KB <= 0) {
    return(NA_real_)
  }
  
  # If B is rejected throughout the sample space,
  # A must also be rejected throughout the sample space.
  if (KB >= K_B_max) {
    return(K_A_max)
  }
  
  K_A_star(
    KB,
    n,
    sigma2,
    theta0
  )
}

# Grid for K_B
grid_KB <- seq(
  1e-4,
  K_B_max,
  length.out = 500
)

# Corresponding coherent boundary for K_A
grid_KA <- sapply(
  grid_KB,
  K_A_star_fun
)

tab_KA_KB <- data.frame(
  K_B = grid_KB,
  K_A_star = grid_KA
)

pdf(
  "relacao_KAstar_em_funcao_de_KB.pdf",
  width = 10,
  height = 4.5
)

par(
  mar = c(4, 4.5, 1, 13)
)

plot(
  grid_KB,
  grid_KA,
  type = "n",
  xlim = c(0, K_B_max),
  ylim = c(0, K_A_max),
  xlab = expression(K[B]),
  ylab = expression(K[A]^"*"),
  cex.lab = 1.1
)

# Incoherent region: below K_A^*(K_B)
polygon(
  c(
    grid_KB,
    rev(grid_KB)
  ),
  c(
    rep(
      0,
      length(grid_KB)
    ),
    rev(grid_KA)
  ),
  col = adjustcolor(
    "red",
    alpha.f = 0.12
  ),
  border = NA
)

# Coherent region: above K_A^*(K_B)
polygon(
  c(
    grid_KB,
    rev(grid_KB)
  ),
  c(
    grid_KA,
    rep(
      K_A_max,
      length(grid_KB)
    )
  ),
  col = adjustcolor(
    "blue",
    alpha.f = 0.08
  ),
  border = NA
)

# Coherence boundary
lines(
  grid_KB,
  grid_KA,
  col = "black",
  lwd = 2
)

# Common-threshold line K_A = K_B
abline(
  a = 0,
  b = 1,
  col = "gray40",
  lty = 2,
  lwd = 1.5
)

# Adaptive K_B
abline(
  v = K_B,
  col = "blue",
  lty = 3,
  lwd = 1
)

# Corresponding K_A^*
abline(
  h = KA_star,
  col = "red",
  lty = 3,
  lwd = 1
)

# Adaptive pair
points(
  K_B,
  KA_star,
  pch = 19,
  cex = 1
)

# Common-threshold pair K_A = K_B = 1/2
points(
  K_common,
  K_common,
  pch = 4,
  cex = 1.2,
  lwd = 2
)

# Region labels
text(
  0.58,
  1.18,
  "Coherent",
  col = "darkblue",
  cex = 1.1
)

text(
  0.52,
  0.23,
  "Incoherent",
  col = "darkred",
  cex = 1.1
)

KB_lab <- formatC(
  K_B,
  format = "f",
  digits = 3
)

KA_lab <- formatC(
  KA_star,
  format = "f",
  digits = 3
)

legend(
  "right",
  inset = c(-0.3, 0),
  xpd = NA,
  legend = c(
    expression(
      K[A]^"*" * "(" * K[B] * ")"
    ),
    expression(
      K[A] == K[B]
    ),
    bquote(
      K[B] == .(KB_lab)
    ),
    bquote(
      K[A]^"*" *
        "(" * .(KB_lab) * ")" ==
        .(KA_lab)
    ),
    expression(
      paste(
        K[A],
        " = ",
        K[B],
        " = ",
        frac(1, 2)
      )
    ),
    "Coherent region",
    "Incoherent region"
  ),
  col = c(
    "black",
    "gray40",
    "blue",
    "red",
    "black",
    adjustcolor(
      "blue",
      alpha.f = 0.25
    ),
    adjustcolor(
      "red",
      alpha.f = 0.25
    )
  ),
  lty = c(
    1,
    2,
    3,
    3,
    NA,
    NA,
    NA
  ),
  lwd = c(
    2,
    1.5,
    1,
    1,
    NA,
    NA,
    NA
  ),
  pch = c(
    NA,
    NA,
    NA,
    NA,
    4,
    15,
    15
  ),
  pt.cex = c(
    NA,
    NA,
    NA,
    NA,
    1.2,
    1.2,
    1.2
  ),
  cex = 0.9,
  bty = "n"
)

dev.off()
