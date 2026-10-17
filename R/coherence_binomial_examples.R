# Binomial model: Examples 1 and 4
# X | theta ~ Binomial(n, theta), theta ~ Beta(1,1)

n <- 10
x <- 0:n

# Bayes factor for H: theta <= theta0
BF_interval <- function(x, theta0, n) {
  I <- pbeta(theta0, x + 1, n - x + 1)
  I_c <- pbeta(theta0, x + 1, n - x + 1, lower.tail = FALSE)
  (1 - theta0) / theta0 * I / I_c
}

# Bayes factor for the sharp hypothesis H: theta = theta0
BF_sharp <- function(x, theta0, n) {
  (n + 1) * dbinom(x, size = n, prob = theta0)
}

# Rejection and non-rejection regions
RR <- function(x, BF, K) x[BF <= K]
NR <- function(x, BF, K) x[BF > K]

# Rejection regions and coherence
coherence_check <- function(x, BF_A, BF_B, K_A, K_B) {
  RR_A <- RR(x, BF_A, K_A)
  RR_B <- RR(x, BF_B, K_B)
  list(RR_A = RR_A, RR_B = RR_B,
       coherent = all(RR_B %in% RR_A))
}

# Sufficient condition in Proposition 2
sufficient_condition <- function(K_A, K_B, pi_A, pi_B) {
  K_A * pi_A / (1 - pi_A) >= K_B * pi_B / (1 - pi_B)
}

# Loss-induced threshold
loss_threshold <- function(pi_H, c_H) {
  (1 - pi_H) / (c_H * pi_H)
}


# ============================================================
# Example 1: sufficient condition for coherence
# A: theta <= 0.3, B: theta <= 0.5
# ============================================================

theta_A <- pi_A <- 0.3
theta_B <- pi_B <- 0.5

BF_A <- BF_interval(x, theta_A, n)
BF_B <- BF_interval(x, theta_B, n)

# Bayes factors reported in Supplementary Table 1
data.frame(x, BF_A, BF_B)

# K_A = 3, K_B = 1: sufficient condition holds
K_A <- 3; K_B <- 1
sufficient_condition(K_A, K_B, pi_A, pi_B)
coherence_check(x, BF_A, BF_B, K_A, K_B)

# K_A = 1/2, K_B = 1: condition fails, but compatible decision holds
K_A <- 1/2; K_B <- 1
sufficient_condition(K_A, K_B, pi_A, pi_B)
coherence_check(x, BF_A, BF_B, K_A, K_B)

# Loss-induced thresholds: c_A = 1, c_B = 2
c_A <- 1; c_B <- 2
K_A <- loss_threshold(pi_A, c_A)
K_B <- loss_threshold(pi_B, c_B)

c(K_A = K_A, K_B = K_B)
sufficient_condition(K_A, K_B, pi_A, pi_B)
coherence_check(x, BF_A, BF_B, K_A, K_B)


# ============================================================
# Example 4: attained infimum
# A: theta = 1/2, B: theta <= 1/2
# ============================================================

theta_A <- theta_B <- 1/2

BF_A <- BF_sharp(x, theta_A, n)
BF_B <- BF_interval(x, theta_B, n)

# Bayes factors used in Example 4
data.frame(x, BF_A, BF_B)

# Fixed K_A = 1
K_A <- 1
RR_A <- RR(x, BF_A, K_A)
NR_A <- NR(x, BF_A, K_A)

list(RR_A = RR_A, NR_A = NR_A)

# K_B^*(1): minimum of BF_B over NR_A(1)
BF_B_NR_A <- BF_B[x %in% NR_A]
K_B_star <- min(BF_B_NR_A)
x_star <- NR_A[which.min(BF_B_NR_A)]

c(x_star = x_star, K_B_star = K_B_star)

# Boundary value K_B = BF_B(7) = K_B^*(1): not coherent
K_B <- K_B_star
coherence_check(x, BF_A, BF_B, K_A, K_B)

# Smaller threshold K_B = BF_B(8): coherent
K_B <- BF_B[x == 8]
coherence_check(x, BF_A, BF_B, K_A, K_B)
