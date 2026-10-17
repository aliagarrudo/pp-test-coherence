# ============================================================================
# Supplementary R code: anorexia application
# ============================================================================
# Data: MASS::anorexia
# Response: weight change = Postwt - Prewt
# Groups: CBT (group 1) and FT (group 2)
#
# Hypotheses:
#   A: mu_CBT = mu_FT = 0        <=> eta = 0 and delta = 0
#   B: mu_CBT = mu_FT            <=> delta = 0
#
# Parameterisation:
#   eta   = (n1 * mu_CBT + n2 * mu_FT) / n
#   delta = mu_FT - mu_CBT
#
# Decision rule:
#   reject H when BF_H <= K_H
#
# The response is divided by the pooled residual standard deviation estimated
# from the same data and is not centred. The resulting residual variance is
# then treated as known and equal to one (plug-in analysis).
#
# This script reproduces the numerical calculations, the rejection-region
# figure used in the paper, and the optional interactive 3D visualisation.
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Packages and output directory
# ----------------------------------------------------------------------------

required_packages <- c("MASS", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(MASS)
library(ggplot2)

output_dir <- "anorexia_application_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ----------------------------------------------------------------------------
# 2. Data preparation and plug-in residual scale
# ----------------------------------------------------------------------------

data("anorexia", package = "MASS")

dat <- subset(anorexia, Treat %in% c("CBT", "FT"))
dat$Treat <- factor(dat$Treat, levels = c("CBT", "FT"))
dat$change <- dat$Postwt - dat$Prewt

# Pooled within-group residual standard deviation on the original scale.
sigma_hat_original <- sigma(lm(change ~ Treat, data = dat))

# Rescale without centring.
dat$Y <- dat$change / sigma_hat_original

# By construction, the plug-in residual variance on this scale is one.
sigma2 <- 1

Y_CBT <- dat$Y[dat$Treat == "CBT"]
Y_FT  <- dat$Y[dat$Treat == "FT"]

n1 <- length(Y_CBT)
n2 <- length(Y_FT)
n  <- n1 + n2
q  <- n1 * n2 / n

mean_CBT <- mean(Y_CBT)
mean_FT  <- mean(Y_FT)

# Observed sufficient statistics.
ybar0 <- (n1 * mean_CBT + n2 * mean_FT) / n
d0    <- mean_FT - mean_CBT

# Group means on the original weight-change scale.
mean_change_CBT <- mean(dat$change[dat$Treat == "CBT"])
mean_change_FT  <- mean(dat$change[dat$Treat == "FT"])


# ----------------------------------------------------------------------------
# 3. Prior specification and Bayes factor functions
# ----------------------------------------------------------------------------

m_eta      <- 0
m_delta    <- 0
tau_eta2   <- 1
tau_delta2 <- 1

# Quantities appearing in Example 5.
M_eta <- sqrt((sigma2 + n * tau_eta2) / sigma2) *
  exp(m_eta^2 / (2 * tau_eta2))

M_B <- sqrt((sigma2 + q * tau_delta2) / sigma2) *
  exp(m_delta^2 / (2 * tau_delta2))

a_eta <- n^2 * tau_eta2 /
  (2 * sigma2 * (sigma2 + n * tau_eta2))

a_delta <- q^2 * tau_delta2 /
  (2 * sigma2 * (sigma2 + q * tau_delta2))

y_A_star <- -(sigma2 * m_eta) / (n * tau_eta2)
d_B_star <- -(sigma2 * m_delta) / (q * tau_delta2)

M_A <- M_eta * M_B

G_eta <- function(ybar) {
  M_eta * exp(-a_eta * (ybar - y_A_star)^2)
}

BF_B <- function(d) {
  M_B * exp(-a_delta * (d - d_B_star)^2)
}

BF_A <- function(ybar, d) {
  G_eta(ybar) * BF_B(d)
}


# ----------------------------------------------------------------------------
# 4. Observed Bayes factors and threshold calibration
# ----------------------------------------------------------------------------

BF_A_obs <- BF_A(ybar0, d0)
BF_B_obs <- BF_B(d0)

K_B <- 1
K_A_uncalibrated <- 1

if (!(K_B > 0 && K_B < M_B)) {
  stop("This application assumes 0 < K_B < M_B.")
}

# Exact calibration from Example 5.
K_A_star <- M_eta * K_B

# Boundary of NR_B(K_B):
#   |d - d_B_star| < c_B.
c_B <- sqrt(log(M_B / K_B) / a_delta)
d_B_boundary <- d_B_star + c(-c_B, c_B)

# Points at which BF_A reaches K_A_star on the boundary BF_B = K_B.
tangency_points <- data.frame(
  ybar = rep(y_A_star, 2),
  d = d_B_boundary
)

# Numerical checks of the calibration.
stopifnot(
  isTRUE(all.equal(BF_B(d_B_boundary), rep(K_B, 2), tolerance = 1e-10)),
  isTRUE(all.equal(
    BF_A(tangency_points$ybar, tangency_points$d),
    rep(K_A_star, 2),
    tolerance = 1e-10
  ))
)


# ----------------------------------------------------------------------------
# 5. Observed and illustrative decisions
# ----------------------------------------------------------------------------

reject <- function(BF, K) BF <= K

# Illustrative sample point used in the paper.
ybar_inc <- 0
d_inc <- 0.6

BF_A_inc <- BF_A(ybar_inc, d_inc)
BF_B_inc <- BF_B(d_inc)

decision_table <- data.frame(
  point = c("Observed", "Illustrative"),
  ybar = c(ybar0, ybar_inc),
  d = c(d0, d_inc),
  BF_A = c(BF_A_obs, BF_A_inc),
  BF_B = c(BF_B_obs, BF_B_inc),
  reject_A_K1 = c(
    reject(BF_A_obs, K_A_uncalibrated),
    reject(BF_A_inc, K_A_uncalibrated)
  ),
  reject_B_K1 = c(
    reject(BF_B_obs, K_B),
    reject(BF_B_inc, K_B)
  ),
  reject_A_calibrated = c(
    reject(BF_A_obs, K_A_star),
    reject(BF_A_inc, K_A_star)
  )
)

cat("\n--- Data and sufficient statistics ---\n")
cat("n1 (CBT) =", n1, "\n")
cat("n2 (FT)  =", n2, "\n")
cat("n         =", n, "\n")
cat("q         =", format(q, digits = 8), "\n")
cat("Residual SD on original scale =", format(sigma_hat_original, digits = 8), "\n")
cat("Mean weight change, CBT =", format(mean_change_CBT, digits = 8), "lb\n")
cat("Mean weight change, FT  =", format(mean_change_FT, digits = 8), "lb\n")
cat("Observed ybar =", format(ybar0, digits = 8), "\n")
cat("Observed d    =", format(d0, digits = 8), "\n")

cat("\n--- Bayes factor constants ---\n")
cat("M_eta   =", format(M_eta, digits = 8), "\n")
cat("M_B     =", format(M_B, digits = 8), "\n")
cat("a_eta   =", format(a_eta, digits = 8), "\n")
cat("a_delta =", format(a_delta, digits = 8), "\n")

cat("\n--- Observed Bayes factors ---\n")
cat("BF_A(ybar0, d0) =", format(BF_A_obs, digits = 8), "\n")
cat("BF_B(d0)        =", format(BF_B_obs, digits = 8), "\n")

cat("\n--- Calibration ---\n")
cat("K_B      =", K_B, "\n")
cat("K_A      =", K_A_uncalibrated, "(uncalibrated)\n")
cat("K_A_star =", format(K_A_star, digits = 8), "(calibrated)\n")

cat("\n--- Decisions ---\n")
print(decision_table, row.names = FALSE)


# ----------------------------------------------------------------------------
# 6. Figure: rejection regions in the (ybar, d) plane
# ----------------------------------------------------------------------------

ellipse_BF_A <- function(K, type, n_points = 800) {
  if (!(K > 0 && K < M_A)) {
    stop("The level K must satisfy 0 < K < max(BF_A).")
  }

  angle <- seq(0, 2 * pi, length.out = n_points)
  radius_y <- sqrt(log(M_A / K) / a_eta)
  radius_d <- sqrt(log(M_A / K) / a_delta)

  data.frame(
    ybar = y_A_star + radius_y * cos(angle),
    d = d_B_star + radius_d * sin(angle),
    K_A = K,
    type = type
  )
}

ellipse_uncalibrated <- ellipse_BF_A(
  K_A_uncalibrated,
  "uncalibrated"
)

ellipse_calibrated <- ellipse_BF_A(
  K_A_star,
  "calibrated"
)

ellipses <- rbind(ellipse_uncalibrated, ellipse_calibrated)

max_y <- max(abs(c(ellipses$ybar, ybar0)), 0.80)
max_d <- max(abs(c(ellipses$d, d0, d_B_boundary)), 0.90)

y_limit <- 1.22 * max_y + 0.15
d_limit <- 1.10 * max_d + 0.08

grid <- expand.grid(
  ybar = seq(-y_limit, y_limit, length.out = 450),
  d = seq(-d_limit, d_limit, length.out = 450)
)

grid$BF_A <- BF_A(grid$ybar, grid$d)

# Calibrated rejection region of A.
grid$alpha_A <- ifelse(grid$BF_A <= K_A_star, 0.28, 0)

# Incoherent region under K_A = K_B = 1:
# RR_B(1) intersect NR_A(1).
grid$alpha_incoherent <- ifelse(
  abs(grid$d - d_B_star) >= c_B &
    grid$BF_A > K_A_uncalibrated,
  0.68,
  0
)

observed_point <- data.frame(ybar = ybar0, d = d0)

observed_label_x <- min(
  ybar0 + 0.055 * y_limit,
  0.88 * y_limit
)
observed_label_y <- d0 + 0.015 * d_limit

rejection_regions_plot <- ggplot() +

  # Calibrated RR_A.
  geom_raster(
    data = grid,
    aes(x = ybar, y = d, alpha = alpha_A),
    fill = "red3",
    interpolate = FALSE
  ) +

  # RR_B.
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = d_B_star - c_B,
    fill = "royalblue3", alpha = 0.28
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = d_B_star + c_B, ymax = Inf,
    fill = "royalblue3", alpha = 0.28
  ) +

  # Incoherent region for the uncalibrated thresholds.
  geom_raster(
    data = grid,
    aes(x = ybar, y = d, alpha = alpha_incoherent),
    fill = "grey55",
    interpolate = FALSE
  ) +
  scale_alpha_identity() +

  # Boundary BF_B = K_B.
  geom_hline(
    yintercept = d_B_boundary,
    linetype = "dashed",
    linewidth = 0.9,
    colour = "royalblue4"
  ) +

  # Level curves BF_A = K_A.
  geom_path(
    data = ellipses,
    aes(
      x = ybar,
      y = d,
      group = type,
      linetype = type
    ),
    colour = "grey25",
    linewidth = 1
  ) +
  scale_linetype_manual(
    name = expression("Level curves of " * A),
    breaks = c("uncalibrated", "calibrated"),
    labels = c(
      expression("uncalibrated " * K[A]),
      expression(K[A]^"*")
    ),
    values = c(
      uncalibrated = "dotted",
      calibrated = "solid"
    )
  ) +

  # Calibration boundary points.
  geom_point(
    data = tangency_points,
    aes(x = ybar, y = d),
    shape = 21,
    size = 4.2,
    stroke = 1,
    fill = "gold",
    colour = "black"
  ) +

  # Observed sample point.
  geom_point(
    data = observed_point,
    aes(x = ybar, y = d),
    shape = 23,
    size = 5,
    stroke = 1.2,
    fill = "white",
    colour = "black"
  ) +
  annotate(
    "text",
    x = observed_label_x,
    y = observed_label_y,
    label = "Observed data",
    hjust = 0,
    vjust = -0.15,
    size = 3.8
  ) +

  # Illustrative incoherent sample point under K_A = K_B = 1.
  annotate(
    "point",
    x = ybar_inc,
    y = d_inc,
    shape = 4,
    size = 3,
    stroke = 1.6,
    colour = "red3"
  ) +

  # Reference axes.
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.5) +

  # Region labels.
  annotate(
    "text",
    x = 0.80 * y_limit,
    y = 0.82 * d_limit,
    label = "RR[B]",
    parse = TRUE,
    colour = "royalblue4",
    fontface = "bold",
    size = 4.8
  ) +
  annotate(
    "text",
    x = 0.80 * y_limit,
    y = -0.82 * d_limit,
    label = "RR[B]",
    parse = TRUE,
    colour = "royalblue4",
    fontface = "bold",
    size = 4.8
  ) +
  annotate(
    "text",
    x = -0.58 * y_limit,
    y = 0.12 * d_limit,
    label = "RR[A]~'for'~K[A]^'*'",
    parse = TRUE,
    colour = "red4",
    fontface = "bold",
    size = 4.2
  ) +

  coord_cartesian(
    xlim = c(-y_limit, y_limit),
    ylim = c(-d_limit, d_limit),
    expand = FALSE
  ) +

  labs(
    x = expression(bar(Y)),
    y = expression(D == bar(Y)[FT] - bar(Y)[CBT]),
    subtitle = bquote(
      paste(
        K[B], " = ", .(formatC(K_B, digits = 3, format = "f")),
        "; ",
        K[A]^"*", " = ",
        .(formatC(K_A_star, digits = 3, format = "f"))
      )
    ),
    caption = expression(
      atop(
        paste(
          "Red: ", RR[A], " calibrated; Blue: ", RR[B],
          "; Purple: overlap."
        ),
        paste(
          "Gray: ", RR[B] %intersection% NR[A],
          " under the uncalibrated thresholds; red cross: illustrative incoherent sample point."
        )
      )
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(rejection_regions_plot)

ggsave(
  file.path(output_dir, "anorexia_rejection_regions.pdf"),
  rejection_regions_plot,
  width = 7,
  height = 5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(output_dir, "anorexia_rejection_regions.png"),
  rejection_regions_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ----------------------------------------------------------------------------
# 7. Interactive 3D Bayes factor surfaces
# ----------------------------------------------------------------------------
# This section reproduces the interactive visualisation linked in the paper.

if (!requireNamespace("plotly", quietly = TRUE)) {
  stop(
    "Package 'plotly' is required. Install it with install.packages('plotly')."
  )
}

if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'htmlwidgets' is required. ",
      "Install it with install.packages('htmlwidgets')."
    )
  )
}


create_interactive_bayes_factor_surfaces <- function(
    K_A_uncalibrated,
    K_A_calibrated,
    K_B,
    n_grid = 120
) {
  
  # --------------------------------------------------------------------------
  # Model quantities
  # --------------------------------------------------------------------------
  
  center_Y <- y_A_star
  center_D <- d_B_star
  
  # --------------------------------------------------------------------------
  # Plot ranges
  # --------------------------------------------------------------------------
  
  radius_D_uncalibrated <- sqrt(
    log(M_A / K_A_uncalibrated) / a_delta
  )
  
  radius_Y_uncalibrated <- sqrt(
    log(M_A / K_A_uncalibrated) / a_eta
  )
  
  radius_D_calibrated <- sqrt(
    log(M_A / K_A_calibrated) / a_delta
  )
  
  radius_Y_calibrated <- sqrt(
    log(M_A / K_A_calibrated) / a_eta
  )
  
  max_Y <- max(
    abs(center_Y - radius_Y_uncalibrated),
    abs(center_Y + radius_Y_uncalibrated),
    abs(center_Y - radius_Y_calibrated),
    abs(center_Y + radius_Y_calibrated),
    abs(ybar0),
    0.80
  )
  
  max_D <- max(
    abs(center_D - radius_D_uncalibrated),
    abs(center_D + radius_D_uncalibrated),
    abs(center_D - radius_D_calibrated),
    abs(center_D + radius_D_calibrated),
    abs(center_D - c_B),
    abs(center_D + c_B),
    abs(d0),
    0.90
  )
  
  Y_limit <- 1.20 * max_Y + 0.15
  D_limit <- 1.15 * max_D + 0.12
  
  # --------------------------------------------------------------------------
  # Grid
  # --------------------------------------------------------------------------
  
  Y_seq <- seq(
    -Y_limit,
    Y_limit,
    length.out = n_grid
  )
  
  D_seq <- seq(
    -D_limit,
    D_limit,
    length.out = n_grid
  )
  
  Z_A <- outer(
    D_seq,
    Y_seq,
    FUN = function(D, Ybar) {
      BF_A(Ybar, D)
    }
  )
  
  Z_B <- outer(
    D_seq,
    Y_seq,
    FUN = function(D, Ybar) {
      BF_B(D)
    }
  )
  
  # --------------------------------------------------------------------------
  # Threshold planes
  # --------------------------------------------------------------------------
  
  plane_A_uncalibrated <- matrix(
    K_A_uncalibrated,
    nrow = length(D_seq),
    ncol = length(Y_seq)
  )
  
  plane_A_calibrated <- matrix(
    K_A_calibrated,
    nrow = length(D_seq),
    ncol = length(Y_seq)
  )
  
  plane_B <- matrix(
    K_B,
    nrow = length(D_seq),
    ncol = length(Y_seq)
  )
  
  # --------------------------------------------------------------------------
  # Explicit level curves
  # --------------------------------------------------------------------------
  
  angle <- seq(
    0,
    2 * pi,
    length.out = 600
  )
  
  # BF_A = K_A uncalibrated
  level_A_uncalibrated_Y <- center_Y +
    radius_Y_uncalibrated * cos(angle)
  
  level_A_uncalibrated_D <- center_D +
    radius_D_uncalibrated * sin(angle)
  
  level_A_uncalibrated_Z <- rep(
    K_A_uncalibrated,
    length(angle)
  )
  
  # BF_A = K_A calibrated
  level_A_calibrated_Y <- center_Y +
    radius_Y_calibrated * cos(angle)
  
  level_A_calibrated_D <- center_D +
    radius_D_calibrated * sin(angle)
  
  level_A_calibrated_Z <- rep(
    K_A_calibrated,
    length(angle)
  )
  
  # BF_B = K_B
  level_B_lower_D <- rep(
    center_D - c_B,
    length(Y_seq)
  )
  
  level_B_upper_D <- rep(
    center_D + c_B,
    length(Y_seq)
  )
  
  level_B_Z <- rep(
    K_B,
    length(Y_seq)
  )
  
  # --------------------------------------------------------------------------
  # Observed and illustrative sample points
  # --------------------------------------------------------------------------
  
  observed_Y <- ybar0
  observed_D <- d0
  
  observed_BF_A <- BF_A(
    observed_Y,
    observed_D
  )
  
  observed_BF_B <- BF_B(
    observed_D
  )
  
  illustrative_Y <- ybar_inc
  illustrative_D <- d_inc
  
  illustrative_BF_A <- BF_A(
    illustrative_Y,
    illustrative_D
  )
  
  illustrative_BF_B <- BF_B(
    illustrative_D
  )
  
  # --------------------------------------------------------------------------
  # Calibration boundary points
  # --------------------------------------------------------------------------
  
  tangency_3d <- data.frame(
    Ybar = tangency_points$ybar,
    D = tangency_points$d
  )
  
  tangency_3d$BF_A <- BF_A(
    tangency_3d$Ybar,
    tangency_3d$D
  )
  
  tangency_3d$BF_B <- BF_B(
    tangency_3d$D
  )
  
  # ==========================================================================
  # Single Plotly object with two independent 3D scenes
  # ==========================================================================
  
  p <- plotly::plot_ly()
  
  # ==========================================================================
  # LEFT PANEL: BF_A(Ybar, D)
  # ==========================================================================
  
  p <- p |>
    
    # BF_A surface
    plotly::add_surface(
      x = Y_seq,
      y = D_seq,
      z = Z_A,
      scene = "scene",
      
      colorscale = list(
        c(0.00, "#FEE5D9"),
        c(0.25, "#FCAE91"),
        c(0.50, "#FB6A4A"),
        c(0.75, "#DE2D26"),
        c(1.00, "#7F0000")
      ),
      
      opacity = 0.92,
      showscale = TRUE,
      
      colorbar = list(
        title = list(
          text = "BF<sub>A</sub>"
        ),
        x = 0.465,
        y = 0.50,
        len = 0.58,
        thickness = 12
      ),
      
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "<b>BF<sub>A</sub></b><br>",
        "Y&#772; = %{x:.3f}<br>",
        "D = %{y:.3f}<br>",
        "BF<sub>A</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
    
    # Uncalibrated threshold plane for A
    plotly::add_surface(
      x = Y_seq,
      y = D_seq,
      z = plane_A_uncalibrated,
      scene = "scene",
      
      colorscale = list(
        c(0, "#555555"),
        c(1, "#555555")
      ),
      
      opacity = 0.10,
      showscale = FALSE,
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "Uncalibrated K<sub>A</sub> = ",
        formatC(
          K_A_uncalibrated,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # Calibrated threshold plane for A
    plotly::add_surface(
      x = Y_seq,
      y = D_seq,
      z = plane_A_calibrated,
      scene = "scene",
      
      colorscale = list(
        c(0, "#555555"),
        c(1, "#555555")
      ),
      
      opacity = 0.13,
      showscale = FALSE,
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "K<sub>A</sub><sup>min</sup> = ",
        formatC(
          K_A_calibrated,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # BF_A = K_A uncalibrated
    plotly::add_trace(
      x = level_A_uncalibrated_Y,
      y = level_A_uncalibrated_D,
      z = level_A_uncalibrated_Z,
      scene = "scene",
      
      type = "scatter3d",
      mode = "lines",
      
      line = list(
        color = "#FF8C00",
        width = 8
      ),
      
      name = "Uncalibrated K<sub>A</sub>",
      
      hovertemplate = paste0(
        "<b>BF<sub>A</sub> = K<sub>A</sub></b><br>",
        "uncalibrated K<sub>A</sub> = ",
        formatC(
          K_A_uncalibrated,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # BF_A = calibrated K_A
    plotly::add_trace(
      x = level_A_calibrated_Y,
      y = level_A_calibrated_D,
      z = level_A_calibrated_Z,
      scene = "scene",
      
      type = "scatter3d",
      mode = "lines",
      
      line = list(
        color = "#00A651",
        width = 9
      ),
      
      name = "K<sub>A</sub><sup>min</sup>",
      
      hovertemplate = paste0(
        "<b>BF<sub>A</sub> = K<sub>A</sub><sup>min</sup></b><br>",
        "K<sub>A</sub><sup>min</sup> = ",
        formatC(
          K_A_calibrated,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # Calibration boundary points
    plotly::add_trace(
      x = tangency_3d$Ybar,
      y = tangency_3d$D,
      z = tangency_3d$BF_A,
      scene = "scene",
      
      type = "scatter3d",
      mode = "markers",
      
      marker = list(
        size = 5,
        color = "#FFD600",
        line = list(
          color = "black",
          width = 1
        )
      ),
      
      name = "Tangency points",
      
      hovertemplate = paste0(
        "<b>Tangency point</b><br>",
        "Y&#772; = %{x:.3f}<br>",
        "D = %{y:.3f}<br>",
        "BF<sub>A</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
    
    # Observed sample point
    plotly::add_trace(
      x = observed_Y,
      y = observed_D,
      z = observed_BF_A,
      scene = "scene",
      
      type = "scatter3d",
      mode = "markers",
      
      marker = list(
        size = 5,
        color = "#00C853",
        line = list(
          color = "black",
          width = 1
        )
      ),
      
      name = "Observed data",
      
      hovertemplate = paste0(
        "<b>Observed data</b><br>",
        "y&#772; = %{x:.4f}<br>",
        "d = %{y:.4f}<br>",
        "BF<sub>A</sub> = %{z:.6f}",
        "<extra></extra>"
      )
    ) |>
    
    # Illustrative incoherent sample point
    plotly::add_trace(
      x = illustrative_Y,
      y = illustrative_D,
      z = illustrative_BF_A,
      scene = "scene",
      
      type = "scatter3d",
      mode = "markers",
      
      marker = list(
        size = 5,
        color = "#FF8C00",
        symbol = "x",
        line = list(
          color = "black",
          width = 1
        )
      ),
      
      name = "Illustrative sample point",
      
      hovertemplate = paste0(
        "<b>Illustrative sample point</b><br>",
        "y&#772; = %{x:.3f}<br>",
        "d = %{y:.3f}<br>",
        "BF<sub>A</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    )
  
  # ==========================================================================
  # RIGHT PANEL: BF_B(Ybar, D)
  # ==========================================================================
  
  p <- p |>
    
    # BF_B surface
    plotly::add_surface(
      x = Y_seq,
      y = D_seq,
      z = Z_B,
      scene = "scene2",
      
      colorscale = list(
        c(0.00, "#EFF3FF"),
        c(0.25, "#BDD7E7"),
        c(0.50, "#6BAED6"),
        c(0.75, "#3182BD"),
        c(1.00, "#08519C")
      ),
      
      opacity = 0.80,
      showscale = TRUE,
      
      colorbar = list(
        title = list(
          text = "BF<sub>B</sub>"
        ),
        x = 1.00,
        y = 0.50,
        len = 0.58,
        thickness = 12
      ),
      
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "<b>BF<sub>B</sub></b><br>",
        "Y&#772; = %{x:.3f}<br>",
        "D = %{y:.3f}<br>",
        "BF<sub>B</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
    
    # Threshold plane for B
    plotly::add_surface(
      x = Y_seq,
      y = D_seq,
      z = plane_B,
      scene = "scene2",
      
      colorscale = list(
        c(0, "#555555"),
        c(1, "#555555")
      ),
      
      opacity = 0.12,
      showscale = FALSE,
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "K<sub>B</sub> = ",
        formatC(
          K_B,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # Lower BF_B = K_B curve
    plotly::add_trace(
      x = Y_seq,
      y = level_B_lower_D,
      z = level_B_Z,
      scene = "scene2",
      
      type = "scatter3d",
      mode = "lines",
      
      line = list(
        color = "black",
        width = 5
      ),
      
      name = "K<sub>B</sub>",
      
      hovertemplate = paste0(
        "<b>BF<sub>B</sub> = K<sub>B</sub></b><br>",
        "K<sub>B</sub> = ",
        formatC(
          K_B,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # Upper BF_B = K_B curve
    plotly::add_trace(
      x = Y_seq,
      y = level_B_upper_D,
      z = level_B_Z,
      scene = "scene2",
      
      type = "scatter3d",
      mode = "lines",
      
      line = list(
        color = "black",
        width = 5
      ),
      
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "<b>BF<sub>B</sub> = K<sub>B</sub></b><br>",
        "K<sub>B</sub> = ",
        formatC(
          K_B,
          digits = 3,
          format = "f"
        ),
        "<extra></extra>"
      )
    ) |>
    
    # Boundary points associated with the calibration
    plotly::add_trace(
      x = tangency_3d$Ybar,
      y = tangency_3d$D,
      z = tangency_3d$BF_B,
      scene = "scene2",
      
      type = "scatter3d",
      mode = "markers",
      
      marker = list(
        size = 5,
        color = "#FFD600",
        line = list(
          color = "black",
          width = 1
        )
      ),
      
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "<b>Boundary point</b><br>",
        "Y&#772; = %{x:.3f}<br>",
        "D = %{y:.3f}<br>",
        "BF<sub>B</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
    
    # Observed sample point
    plotly::add_trace(
      x = observed_Y,
      y = observed_D,
      z = observed_BF_B,
      scene = "scene2",
      
      type = "scatter3d",
      mode = "markers",
      
      marker = list(
        size = 5,
        color = "#00C853",
        line = list(
          color = "black",
          width = 1
        )
      ),
      
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "<b>Observed data</b><br>",
        "y&#772; = %{x:.4f}<br>",
        "d = %{y:.4f}<br>",
        "BF<sub>B</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
    
    # Illustrative sample point
    plotly::add_trace(
      x = illustrative_Y,
      y = illustrative_D,
      z = illustrative_BF_B,
      scene = "scene2",
      
      type = "scatter3d",
      mode = "markers",
      
      marker = list(
        size = 5,
        color = "#FF8C00",
        symbol = "x",
        line = list(
          color = "black",
          width = 1
        )
      ),
      
      showlegend = FALSE,
      
      hovertemplate = paste0(
        "<b>Illustrative sample point</b><br>",
        "y&#772; = %{x:.3f}<br>",
        "d = %{y:.3f}<br>",
        "BF<sub>B</sub> = %{z:.4f}",
        "<extra></extra>"
      )
    )
  
  # ==========================================================================
  # Layout
  # ==========================================================================
  
  p <- p |>
    plotly::layout(
      
      # Left scene: BF_A
      scene = list(
        
        domain = list(
          x = c(0.00, 0.46),
          y = c(0.08, 1.00)
        ),
        
        xaxis = list(
          title = list(
            text = "\u0232"
          ),
          range = c(
            -Y_limit,
            Y_limit
          )
        ),
        
        yaxis = list(
          title = list(
            text = "D"
          ),
          range = c(
            -D_limit,
            D_limit
          )
        ),
        
        zaxis = list(
          title = list(
            text = "BF<sub>A</sub>"
          ),
          range = c(
            0,
            1.15 * M_A
          )
        ),
        
        aspectmode = "manual",
        
        aspectratio = list(
          x = 1.00,
          y = 1.00,
          z = 0.90
        ),
        
        camera = list(
          eye = list(
            x = 1.55,
            y = 1.45,
            z = 1.15
          )
        )
      ),
      
      # Right scene: BF_B
      scene2 = list(
        
        domain = list(
          x = c(0.54, 1.00),
          y = c(0.08, 1.00)
        ),
        
        xaxis = list(
          title = list(
            text = "\u0232"
          ),
          range = c(
            -Y_limit,
            Y_limit
          )
        ),
        
        yaxis = list(
          title = list(
            text = "D"
          ),
          range = c(
            -D_limit,
            D_limit
          )
        ),
        
        zaxis = list(
          title = list(
            text = "BF<sub>B</sub>"
          ),
          range = c(
            0,
            1.25 * M_B
          )
        ),
        
        aspectmode = "manual",
        
        aspectratio = list(
          x = 1.00,
          y = 1.00,
          z = 0.95
        ),
        
        camera = list(
          eye = list(
            x = 1.55,
            y = 1.45,
            z = 1.15
          )
        )
      ),
      
      # Panel labels
      annotations = list(
        
        list(
          text = "<b>BF<sub>A</sub>(Y&#772;, D)</b>",
          x = 0.22,
          y = 1.03,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(
            size = 16
          )
        ),
        
        list(
          text = "<b>BF<sub>B</sub>(Y&#772;, D)</b>",
          x = 0.78,
          y = 1.03,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(
            size = 16
          )
        )
      ),
      
      # Legend
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = -0.03,
        yanchor = "top"
      ),
      
      margin = list(
        l = 20,
        r = 70,
        b = 75,
        t = 80
      )
    )
  
  return(p)
}


# ----------------------------------------------------------------------------
# Create the interactive figure
# ----------------------------------------------------------------------------

interactive_bayes_factor_surfaces <-
  create_interactive_bayes_factor_surfaces(
    K_A_uncalibrated = K_A_uncalibrated,
    K_A_calibrated = K_A_star,
    K_B = K_B,
    n_grid = 120
  )

print(interactive_bayes_factor_surfaces)


# ----------------------------------------------------------------------------
# Save the interactive figure as a standalone HTML file
# ----------------------------------------------------------------------------

htmlwidgets::saveWidget(
  widget = interactive_bayes_factor_surfaces,
  file = file.path(
    output_dir,
    "interactive_bayes_factor_surfaces.html"
  ),
  selfcontained = TRUE
)
