library(ggplot2)
library(dplyr)
library(patchwork)
library(sdmTMB)

set.seed(2)

## Simulate Grid & Spatial Mesh
n <- 100
grid <- expand.grid(
  x = seq(0, 1, length.out = n),
  y = seq(0, 1, length.out = n)
)

# Build a cutoff mesh for spatial field simulation
mesh <- make_mesh(grid, xy_cols = c("x", "y"), cutoff = 0.05)

## Simulate density fields
sim_density <- sdmTMB_simulate(
  formula = ~1,
  data = grid,
  mesh = mesh,
  family = gaussian(),
  range = 0.35,
  sigma_O = 1.0,
  seed = 103,
  B = 0,
  phi = 0.01
)

sim_utility <- sdmTMB_simulate(
  formula = ~1,
  data = grid,
  mesh = mesh,
  family = gaussian(),
  range = 0.25,
  sigma_O = 1.0,
  seed = 202,
  B = 0,
  phi = 0.01
)

# Scale and transform simulated field values
grid$density <- exp(sim_density$omega_s + 1)
grid$utility <- as.numeric(scale(sim_utility$omega_s))
grid$log_density <- as.numeric(scale(log(grid$density)))

# This is all simplified for visualization purposes from
# what we're actually fitting
grid$prob_conn <- exp(grid$log_density)
grid$prob_conn <- grid$prob_conn / sum(grid$prob_conn)

grid$prob_dc <- exp(grid$utility)
grid$prob_dc <- grid$prob_dc / sum(grid$prob_dc)

nstations <- 500

conn_pts <- grid[sample(nrow(grid), nstations, replace = TRUE, prob = grid$prob_conn), ]

dc_pts <- grid[sample(nrow(grid), nstations, replace = TRUE, prob = grid$prob_dc), ]


joint_theme <-
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 11,
      margin = margin(b = 3)
    ),
    plot.margin = margin(1, 1, 1, 1)
  )

lim <- max(abs(c(grid$utility, grid$log_density)))

# Top Panel
p_density <-
  ggplot(grid, aes(x, y)) +
  geom_raster(aes(fill = density)) +
  coord_equal(expand = FALSE) +
  scale_fill_viridis_c(begin = 0.8, end = 0.0) +
  labs(title = "True density") +
  joint_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

# Joint Model Panels
p_conn_driver <-
  ggplot(grid, aes(x, y)) +
  geom_raster(aes(fill = log_density)) +
  coord_equal(expand = FALSE) +
  scale_fill_gradient2(
    low = "#2166ac", # more saturated blue and red
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-lim, lim)
  ) +
  labs(title = "Sampling field: joint") +
  joint_theme

p_conn_samples <-
  ggplot(grid, aes(x, y)) +
  geom_raster(aes(fill = density), alpha = 0.7) +
  geom_point(
    data = conn_pts,
    colour = "black",
    alpha = 0.20,
    # size = 0.20
  ) +
  coord_equal(expand = FALSE) +
  scale_fill_viridis_c(begin = 0.8, end = 0.0) +
  labs(title = "Realized sampling: joint") +
  joint_theme

# Discrete-Choice Model Panels
p_dc_driver <-
  ggplot(grid, aes(x, y)) +
  geom_raster(aes(fill = utility)) +
  coord_equal(expand = FALSE) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-lim, lim)
  ) +
  labs(title = "Sampling field: DC") +
  joint_theme

p_dc_samples <-
  ggplot(grid, aes(x, y)) +
  geom_raster(aes(fill = density), alpha = 0.7) +
  geom_point(
    data = dc_pts,
    colour = "black",
    alpha = 0.20,
    # size = 0.20
  ) +
  coord_equal(expand = FALSE) +
  scale_fill_viridis_c(begin = 0.8, end = 0.0) +
  labs(title = "Realized sampling: DC") +
  joint_theme


col_conn <- (p_conn_driver / p_conn_samples) +
  plot_annotation(title = "Joint areal model") &
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12, margin = margin(b = 4)))

col_dc <- (p_dc_driver / p_dc_samples) +
  plot_annotation(title = "Discrete-choice model") &
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12, margin = margin(b = 4)))

top_row <- (plot_spacer() | p_density | plot_spacer()) +
  plot_layout(widths = c(1, 2, 1))

bottom_row <- col_conn | col_dc

figure1 <- (top_row / bottom_row) +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 13)
  )

figure1

ggsave(plot = figure1, filename = "simulation_paper/figures/Figure_01.png", width = 7, height = 7)
