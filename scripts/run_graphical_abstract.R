library(ggplot2)


set.seed(42)
# --- Condition 1 cells ---
n <- 150
theta1 <- runif(n, 0, 2 * pi)
r1 <- sqrt(runif(n, 0, 1)) * 0.5
# x1 <- r1 * cos(theta1) - 0.6
x1 <- r1 * cos(theta1) # mixed
y1 <- r1 * sin(theta1) * 0.8

# --- Condition 2 cells ---
theta2 <- runif(n, 0, 2 * pi)
r2 <- sqrt(runif(n, 0, 1)) * 0.5
# x2 <- r2 * cos(theta2) + 0.6
x2 <- r2 * cos(theta2) # mixed
y2 <- r2 * sin(theta2) * 0.8

df <- data.frame(
  x = c(x1, x2),
  y = c(y1, y2),
  condition = rep(c("Condition 1", "Condition 2"), each = n)
)

p = ggplot(df, aes(x, y, color = condition)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("Condition 1" = "blue",
                                "Condition 2" = "red")) +
  theme_void() +
  theme(legend.position = "bottom") +
  coord_equal()
p


# --- Assign module membership ---
# Condition 1: localized modules
# Module A = cells near (-0.75, 0.15), Module B = cells near (-0.40, -0.20)
dist_A1 <- sqrt((x1 - (-0.75))^2 + (y1 - 0.15)^2)
dist_B1 <- sqrt((x1 - (-0.40))^2 + (y1 - (-0.20))^2)

module1 <- rep("None", n)
module1[dist_A1 < 0.18] <- "Module A"
module1[dist_B1 < 0.18] <- "Module B"

# Condition 2: dispersed — randomly sprinkle a few across the whole cloud
module2 <- rep("None", n)
set.seed(99)
module2[sample(which(module2 == "None"), 12)] <- "Module A"
module2[sample(which(module2 == "None"), 12)] <- "Module B"

df <- data.frame(
  x = c(x1, x2),
  y = c(y1, y2),
  condition = rep(c("Condition 1", "Condition 2"), each = n),
  module = c(module1, module2)
)

# --- Alpha: bright in Cond1, very faint in Cond2 ---
df$alpha <- ifelse(df$module == "None", 0.4,
                   ifelse(df$condition == "Condition 1", 0.9, 0.25))

# --- Size: module cells slightly larger ---
df$pt_size <- ifelse(df$module == "None", 1.2, 2.2)

# --- Plot order: None first so modules draw on top ---
df <- df[order(df$module == "None", decreasing = TRUE), ]

p <- ggplot(df, aes(x, y)) +
  geom_point(aes(color = module, alpha = alpha), size = 1.5) +
  scale_color_manual(
    values = c("None"     = "#D0D0D0",
               "Module A" = "#FF9800",
               "Module B" = "#4CAF50"),
    breaks = c("Module A", "Module B", "None"),
    labels = c("Module A", "Module B", "No expression")
  ) +
  scale_alpha_identity() +
  facet_wrap(~ condition) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank()
  ) +
  coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

p



# correlation simulation
library(ggplot2)
set.seed(42)
n <- 200

# Condition A (blue) - highly correlated
x_a <- rnorm(n)
y_a <- x_a + rnorm(n, sd = 0.3)

# Condition B (red) - no correlation
x_b <- rnorm(n)
y_b <- rnorm(n)

df <- data.frame(
  x = c(x_a, x_b),
  y = c(y_a, y_b),
  condition = rep(c("Condition A", "Condition B"), each = n)
)

# Compute R labels per condition
r_labels <- aggregate(cbind(x, y) ~ condition, df, function(v) v)
labels_df <- data.frame(
  condition = c("Condition A", "Condition B"),
  label = c(
    paste0("R = ", round(cor(x_a, y_a), 2)),
    paste0("R = ", round(cor(x_b, y_b), 2))
  )
)

ggplot(df, aes(x, y, color = condition)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text(data = labels_df, aes(x = Inf, y = -Inf, label = label),
            hjust = 1.1, vjust = -0.5, size = 6, show.legend = FALSE) +
  scale_color_manual(values = c("Condition A" = "blue", "Condition B" = "red")) +
  facet_wrap(~condition) +
  labs(x = "Gene X expression", y = "Gene Y expression") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# boxplot
library(ggplot2)
set.seed(42)
n <- 100

df <- data.frame(
  value = c(rnorm(n, mean = 8, sd = 1), rnorm(n, mean = 2, sd = 1)),
  condition = rep(c("Condition 1", "Condition 2"), each = n)
)

ggplot(df, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  scale_fill_manual(values = c("Condition 1" = "blue", "Condition 2" = "red")) +
  labs(x = NULL, y = "Expression") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))
