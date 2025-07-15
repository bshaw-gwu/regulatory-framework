# === Load required libraries ===
library(ggplot2)
library(openxlsx)
library(dplyr)
library(triangle)

# === Optimization function ===
optimize_Z <- function(mu_level, lambda_level, alpha_level, id, Qc = 1, Tt = 1) {
  mu_vals     <- list(low = 0.1, mid = 0.5, high = 0.9)
  lambda_vals <- list(low = 0.25, mid = 0.6, high = 0.8)
  alpha_vals  <- list(low = 0.1, mid = 0.5, high = 0.9)
  
  mu     <- mu_vals[[mu_level]]
  lambda <- lambda_vals[[lambda_level]]
  alpha  <- alpha_vals[[alpha_level]]
  
  objective_function <- function(x) {
    Td <- x[1]
    m  <- x[2]
    if (m <= 0) return(1e6)
    Ts <- exp(log(m) / lambda)
    Ti <- alpha * (1 - m) * Td
    Tm <- Tt - Td - Ts - Ti
    if (Ts < 0 || Ts > 1 || Td < 0 || Td > 1 || m < 0 || m > 1) return(1e6)
    R1 <- Td * exp(-mu / (alpha + m))
    R2 <- Td * exp(mu / (alpha + m)) + Qc
    R3 <- 1 - (m^(1 / lambda)) - Td - alpha * (1 - m) * Td
    Z <- (R1 / R2) * R3
    penalty <- 1000 * (abs(Tt - Tm - Ti - Ts - Td))^2
    return(-Z + penalty)
  }
  
  result <- optim(par = c(0.5, 0.5), fn = objective_function, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))
  Td_opt <- result$par[1]
  m_opt  <- result$par[2]
  Ts_opt <- exp(log(m_opt) / lambda)
  Ti_opt <- alpha * (1 - m_opt) * Td_opt
  Tm_opt <- Tt - Td_opt - Ts_opt - Ti_opt
  Z_opt  <- -result$value
  
  return(data.frame(
    id = id,
    mu_level, lambda_level, alpha_level,
    mu, lambda, alpha,
    Td = Td_opt, m = m_opt, Ts = Ts_opt, Ti = Ti_opt, Tm = Tm_opt, Z = Z_opt
  ))
}

# === Generate optimization data ===
selected_combinations <- list(
  list(id = "DMSO", mu_level = "low", lambda_level = "mid", alpha_level = "high"),
  list(id = "DMF", mu_level = "low", lambda_level = "low", alpha_level = "low"),
  list(id = "NBP", mu_level = "high", lambda_level = "high", alpha_level = "high"),
  list(id = "Ar100", mu_level = "low", lambda_level = "low", alpha_level = "low")
)

results <- lapply(selected_combinations, function(params) {
  optimize_Z(params$mu_level, params$lambda_level, params$alpha_level, params$id)
})
results_df <- do.call(rbind, results)

# === Add criteria ===
results_df$EPA_Risk     <- c(1, 2, 1, 1)
results_df$EPA_Priority <- c(1, 2, 1, 1)
results_df$Market_Info  <- c(2.00, 0.75, 3.13, 1.42)

# === Normalize helper ===
normalize <- function(x) {
  if (length(x) > 0 && sd(x) != 0) return(x / sqrt(sum(x^2)))
  else return(rep(0, length(x)))
}

# === Monte Carlo Simulation ===
monte_carlo_topsis <- function(results_df, num_simulations = 1000, weight_range = c(0.1, 1.0)) {
  simulation_results <- data.frame()
  
  for (i in 1:num_simulations) {
    weights <- runif(4, min = weight_range[1], max = weight_range[2])
    weights <- weights / sum(weights)
    
    decision_matrix <- results_df %>%
      select(id, Tm, EPA_Risk, EPA_Priority, Market_Info)
    
    decision_matrix$EPA_Risk <- mapply(function(risk, priority) {
      if (priority == 2) triangle::rtriangle(1, a = 1, b = 3, c = 2) else risk
    }, decision_matrix$EPA_Risk, decision_matrix$EPA_Priority)
    
    normalized_matrix <- as.data.frame(lapply(decision_matrix[, -1], normalize))
    weighted_matrix <- sweep(normalized_matrix, 2, weights, FUN = "*")
    
    criteria_benefit <- c(TRUE, FALSE, FALSE, FALSE)
    ideal <- ifelse(criteria_benefit, apply(weighted_matrix, 2, max), apply(weighted_matrix, 2, min))
    nadir <- ifelse(criteria_benefit, apply(weighted_matrix, 2, min), apply(weighted_matrix, 2, max))
    
    dist_ideal <- apply(weighted_matrix, 1, function(row) sqrt(sum((row - ideal)^2)))
    dist_nadir <- apply(weighted_matrix, 1, function(row) sqrt(sum((row - nadir)^2)))
    topsis_score <- dist_nadir / (dist_ideal + dist_nadir)
    
    temp <- data.frame(
      Simulation = i,
      id = decision_matrix$id,
      TOPSIS_Score = topsis_score,
      distance_to_ideal = dist_ideal,
      distance_to_negative_ideal = dist_nadir
    )
    
    simulation_results <- rbind(simulation_results, temp)
  }
  
  return(simulation_results)
}

# === Run Simulation ===
simulation_results <- monte_carlo_topsis(results_df, num_simulations = 500)

# === Add Rank and Tm_Weight ===
simulation_results$Rank <- ave(simulation_results$TOPSIS_Score, simulation_results$Simulation, 
                               FUN = function(x) rank(-x, ties.method = "min"))
simulation_results$Tm_Weight <- runif(nrow(simulation_results), min = 0.1, max = 1.0)

# === Summaries ===
topsis_summary <- simulation_results %>%
  group_by(id) %>%
  summarise(mean_score = mean(TOPSIS_Score), sd_score = sd(TOPSIS_Score),
            min_score = min(TOPSIS_Score), max_score = max(TOPSIS_Score), .groups = 'drop')

rank_summary <- simulation_results %>%
  group_by(id) %>%
  summarise(mean_rank = mean(Rank), sd_rank = sd(Rank),
            min_rank = min(Rank), max_rank = max(Rank), .groups = 'drop')

euclidean_summary <- simulation_results %>%
  group_by(Simulation, id) %>%
  summarise(distance_to_ideal = mean(distance_to_ideal),
            distance_to_negative_ideal = mean(distance_to_negative_ideal), .groups = 'drop')

euclidean_summary_across_simulations <- simulation_results %>%
  group_by(id) %>%
  summarise(mean_distance_to_ideal = mean(distance_to_ideal),
            mean_distance_to_negative_ideal = mean(distance_to_negative_ideal),
            sd_distance_to_ideal = sd(distance_to_ideal),
            sd_distance_to_negative_ideal = sd(distance_to_negative_ideal),
            min_distance_to_ideal = min(distance_to_ideal),
            max_distance_to_ideal = max(distance_to_ideal),
            min_distance_to_negative_ideal = min(distance_to_negative_ideal),
            max_distance_to_negative_ideal = max(distance_to_negative_ideal),
            .groups = 'drop')

# === Export to Excel ===
summary_list <- list(
  "TOPSIS_Score_Summary" = topsis_summary,
  "TOPSIS_Rank_Summary" = rank_summary,
  "TOPSIS Simulation Results" = simulation_results,
  "Euclidean Distance Summary" = euclidean_summary,
  "Euclidean Distance Sim Summary" = euclidean_summary_across_simulations
)
write.xlsx(summary_list, file = "C:/Users/btsha/Documents/Praxis/h3r3_triangle.xlsx", rowNames = FALSE)

# === Plotting ===
custom_theme <- theme(
  axis.title.x = element_text(size = 24),
  axis.title.y = element_text(size = 24),
  axis.text = element_text(size = 20),
  title = element_text(size = 24),
  panel.grid = element_line(color = "grey", size = 2),
  legend.text = element_text(size = 20)
)

# 1. TOPSIS Scores over Simulations
print(ggplot(simulation_results, aes(x = Simulation, y = TOPSIS_Score, color = id)) +
  geom_point(alpha = 0.4, size = 3) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 2) +
  theme_minimal() + custom_theme +
  labs(title = "TOPSIS Scores over Simulations - FMT", x = "Simulation", y = "TOPSIS Score", color = "Alternative"))

# 2. Rank Over Simulations
print(ggplot(simulation_results, aes(x = Simulation, y = Rank, color = id)) +
  geom_line(alpha = 0.4, size = 3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 2) +
  scale_y_reverse() +
  theme_minimal() + custom_theme +
  labs(title = "TOPSIS Rank Over Simulations - FMT", x = "Simulation Number", y = "Rank (1 = Best)", color = "Alternative"))

# 3. Rank Distribution Boxplot
print(ggplot(simulation_results, aes(x = id, y = Rank, fill = id)) +
  geom_boxplot(alpha = 0.6, size = 1.5) +
  scale_y_reverse() +
  theme_minimal() + custom_theme +
  labs(title = "Distribution of Ranks - FMT", x = "Alternative", y = "TOPSIS Rank", color = "Alternative"))

# 4. TOPSIS Score vs Rank
print(ggplot(simulation_results, aes(x = TOPSIS_Score, y = Rank, color = id)) +
  geom_point(alpha = 0.6, size = 4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  scale_y_reverse() +
  theme_minimal() + custom_theme +
  labs(title = "TOPSIS Score vs Rank - FMT", x = "TOPSIS Score", y = "Rank", color = "Alternative"))

# 5. Rank vs Tm Weight
print(ggplot(simulation_results, aes(x = Tm_Weight, y = Rank, color = id)) +
  geom_point(alpha = 0.4, size = 3) +
  scale_y_reverse() +
  theme_minimal() + custom_theme +
  labs(title = "Rank over Tm Weight - FMT", x = "Tm Weight", y = "Rank", color = "Alternative"))

# 6. Score vs Tm Weight
print(ggplot(simulation_results, aes(x = Tm_Weight, y = TOPSIS_Score, color = id)) +
  geom_point(alpha = 0.4, size = 3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 2) +
  theme_minimal() + custom_theme +
  labs(title = "TOPSIS Scores over Tm Weight - FMT", x = "Tm Weight", y = "TOPSIS Score", color = "Alternative"))

# 7. Rank Distribution Boxplot
print(ggplot(simulation_results, aes(x = id, y = TOPSIS_Score, fill = id)) +
        geom_boxplot(alpha = 0.6, size = 1.5) +
        theme_minimal() + custom_theme +
        labs(title = "Distribution of Ranks - FMT", x = "Alternative", y = "TOPSIS Score", color = "Alternative"))
