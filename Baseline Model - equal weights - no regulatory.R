# === Load required libraries ===
library(ggplot2)
library(openxlsx)
library(dplyr)

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
  
  result <- optim(
    par = c(0.5, 0.5),
    fn = objective_function,
    method = "L-BFGS-B",
    lower = c(0, 0),
    upper = c(1, 1)
  )
  
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
  list(id = "DMSO", mu_level = "low",  lambda_level = "mid", alpha_level = "high"),
  list(id = "DMF", mu_level = "low",  lambda_level = "low", alpha_level = "low"),
  list(id = "NBP", mu_level = "high", lambda_level = "high", alpha_level = "high"),
  list(id = "Ar100", mu_level = "low",  lambda_level = "low", alpha_level = "low")
)

results <- lapply(selected_combinations, function(params) {
  optimize_Z(
    mu_level = params$mu_level,
    lambda_level = params$lambda_level,
    alpha_level = params$alpha_level,
    id = params$id
  )
})
results_df <- do.call(rbind, results)

# === Add example criteria (without EPA_Priority) ===
results_df$Market_Info  <- c(2.00, 0.75, 3.13, 1.42)   # Higher is worse (cost criterion)

# === Normalize helper function ===
normalize <- function(x) {
  if (length(x) > 0 && sd(x) != 0) {
    return(x / sqrt(sum(x^2)))
  } else {
    return(rep(0, length(x)))  # Return zero vector for degenerate cases
  }
}

# === Monte Carlo Simulation function ===
monte_carlo_topsis <- function(results_df, num_simulations = 1000) {
  simulation_results <- data.frame()
  euclidean_summary <- data.frame()  # To store Euclidean Distance summaries
  
  for (i in 1:num_simulations) {
    # Assign equal weights to each criterion
    weight_tm <- 0.5
    weight_market <- 0.5
    
    weights <- c(weight_tm, weight_market)
    criteria_benefit <- c(TRUE, FALSE)  # Tm is a benefit, Market_Info is a cost
    
    # Remove EPA_Priority from decision matrix
    decision_matrix <- results_df %>%
      select(Tm, Market_Info)  # Now only select the remaining criteria
    
    normalized_matrix <- as.data.frame(lapply(decision_matrix, normalize))
    weighted_matrix <- sweep(normalized_matrix, 2, weights, FUN = "*")
    
    ideal_solution <- ifelse(criteria_benefit, apply(weighted_matrix, 2, max), apply(weighted_matrix, 2, min))
    negative_ideal <- ifelse(criteria_benefit, apply(weighted_matrix, 2, min), apply(weighted_matrix, 2, max))
    
    distance_to_ideal <- apply(weighted_matrix, 1, function(row) sqrt(sum((row - ideal_solution)^2)))
    distance_to_negative <- apply(weighted_matrix, 1, function(row) sqrt(sum((row - negative_ideal)^2)))
    
    topsis_score <- distance_to_negative / (distance_to_negative + distance_to_ideal)
    
    # Collect Euclidean Distance data
    euclidean_summary <- rbind(euclidean_summary, data.frame(
      Simulation = i,
      Euclidean_Distance_to_Ideal = mean(distance_to_ideal),
      Euclidean_Distance_to_Negative = mean(distance_to_negative)
    ))
    
    temp_results <- data.frame(
      Simulation = i,
      Weight_Tm = weight_tm,
      id = results_df$id,
      TOPSIS_Score = topsis_score
    )
    
    simulation_results <- rbind(simulation_results, temp_results)
  }
  
  return(list(simulation_results = simulation_results, euclidean_summary = euclidean_summary))
}

# === Run the simulation ===
simulation_data <- monte_carlo_topsis(results_df, num_simulations = 500)
simulation_results <- simulation_data$simulation_results
euclidean_summary <- simulation_data$euclidean_summary

# === View results ===
head(simulation_results)
head(euclidean_summary)

# === Add code to calculate TOPSIS Summary Statistics ===
topsis_summary <- simulation_results %>%
  group_by(id) %>%
  summarise(
    mean_score = mean(TOPSIS_Score),
    sd_score = sd(TOPSIS_Score),
    min_score = min(TOPSIS_Score),
    max_score = max(TOPSIS_Score),
    range_score = max(TOPSIS_Score) - min(TOPSIS_Score)
  )

print("TOPSIS Summary Statistics:")
print(topsis_summary)

# === Calculate TOPSIS Rank ===
simulation_results$Rank <- ave(simulation_results$TOPSIS_Score, simulation_results$Simulation, FUN = function(x) rank(-x, ties.method = "min"))

# === TOPSIS Rank Summary Statistics ===
rank_summary <- simulation_results %>%
  group_by(id) %>%
  summarise(
    mean_rank = mean(Rank),
    sd_rank = sd(Rank),
    min_rank = min(Rank),
    max_rank = max(Rank),
    range_rank = max(Rank) - min(Rank)
  )

print("TOPSIS Rank Summary Statistics:")
print(rank_summary)

# === Euclidean Distance Summary across all simulations ===
euclidean_summary_across_simulations <- euclidean_summary %>%
  summarise(
    mean_distance_to_ideal = mean(Euclidean_Distance_to_Ideal),
    sd_distance_to_ideal = sd(Euclidean_Distance_to_Ideal),
    mean_distance_to_negative = mean(Euclidean_Distance_to_Negative),
    sd_distance_to_negative = sd(Euclidean_Distance_to_Negative)
  )

print("Euclidean Distance Summary Across Simulations:")
print(euclidean_summary_across_simulations)

# === Save to Excel ===
summary_list <- list(
  "TOPSIS_Score_Summary" = topsis_summary,
  "TOPSIS_Rank_Summary" = rank_summary, 
  "TOPSIS Simulation Results" = simulation_results,
  "Euclidean Distance Summary" = euclidean_summary,
  "Euclidean Distance Sim Summary" = euclidean_summary_across_simulations
)

# Append the summary statistics to the same Excel file
write.xlsx(summary_list, file = "C:/Users/btsha/Documents/Praxis/h2r1-baseline-ew.xlsx", rowNames = FALSE, append = TRUE)


# === Plotting (example) ===
custom_theme <- theme(
  axis.title.x = element_text(size = 24),
  axis.title.y = element_text(size = 24),
  axis.text = element_text(size = 20),
  title = element_text(size = 24),
  panel.grid = element_line(color = "grey", size = 2),
  legend.text = element_text(size = 20)
)

# 1. Plot TOPSIS Scores over Simulations
print(ggplot(simulation_results, aes(x = Simulation, y = TOPSIS_Score, color = id)) +
        geom_point(alpha = 0.4, size = 3) +
        geom_smooth(method = "loess", se = TRUE, linewidth = 2) +
        theme_minimal() +
        labs(
          title = "TOPSIS Scores over Simulations (BASELINE - equal weights)",
          x = "Simulation",
          y = "TOPSIS Score",
          color = "Alternative"
        ) + custom_theme)

# 2. Plot TOPSIS Rank Over Simulations
simulation_results$Rank <- ave(simulation_results$TOPSIS_Score, simulation_results$Simulation, FUN = function(x) rank(-x, ties.method = "min"))

print(ggplot(simulation_results, aes(x = Simulation, y = Rank, color = id)) +
        geom_line(alpha = 0.4, size = 3) +
        geom_smooth(method = "loess", se = FALSE, linewidth = 2) +
        scale_y_reverse(breaks = 1:n_distinct(simulation_results$id)) +  # Rank 1 is best
        theme_minimal() +
        labs(
          title = "TOPSIS Rank Over Simulations by Alternative (BASELINE - equal weights)",
          x = "Simulation Number",
          y = "TOPSIS Rank (1 = Best)",
          color = "Alternative"
        ) + custom_theme)

# 3. Rank Distribution Boxplot
print(ggplot(simulation_results, aes(x = id, y = Rank, fill = id)) +
        geom_boxplot(alpha = 0.6, size = 1.5) +
        scale_y_reverse() +
        theme_minimal() +
        labs(
          title = "Distribution of TOPSIS Ranks Across Simulations (BASELINE - equal weights)",
          x = "Alternative",
          y = "TOPSIS Rank (1 = Best)",
          color = "Alternative"
        ) + custom_theme)

# 4. Scatter Plot: TOPSIS Score vs Rank
print(ggplot(simulation_results, aes(x = TOPSIS_Score, y = Rank, color = id)) +
        geom_point(alpha = 0.6, size = 14) +
        geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
        scale_y_reverse() +  # Rank 1 is the best
        theme_minimal() +
        labs(
          title = "TOPSIS Score vs Rank for Each Alternative (BASELINE - equal weights)",
          x = "TOPSIS Score",
          y = "TOPSIS Rank (1 = Best)",
          color = "Alternative"
        ) + custom_theme)

# 5. Plot Tm vs Rank
simulation_results$Tm_Weight <- runif(nrow(simulation_results), min = 0.1, max = 1.0)

print(ggplot(simulation_results, aes(x = Tm_Weight, y = Rank, color = id)) +
        geom_point(aes(color = factor(id)), alpha = 0.4, size = 3) +
        scale_y_reverse() +  # Rank 1 is the best
        theme_minimal() +
        labs(
          title = "Rank over Tm weight (BASELINE - equal weights)",
          x = "Tm weight",
          y = "Rank",
          color = "Alternative"
        ) + custom_theme)

# 6. Plot Tm vs TOPSIS Score
print(ggplot(simulation_results, aes(x = Tm_Weight, y = TOPSIS_Score, color = id)) +
        geom_point(aes(color = factor(id)), alpha = 0.4, size = 3) +
        theme_minimal() +
        labs(
          title = "TOPSIS Scores over Tm weight (BASELINE - equal weights)",
          x = "Tm weight",
          y = "TOPSIS Score",
          color = "Alternative"
        ) + custom_theme)
