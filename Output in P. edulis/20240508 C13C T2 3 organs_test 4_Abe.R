### Step 1. Input experimental data
{
data <- read.csv("C:/Users/Sorbus/Downloads/work log on task 2/T2exp.csv", header = TRUE, sep = ',')

# Print the first few rows of the data
print(head(data))

LmT <- 0.327209279
BmT <- 0.976583259
CmT <- 17.8874446
BsmT <- 2.414524679
C1L <- 2.313385459
k1L <- 1.489827803
C2L <- 0.132945945
C1B <- 0.221346507
k1B <- 0.46786613
C2B <- 0.024940457
#C1C <- 1.001632481 here we did not sample after 13C pulse labelling due to concern treatment effect of digging holes in the culm until D14
#k1C <- 0.415758923
C2C <- 0.000967317


Total13C = 8.00462335 + 2.405197254
}

### Step 2. operate the equilibrium approach
{
## Leaf
{
# Define function to calculate RMSE
leaf_rmse <- function(params) {
  kRL <- params[1]
  kLimsyn <- params[2]
  aAcRL <- params[3]
  lambda_AcRL <- params[4]
  
  t_values <- seq(0, 30)  # t from 0 to 30
  
  # Calculate values from leaves
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QLnetchange_theoretical <- numeric(length(t_values))
  QLnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QLnetchange_theoretical[i] <- QL_theoretical[i] - QL_theoretical[i-1]
  }
  QLim_theoretical <- (C2L * 10 * LmT) * (1 - exp(-kLimsyn * t_values))
  QLim_theoretical <- ifelse(QLim_theoretical > 0, QLim_theoretical, 0)
  QLm_theoretical <- QL_theoretical - QLim_theoretical
  QLm_theoretical <- ifelse(QLm_theoretical > 0, QLm_theoretical, 0)
  
  RL_experimental <- data$RL_experimental
  AcRL_experimental <- data$AcRL_experimental
  
  RL_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    RL_theoretical[i] <- QLm_theoretical[i-1] * kRL
  }
  
  AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))
  
  # Constraint for leaf
  QL_constraint <- QLim_theoretical + QLm_theoretical
  
  rmse_value <- sqrt(mean((QL_constraint - QL_theoretical)^2)) + sqrt(mean((AcRL_experimental - AcRL_theoretical)^2)) + sqrt(mean((RL_experimental - RL_theoretical)^2))
  
  return(rmse_value)
}

# Initial guess
initial_guess_leaf <- c(0.01, 0.01, 0.01, 0.01)

# Minimize RMSE using 'optim' function
result_leaf <- optim(par = initial_guess_leaf, fn = leaf_rmse, method = "L-BFGS-B", control = list(maxit = 10000))

# Optimal values
kRL <- result_leaf$par[1]
kLimsyn <- result_leaf$par[2]
aAcRL <- result_leaf$par[3]
lambda_AcRL <- result_leaf$par[4]

# Print optimal values
cat("Optimal values:\n")
cat("kRL:", kRL, "\n")
cat("kLimsyn:", kLimsyn, "\n")
cat("aAcRL:", aAcRL, "\n")
cat("lambda_AcRL:", lambda_AcRL, "\n")
}

## Branch
{
# Define function to calculate RMSE 
branch_rmse <- function(params) {
  kRB <- params[1]
  kBimsyn <- params[2]
  aAcRB <- params[3]
  lambda_AcRB <- params[4]
  
  t_values <- seq(0, 30)  # t from 0 to 30
  
  # Calculate values from branch
  QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  QBnetchange_theoretical <- numeric(length(t_values))
  QBnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QBnetchange_theoretical[i] <- QB_theoretical[i] - QB_theoretical[i-1]
  }
  QBim_theoretical <- (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))
  QBim_theoretical <- ifelse(QBim_theoretical > 0, QBim_theoretical, 0)
  QBm_theoretical <- QB_theoretical - QBim_theoretical
  QBm_theoretical <- ifelse(QBm_theoretical > 0, QBm_theoretical, 0)
  
  RB_experimental <- data$RB_experimental
  AcRB_experimental <- data$AcRB_experimental
  
  RB_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    RB_theoretical[i] <- QBm_theoretical[i-1] * kRB
  }
  
  AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))
  
  # Constraint for branch
  QB_constraint <- QBim_theoretical + QBm_theoretical
  
  rmse_value <- sqrt(mean((QB_constraint - QB_theoretical)^2)) + sqrt(mean((AcRB_experimental - AcRB_theoretical)^2)) + sqrt(mean((RB_experimental - RB_theoretical)^2))
  
  return(rmse_value)
}

# Initial guess
initial_guess_branch <- c(0.01, 0.01, 0.01, 0.01)

# Minimize RMSE
result_branch <- optim(par = initial_guess_branch, fn = branch_rmse, method = "L-BFGS-B", control = list(maxit = 10000))

# Optimal values
kRB <- result_branch$par[1]
kBimsyn <- result_branch$par[2]
aAcRB <- result_branch$par[3]
lambda_AcRB <- result_branch$par[4]

# Print optimal values
cat("Optimal values:\n")
cat("kRB:", kRB, "\n")
cat("kBimsyn:", kBimsyn, "\n")
cat("aAcRB:", aAcRB, "\n")
cat("lambda_AcRB:", lambda_AcRB, "\n")
}

## Culm 
{
# Define function to calculate RMSE
culm_rmse <- function(params) {
  k1C <- params[1]
  C1C <- params[2]
  kRC <- params[3]
  kCimsyn <- params[4]
  aAcRC <- params[5]
  lambda_AcRC <- params[6]
  kLBCExsyn <- params[7]
  #kBCsyn <- params[8]
  #kExsym <- params[9]
  #we initially found mature culms transfer carbon in a quick way, so we assumed that mobile pool as all in one k (kLBCExsyn)
  
  t_values <- seq(0, 30)  # t from 0 to 30

  # Relatives from leaves
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QLnetchange_theoretical <- numeric(length(t_values))
  QLnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QLnetchange_theoretical[i] <- QL_theoretical[i] - QL_theoretical[i-1]
  }
  QLim_theoretical <- (C2L * 10 * LmT) * (1 - exp(-kLimsyn * t_values))
  QLim_theoretical <- ifelse(QLim_theoretical > 0, QLim_theoretical, 0)
  QLm_theoretical <- QL_theoretical - QLim_theoretical
  QLm_theoretical <- ifelse(QLm_theoretical > 0, QLm_theoretical, 0)
  
  RL_experimental <- data$RL_experimental
  AcRL_experimental <- data$AcRL_experimental
  
  RL_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    RL_theoretical[i] <- QLm_theoretical[i-1] * kRL
  }
  
  AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))
  
    
  # Relatives from branch
  QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  QBnetchange_theoretical <- numeric(length(t_values))
  QBnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QBnetchange_theoretical[i] <- QB_theoretical[i] - QB_theoretical[i-1]
  }
  QBim_theoretical <- (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))
  QBim_theoretical <- ifelse(QBim_theoretical > 0, QBim_theoretical, 0)
  QBm_theoretical <- QB_theoretical - QBim_theoretical
  QBm_theoretical <- ifelse(QBm_theoretical > 0, QBm_theoretical, 0)
  
  # Relatives from culm
  RC_experimental <- data$RC_experimental
  AcRC_experimental <- data$AcRC_experimental
  
  QCm_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    QCm_theoretical[i] <- kLBCExsyn * QBm_theoretical[i-1] - kRC * QCm_theoretical[i-1] - kLBCExsyn * QCm_theoretical[i-1] - kCimsyn * QCm_theoretical[i-1]
  }
  QCm_theoretical <- ifelse(QCm_theoretical > 0, QCm_theoretical, 0)

  #RC_theoretical <- numeric(length(t_values))
  #for (i in 2:length(t_values)) {
  #RC_theoretical <- kRC * QCm_theoretical[i-1]
  #}
  #RC_theoretical <- ifelse(RC_theoretical > 0, RC_theoretical, 0)
  #RC is irregular, no theoretical but experimental only
  
  AcRC_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    AcRC_theoretical[i] <- RC_experimental[i-1] + RC_experimental[i]
  }
  AcRC_theoretical <- aAcRC * (1 - exp(-lambda_AcRC * t_values))       
  
    
  QCim_theoretical <- (C2C * 10 * CmT) * (1 - exp(-kCimsyn * t_values))
  QCim_theoretical <- ifelse(QCim_theoretical > 0, QCim_theoretical, 0)
  
  QC_theoretical <- QCm_theoretical + QCim_theoretical
  
  QCnetchange_theoretical <- numeric(length(t_values))
  QCnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QCnetchange_theoretical[i] <- (CmT * (C1C * exp(-k1C * (t_values[i])) + C2C) * 10) - (CmT * (C1C * exp(-k1C * t_values[i-1]) + C2C) * 10)
    if (QCnetchange_theoretical[i] < 0) {
      QCnetchange_theoretical[i] <- (C2C * 10 * CmT) * (1 - exp(-k1C * t_values[i]))
    }
  }
  
  Ex_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    Ex_theoretical[i] <- kLBCExsyn * QCm_theoretical[i]
  }
  
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QLnetchange_theoretical <- numeric(length(t_values))
  QLnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QLnetchange_theoretical[i] <- QL_theoretical[i] - QL_theoretical[i-1]
  }
  
  
  RL_experimental <- data$RL_experimental
  RB_experimental <- data$RB_experimental
  RC_experimental <- data$RC_experimental
  
  if (any(QCnetchange_theoretical > 0)) {
    netdecrease_constraint <- -QLnetchange_theoretical - QBnetchange_theoretical
  } else {
    netdecrease_constraint <- -QLnetchange_theoretical - QBnetchange_theoretical - QCnetchange_theoretical
  }
  
  if (any(QCnetchange_theoretical > 0)) {
    netincrease_constraint <- QCnetchange_theoretical + Ex_theoretical + RL_experimental + RB_experimental + RC_experimental
  } else {
    netincrease_constraint <- Ex_theoretical + RL_experimental + RB_experimental + RC_experimental
  }
  
  rmse_value <- sqrt(mean((AcRC_experimental - AcRC_theoretical)^2))
  
  return(rmse_value)
}

initial_guess_culm <- c(0.415758923, 1.001632481, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

# Constraints
constraints <- list(
  function(x) x[1], 
  function(x) x[2],  
  function(x) x[3], 
  function(x) x[4], 
  function(x) x[5], 
  function(x) x[6], 
  function(x) x[7],  
  function(x) x[8],
  function(x) x[9]
)

# Minimize RMSE
result_culm <- optim(par = initial_guess_culm, fn = culm_rmse, method = "L-BFGS-B", control = list(maxit = 10000))

# Optimal values
k1C <- result_culm$par[1]
C1C <- result_culm$par[2]
kRC <- result_culm$par[3]
kCimsyn <- result_culm$par[4]
aAcRC <- result_culm$par[5]
lambda_AcRC <- result_culm$par[6]
kLBCExsyn <- result_culm$par[7]
#kBCsyn <- result_culm$par[8]
#kExsym <- result_culm$par[9]
#note: here we use the real flux amount and it was irregular and cannot get any optimal value 


# Print optimal values
cat("Optimal values:\n")
cat("k1C:", k1C, "\n")
cat("C1C:", C1C, "\n")
cat("kRC:", kRC, "\n")
cat("kCimsyn:", kCimsyn, "\n")
cat("aAcRC:", aAcRC, "\n")
cat("lambda_AcRC:", lambda_AcRC, "\n")
cat("kLBCExsyn:", kLBCExsyn, "\n")
}

## Export
{
#install.packages("pracma")
#library("pracma")
#install pracma if R require it

# Define function to calculate RMSE
exp_rmse <- function(params) {
  aAcEx <- params[1]
  lambda_AcEx <- params[2]
  kLBCExsyn <- params[3]
  
  t_values <- seq(0, 30)  # t from 0 to 30

  # Relatives from leaves
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QLnetchange_theoretical <- numeric(length(t_values))
  QLnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QLnetchange_theoretical[i] <- QL_theoretical[i] - QL_theoretical[i-1]
  }
  QLim_theoretical <- (C2L * 10 * LmT) * (1 - exp(-kLimsyn * t_values))
  QLim_theoretical <- ifelse(QLim_theoretical > 0, QLim_theoretical, 0)
  QLm_theoretical <- QL_theoretical - QLim_theoretical
  QLm_theoretical <- ifelse(QLm_theoretical > 0, QLm_theoretical, 0)
  
  RL_experimental <- data$RL_experimental
  AcRL_experimental <- data$AcRL_experimental
  
  RL_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    RL_theoretical[i] <- QLm_theoretical[i-1] * kRL
  }
  
  AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))
  
    
  # Relatives from branch
  QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  QBnetchange_theoretical <- c(0, diff(QB_theoretical))
  QBim_theoretical <- (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))
  QBim_theoretical <- ifelse(QBim_theoretical > 0, QBim_theoretical, 0)
  QBm_theoretical <- QB_theoretical - QBim_theoretical
  QBm_theoretical <- ifelse(QBm_theoretical > 0, QBm_theoretical, 0)
  
  RB_experimental <- data$RB_experimental
  AcRB_experimental <- data$AcRB_experimental
  
  RB_theoretical <- c(0, diff(QBm_theoretical) * kRB)
  AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))
  
  # Calculate values from culm
  
  RC_experimental <- data$RC_experimental
  AcRC_experimental <- data$AcRC_experimental
  
  QCm_theoretical <- c(0, diff(RC_experimental) / kRC)
  
  AcRC_theoretical <- aAcRC * (1 - exp(-lambda_AcRC * t_values))       
  
  QCm_theoretical <- c(0, kLBCExsyn * QBm_theoretical[-length(QBm_theoretical)] - kRC * QCm_theoretical[-length(QCm_theoretical)] - kLBCExsyn * QCm_theoretical[-length(QCm_theoretical)] - kCimsyn * QCm_theoretical[-length(QCm_theoretical)])
  QCm_theoretical <- ifelse(QCm_theoretical > 0, QCm_theoretical, 0)
  
  QCim_theoretical <- c(0, (C2C * 10 * CmT) * (1 - exp(-kCimsyn * t_values[-length(t_values)])))
  QCim_theoretical <- ifelse(QCim_theoretical > 0, QCim_theoretical, 0)
  
  QC_theoretical <- QCm_theoretical + QCim_theoretical
  
  QCnetchange_theoretical <- c(0, diff(QC_theoretical))
  
  Ex_theoretical <- c(0, diff(kLBCExsyn * QCm_theoretical))
  
  AcEx_theoretical <- c(0, cumsum(Ex_theoretical))
  AcEx_theoretical <- aAcEx * (1 - exp(-lambda_AcEx * t_values))
  
  Ex_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    Ex_theoretical[i+1] <- AcEx_theoretical[i+1] - AcEx_theoretical[i] 
  }
  
  
  RL_experimental <- data$RL_experimental
  RB_experimental <- data$RB_experimental
  RC_experimental <- data$RC_experimental
  AcRL_experimental <- data$AcRL_experimental
  AcRB_experimental <- data$AcRB_experimental
  AcRC_experimental <- data$AcRC_experimental
  
  RL_theoretical <- RL_experimental
  
  AcRL_theoretical <- c(0, cumsum(RL_theoretical))
  AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))       
  
  AcRB_theoretical <- c(0, cumsum(RB_theoretical))
  AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))   
  
  QL_constraint <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QB_constraint <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  
  QLnetchange_theoretical <- c(0, diff(QL_theoretical))
  QBnetchange_theoretical <- c(0, diff(QB_theoretical))
  
  # General constraints
  Total_theoretical <- (QL_constraint + QB_constraint + QC_theoretical + AcEx_theoretical + AcRL_experimental + AcRB_experimental + AcRC_experimental)
  Total_constraint <- Total13C
  
  netdecrease_constraint <- ifelse(any(QCnetchange_theoretical > 0), -QLnetchange_theoretical - QBnetchange_theoretical, -QLnetchange_theoretical - QBnetchange_theoretical - QCnetchange_theoretical)
  
  netincrease_constraint <- ifelse(any(QCnetchange_theoretical > 0), QCnetchange_theoretical + Ex_theoretical + RL_experimental + RB_experimental + RC_experimental, Ex_theoretical + RL_experimental + RB_experimental + RC_experimental)
  
  rmse_value <- sqrt(mean((netdecrease_constraint - netincrease_constraint)^2)) + sqrt(mean((Total_constraint - Total_theoretical)^2))

  
  if (!is.finite(rmse_value)) {
    rmse_value <- Inf  # Set to infinity if RMSE value is not finite
  }
    
  return(rmse_value)
}

# Set initial guess for parameters
initial_guess_exp <- c(0.01, 0.01, 0.01)

# Minimize RMSE
result_exp <- optim(par = initial_guess_exp, fn = exp_rmse, method = "L-BFGS-B", control = list(maxit = 80000))

# Optimal values
aAcEx <- result_exp$par[1]
lambda_AcEx <- result_exp$par[2]
kLBCExsyn <- result_exp$par[3]

# Print optimal values
cat("Optimal values:\n")
cat("aAcEx:", aAcEx, "\n")
cat("lambda_AcEx:", lambda_AcEx, "\n")
cat("kLBCExsyn:", kLBCExsyn, "\n")
}

## Final optimal for all synthesised kinetic inside the mature bamboo
{
#install.packages("pracma")
#library("pracma")
#install pracma if R require it

# Define function to calculate RMSE
rmse <- function(params) {
  kLBCExsyn <- params[1]
  kCimsyn <- params[2]
  k1C <- params[3]
  
  t_values <- seq(0, 30)  # t from 0 to 30
  
  # Relatives from leaves
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QLnetchange_theoretical <- numeric(length(t_values))
  QLnetchange_theoretical[1] <- 0
  for (i in 2:length(t_values)) {
    QLnetchange_theoretical[i] <- QL_theoretical[i] - QL_theoretical[i-1]
  }
  QLim_theoretical <- (C2L * 10 * LmT) * (1 - exp(-kLimsyn * t_values))
  QLim_theoretical <- ifelse(QLim_theoretical > 0, QLim_theoretical, 0)
  QLm_theoretical <- QL_theoretical - QLim_theoretical
  QLm_theoretical <- ifelse(QLm_theoretical > 0, QLm_theoretical, 0)
  
  RL_experimental <- data$RL_experimental
  AcRL_experimental <- data$AcRL_experimental
  
  RL_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    RL_theoretical[i] <- QLm_theoretical[i-1] * kRL
  }
  
  AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))
  
  
  # Relatives from branch
  QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  QBnetchange_theoretical <- c(0, diff(QB_theoretical))
  QBim_theoretical <- (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))
  QBim_theoretical <- ifelse(QBim_theoretical > 0, QBim_theoretical, 0)
  QBm_theoretical <- QB_theoretical - QBim_theoretical
  QBm_theoretical <- ifelse(QBm_theoretical > 0, QBm_theoretical, 0)
  
  RB_experimental <- data$RB_experimental
  AcRB_experimental <- data$AcRB_experimental
  
  RB_theoretical <- c(0, diff(QBm_theoretical) * kRB)
  AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))
  
  # Calculate values from culm
  
  RC_experimental <- data$RC_experimental
  AcRC_experimental <- data$AcRC_experimental
  
  QCm_theoretical <- c(0, diff(RC_experimental) / kRC)
  
  AcRC_theoretical <- aAcRC * (1 - exp(-lambda_AcRC * t_values))       
  
  QCm_theoretical <- c(0, kLBCExsyn * QBm_theoretical[-length(QBm_theoretical)] - kRC * QCm_theoretical[-length(QCm_theoretical)] - kLBCExsyn * QCm_theoretical[-length(QCm_theoretical)] - kCimsyn * QCm_theoretical[-length(QCm_theoretical)])
  QCm_theoretical <- ifelse(QCm_theoretical > 0, QCm_theoretical, 0)
  
  QCim_theoretical <- c(0, (C2C * 10 * CmT) * (1 - exp(-kCimsyn * t_values[-length(t_values)])))
  QCim_theoretical <- ifelse(QCim_theoretical > 0, QCim_theoretical, 0)
  
  QC_theoretical <- QCm_theoretical + QCim_theoretical
  
  QCnetchange_theoretical <- c(0, diff(QC_theoretical))
  
  Ex_theoretical <- c(0, diff(kLBCExsyn * QCm_theoretical))
  
  AcEx_theoretical <- c(0, cumsum(Ex_theoretical))
  AcEx_theoretical <- aAcEx * (1 - exp(-lambda_AcEx * t_values))
  
  Ex_theoretical <- numeric(length(t_values))
  for (i in 2:length(t_values)) {
    Ex_theoretical[i+1] <- AcEx_theoretical[i+1] - AcEx_theoretical[i] 
  }
  
  
  RL_experimental <- data$RL_experimental
  RB_experimental <- data$RB_experimental
  RC_experimental <- data$RC_experimental
  AcRL_experimental <- data$AcRL_experimental
  AcRB_experimental <- data$AcRB_experimental
  AcRC_experimental <- data$AcRC_experimental
  
  RL_theoretical <- RL_experimental
  
  AcRL_theoretical <- c(0, cumsum(RL_theoretical))
  AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))       
  
  AcRB_theoretical <- c(0, cumsum(RB_theoretical))
  AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))   
  
  QL_constraint <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QB_constraint <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
  QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
  
  QLnetchange_theoretical <- c(0, diff(QL_theoretical))
  QBnetchange_theoretical <- c(0, diff(QB_theoretical))
  
  # General constraints
  Total_theoretical <- (QL_constraint + QB_constraint + QC_theoretical + AcEx_theoretical + AcRL_experimental + AcRB_experimental + AcRC_experimental)
  Total_constraint <- Total13C
  
  netdecrease_constraint <- ifelse(any(QCnetchange_theoretical > 0), -QLnetchange_theoretical - QBnetchange_theoretical, -QLnetchange_theoretical - QBnetchange_theoretical - QCnetchange_theoretical)
  
  netincrease_constraint <- ifelse(any(QCnetchange_theoretical > 0), QCnetchange_theoretical + Ex_theoretical + RL_experimental + RB_experimental + RC_experimental, Ex_theoretical + RL_experimental + RB_experimental + RC_experimental)
  
  rmse_value <- sqrt(mean((netdecrease_constraint - netincrease_constraint)^2)) + sqrt(mean((Total_constraint - Total_theoretical)^2))
  
  
  if (!is.finite(rmse_value)) {
    rmse_value <- Inf  # Set to infinity if RMSE value is not finite
  }
  
  return(rmse_value)
}

# Set initial guess for parameters
initial_guess <- c(0.01, 0.01, 0.415758923, 0.01)

# Minimize RMSE
result <- optim(par = initial_guess, fn = rmse, method = "L-BFGS-B", control = list(maxit = 80000))

# Optimal values
kLBCExsyn <- result$par[1]
kCimsyn <- result$par[2]
k1C <- result$par[3]

# Print optimal values
cat("Optimal values:\n")
cat("kLBCExsyn:", kLBCExsyn, "\n")
cat("kCimsyn:", kCimsyn, "\n")
cat("k1C:", k1C, "\n")
}

## Final fitting for kRC 
{
  # Define function to calculate RMSE
  kRC_rmse <- function(params) {
    kRC <- params[1]

    t_values <- seq(0, 30)  # t from 0 to 30
    
    # Relatives from leaves
    QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
    QLnetchange_theoretical <- numeric(length(t_values))
    QLnetchange_theoretical[1] <- 0
    for (i in 2:length(t_values)) {
      QLnetchange_theoretical[i] <- QL_theoretical[i] - QL_theoretical[i-1]
    }
    QLim_theoretical <- (C2L * 10 * LmT) * (1 - exp(-kLimsyn * t_values))
    QLim_theoretical <- ifelse(QLim_theoretical > 0, QLim_theoretical, 0)
    QLm_theoretical <- QL_theoretical - QLim_theoretical
    QLm_theoretical <- ifelse(QLm_theoretical > 0, QLm_theoretical, 0)
    
    RL_experimental <- data$RL_experimental
    AcRL_experimental <- data$AcRL_experimental
    
    RL_theoretical <- numeric(length(t_values))
    for (i in 2:length(t_values)) {
      RL_theoretical[i] <- QLm_theoretical[i-1] * kRL
    }
    
    AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))
    
    
    # Relatives from branch
    QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
    QBnetchange_theoretical <- c(0, diff(QB_theoretical))
    QBim_theoretical <- (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))
    QBim_theoretical <- ifelse(QBim_theoretical > 0, QBim_theoretical, 0)
    QBm_theoretical <- QB_theoretical - QBim_theoretical
    QBm_theoretical <- ifelse(QBm_theoretical > 0, QBm_theoretical, 0)
    
    RB_experimental <- data$RB_experimental
    AcRB_experimental <- data$AcRB_experimental
    
    RB_theoretical <- c(0, diff(QBm_theoretical) * kRB)
    AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))
    
    # Calculate values from culm
    
    RC_experimental <- data$RC_experimental
    AcRC_experimental <- data$AcRC_experimental
    
    QCm_theoretical <- c(0, diff(RC_experimental) / kRC)
    
    AcRC_theoretical <- aAcRC * (1 - exp(-lambda_AcRC * t_values))       
    
    QCm_theoretical <- c(0, kLBCExsyn * QBm_theoretical[-length(QBm_theoretical)] - kRC * QCm_theoretical[-length(QCm_theoretical)] - kLBCExsyn * QCm_theoretical[-length(QCm_theoretical)] - kCimsyn * QCm_theoretical[-length(QCm_theoretical)])
    QCm_theoretical <- ifelse(QCm_theoretical > 0, QCm_theoretical, 0)
    
    QCim_theoretical <- c(0, (C2C * 10 * CmT) * (1 - exp(-kCimsyn * t_values[-length(t_values)])))
    QCim_theoretical <- ifelse(QCim_theoretical > 0, QCim_theoretical, 0)
    
    QC_theoretical <- QCm_theoretical + QCim_theoretical
    
    QCnetchange_theoretical <- c(0, diff(QC_theoretical))
    
    Ex_theoretical <- c(0, diff(kLBCExsyn * QCm_theoretical))
    
    AcEx_theoretical <- c(0, cumsum(Ex_theoretical))
    AcEx_theoretical <- aAcEx * (1 - exp(-lambda_AcEx * t_values))
    
    Ex_theoretical <- numeric(length(t_values))
    for (i in 2:length(t_values)) {
      Ex_theoretical[i+1] <- AcEx_theoretical[i+1] - AcEx_theoretical[i] 
    }
    
    
    RL_experimental <- data$RL_experimental
    RB_experimental <- data$RB_experimental
    RC_experimental <- data$RC_experimental
    AcRL_experimental <- data$AcRL_experimental
    AcRB_experimental <- data$AcRB_experimental
    AcRC_experimental <- data$AcRC_experimental
    
    RL_theoretical <- RL_experimental
    
    AcRL_theoretical <- c(0, cumsum(RL_theoretical))
    AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))       
    
    AcRB_theoretical <- c(0, cumsum(RB_theoretical))
    AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))   
    
    QL_constraint <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
    QB_constraint <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
    QL_theoretical <- (LmT * (C1L * exp(-k1L * t_values) + C2L) * 10)
    QB_theoretical <- (BmT * (C1B * exp(-k1B * t_values) + C2B) * 10)
    
    QLnetchange_theoretical <- c(0, diff(QL_theoretical))
    QBnetchange_theoretical <- c(0, diff(QB_theoretical))
    
    # General constraints
    Total_theoretical <- (QL_constraint + QB_constraint + QC_theoretical + AcEx_theoretical + AcRL_experimental + AcRB_experimental + AcRC_experimental)
    Total_constraint <- Total13C
    
    netdecrease_constraint <- ifelse(any(QCnetchange_theoretical > 0), -QLnetchange_theoretical - QBnetchange_theoretical, -QLnetchange_theoretical - QBnetchange_theoretical - QCnetchange_theoretical)
    
    netincrease_constraint <- ifelse(any(QCnetchange_theoretical > 0), QCnetchange_theoretical + Ex_theoretical + RL_experimental + RB_experimental + RC_experimental, Ex_theoretical + RL_experimental + RB_experimental + RC_experimental)
    
    
    
    rmse_value <- sqrt(mean((netdecrease_constraint - netincrease_constraint)^2)) + sqrt(mean((Total_constraint - Total_theoretical)^2))
    
    
    if (!is.finite(rmse_value)) {
      rmse_value <- Inf  # Set to infinity if RMSE value is not finite
    }
    
    return(rmse_value)
  }
  
  # Set initial guess for parameters
  initial_guess_kRC <- c(0.01)
  
  # Minimize RMSE
  result_kRC <- optim(par = initial_guess_kRC, fn = kRC_rmse, method = "L-BFGS-B", control = list(maxit = 80000))
  
  # Optimal values
  kRC <- result$par[1]

  
  # Print optimal values
  cat("Optimal values:\n")
  cat("kRC:", kRC, "\n")
}
}

### Step 3. Calculate the optimised output
{
# Define parameters
t_values <- 0:30

# Calculate theoretical values
QL_theoretical <- LmT * (C1L * exp(-k1L * t_values) + C2L) * 10
QB_theoretical <- BmT * (C1B * exp(-k1B * t_values) + C2B) * 10
QBm_theoretical <- QB_theoretical - (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))

QCm_theoretical <- numeric(length(t_values))
QCm_theoretical[1] <- 0
for (i in 2:length(t_values)) {
  QCm_theoretical[i] <- kLBCExsyn * QBm_theoretical[i-1] - kRC * QCm_theoretical[i-1] - kLBCExsyn * QCm_theoretical[i-1] - kCimsyn * QCm_theoretical[i-1]
}
QCm_theoretical <- ifelse(QCm_theoretical > 0, QCm_theoretical, 0)

QLim_theoretical <- (C2L * 10 * LmT) * (1 - exp(-kLimsyn * t_values))
QBim_theoretical <- (C2B * 10 * BmT) * (1 - exp(-kBimsyn * t_values))

QCim_theoretical <- numeric(length(t_values))
QCim_theoretical[1] <- 0
for (i in 2:length(t_values)) {
  QCim_theoretical[i] <- (C2C * 10 * CmT) * (1 - exp(-kCimsyn * t_values[i]))
}
QCim_theoretical <- ifelse(QCim_theoretical > 0, QCim_theoretical, 0)

QLm_theoretical <- QL_theoretical - QLim_theoretical
QBm_theoretical <- QB_theoretical - QBim_theoretical
QC_theoretical <- QCm_theoretical + QCim_theoretical

# Calculate net change
QLnetchange_theoretical <- c(0, diff(QL_theoretical))
QBnetchange_theoretical <- c(0, diff(QB_theoretical))
QCnetchange_theoretical <- numeric(length(t_values))
QCnetchange_theoretical[1] <- 0
for (i in 2:length(t_values)) {
  QCnetchange_theoretical[i] <- (CmT * (C1C * exp(-k1C * (t_values[i])) + C2C) * 10) - (CmT * (C1C * exp(-k1C * t_values[i-1]) + C2C) * 10)
  if (QCnetchange_theoretical[i] < 0) {
    QCnetchange_theoretical[i] <- (C2C * 10 * CmT) * (1 - exp(-k1C * t_values[i]))
  }


# Calculate accumulation
#AcRL_theoretical <- aAcRL * (1 - exp(-lambda_AcRL * t_values))
#AcRB_theoretical <- aAcRB * (1 - exp(-lambda_AcRB * t_values))
#AcRC_theoretical <- aAcRC * (1 - exp(-lambda_AcRC * t_values))
AcRL_theoretical <- data$AcRL_experimental
AcRB_theoretical <- data$AcRB_experimental
AcRC_theoretical <- data$AcRC_experimental
AcEx_theoretical <- aAcEx * (1 - exp(-lambda_AcEx * t_values))

#RL_theoretical <- numeric(length(t_values))
#RB_theoretical <- numeric(length(t_values))
#RC_theoretical <- numeric(length(t_values))
RL_theoretical <- data$RL_experimental
RB_theoretical <- data$RB_experimental
RC_experimental <- data$RC_experimental
Ex_theoretical <- numeric(length(t_values))
for (i in 2:length(t_values)) {
  Ex_theoretical[i] <- AcEx_theoretical[i] - AcEx_theoretical[i-1] 
}


# Calculate total theoretical value
Total13C_theoretical <- QLm_theoretical + QBm_theoretical + QCm_theoretical +
  QLim_theoretical + QBim_theoretical + QCim_theoretical +
  AcEx_theoretical + AcRL_theoretical + AcRB_theoretical + AcRC_theoretical
}

# Create a data frame to store the results
df <- data.frame(
  t = t_values,
  QL = QL_theoretical,
  QB = QB_theoretical,
  QC = QC_theoretical,
  QLim = QLim_theoretical,
  QBim = QBim_theoretical,
  QCim = QCim_theoretical,
  QLm = QLm_theoretical,
  QBm = QBm_theoretical,
  QCm = QCm_theoretical,
  QLnetchange = QLnetchange_theoretical,
  QBnetchange = QBnetchange_theoretical,
  QCnetchange = QCnetchange_theoretical,
  RL = RL_theoretical,
  RB = RB_theoretical,
  RC = RC_experimental,
  Ex = Ex_theoretical,
  AcRL = AcRL_theoretical,
  AcRB = AcRB_theoretical,
  AcRC = AcRC_theoretical,
  AcEx = AcEx_theoretical,
  Total13C_theoretical = Total13C_theoretical
)

# Save the data frame to a CSV file
write.csv(df, file = "theoretical_valuesfittoexperimentsT2_R_4_Abe.csv", row.names = FALSE)
}

### Step 4. Visualisation
{
#install.packages("ggplot2")
library(ggplot2)

# Plotting
p <- ggplot(df, aes(x = t)) +
  geom_line(aes(y = QL, color = "QL")) +
  geom_line(aes(y = QB, color = "QB")) +
  geom_line(aes(y = QC, color = "QC")) +
  geom_line(aes(y = QLim, color = "QLim")) +
  geom_line(aes(y = QBim, color = "QBim")) +
  geom_line(aes(y = QCim, color = "QCim")) +
  geom_line(aes(y = QLm, color = "QLm")) +
  geom_line(aes(y = QBm, color = "QBm")) +
  geom_line(aes(y = QCm, color = "QCm")) +
  geom_line(aes(y = Ex, color = "Ex")) +
  geom_line(aes(y = AcEx, color = "AcEx")) +
  geom_line(aes(y = AcRL, color = "AcRL")) +
  geom_line(aes(y = AcRB, color = "AcRB")) +
  geom_line(aes(y = AcRC, color = "AcRC")) +
  geom_line(aes(y = Total13C_theoretical, color = "Total13C_theoretical")) +
  labs(x = "Day after labelling", y = "13C (g)") +
  scale_color_manual(name = "Variables",
                     values = c(QL = "darkgreen", QB = "green", QC = "lightgreen", QLim = "darkblue", 
                                QBim = "blue", QCim = "lightblue", QLm = "darkred", QBm = "red", 
                                QCm = "brown", Ex = "orange", AcEx = "darkorange", AcRLexp = "darkgray", 
                                AcRB = "gray", AcRC = "lightgray", 
                                Total13C_theoretical = "black")) +
  theme_minimal()

# Save as JPEG with 600 DPI
ggsave("C:/Users/Sorbus/Downloads/work log on task 2/T2_3organs_plot_R4_Abe.jpg", p, dpi = 600, width = 10, height = 6)
}
