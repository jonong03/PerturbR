# =============================================================================
# Shiny App: 3D Visualization for Mapping Correlation to R-squared (Use case)
# =============================================================================
# Required packages:
#   shiny, plotly, mvtnorm, rethinking
# =============================================================================

library(shiny)
library(plotly)
library(mvtnorm)    # rmvnorm
library(rethinking) # rlkjcorr

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
calc_psi <- function(R, n = 1) {
  R <- as.matrix(R)
  p <- nrow(R)
  
  if (ncol(R) != p) { stop("R must be a square matrix.") }
  
  idx <- list()
  k <- 1
  
  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      idx[[k]] <- c(i, j)
      k <- k + 1
    }
  }
  
  k <- length(idx)
  psi <- diag(k)
  
  for (i in seq_len(k)) {
    e <- idx[[i]][1]
    f <- idx[[i]][2]
    
    psi[i, i] <- (1 - R[e, f]^2)^2
    
    if (i > 1) {
      for (j in seq_len(i - 1)) {
        g <- idx[[j]][1]
        h <- idx[[j]][2]
        
        tmp <- c(
          (R[e, g] - R[e, f] * R[f, g]) *
            (R[f, h] - R[f, g] * R[g, h]),
          
          (R[e, h] - R[e, g] * R[g, h]) *
            (R[f, g] - R[f, e] * R[e, g]),
          
          (R[e, g] - R[e, h] * R[h, g]) *
            (R[f, h] - R[f, e] * R[e, h]),
          
          (R[e, h] - R[e, f] * R[f, h]) *
            (R[f, g] - R[f, h] * R[h, g])
        )
        
        psi[i, j] <- 0.5 * sum(tmp)
        psi[j, i] <- psi[i, j]
      }
    }
  }
  
  return(psi/n)
}
rwaldcloud<- function(R0, Psi, n= 1e5, B=1, tol=1e-5, verbose=FALSE){
  # Psi must be defined at sample size = 1
  # If Psi is not defined, it will be estimated by calc_psi(R0, n= n)
  
  d<- ncol(R0)
  ps<- d*(d-1)/2
  out <- array(NA_real_, dim = c(d, d, B))
  
  if(missing(Psi)){
    Psi <- calc_psi(R0, n = n)
  }
  
  # Eigen Decomposition and check if R0 (and Psi0) can be decomposed
  eig<- eigen(Psi, symmetric = TRUE)
  if (any(eig$values < -tol)) {
    stop("Psi0 is not positive semidefinite.")
  }
  
  A <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = ps)
  
  r0 <- R0[lower.tri(R0)]
  b <- 0
  iter <- 0
  non.psd.count <- 0
  max.iter <- 100
  
  while (b < B && iter < max.iter) {
    iter <- iter + 1
    
    z <- matrix(rnorm(ps), ncol = 1)
    rhatb <- as.vector(r0 + A %*% z)
    
    Rhatb <- diag(d)
    Rhatb[lower.tri(Rhatb)] <- rhatb
    Rhatb[upper.tri(Rhatb)] <- t(Rhatb)[upper.tri(Rhatb)]
    
    evals <- eigen(Rhatb, symmetric = TRUE, only.values = TRUE)$values
    
    if (all(evals > tol)) {
      b <- b + 1
      out[, , b] <- Rhatb
    } else {
      non.psd.count <- non.psd.count + 1
    }
  }
  
  if (verbose) {
    cat("Accepted:", b, "; Rejected:", non.psd.count,
        "; Total iterations:", iter, "\n")
  }
  
  if (b < B) {
    stop(sprintf("Could not generate %d positive-definite samples within max.iter = %d.", B, max.iter))
  }
  
  #attr(out, "accepted") <- b
  #attr(out, "rejected") <- non.psd.count
  #attr(out, "iterations") <- iter
  
  return(out)
}
wald.test <- function(Rpop, Rsample, alpha = 0.05, asy.n = 100000, fisherz = FALSE) {
  # Rpop     : Target (true) correlation matrix, this will be treated as the center and used to compute variance
  # Rsample  : Estimated correlation matrix
  # alpha    : Type I error rate
  # indep    : If TRUE, use only marginal variances (diagonal Psi)
  # asy.n    : Sample size for asymptotic covariance
  # fisherz  : If TRUE, use Fisher-z transformation
  # pd_check : If TRUE, reject (return FALSE) when Rsample is not PD
  # pd_tol   : PD tolerance for smallest eigenvalue
  
  asyCov.z <- function(R, asy.n= 10000) {
    # R: untransformed correlation matrix
    
    rvec <- R[lower.tri(R)]
    cov <- calc_psi(R, n= asy.n)
    
    #denominator: 1-rho^2
    #denominator[i,j] = (1-r_i^2)*(1-r_j^2)
    weights <- 1 - rvec^2
    denom_matrix <- outer(weights, weights)
    
    #Steiger defines psi = N * sigma
    psi_matrix <- asy.n * cov
    
    #Transformation: Eq10 and Eq11
    c_matrix <- psi_matrix / denom_matrix
    zcov_matrix <- c_matrix / (asy.n-3)
    diag(zcov_matrix)<- 1/ (asy.n-3)
    
    return(zcov_matrix)  
  }
  
  ## --- vectorize correlations ---
  Rpop_vec    <- Rpop[lower.tri(Rpop)]
  Rsample_vec <- Rsample[lower.tri(Rsample)]
  ps <- length(Rpop_vec)  # dimension = p*(p-1)/2
  
  ## --- Fisher z transform if requested ---
  if (fisherz) {
    center_vec <- atanh(Rpop_vec)
    sample_vec <- atanh(Rsample_vec)
    
    Psi <- asyCov.z(Rpop, asy.n = asy.n)     # asymptotic covariance in z-space
  } else {
    center_vec <- Rpop_vec
    sample_vec <- Rsample_vec
    
    Psi <- calc_psi(Rpop, n = asy.n)  # asymptotic covariance in r-space
  }
  
  Psi0 <- Psi* asy.n   # covariance when sample size = 1, use this to study eigenstructure
  
  if(qr(Psi)$rank != ps) { return(c(T=NA, df=NA, cricval=NA, pval=NA, reject= NA)) }
  
  ## --- Mahalanobis distance ---
  distance = mahalanobis(x= sample_vec,center = center_vec,cov = Psi)
  df= ps
  
  ## --- chi-square cutoff ---
  cricval<- qchisq(1 - alpha, df = df)
  pval<- as.numeric(pchisq(distance, df= df, lower.tail = F))
  reject= as.logical(distance >= cricval)
  out = data.frame(T=distance, df= df, cricval= as.numeric(cricval), pval= pval, reject= reject) 
  
  ## --- ellipsoid membership ---
  return(out)
}

gen_AR1 <- function(p, rho = 0.6) {
  Rxx <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) Rxx[i, j] <- 1
      if (i != j) Rxx[i, j] <- rho^(abs(i - j))
    }
  }
  return(Rxx)
}

computeRsq <- function(M) {
  p   <- ncol(M)
  Rxy <- as.matrix(M[2:p, 1], ncol = 1)
  Rxx <- M[2:p, 2:p]
  Rsq <- t(Rxy) %*% solve(Rxx) %*% Rxy
  return(Rsq)
}

make_mcmc<- function(S, asy.n, NCHAIN=1000, init= S, alpha=0.05, 
                     scale.init= 1, adaptive= TRUE, M= 100, acceptance.target = 0.44, 
                     seed= 123){
  p<- ncol(S)    #p: dimension of S (correlation matrix)
  d<- p*(p-1)/2  #d: dimension of correlation space
  accept <- 0
  q<- 0         # to trace scale_factor
  
  trace_scale<- array(NA_real_, dim=c(NCHAIN))
  scale_factor<- scale.init
  set.seed(seed)
  
  RCHAIN<- array(NA, dim= c(p,p,NCHAIN))
  RCHAIN[,,1] <- init
  
  for (i in 2:NCHAIN){
    
    Rcurrent<- RCHAIN[,,i-1]
    rcurrent<- Rcurrent[lower.tri(Rcurrent)]
    
    # Step 1: Draw R* (named as Rs)
    Psi_Rcurrent<- calc_psi(Rcurrent, n = asy.n)
    
    Rs<- rwaldcloud(Rcurrent, Psi= scale_factor*Psi_Rcurrent, B = 1) %>% drop()
    rs<- Rs[lower.tri(Rs)]
    
    Psi_Rs<- calc_psi(Rs, n= asy.n)
    # Step 2: Acceptance Ratio
    partA<- dmvnorm(x= rcurrent, mean = rs, sigma= Psi_Rs, log= T)
    partB<- dmvnorm(x= rs, mean= rcurrent, sigma= Psi_Rcurrent, log=T)
    logMH = partA - partB 
    
    partC<- wald.test(Rpop=Rs, Rsample= S, alpha= alpha, asy.n= asy.n)  # Rstar is at the center
    IR   <- partC$pval > alpha   # accept Rs is a potential R0 for S
    
    # Step 3: accept or reject
    accept_proposal <- !is.na(partC[4]) && IR && (logMH > log(runif(1)))
    if (accept_proposal) {
      RCHAIN[, , i] <- Rs
      accept <- accept + 1
    } else {
      RCHAIN[, , i] <- Rcurrent
    }
    
    a_rate = accept/i
    cat("Acceptance Rate: ", round(a_rate,3),"\n")
    
    if(adaptive==TRUE){
      # update scale parameter every M steps
      ### If acceptance rate > 44% (acceptance.target parameter) 
      ############################ (step size is too small, then update scale to min(acceptance, 0.75)/ 0.44)
      ### If acceptance rate < 44% (step size is too big, then update scale to max(acceptance, 0.20)/ 0.44)
      ##### M: adaptive block
      if(i/M == floor(i/M)){
        a_rate=ifelse(a_rate>.75,.75,ifelse(a_rate<.1,.1,a_rate))
        scale_factor=scale_factor*a_rate/acceptance.target
        # cat("Scale Factor: ", scale_factor, "\n")
      }
    }
    trace_scale[i]<- scale_factor
  }
  attr(RCHAIN, "scale_trace") <- trace_scale
  
  return(RCHAIN)
}

simdriver <- function(Rxx, beta, k = 0.85, ss = 100, nchain = 1000,
                      alpha = 0.05, adaptive = FALSE,
                      acceptance.target = 0.44) {
  p      <- ncol(Rxx)
  S      <- t(beta) %*% Rxx %*% beta
  sigma2 <- S * (1 - k) / k

  X   <- rmvnorm(n = ss, mean = rep(0, p), sigma = Rxx)
  eps <- rnorm(n = ss, mean = 0, sd = sqrt(sigma2))
  Y   <- X %*% beta + eps
  dt  <- cbind(Y, X)

  Rhat       <- cor(cbind(Y, X))
  khat       <- computeRsq(Rhat)

  RCHAIN     <- make_mcmc(S = Rhat, asy.n = ss, NCHAIN = nchain,
                           alpha = alpha, adaptive = adaptive,
                           acceptance.target = acceptance.target)
  khat.chain <- sapply(seq_len(nchain), function(i) computeRsq(RCHAIN[, , i]))

  RBOOT <- array(NA_real_, dim = c(p + 1, p + 1, nchain))
  for (i in seq_len(nchain)) {
    boot.id      <- sample(ss, ss, replace = TRUE)
    RBOOT[, , i] <- cor(dt[boot.id, ])
  }
  khat.boot <- sapply(seq_len(nchain), function(i) computeRsq(RBOOT[, , i]))

  list(
    Rhat       = Rhat,
    mcmc       = RCHAIN,
    boot       = RBOOT,
    khat.point = khat,
    khat.chain = khat.chain,
    khat.boot  = khat.boot
  )
}

# helper: format a 3x3 matrix as an HTML table with y, x1, x2 row/col headers
mat_html <- function(M, label) {
  hdrs  <- c("y", "x\u2081", "x\u2082")   # y, x1, x2 with subscript unicode
  nr    <- nrow(M)
  nc    <- ncol(M)

  # header row
  header_cells <- c(
    list(tags$th(style = "padding:3px 8px; background:#eeeeee;", "")),  # corner
    lapply(hdrs, function(h)
      tags$th(style = "padding:3px 8px; text-align:center; background:#eeeeee;
                        font-weight:bold; border-bottom:1px solid #aaaaaa;", h)
    )
  )
  header_row <- tags$tr(header_cells)

  # data rows
  data_rows <- lapply(seq_len(nr), function(i) {
    row_header <- tags$th(
      style = "padding:3px 8px; text-align:center; background:#eeeeee;
               font-weight:bold; border-right:1px solid #aaaaaa;",
      hdrs[i]
    )
    cells <- lapply(seq_len(nc), function(j) {
      tags$td(
        style = "padding:3px 10px; text-align:center;",
        formatC(round(M[i, j], 3), format = "f", digits = 3)
      )
    })
    tags$tr(c(list(row_header), cells))
  })

  tags$div(
    style = "display:inline-block; margin-right:36px; vertical-align:top;",
    tags$div(
      style = "font-weight:bold; font-size:1rem; margin-bottom:6px;",
      HTML(label)
    ),
    tags$table(
      style = "border-collapse:collapse; font-size:0.95rem; font-family:monospace;
               border:1px solid #cccccc; background:#ffffff;",
      tags$thead(header_row),
      tags$tbody(data_rows)
    )
  )
}

# helper: R2 range display block
r2_range_html <- function(label, lo, hi) {
  tags$div(
    style = "display:inline-block; vertical-align:top; margin-right:36px;
             background:#f9f9f9; border:1px solid #cccccc; border-radius:5px;
             padding:10px 16px; min-width:160px;",
    tags$div(style = "font-weight:bold; font-size:1rem; margin-bottom:6px;",
             HTML(label)),
    tags$div(style = "font-family:monospace; font-size:0.95rem; color:#333333;",
      tags$span(style="color:#777;", "min: "),
      tags$span(formatC(round(lo, 3), format = "f", digits = 3)),
      tags$br(),
      tags$span(style="color:#777;", "max: "),
      tags$span(formatC(round(hi, 3), format = "f", digits = 3))
    )
  )
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

ui <- fluidPage(

  tags$head(
    tags$style(HTML("

      /* ── Base ── */
      body {
        background: #ffffff;
        color: #1a1a1a;
        font-family: Georgia, 'Times New Roman', serif;
        font-size: 17px;
      }

      /* ── Layout shell: sidebar shrunk by 30% from 320px = 224px ── */
      .app-shell {
        display: grid;
        grid-template-columns: 224px 1fr;
        grid-template-rows: auto 1fr;
        min-height: 100vh;
      }

      /* ── Header ── */
      .app-header {
        grid-column: 1 / -1;
        padding: 20px 30px 16px;
        border-bottom: 2px solid #1a1a1a;
        background: #ffffff;
      }
      .app-header h1 {
        font-size: 1.5rem;
        font-weight: bold;
        color: #1a1a1a;
        margin: 0 0 2px 0;
      }
      .app-header .subtitle {
        font-size: 1rem;
        color: #555555;
      }

      /* ── Sidebar ── */
      .sidebar-panel {
        grid-column: 1;
        padding: 18px 14px 28px;
        border-right: 1px solid #cccccc;
        background: #f9f9f9;
        overflow-y: auto;
      }

      .param-section { margin-bottom: 22px; }

      .param-section-label {
        font-size: 1.1rem;
        font-weight: bold;
        letter-spacing: 0.05em;
        text-transform: uppercase;
        color: #222222;
        margin-bottom: 10px;
        padding-bottom: 5px;
        border-bottom: 2px solid #aaaaaa;
      }

      /* Shiny input overrides */
      .form-group { margin-bottom: 10px; }

      label {
        font-size: 1.05rem;
        color: #222222;
        display: block;
        margin-bottom: 3px;
      }

      input[type='number'],
      input[type='text'] {
        background: #ffffff !important;
        border: 1px solid #aaaaaa !important;
        border-radius: 4px !important;
        color: #1a1a1a !important;
        font-size: 1.05rem !important;
        padding: 5px 8px !important;
        width: 100% !important;
      }
      input[type='number']:focus,
      input[type='text']:focus {
        border-color: #333333 !important;
        outline: none !important;
        box-shadow: 0 0 0 2px rgba(0,0,0,0.12) !important;
      }

      /* Run button */
      #run_btn {
        width: 100%;
        margin-top: 4px;
        padding: 10px 0;
        background: #1a1a1a;
        border: none;
        border-radius: 4px;
        color: #ffffff;
        font-size: 1.1rem;
        font-weight: bold;
        letter-spacing: 0.05em;
        cursor: pointer;
        transition: background 0.2s;
      }
      #run_btn:hover  { background: #444444; }
      #run_btn:active { background: #000000; }

      /* ── Main panel ── */
      .main-panel {
        grid-column: 2;
        padding: 22px 26px;
        display: flex;
        flex-direction: column;
        gap: 14px;
        background: #ffffff;
      }

      /* Status bar */
      .status-bar {
        font-size: 1.05rem;
        color: #777777;
        min-height: 20px;
      }
      .status-bar.done { color: #2a7a2a; font-weight: bold; }

      /* Output info row: matrices + R2 ranges */
      .output-info-row {
        display: flex;
        align-items: flex-start;
        gap: 0px;
        flex-wrap: wrap;
        padding: 8px 0 4px 0;
      }

      /* Plot container */
      .plot-container {
        flex: 1;
        border: 1px solid #dddddd;
        border-radius: 6px;
        padding: 4px;
        min-height: 520px;
        background: #ffffff;
      }

      /* Info chips */
      .info-chips {
        display: flex;
        gap: 10px;
        flex-wrap: wrap;
      }
      .chip {
        font-size: 1.05rem;
        padding: 5px 14px;
        border-radius: 20px;
        background: #f0f0f0;
        border: 1px solid #cccccc;
        color: #333333;
      }
      .chip span { font-weight: bold; color: #1a1a1a; }

    "))
  ),

  div(class = "app-shell",

    # Header
    div(class = "app-header",
      h1("3D Correlation Space vs Rsquared"),
      span(class = "subtitle", "\u03B8 Space \u00B7 MCMC \u00B7 Bootstrap \u00B7 Full Space")
    ),

    # ── Sidebar ──────────────────────────────────────────────────────────────
    div(class = "sidebar-panel",

      # 1. Population Parameters
      div(class = "param-section",
        div(class = "param-section-label", "Population Parameters"),
        numericInput("rho_j",
          label = HTML("\u03C1<sub>j</sub> &mdash; AR(1) autocorrelation"),
          value = 0.5, min = -0.99, max = 0.99, step = 0.05),
        numericInput("beta1",
          label = HTML("\u03B2<sub>1</sub>"),
          value = 0.6, step = 0.05),
        numericInput("beta2",
          label = HTML("\u03B2<sub>2</sub>"),
          value = 0.6, step = 0.05),
        numericInput("ktrue",
          label = HTML("R\u00B2 &mdash; true population R-squared"),
          value = 0.7, min = 0.01, max = 0.99, step = 0.05)
      ),

      # 2. Sampling
      div(class = "param-section",
        div(class = "param-section-label", "Sampling"),
        numericInput("ss_j",
          label = "Sample size",
          value = 30, min = 10, max = 5000, step = 10)
      ),

      # 3. Perturbation
      div(class = "param-section",
        div(class = "param-section-label", "Perturbation"),
        numericInput("nchain",
          label = "Number of perturbations",
          value = 2000, min = 100, max = 50000, step = 500),
        numericInput("alpha",
          label = HTML("Type I error for MCMC (\u03B1)"),
          value = 0.05, min = 0.001, max = 0.2, step = 0.005)
      ),

      # Run button
      div(class = "param-section",
        actionButton("run_btn", "\u25B6  Run Simulation")
      )

    ), # end sidebar

    # ── Main panel ────────────────────────────────────────────────────────────
    div(class = "main-panel",

      uiOutput("status_ui"),

      # Matrices + R2 ranges above the plot
      uiOutput("output_info_ui"),

      div(class = "info-chips",
        uiOutput("chips_ui")
      ),

      div(class = "plot-container",
        plotlyOutput("scatter3d", height = "100%")
      )

    ) # end main panel

  ) # end app-shell
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {

  sim_result <- eventReactive(input$run_btn, {
    withProgress(message = "Running simulation...", value = 0, {

      nvar   <- 2
      rho_j  <- input$rho_j
      beta_j <- c(input$beta1, input$beta2)
      ss_j   <- input$ss_j
      nchain <- input$nchain
      ktrue  <- input$ktrue
      alpha  <- input$alpha

      incProgress(0.1, detail = "Building AR(1) matrix")
      M_j <- gen_AR1(rho = rho_j, p = nvar)

      incProgress(0.1, detail = "Computing population R")
      Rxx    <- M_j
      beta   <- beta_j
      p      <- nvar
      S      <- t(beta) %*% Rxx %*% beta
      sigma2 <- S * (1 - ktrue) / ktrue
      num    <- t(beta) %*% Rxx
      denom  <- sqrt((t(beta) %*% Rxx %*% beta) + sigma2)
      r_yx   <- as.vector(num) / c(denom)
      Rpop   <- diag(p + 1)
      Rpop[1, 1] <- 1
      Rpop[2:(p+1), 2:(p+1)] <- Rxx
      Rpop[1, 2:(p+1)] <- r_yx
      Rpop[2:(p+1), 1] <- r_yx

      incProgress(0.2, detail = "Simulating data & MCMC")
      runs_j <- simdriver(
        Rxx    = M_j,
        beta   = beta_j,
        k      = ktrue,
        alpha  = alpha,
        ss     = ss_j,
        nchain = nchain,
        adaptive = FALSE,
        acceptance.target = 0.44
      )

      Rhat       <- runs_j$Rhat
      Rmcmc      <- runs_j$mcmc
      Rboot      <- runs_j$boot
      khat.chain <- runs_j$khat.chain
      khat.boot  <- runs_j$khat.boot

      # R2 ranges
      r2_mcmc_lo <- min(khat.chain,  na.rm = TRUE)
      r2_mcmc_hi <- max(khat.chain,  na.rm = TRUE)
      r2_boot_ci <- quantile(khat.boot, probs = c(alpha / 2, 1 - alpha / 2),
                             na.rm = TRUE)
      r2_boot_lo <- r2_boot_ci[1]
      r2_boot_hi <- r2_boot_ci[2]

      incProgress(0.3, detail = "Sampling full correlation space")
      fullspace <- rlkjcorr(n = 10000, K = 3, eta = 1)

      incProgress(0.3, detail = "Done")

      list(
        Rpop       = Rpop,
        Rhat       = Rhat,
        Rmcmc      = Rmcmc,
        Rboot      = Rboot,
        khat.chain = khat.chain,
        khat.boot  = khat.boot,
        fullspace  = fullspace,
        nchain     = nchain,
        ktrue      = ktrue,
        khat.point = runs_j$khat.point,
        r2_mcmc_lo = r2_mcmc_lo,
        r2_mcmc_hi = r2_mcmc_hi,
        r2_boot_lo = r2_boot_lo,
        r2_boot_hi = r2_boot_hi,
        alpha      = alpha
      )
    })
  })

  # Status
  output$status_ui <- renderUI({
    if (input$run_btn == 0) {
      div(class = "status-bar", "Awaiting simulation run...")
    } else {
      req(sim_result())
      div(class = "status-bar done", "\u2713  Simulation complete")
    }
  })

  # Output info: theta_pop matrix, theta_hat matrix, R2 ranges
  output$output_info_ui <- renderUI({
    req(sim_result())
    res <- sim_result()

    alpha_pct_lo <- round(res$alpha / 2 * 100, 1)
    alpha_pct_hi <- round((1 - res$alpha / 2) * 100, 1)

    div(class = "output-info-row",
      # theta_pop matrix
      mat_html(res$Rpop,
        "\u03B8<sub>pop</sub> &mdash; Population Correlation"),

      # theta_hat matrix
      mat_html(res$Rhat,
        "\u03B8&#x0302; &mdash; Sample Correlation"),

      # R2 range from MCMC (min / max)
      r2_range_html(
        "R\u00B2 range &mdash; MCMC<br><small style='font-weight:normal;color:#777;'>[min, max]</small>",
        res$r2_mcmc_lo, res$r2_mcmc_hi
      ),

      # R2 range from Bootstrap (alpha/2 and 1-alpha/2 quantiles)
      r2_range_html(
        paste0("R\u00B2 CI &mdash; Bootstrap<br><small style='font-weight:normal;color:#777;'>[",
               alpha_pct_lo, "%, ", alpha_pct_hi, "%]</small>"),
        res$r2_boot_lo, res$r2_boot_hi
      )
    )
  })

  # Info chips
  output$chips_ui <- renderUI({
    req(sim_result())
    res <- sim_result()
    tagList(
      div(class = "chip", HTML("R\u00B2 true: "),  tags$span(round(res$ktrue, 3))),
      div(class = "chip", HTML("R\u00B2 hat: "),   tags$span(round(res$khat.point, 3))),
      div(class = "chip", "n = ",                   tags$span(input$ss_j)),
      div(class = "chip", "perturbations = ",       tags$span(res$nchain))
    )
  })

  # 3D Plotly
  output$scatter3d <- renderPlotly({
    req(sim_result())
    res <- sim_result()

    Rpop       <- res$Rpop
    Rhat       <- res$Rhat
    Rmcmc      <- res$Rmcmc
    Rboot      <- res$Rboot
    khat.chain <- res$khat.chain
    khat.boot  <- res$khat.boot
    fullspace  <- res$fullspace

    # Combined range for shared colorbar
    all_khat <- c(khat.chain, khat.boot)
    khat_min <- min(all_khat, na.rm = TRUE)
    khat_max <- max(all_khat, na.rm = TRUE)

    # Viridis colorscale
    viridis_scale <- list(
      list(0,     "#440154"), list(0.111, "#482878"),
      list(0.222, "#3e4989"), list(0.333, "#31688e"),
      list(0.444, "#26828e"), list(0.556, "#1f9e89"),
      list(0.667, "#35b779"), list(0.778, "#6ece58"),
      list(0.889, "#b5de2b"), list(1,     "#fde725")
    )

    plot_ly(
      x = Rpop[1, 2], y = Rpop[1, 3], z = Rpop[2, 3],
      type = "scatter3d", mode = "markers",
      marker = list(color = "black", size = 8),
      name  = "\u03B8 (population)"
    ) %>%

    add_markers(
      x = Rhat[1, 2], y = Rhat[1, 3], z = Rhat[2, 3],
      marker = list(color = "red", size = 8),
      name = "\u03B8\u0302 (sample)"
    ) %>%

    add_markers(
      x = ~Rmcmc[1, 2, ], y = ~Rmcmc[1, 3, ], z = ~Rmcmc[2, 3, ],
      marker = list(
        color      = khat.chain,
        colorscale = viridis_scale,
        cmin       = khat_min,
        cmax       = khat_max,
        opacity    = 0.3,
        size       = 3,
        showscale  = FALSE
      ),
      name = "\u03B8\u0302 MCMC"
    ) %>%

    add_markers(
      x = ~Rboot[1, 2, ], y = ~Rboot[1, 3, ], z = ~Rboot[2, 3, ],
      marker = list(
        color      = khat.boot,
        colorscale = viridis_scale,
        cmin       = khat_min,
        cmax       = khat_max,
        opacity    = 0.3,
        size       = 3,
        showscale  = TRUE,
        colorbar   = list(
          title        = list(text = "R\u00B2 (R-squared)", side = "right",
                              font = list(size = 13)),
          tickfont     = list(size = 12),
          len          = 0.6,
          thickness    = 18,
          outlinecolor = "#cccccc",
          outlinewidth = 1
        )
      ),
      name = "\u03B8\u0302 Bootstrap"
    ) %>%

    add_markers(
      x = ~fullspace[, 1, 2],
      y = ~fullspace[, 1, 3],
      z = ~fullspace[, 2, 3],
      marker = list(color = "gray", opacity = 0.2, size = 2),
      name   = "Full space"
    ) %>%

    layout(
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      font  = list(family = "Georgia, serif", color = "#1a1a1a", size = 13),
      scene = list(
        bgcolor = "white",
        xaxis = list(title = "r(Y, X1)",  gridcolor = "#dddddd"),
        yaxis = list(title = "r(Y, X2)",  gridcolor = "#dddddd"),
        zaxis = list(title = "r(X1, X2)", gridcolor = "#dddddd")
      ),
      legend = list(
        bgcolor     = "rgba(255,255,255,0.9)",
        bordercolor = "#cccccc",
        borderwidth = 1,
        font        = list(color = "#1a1a1a", size = 13)
      ),
      margin = list(l = 0, r = 0, t = 10, b = 0)
    )
  })

}

# -----------------------------------------------------------------------------
# Launch
# -----------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
