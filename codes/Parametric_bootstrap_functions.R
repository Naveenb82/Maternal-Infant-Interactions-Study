#######Claude Sonnet 4.6 was used for the development of these wrapper functions
#######These functions implement parametric bootstrapping fixed effect confidence intervals
#######using lme4::bootmer(). The generated functions were reviewed by the authors.

boot_fixef_ci <- function(model,
                          terms = NULL,
                          nsim = 1000,
                          level = 0.95,
                          type = "parametric",
                          use_u = TRUE,
                          parallel = c("no", "multicore", "snow"),
                          ncpus = 24,
                          seed = NULL) {
  
  parallel <- match.arg(parallel)
  if (!is.null(seed)) set.seed(seed)
  
  fe <- fixef(model)
  term_names <- names(fe)
  
  
  if (is.null(terms)) terms <- term_names
  
  
  missing_terms <- setdiff(terms, term_names)
  if (length(missing_terms) > 0) {
    stop("These terms are not in fixef(model): ",
         paste(missing_terms, collapse = ", "),
         "\nAvailable terms: ", paste(term_names, collapse = ", "))
  }
  
  
  fun <- function(fit) {
    fx <- fixef(fit)
    as.numeric(fx[terms])
  }
  
  
  bb <- bootMer(
    model,
    FUN = fun,
    nsim = nsim,
    type = type,
    use.u = use_u,
    parallel = parallel,
    ncpus = ncpus
  )
  
  
  alpha <- (1 - level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  
  boot_mat <- bb$t
  if (is.null(dim(boot_mat))) boot_mat <- matrix(boot_mat, ncol = length(terms))
  colnames(boot_mat) <- terms
  
  q <- t(apply(boot_mat, 2, quantile, probs = probs, na.rm = TRUE))
  colnames(q) <- c("lo", "mid", "hi")
  
  est <- fe[terms]
  
  out <- data.frame(
    term = terms,
    estimate = as.numeric(est),
    lo = q[, "lo"],
    mid = q[, "mid"],
    hi = q[, "hi"],
    row.names = NULL
  )
  
  
  out$OR_est <- exp(out$estimate)
  out$OR_lo  <- exp(out$lo)
  out$OR_mid <- exp(out$mid)
  out$OR_hi  <- exp(out$hi)
  
  
  attr(out, "boot_draws") <- boot_mat
  attr(out, "bootMer_object") <- bb
  out
}

###############################

boot_glmm_marginal_plot <- function(
    mod,
    data,
    age_var,
    group_vars,              
    facet_var = NULL,        
    color_var = NULL,        
    trim = c(.05, .95),      
    trim_by_group = FALSE,   
    condition = NULL,
    n_grid = 60,             
    n_re = 2000,             
    nsim = 500,              
    seed = 123,
    parallel = TRUE,
    ncpus = max(1, parallel::detectCores() - 1)
) {
  stopifnot(inherits(mod, "glmerMod"))
  stopifnot(length(group_vars) >= 1)
  if (is.null(color_var)) color_var <- group_vars[1]
  
  
  get_levels <- function(v) {
    if (is.factor(data[[v]])) levels(data[[v]]) else sort(unique(data[[v]]))
  }
  
  
  if (!trim_by_group) {
    
    age_grid <- sort(unique(data[[age_var]]))
    
    grid_list <- c(
      setNames(list(age_grid), age_var),
      lapply(group_vars, get_levels)
    )
    names(grid_list)[-1] <- group_vars
    
    grid <- expand.grid(
      grid_list,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    if(!is.null(condition)) {
      for(v in names(condition)) {
        grid[[v]] <- condition[[v]]
      }
    }
    
  } else {
    
    dsub <- data[, c(age_var, group_vars), drop = FALSE]
    
    
    for (v in group_vars) {
      if (is.factor(data[[v]])) {
        dsub[[v]] <- factor(dsub[[v]], levels = levels(data[[v]]))
      }
    }
    
    split_list <- split(dsub, interaction(dsub[group_vars], drop = TRUE))
    grid_list <- lapply(split_list, function(df) {
      if (nrow(df) < 5) return(NULL)
      rng <- quantile(df[[age_var]], probs = trim, na.rm = TRUE)
      age_grid <- seq(rng[1], rng[2], length.out = n_grid)
      g <- data.frame(setNames(list(age_grid), age_var))
      for (v in group_vars) g[[v]] <- df[[v]][1]
      g
    })
    grid <- do.call(rbind, grid_list)
    rownames(grid) <- NULL
  }
  
  
  mf <- model.frame(mod)
  all_vars <- all.vars(formula(mod))  
  resp <- all_vars[1]
  preds <- setdiff(all_vars, resp)
  
  missing_in_grid <- setdiff(preds, c(age_var, group_vars))
  
  re_terms <- names(lme4::getME(mod, "flist"))
  missing_in_grid <- setdiff(missing_in_grid, re_terms)
  
  for (v in missing_in_grid) {
    if (!v %in% names(mf)) next
    if (is.factor(mf[[v]])) {
      grid[[v]] <- factor(levels(mf[[v]])[1], levels = levels(mf[[v]]))
    } else if (is.logical(mf[[v]])) {
      grid[[v]] <- FALSE
    } else {
      grid[[v]] <- median(mf[[v]], na.rm = TRUE)
    }
  }
  
  
  for (v in names(grid)) {
    if (v %in% names(mf) && is.factor(mf[[v]])) {
      grid[[v]] <- factor(grid[[v]], levels = levels(mf[[v]]))
    }
  }
  
  
  boot_fun <- function(fit) {
    X <- model.matrix(delete.response(terms(fit)), grid)
    eta <- as.vector(X %*% lme4::fixef(fit))
    
    
    vc <- lme4::VarCorr(fit)
    
    sd_b <- as.numeric(attr(vc[[1]], "stddev"))[1]
    
    b_sim <- rnorm(n_re, 0, sd_b)
    vapply(eta, function(e) mean(plogis(e + b_sim)), numeric(1))
  }
  
  
  set.seed(seed)
  
  if (parallel) {
    RNGkind("L'Ecuyer-CMRG")
    cl <- parallel::makeCluster(ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("grid", "n_re"), envir = environment())
    parallel::clusterEvalQ(cl, library(lme4))
    
    boot_res <- lme4::bootMer(
      mod, FUN = boot_fun, nsim = nsim, type = "parametric", use.u = FALSE,
      parallel = "snow", ncpus = ncpus, cl = cl
    )
  } else {
    boot_res <- lme4::bootMer(
      mod, FUN = boot_fun, nsim = nsim, type = "parametric", use.u = FALSE
    )
  }
  
  
  ci <- apply(boot_res$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  plot_df <- cbind(
    grid,
    estimate  = colMeans(boot_res$t, na.rm = TRUE),
    conf.low  = ci[1, ],
    conf.high = ci[2, ]
  )
  
  
  p <- ggplot(
    plot_df,
    aes_string(x = age_var, y = "estimate",
               color = color_var, fill = color_var,
               group = color_var)
  ) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
    geom_line(linewidth = 1.1) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', size = 15, hjust = 0.5),
          axis.text.x = element_text(face = 'bold', size = 15, 
                                     hjust = 0.5),
          axis.text.y = element_text(face = 'bold', size = 15),
          axis.title = element_text(face = 'bold', size = 17),
          panel.grid = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 15, face = 'bold'),
          legend.key.size = unit(0.2, 'in'),
          legend.title = element_text(size = 15, face = 'bold'),
          panel.spacing = unit(2, "lines"))
  
  coord_cartesian(ylim = c(0, 1))
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }
  
  list(plot = p, plot_df = plot_df, boot = boot_res, grid = grid)
}

########################################################

boot_glmm_pointrange <- function(
    mod,
    data,
    focal_var,                   
    condition       = NULL,      
    n_re            = 1500,      
    nsim            = 1000,      
    seed            = 123,
    parallel        = TRUE,
    ncpus           = max(1, parallel::detectCores() - 1)
) {
  stopifnot(inherits(mod, "glmerMod"))
  
  
  focal_levels <- if (is.factor(data[[focal_var]])) {
    levels(data[[focal_var]])
  } else {
    sort(unique(data[[focal_var]]))
  }
  
  
  grid <- data.frame(
    setNames(list(factor(focal_levels, levels = focal_levels)), focal_var),
    stringsAsFactors = FALSE
  )
  
  
  mf        <- model.frame(mod)
  all_vars  <- all.vars(formula(mod))
  resp      <- all_vars[1]
  preds     <- setdiff(all_vars, resp)
  re_terms  <- names(lme4::getME(mod, "flist"))
  miss_vars <- setdiff(preds, c(focal_var, re_terms))
  
  
  for (v in miss_vars) {
    if (!v %in% names(mf)) next
    
    if (!is.null(condition) && v %in% names(condition)) {
      # User explicitly specified this value
      val <- condition[[v]]
      if (is.factor(mf[[v]])) {
        grid[[v]] <- factor(val, levels = levels(mf[[v]]))
      } else {
        grid[[v]] <- val
      }
      message(sprintf("  [condition] '%s' -> %s", v, as.character(val)))
      
    } else if (is.factor(mf[[v]])) {
      ref <- levels(mf[[v]])[1]
      grid[[v]] <- factor(ref, levels = levels(mf[[v]]))
      message(sprintf("  [default]   '%s' -> reference level: '%s'", v, ref))
      
    } else if (is.logical(mf[[v]])) {
      grid[[v]] <- FALSE
      message(sprintf("  [default]   '%s' -> FALSE", v))
      
    } else {
      # numeric: use condition value if supplied, else 0
      grid[[v]] <- 0
      message(sprintf("  [default]   '%s' -> 0", v))
    }
  }
  
  
  for (v in names(grid)) {
    if (v %in% names(mf) && is.factor(mf[[v]])) {
      grid[[v]] <- factor(grid[[v]], levels = levels(mf[[v]]))
    }
  }
  
  message(sprintf("\nPrediction grid (%d rows):", nrow(grid)))
  print(grid)
  
  
  boot_fun <- function(fit) {
    X   <- model.matrix(delete.response(terms(fit)), grid)
    eta <- as.vector(X %*% lme4::fixef(fit))
    
    vc    <- lme4::VarCorr(fit)
    sd_b  <- as.numeric(attr(vc[[1]], "stddev"))[1]
    b_sim <- rnorm(n_re, 0, sd_b)
    
    vapply(eta, function(e) mean(plogis(e + b_sim)), numeric(1))
  }
  
  
  set.seed(seed)
  
  if (parallel) {
    RNGkind("L'Ecuyer-CMRG")
    cl <- parallel::makeCluster(ncpus)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("grid", "n_re"), envir = environment())
    parallel::clusterEvalQ(cl, library(lme4))
    
    boot_res <- lme4::bootMer(
      mod, FUN = boot_fun, nsim = nsim,
      type = "parametric", use.u = FALSE,
      parallel = "snow", ncpus = ncpus, cl = cl
    )
  } else {
    boot_res <- lme4::bootMer(
      mod, FUN = boot_fun, nsim = nsim,
      type = "parametric", use.u = FALSE
    )
  }
  
  
  ci <- apply(boot_res$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  
  plot_df <- data.frame(
    x         = factor(focal_levels, levels = focal_levels),
    predicted = colMeans(boot_res$t, na.rm = TRUE),
    conf.low  = ci[1, ],
    conf.high = ci[2, ],
    row.names = NULL
  )
  
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_pointrange(
      data    = plot_df,
      mapping = ggplot2::aes(
        x     = x,
        y     = predicted,
        ymin  = conf.low,
        ymax  = conf.high,
        group = x,
        color = x,
        fill  = x
      ),
      lwd   = 0.9,
      size  = 2,
      shape = 16
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0.2)) +
    ggplot2::xlab(focal_var) +
    ggplot2::ylab("Predicted probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5),
      axis.text.y = ggplot2::element_text(face = "bold", size = 15),
      axis.title  = ggplot2::element_text(face = "bold", size = 17),
      panel.grid  = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  list(plot_df = plot_df, gg = p, boot = boot_res, grid = grid)
}


########################################################################