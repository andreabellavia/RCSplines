# Last update: 2025-04-09
# Andrea Discacciati & Andrea Bellavia

# Function to derive data for OR and HR
data_or_hr <- function(model,                  # required
                       covariate,              # required
                       referent.value,         # required
                       data,                   # optional, default = model$call$data
                       conf.level = 0.95,      # optional, default = 0.95
                       length.grid.x = 100){   # optional, default = 100
  
  coef.glmgee <- function(x) { 
    setNames(as.vector(x$coefficients), rownames(x$coefficients))
  }
  
  if (isTRUE(base::missing(data))) 
    data <- eval(model$call$data, parent.frame())
  x <- eval(substitute(covariate), data, parent.frame())
  name.x <- deparse(substitute(covariate))
  grid.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = length.grid.x)
  
  coefs.regex <- paste0("rcs\\(", 
                        name.x, 
                        ".*\\)",
                        name.x,
                        "'*$")
  coefs.rcs <- grep(coefs.regex, names(coef(model)))
  names(coefs.rcs) <- names(coef(model))[coefs.rcs]
  n.knots <- length(coefs.rcs)+1L
  
  # store and print p-values
  p.overall <- base::format.pval(aod::wald.test(vcov(model), coef(model), 
                                                Terms = coefs.rcs)$result$chi2["P"],
                                 eps = 0.0001, digits = 3, scientific = FALSE)
  p.nonlinearity <- base::format.pval(aod::wald.test(vcov(model), coef(model), 
                                                     Terms = coefs.rcs[-1])$result$chi2["P"],
                                      eps = 0.0001, digits = 3, scientific = FALSE)
  cat(sprintf("\nP-value for overall association: %s\nP-value for non-linearity: %s\n\n", 
              p.overall, p.nonlinearity))
  
  # data.frame
  knots.location <- Hmisc::rcspline.eval(x,
                                         nk = n.knots,
                                         knots.only = TRUE)
  
  grid.rcs <- Hmisc::rcspline.eval(grid.x,
                                   knots = knots.location,
                                   inclx = TRUE)
  ref.rcs <-   Hmisc::rcspline.eval(referent.value,
                                    knots = knots.location,
                                    inclx = TRUE)
  colnames(ref.rcs) <- colnames(grid.rcs) <- names(coefs.rcs)
  
  contrast.matrix <- sweep(grid.rcs, 2, ref.rcs, FUN = "-")
  
  log.or <- Epi::ci.lin(model,
                        ctr.mat = contrast.matrix,
                        subset = coefs.rcs,
                        Exp = TRUE,
                        alpha = 1-conf.level)
  
  data.out <- as.data.frame(
    cbind(grid.x,
          log.or)
  )
  colnames(data.out) <- c(name.x, colnames(log.or))
  attr(data.out, "knots.location") <- knots.location
  attr(data.out, "p.overall") <- p.overall
  attr(data.out, "p.nonlinearity") <- p.nonlinearity
  
  data.out
}

# Function to derive predicted event probabilities
data_probability <- function(model,                   # required
                             covariate,               # required
                             time,                    # required with coxph models
                             newdata,                 # required if any extra covariate
                             data,                    # optional, default = model$call$data
                             conf.level = 0.95,       # optional, default = 0.95
                             length.grid.x = 100){    # optional, default = 100
  
  if (isTRUE(!inherits(model, "glm") & !inherits(model, "glmgee") & !inherits(model, "coxph")))
    stop("Objects of class glm, glmgee, or coxph are supported")
  
  isCOX <- inherits(model, "coxph")
  
  if (isTRUE(base::missing(data))) 
    data <- eval(model$call$data, parent.frame())
  x <- eval(substitute(covariate), data, parent.frame())
  name.x <- deparse(substitute(covariate))
  grid.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = length.grid.x)
  newdata2 <- as.data.frame(grid.x)
  colnames(newdata2) <- name.x
  
  if (isFALSE(base::missing(newdata))) { # if newdata is not missing, full join grid.x with newdata
    newdata2 <- merge(newdata2, newdata)
    colnames(newdata2) <- c(name.x, colnames(newdata))
  }
  
  if (isFALSE(isCOX)) { # logit
    ci.pred <- function (obj, newdata, alpha, varest = "robust") { # based on Epi::ci.pred
      predict.glmgee <- function(...) {
        pred.mat <- glmtoolbox:::predict.glmgee(...)
        asplit(pred.mat, 2)
      }
      zz <- predict(obj, newdata = newdata, se.fit = TRUE, 
                    type = "link", varest = varest)
      zz <- cbind(zz$fit, zz$se.fit) %*% Epi::ci.mat(alpha = alpha)
      return(obj$family$linkinv(zz))
    }
    
    prob <- ci.pred(model,
                    newdata2,
                    alpha = 1-conf.level)
    
    data.out <- as.data.frame(
      cbind(newdata2,
            prob)
    )
    colnames(data.out) <- c(colnames(newdata2), colnames(prob))
  } else { # cox
    as.data.frame.summary.survfit <- function(x) {
      alpha <- 1-x$conf.int
      df.out <- data.frame(
        1-as.vector(t(x$surv)),
        1-as.vector(t(x$lower)),
        1-as.vector(t(x$upper))
      )
      colnames.ci.mat <- c("Estimate", # from Epi:::ci.mat
                           paste(formatC(100 * alpha/2, format = "f", digits = 1), "%", sep = ""), 
                           paste(formatC(100 * (1-alpha/2), format = "f", digits = 1), "%", sep = ""))
      colnames(df.out) <- colnames.ci.mat
      
      df.out
    }
    
    prob <- survival:::summary.survfit(
      survival::survfit(model, 
                        newdata2,
                        se.fit = TRUE, 
                        conf.type = "log-log",
                        conf.int = conf.level),
      time = time)
    
    data.out <- cbind(
      time = rep(time, each = nrow(newdata2)),
      newdata2,
      as.data.frame(prob)
    )
  }
  
  data.out
}
