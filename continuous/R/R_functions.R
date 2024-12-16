# last update: 2024-12-06
# Andrea D. & Andrea B.

# Function to derive data for OR and HR
data_or_hr <- function(model,
                       covariate,
                       data,
                       n.knots,
                       referent.value){
  
  x <- eval(substitute(covariate), data, parent.frame())
  name.x <- deparse(substitute(covariate))
  grid.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 500)
  
  coefs.regex <- paste0("rcs\\(", 
                        name.x, 
                        ".*\\)",
                        name.x,
                        "'*$")
  coefs.rcs <- grep(coefs.regex, names(coef(model)))
  names(coefs.rcs) <- names(coef(model))[coefs.rcs]

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
                        Exp = TRUE)
  
  data.out <- as.data.frame(
    cbind(grid.x,
          log.or)
  )
  colnames(data.out) <- c(name.x, colnames(log.or))
  attr(data.out, "knots.location") <- knots.location
  
  p.overall <- base::format.pval(aod::wald.test(vcov(model), coef(model), 
                                                Terms = coefs.rcs)$result$chi2["P"],
                                 eps = 0.0001, digits = 3, scientific = FALSE)
  p.nonlinearity <- base::format.pval(aod::wald.test(vcov(model), coef(model), 
                                                     Terms = coefs.rcs[-1])$result$chi2["P"],
                                      eps = 0.0001, digits = 3, scientific = FALSE)
  cat(sprintf("\nP-value for overall association: %s\nP-value for non-linearity: %s\n\n", 
              p.overall, p.nonlinearity))
  
  data.out
}

# Function to derive model prediction
data_probability <- function(model,
                             covariate,
                             data,
                             time){
  
  isCOX <- inherits(model, "coxph")
  
  x <- eval(substitute(covariate), data, parent.frame())
  name.x <- deparse(substitute(covariate))
  grid.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 500)
  
  newdata <- as.data.frame(grid.x)
  colnames(newdata) <- name.x
  
  if (!isCOX) {
    prob <- Epi::ci.pred(model,
                         newdata)
    
    data.out <- as.data.frame(
      cbind(grid.x,
            prob)
    )
    colnames(data.out) <- c(name.x, colnames(prob))
  } else {
    prob <- survival:::summary.survfit(
      survival::survfit(model, 
                        newdata,
                        se.fit = TRUE, 
                        conf.type = "log-log"),
      time = time)
    
    data.out <- data.frame(
      rep(time, each = length(grid.x)),
      grid.x,
      1-as.vector(t(prob$surv)),
      1-as.vector(t(prob$lower)),
      1-as.vector(t(prob$upper))
    )
    colnames(data.out) <- c("time", name.x, "Estimate", "2.5%", "97.5%")
  }
  
  data.out
}
