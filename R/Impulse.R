
#NPARAM=6
fit_impulse=function(vec,timepoints,n_iters = 100){
      NPARAM=6
      mat1 = c(0,abs(vec[2:length(vec)] - vec[1:(length(vec) - 1)]) * vec[2:length(vec)]) # abs(v(k) - v(k-1) * v(k))
      mn_beta = 0
      mx_beta = 10
      middle_ind = which(timepoints > min(timepoints) & timepoints < max(timepoints)) #position between min and max time
      min_ind = which(timepoints == min(timepoints))[1] #timepoint min position
      max_ind = which(timepoints == max(timepoints))[1] #timepoint max position
      mx_tm = max(timepoints) #max time
      mn_tm = min(timepoints) #min time
      tmp_ind = which.max(mat1) # peak timepoint position
      peaktime = timepoints[tmp_ind] #peak time point
      peak = min(which(timepoints == peaktime)) 
      beta1 = abs(log(abs((vec[peak] - vec[min_ind])/(peaktime - timepoints[min_ind]))));
    
      ### set beta to 1 if calculcated value is infinite
      if (is.finite(beta1) == FALSE || is.na(beta1)) { beta1 = 1 }

      ### define start value theta based on the gene's data
      orig_theta = c(beta1,
         vec[min_ind],                    # h0
         vec[peak],                       # h1
         vec[max_ind],                    # h2
         (peaktime - timepoints[min_ind])/2,            # t1
         peaktime + (timepoints[max_ind] - peaktime)/2) # t2
      names(orig_theta) = c("beta","h0","h1","h2","t1","t2")
      theta = orig_theta
    
      ### replace theta estimates of 0 by small value, because optimization
      ### function has problems with 0 as input
      rand_num <- runif(length(which(theta == 0)),0,0.0001)
      theta[names(which(theta == 0))] <- rand_num
      names(theta) = c("beta","h0","h1","h2","t1","t2")
    
    toll <- cbind(theta,matrix(start_val,7,3)[1:NPARAM,][,1:3])
    
    # fitting via quasi-Newton method (optim, "BFGS")
    fitout1=apply(toll,2,function(x){unlist(optim(x, two_impulses, x_vec = timepoints,
                                                  y_vec = as.numeric(vec),method = "BFGS" )[c("par","value")])})
    # again fitting via PORT routines (nlminb)
    fitout2=apply(toll,2,function(x){unlist(nlminb(x, two_impulses, x_vec = timepoints,
                                                   y_vec = as.numeric(vec) )[c("par","objective")])})

    fmin_outs <- cbind(fitout1,fitout2)
    
    if (is.null(dim(fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])]))) {
          pvec_and_SSE = fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])]
       } else {
          ### if two or more randomization have the same minimum Sum of
          ### Squared Errors (SSE), choose the first one
          pvec_and_SSE = fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])][,1]
       }
}

################################################################################

#' Impulse model value prediction
#'
#' Calculates impulse model values for given timepoints and predicted
#' impulse parameters.
#' @aliases calc_impulse
#' @param theta numerical vector of impulse parameters with the order
#' beta, h0, h1, h2, t1, t2.
#' @param timepoints numercial vector of time point(s).
#' @return The predicted impulse model values for the given time point(s).
#' @seealso \code{\link{impulse_DE}}, \code{\link{plot_impulse}}.
#' @author Jil Sander
#' @references Chechik, G. and Koller, D. (2009) Timing of Gene Expression
#' Responses to Envi-ronmental Changes. J. Comput. Biol., 16, 279-290.
#' @examples
#' #' theta vector in the order beta, h0, h1, h2, t1, t2
#' theta <- c(9.9, 14.7, 17.0, 16.9, -0.1, 37.0)
#' #' time points
#' timepoints <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' #' calculate impulse values
#' impulse_values <- calc_impulse(theta, timepoints)
#' @export
calc_impulse <- function(theta,timepoints){
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  res = NULL
   for (i in 1:length(timepoints)) {
    res[i] = (1/h1) * (h0 + (h1 - h0) * (1/(1 + exp(-beta1*(timepoints[i] - t1))))) *
        (h2 + (h1 - h2) * (1/(1 + exp(beta1*(timepoints[i] - t2)))))
  }
  res = unlist(res)
  return(res)
}

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++     Calculation of input for Optimization    ++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### specifies the function of squared errors which is minimized to fit
### parameters
two_impulses <- function(theta,x_vec,y_vec){
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
   f = sum((unlist(lapply(x_vec, function(x) {(1/h1) * (h0 + (h1 - h0) *
        (1/(1 + exp(-beta1*(x - t1))))) * (h2 + (h1 - h2) *
        (1/(1 + exp(beta1*(x - t2)))))})) - y_vec)^2)
   return(f)
}
