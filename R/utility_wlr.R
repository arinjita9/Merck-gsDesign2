#  Copyright (c) 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the gsDesign2 program.
#
#  gsDesign2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Create npsurvSS arm object
#'
#' @inheritParams gs_info_ahr
#' @param total_time Total analysis time.
#'
#' @return A list of the two arms.
#'
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Validate if there is only one stratum.
#'    \item Calculate the accrual duration.
#'    \item calculate the accrual intervals.
#'    \item Calculate the accrual parameter as the proportion of enrollment rate*duration.
#'    \item Set cure proportion to zero.
#'    \item set survival intervals and shape.
#'    \item Set fail rate in fail_rate to the Weibull scale parameter for the survival distribution in the arm 0.
#'    \item Set the multiplication of hazard ratio and fail rate to the Weibull scale parameter
#'    for the survival distribution in the arm 1.
#'    \item Set the shape parameter to one as the exponential distribution for
#'    shape parameter for the loss to follow-up distribution
#'    \item Set the scale parameter to one as the scale parameter for the loss to follow-up
#'     distribution since the exponential distribution is supported only
#'    \item Create arm 0 using \code{npsurvSS::create_arm()} using the parameters for arm 0.
#'    \item Create arm 1 using \code{npsurvSS::create_arm()} using the parameters for arm 1.
#'    \item Set the class of the two arms.
#'    \item Return a list of the two arms.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @export
#'
#' @examples
#' enroll_rate <- define_enroll_rate(
#'   duration = c(2, 2, 10),
#'   rate = c(3, 6, 9)
#' )
#'
#' fail_rate <- define_fail_rate(
#'   duration = c(3, 100),
#'   fail_rate = log(2) / c(9, 18),
#'   hr = c(.9, .6),
#'   dropout_rate = .001
#' )
#'
#' gs_create_arm(enroll_rate, fail_rate, ratio = 1)
# Copyright (c) 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
# All rights reserved.

#' Create npsurvSS arm object for stratified or unstratified designs
#'
#' @inheritParams gs_info_ahr
#' @param total_time Total analysis time.
#'
#' @return A list of lists of the two arms per stratum.
#' @export
gs_create_arm <- function(enroll_rate, fail_rate, ratio, total_time = 1e6) {
  strata <- unique(enroll_rate$stratum)
  
  # Map out arms for each stratum
  strata_arms <- lapply(strata, function(str) {
    str_enroll <- enroll_rate[enroll_rate$stratum == str, ]
    str_fail <- fail_rate[fail_rate$stratum == str, ]
    
    accr_time <- sum(str_enroll$duration)
    accr_interval <- cumsum(str_enroll$duration)
    accr_param <- str_enroll$rate * str_enroll$duration / sum(str_enroll$rate * str_enroll$duration)
    
    surv_cure <- 0 
    surv_interval <- c(0, c(cumsum(head(str_fail$duration, -1)), Inf))
    surv_shape <- 1 
    surv_scale0 <- str_fail$fail_rate
    surv_scale1 <- str_fail$hr * str_fail$fail_rate
    loss_shape <- 1 
    loss_scale <- str_fail$dropout_rate[1]
    
    # Control Group
    arm0 <- npsurvSS::create_arm(
      size = 1, accr_time = accr_time, accr_dist = "pieceuni",
      accr_interval = accr_interval, accr_param = accr_param,
      surv_cure = surv_cure, surv_interval = surv_interval,
      surv_shape = surv_shape, surv_scale = surv_scale0,
      loss_shape = loss_shape, loss_scale = loss_scale, total_time = total_time
    )
    
    # Experimental Group
    arm1 <- npsurvSS::create_arm(
      size = ratio, accr_time = accr_time, accr_dist = "pieceuni",
      accr_interval = accr_interval, accr_param = accr_param,
      surv_cure = surv_cure, surv_interval = surv_interval,
      surv_shape = surv_shape, surv_scale = surv_scale1,
      loss_shape = loss_shape, loss_scale = loss_scale, total_time = total_time
    )
    
    class(arm0) <- c("list", "arm")
    class(arm1) <- c("list", "arm")
    
    list(arm0 = arm0, arm1 = arm1)
  })
  
  names(strata_arms) <- strata
  return(strata_arms)
}

#' Stratified Weight Log-Rank Delta (Mean) Computation
#' @export
gs_delta_wlr <- function(strata_arms, tmax = NULL, weight = "logrank", approx = "asymptotic", normalization = FALSE) {
  
  # Calculate delta per stratum and sum them up (weighted by stratum size/proportion if applicable)
  delta_each <- sapply(strata_arms, function(arm_pair) {
    arm0 <- arm_pair$arm0
    arm1 <- arm_pair$arm1
    
    if (is.null(tmax)) {
      tmax <- arm0$total_time
    }
    
    p1 <- arm1$size / (arm0$size + arm1$size)
    p0 <- 1 - p1
    
    if (identical(weight, "logrank")) {
      weight_fun <- wlr_weight_1
    } else {
      weight_fun <- switch(weight$method,
        "fh" = function(x, arm0, arm1) { wlr_weight_fh(x, arm0, arm1, rho = weight$param$rho, gamma = weight$param$gamma) },
        "mb" = function(x, arm0, arm1) { wlr_weight_mb(x, arm0, arm1, tau = weight$param$tau, w_max = weight$param$w_max) }
      )
    }
    
    if (approx == "event_driven") {
      if (sum(arm0$surv_shape != arm1$surv_shape) > 0 || length(unique(arm1$surv_scale / arm0$surv_scale)) > 1) {
        stop("gs_delta_wlr(): Hazard is not proportional over time.", call. = FALSE)
      }
      theta <- c(arm0$surv_shape * log(arm1$surv_scale / arm0$surv_scale))[1]
      nu <- p0 * prob_event(arm0, tmax = tmax) + p1 * prob_event(arm1, tmax = tmax)
      return(theta * p0 * p1 * nu)
      
    } else if (approx == "asymptotic") {
      return(stats::integrate(function(x) {
        term0 <- p0 * prob_risk(arm0, x, tmax)
        term1 <- p1 * prob_risk(arm1, x, tmax)
        term <- (term0 * term1) / (term0 + term1)
        term <- ifelse(is.na(term), 0, term)
        weight_fun(x, arm0, arm1) * term * (npsurvSS::hsurv(x, arm1) - npsurvSS::hsurv(x, arm0))
      }, lower = 0, upper = tmax, rel.tol = 1e-5)$value)
      
    } else if (approx == "generalized_schoenfeld") {
      return(stats::integrate(function(x) {
        log_hr_ratio <- if (normalization) { 1 } else { log(npsurvSS::hsurv(x, arm1) / npsurvSS::hsurv(x, arm0)) }
        weight_fun(x, arm0, arm1) * log_hr_ratio * p0 * prob_risk(arm0, x, tmax) * p1 * prob_risk(arm1, x, tmax) / 
          (p0 * prob_risk(arm0, x, tmax) + p1 * prob_risk(arm1, x, tmax))^2 * (p0 * dens_event(arm0, x, tmax) + p1 * dens_event(arm1, x, tmax))
      }, lower = 0, upper = tmax)$value)
    } else {
      stop("gs_delta_wlr(): Please specify a valid approximation for the mean.", call. = FALSE)
    }
  })
  
  # Total delta is the sum of delta across strata
  return(sum(delta_each))
}

#' Stratified Weight Log-Rank Sigma2 (Variance) Computation
#' @export
gs_sigma2_wlr <- function(strata_arms, tmax = NULL, weight = "logrank", approx = "asymptotic") {
  
  sigma2_each <- sapply(strata_arms, function(arm_pair) {
    arm0 <- arm_pair$arm0
    arm1 <- arm_pair$arm1
    
    if (is.null(tmax)) {
      tmax <- arm0$total_time
    }
    
    p1 <- arm1$size / (arm0$size + arm1$size)
    p0 <- 1 - p1
    
    enroll_duration <- arm0$accr_interval |> diff()
    enroll_total_duration <- arm0$accr_interval |> max()
    n_enroll_piece <- length(enroll_duration)
    enroll_relative_rate <- rep(-10, n_enroll_piece)
    
    for (s in 1:n_enroll_piece) {
      enroll_relative_rate[s] <- arm0$accr_param[s] / arm0$accr_param[n_enroll_piece] * enroll_duration[n_enroll_piece] / enroll_duration[s]
    }
    
    if (identical(weight, "logrank")) {
      weight_fun <- wlr_weight_1
    } else {
      weight_fun <- switch(weight$method,
        "fh" = function(x, arm0, arm1) { wlr_weight_fh(x, arm0, arm1, rho = weight$param$rho, gamma = weight$param$gamma, tau = if("tau" %in% names(weight$param)){weight$param$tau}else{NULL}) },
        "mb" = function(x, arm0, arm1) { wlr_weight_mb(x, arm0, arm1, tau = weight$param$tau, w_max = weight$param$w_max) }
      )
    }
    
    if (approx == "event_driven") {
      nu <- p0 * prob_event(arm0, tmax = tmax) + p1 * prob_event(arm1, tmax = tmax)
      return(p0 * p1 * nu)
      
    } else if (approx %in% c("asymptotic", "generalized_schoenfeld")) {
      if (tmax < enroll_total_duration) {
        arm0$accr_time <- tmax
        arm1$accr_time <- tmax
        arm0$accr_interval <- c(tmax, arm0$accr_interval)[which(c(tmax, arm0$accr_interval) <= tmax)] |> sort()
        arm1$accr_interval <- arm0$accr_interval
        truncated_enroll_duration <- diff(arm0$accr_interval)
        arm0$accr_param <- truncated_enroll_duration * enroll_relative_rate[1:length(truncated_enroll_duration)] / sum(truncated_enroll_duration * enroll_relative_rate[1:length(truncated_enroll_duration)])
        arm1$accr_param <- arm0$accr_param
      }
      
      return(stats::integrate(function(x) {
        weight_fun(x, arm0, arm1)^2 * p0 * prob_risk(arm0, x, tmax) * p1 * prob_risk(arm1, x, tmax) / 
          (p0 * prob_risk(arm0, x, tmax) + p1 * prob_risk(arm1, x, tmax))^2 * (p0 * dens_event(arm0, x, tmax) + p1 * dens_event(arm1, x, tmax))
      }, lower = 0, upper = tmax)$value)
    } else {
      stop("gs_sigma2_wlr(): Please specify a valid approximation for the mean.", call. = FALSE)
    }
  })
  
  # Total variance is the sum of variances across strata
  return(sum(sigma2_each))
}


#' For a subject in the provided arm, calculate the probability he or
#' she is observed to be at risk at `time=teval` after enrollment.
#' @noRd
prob_risk <- function(arm, teval, tmax) {
  if (is.null(tmax)) {
    tmax <- arm$total_time
  }

  npsurvSS::psurv(teval, arm, lower.tail = FALSE) *
    npsurvSS::ploss(teval, arm, lower.tail = FALSE) *
    npsurvSS::paccr(pmin(arm$accr_time, tmax - teval), arm)
}

#' For a subject in the provided arm, calculate the density of event
#' at `time=teval` after enrollment.
#' @noRd
dens_event <- function(arm, teval, tmax = NULL) {
  if (is.null(tmax)) {
    tmax <- arm$total_time
  }

  npsurvSS::dsurv(teval, arm) *
    npsurvSS::ploss(teval, arm, lower.tail = FALSE) *
    npsurvSS::paccr(pmin(arm$accr_time, tmax - teval), arm)
}

#' For a subject in the provided arm, calculate the probability he or
#' she is observed to have experienced an event by `time=teval` after enrollment.
#' @noRd
prob_event <- function(arm, tmin = 0, tmax = arm$total_time) {
  UseMethod("prob_event", arm)
}

#' prob_event for arm of class `arm`
#' @noRd
prob_event.arm <- function(arm, tmin = 0, tmax = arm$total_time) {
  l <- length(tmax)
  if (l == 1) {
    return(stats::integrate(function(x) dens_event(arm, x, tmax = tmax), lower = tmin, upper = tmax)$value)
  } else {
    if (length(tmin) == 1) {
      tmin <- rep(tmin, l)
    }
    return(sapply(seq(l), function(i) prob_event(arm, tmin[i], tmax[i])))
  }
}

#' @noRd
#' @examples
#' enroll_rate <- define_enroll_rate(
#'   duration = c(2, 2, 10),
#'   rate = c(3, 6, 9)
#' )
#'
#' fail_rate <- define_fail_rate(
#'   duration = c(3, 100),
#'   fail_rate = log(2) / c(9, 18),
#'   hr = c(.9, .6),
#'   dropout_rate = .001
#' )
#'
#' x <- gs_create_arm(enroll_rate, fail_rate, ratio = 1)
#' arm0 <- x$arm0
#' arm1 <- x$arm1
#' gs_delta_wlr(arm0, arm1, tmax = 10, weight = "logrank")
#' gs_delta_wlr(arm0, arm1, tmax = 10, weight = list(method = "mb", param = list(tau = Inf, w_max = 2)))

#' @noRd
almost_equal <- function(x, k, tol = .Machine$double.eps^0.5) {
  abs(x - k) < tol
}
