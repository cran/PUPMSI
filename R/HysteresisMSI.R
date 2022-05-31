#' @title Hysteresis Area, Brunauer Classification System
#' @description Hysteresis area evaluation via trapezoidal approximation.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import ggplot2
#' @import stats
#' @return the measurement of hysteresis, classification of isotherms, and graphical visualization for the observed values of moisture sorption isotherms.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples HysteresisMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Caurie, M. (2007) <doi:10.1111/j.1365-2621.2006.01203.x> Hysteresis phenomenon in foods. International Journal of Food Science and Technology, 42(1), 45-49.
#' @references Brunauer, S., et al. (1940) <doi:10.1021/ja01864a025> On a Theory of the van der Waals Adsorption of Gases. Journal of the American Chemical Society, 62(7), 1723-1732.
#' @references  Blahovec J., & Yanniotis S. (2009) <doi:10.1016/j.jfoodeng.2008.08.007> Modified classification of sorption isotherms. J Food Eng. 2009 Mar; 91 (1): 72-77
#' @export

HysteresisMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSI1 <- data.frame(x1, y1, y2)

  Hysteresis <- function(WaterAct, AdsorpM, DesorpM)
  {
    TrapezoidApproxAds1<-function (x1, y1)
    {
    if (missing(y1)) {
      if (length(x1) == 0)
        return(0)
      y1 <- x1
      x1 <- seq(along = x1)
      MSIH1 <- data.frame(x1, y1)
    }
    if (length(x1) == 0 && length(y1) == 0)
      return(0)
    if (!(is.numeric(x1) || is.complex(x1)) || !(is.numeric(y1) ||
                                                 is.complex(y1)))
      stop("Arguments 'x' and 'y' must be real or complex vectors.")
    m <- length(x1)
    if (length(y1) != m)
      stop("Arguments 'x', 'y' must be vectors of the same length.")
    if (m <= 1)
      return(0)
    xp <- c(x1, x1[m:1])
    yp <- c(numeric(m), y1[m:1])
    n <- 2 * m
    p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
    p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
    return(0.5 * (p1 - p2))
  }

    TrapezoidApproxDes1<-function (x1, y2)
    {
    if (missing(y2)) {
      if (length(x1) == 0)
        return(0)
      y2 <- x1
      x1 <- seq(along = x1)
      MSIH1 <- data.frame(x1, y2)
    }
    if (length(x1) == 0 && length(y2) == 0)
      return(0)
    if (!(is.numeric(x1) || is.complex(x1)) || !(is.numeric(y2) ||
                                                 is.complex(y2)))
      stop("Arguments 'x' and 'y' must be real or complex vectors.")
    m <- length(x1)
    if (length(y2) != m)
      stop("Arguments 'x', 'y' must be vectors of the same length.")
    if (m <= 1)
      return(0)
    xp <- c(x1, x1[m:1])
    yp <- c(numeric(m), y2[m:1])
    n <- 2 * m
    p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
    p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
    return(0.5 * (p1 - p2))
  }

    AreaAdsorption1 <- TrapezoidApproxAds1(x1, y1)
    AreaDesorption1 <- TrapezoidApproxDes1(x1, y2)
    AreaHysteresis1 <- abs(AreaDesorption1 - AreaAdsorption1)

    message("HYSTERESIS OF MOISTURE SOPRTION ISOTHERM")
    message("Area of Hysteresis")
    Hyster <- data.frame("Experimental", AreaHysteresis1)
    names(Hyster) <- c("Data","Area")
    print(Hyster, row.names = F, right = F)
  }

  IsoClass <- function(WaterAct, AdsorpM, DesorpM)
  {
    A.eqn2 <- y1 ~ ((x1/(q1 + (w1*x1))) + (x1/(q2 - (w2*x1))))
    D.eqn2 <- y2 ~ ((x1/(q1 + (w1*x1))) + (x1/(q2 - (w2*x1))))

    start2 <- data.frame(q1 = c(0,1), q2 = c(0,1), w1 = c(0, 1), w2 = c(0, 1))

    suppressWarnings(A.fit3 <- nls2::nls2(A.eqn2, data = MSI1, start = start2,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "default"))
    suppressWarnings(D.fit3 <- nls2::nls2(D.eqn2, data = MSI1, start = start2,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "default"))

    A.start2 <- list(q1 = stats::coef(A.fit3)[1],  q2 = stats::coef(A.fit3)[2], w1 = stats::coef(A.fit3)[3], w2 = stats::coef(A.fit3)[4])
    D.start2 <- list(q1 = stats::coef(D.fit3)[1],  q2 = stats::coef(D.fit3)[2], w1 = stats::coef(D.fit3)[3], w2 = stats::coef(D.fit3)[4])

    suppressWarnings(A.fit4 <- nls2::nls2(A.eqn2, data = MSI1, start = A.start2,
                       control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
    suppressWarnings(D.fit4 <- nls2::nls2(D.eqn2, data = MSI1, start = D.start2,
                       control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

    a1 <- stats::coef(A.fit4)[1]
    a2 <- stats::coef(A.fit4)[2]
    b1 <- stats::coef(A.fit4)[3]
    b2 <- stats::coef(A.fit4)[4]

    a3 <- stats::coef(D.fit4)[1]
    a4 <- stats::coef(D.fit4)[2]
    b3 <- stats::coef(D.fit4)[3]
    b4 <- stats::coef(D.fit4)[4]

    D10a <- (((a2)^2 * b1)-((a1)^2 * b2))/(a1 + a2)^2
    X4a <- (b1 - b2)/(a1 + a2)
    X3a <- (b1*b2)/(a1+a2)
    Rfia <- (1 - (X3a * ((2+X4a)/(D10a))))/((1+X4a)^2)

    D10b <- (((a4)^2 * b3)-((a3)^2 * b4))/(a3 + a4)^2
    X4b <- (b3 - b4)/(a3 + a4)
    X3b <- (b3*b4)/(a3+a4)
    Rfib <- (1 - (X3b * ((2+X4b)/(D10b))))/((1+X4b)^2)

    IsoclassAds <- {if (D10a > 0 & Rfia > 0 & X4a > 0) {'Type 1 (Langmuir-like)'}
    else if (D10a > 0 & Rfia < 0 & X4a <= 0.1 & X4a >= -0.1) {'Type II-a (GAB-like)'}
    else if (D10a > 0 & Rfia < 0 & X4a < -0.1) {'Type II-b (GAB-like)'}
    else if (D10a > 0 & Rfia < 0 & X4a > 0.1) {'Type II-c (GAB-like)'}
    else if  (D10a < 0 & Rfia > 0 & X4a < 0) {'Type III (Solution-like)'}
    else {'No data for the Isotherm Class'}
    }

    IsoclassDes <- {if (D10b > 0 & Rfib > 0 & X4b > 0) {'Type 1 (Langmuir-like)'}
    else if (D10b > 0 & Rfib < 0 & X4b <= 0.1 & X4b >= -0.1) {'Type II-a (GAB-like)'}
    else if (D10b > 0 & Rfib < 0 & X4b < -0.1) {'Type II-b (GAB-like)'}
    else if (D10b > 0 & Rfib < 0 & X4b > 0.1) {'Type II-c (GAB-like)'}
    else if  (D10b < 0 & Rfib > 0 & X4b < 0) {'Type III (Solution-like)'}
    else {'No data for the Isotherm Class'}
    }

    message("CLASSIFICATION OF ISOTHERM")
    message("Bruanauer's Classification System based on Blahovec et al. (2009) equation")
    Isotherm <- c("Adsorption", "Desorption")
    BCSys <- c(IsoclassAds, IsoclassDes)
    IsoType <- data.frame(Isotherm, BCSys)
    names(IsoType) <- c("Isotherm","Type")
    print(IsoType, row.names=F, right = F)

}

  HysPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    HysPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_ribbon(ggplot2::aes(x=x1, ymin=y1, ymax=y2, fill = "Area")) +
      ggplot2::geom_line(ggplot2::aes(x=x1, y=y1, colour = "Adsorption"), size = 1) +
      ggplot2::geom_line(ggplot2::aes(x=x1, y=y2, colour = "Desorption"), size = 1) +
      ggplot2::geom_point(ggplot2::aes(x= x1, y=y1, shape  = "Adsorption"), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x=x1, y=y2, shape  = "Desorption"), size = 2) +
      ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", fill = NA),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     legend.key = ggplot2::element_rect(fill = "white", colour = "black"),
                     legend.title = ggplot2::element_text(face = "bold", size = 10),
                     legend.position = "right",
                     legend.box.background = ggplot2::element_rect(),
                     legend.box.margin = ggplot2::margin(1,1,1,1),
                     legend.text = ggplot2::element_text(size = 8)) +
      ggplot2::scale_x_continuous(limits = c(0,1)) +
      ggplot2::scale_y_continuous(limits = c(minlim1,maxlim1)) +
      ggplot2::scale_shape_manual(name = "Isotherm", values = c(15, 19)) +
      ggplot2::scale_linetype_manual(name = "Isotherm", values = c(1,1)) +
      ggplot2::scale_color_manual(name = "Isotherm", values = c("blue", "red")) +
      ggplot2::scale_fill_manual(name = "Hysteresis", values = c("grey60"))

    suppressWarnings(print(HysPlot))

  }

  Hysteresis(WaterAct, AdsorpM, DesorpM)
  IsoClass(WaterAct, AdsorpM, DesorpM)
  HysPlot(WaterAct, AdsorpM, DesorpM)
}


