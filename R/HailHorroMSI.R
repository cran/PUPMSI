#' @title Hailwood-Horrobin (HH) Moisture Sorption Isotherm
#' @description Hailwood-Horrobin (HH) model is an example of multilayer surface sorption model, is suitable for analysis of experimental wood moisture sorption (WMS) isotherms.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Hailwood-Horrobin (HH) Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples HailHorroMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Hailwood, A. J., & Horrobin, S. (1946) <doi:10.1039/TF946420B084> Absorption of water by polymers: Analysis in terms of a simple model. Transactions of the Faraday Society, 42(0), B084-B092.
#' @export

HailHorroMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSIHailHorro1 <- data.frame(x1, y1, y2)

  A.eqn1 <- y1 ~ x1/((C) + (B*x1) - (A*(x1^2)))
  D.eqn1 <- y2 ~ x1/((C) + (B*x1) - (A*(x1^2)))

  start1 <- data.frame(A = c(-100, 100), B = c(-100,100), C = c(-100,100))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSIHailHorro1, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSIHailHorro1, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  A.start1 <- list(A = stats::coef(A.fit1)[1], B = stats::coef(A.fit1)[2], C = stats::coef(A.fit1)[3])
  D.start1 <- list(A = stats::coef(D.fit1)[1], B = stats::coef(D.fit1)[2], C = stats::coef(D.fit1)[3])

  suppressWarnings(A.fit2 <- nls2::nls2(A.eqn1, data = MSIHailHorro1, start = A.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit2 <- nls2::nls2(D.eqn1, data = MSIHailHorro1, start = D.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))


  y3 <- stats::predict(A.fit2)
  y4 <- stats::predict(D.fit2)

  Isotherm <- c("Adsorption", "Desorption")

  message("HAILWOOD-HORROBIN MOISTURE SORPTION MODEL")

  message("Hailwood-Horrobin Adsorption Isotherm Parameters")
  print(summary(A.fit2))

  message("Hailwood-Horrobin Desorption Isotherm Parameters")
  print(summary(D.fit2))

  message("Predicted Values of Moisture Content from Hailwood-Horrobin Moisture Sorption Model")
  HailHorroPredict <- data.frame(x1, y3, y4)
  names(HailHorroPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(HailHorroPredict, row.names = F, right = F)

  HailHorroMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Hailwood-Horrobin Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit2), stats::AIC(D.fit2))
    MSIBIC <- c(stats::BIC(A.fit2), stats::BIC(D.fit2))
    HailHorroCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(HailHorroCriterion) <- c("Isotherm","AIC", "BIC")
    print(HailHorroCriterion, row.names=F, right = F)

    message("Error Analysis for Hailwood-Horrobin Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1 - y3)^2)/(length(y1)-2)), (sqrt(sum(y2 - y4)^2)/(length(y2)-2)))
    HailHorroError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(HailHorroError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(HailHorroError, row.names=F, right = F)
  }

  HailHorroMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    l1 <- stats::coef(A.fit2)[1]
    m1 <- stats::coef(A.fit2)[2]
    n1 <- stats::coef(A.fit2)[3]
    o1 <- sqrt((m1)^2 + (4*l1*n1))

    l2 <- stats::coef(D.fit2)[1]
    m2 <- stats::coef(D.fit2)[2]
    n2 <- stats::coef(D.fit2)[3]
    o2 <- sqrt((m2)^2 + (4*l2*n2))

    Ks1 <- (o1-m1)/(2*l1)
    Kh1 <- (o1+m1)/(o1-m1)
    Mo1 <- 1/o1

    Ks2 <- (o2-m2)/(2*l2)
    Kh2 <- (o2+m2)/(o2-m2)
    Mo2 <- 1/o2

    message("Constants of Hailwood-Horrobin Sorption Model")
    Const1 <- c(Kh1, Kh2)
    Const2 <- c(Ks1, Ks2)
    Const3 <- c(Mo1, Mo2)
    HailHorroConstant <- data.frame(Isotherm, Const1, Const2, Const3)
    names(HailHorroConstant) <- c("Isotherm","Kh", "Ks", "Mo")
    print(HailHorroConstant, row.names = F, right = F)
  }

  HailHorroMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    HailHorroPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Hailwood-Horrobin Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "HailHorro (Adsorption)", color = "HailHorro (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "HailHorro (Desorption)", color = "HailHorro (Desorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
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
      ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.10)) +
      ggplot2::scale_y_continuous(limits = c(minlim1,maxlim1)) +
      ggplot2::scale_shape_manual(name = "Experimental Data", values = c(15, 19)) +
      ggplot2::scale_linetype_manual(name = "Model", values = c(1,1)) +
      ggplot2::scale_color_manual(name = "Model", values = c("blue", "red"))

    suppressWarnings(print(HailHorroPlot))

  }

  HailHorroMSIStats(WaterAct, AdsorpM, DesorpM)
  HailHorroMSIConstant(WaterAct, AdsorpM, DesorpM)
  HailHorroMSIPlot(WaterAct, AdsorpM, DesorpM)
}

