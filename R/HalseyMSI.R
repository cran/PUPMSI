#' @title Halsey Moisture Sorption Isotherm
#' @description Halsey Isotherm is a 2-parameter model which expresses condensation of multilayers at proportionally large distances from the surface considering the assumption that a molecule's potential energy is inversely proportional to the nth power of its distance from the surface.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import minpack.lm
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Halsey Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples HalseyMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Halsey, G. (1948) <doi:10.1063/1.1746689> Physical adsorption on non-uniform surfaces. The Journal of Chemical Physics, 16(10), 931-937.
#' @export

HalseyMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSIHalsey1 <- data.frame(x1, y1, y2)

  start1 <- list(a = 1, r = 1)

  A.eqn1 <- y1 ~ (a/-log(x1))^(1/r)
  D.eqn1 <- y2 ~ (a/-log(x1))^(1/r)

  suppressWarnings(A.fit1 <- minpack.lm::nlsLM(A.eqn1, data = MSIHalsey1, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE)))
  suppressWarnings(D.fit1 <- minpack.lm::nlsLM(D.eqn1, data = MSIHalsey1, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE)))

  A.start1 <- list(a = stats::coef(A.fit1)[1], r = stats::coef(A.fit1)[2])
  D.start1 <- list(a = stats::coef(D.fit1)[1], r = stats::coef(D.fit1)[2])

  suppressWarnings(A.fit2 <- nls2::nls2(A.eqn1, data = MSIHalsey1, start = A.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit2 <- nls2::nls2(D.eqn1, data = MSIHalsey1, start = D.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y3 <- stats::predict(A.fit2)
  y4 <- stats::predict(D.fit2)

  Isotherm <- c("Adsorption", "Desorption")

  message("HALSEY MOISTURE SORPTION MODEL")

  message("Halsey Adsorption Isotherm Parameters")
  print(summary(A.fit2))

  message("Halsey Desorption Isotherm Parameters")
  print(summary(D.fit2))

  message("Predicted Values of Moisture Content from Halsey Moisture Sorption Model")
  HalseyPredict <- data.frame(x1, y3, y4)
  names(HalseyPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(HalseyPredict, row.names = F, right = F)

  HalseyMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Halsey Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit2), stats::AIC(D.fit2))
    MSIBIC <- c(stats::BIC(A.fit2), stats::BIC(D.fit2))
    HalseyCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(HalseyCriterion) <- c("Isotherm","AIC", "BIC")
    print(HalseyCriterion, row.names=F, right = F)

    message("Error Analysis for Halsey Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1 - y3)^2)/(length(y1)-2)), (sqrt(sum(y2 - y4)^2)/(length(y2)-2)))
    HalseyError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(HalseyError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(HalseyError, row.names=F, right = F)
  }

  HalseyMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Halsey Sorption Model")
    Const1 <- c(stats::coef(A.fit2)[1], stats::coef(D.fit2)[1])
    Const2 <- c(stats::coef(A.fit2)[2], stats::coef(D.fit2)[2])
    HalseyConstant <- data.frame(Isotherm, Const1, Const2)
    names(HalseyConstant) <- c("Isotherm","a", "r")
    print(HalseyConstant, row.names = F, right = F)
  }

  HalseyMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    HalseyPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Halsey Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "Halsey (Adsorption)", color = "Halsey (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "Halsey (Desorption)", color = "Halsey (Desorption)"),
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

    suppressWarnings(print(HalseyPlot))

  }

  HalseyMSIStats(WaterAct, AdsorpM, DesorpM)
  HalseyMSIConstant(WaterAct, AdsorpM, DesorpM)
  HalseyMSIPlot(WaterAct, AdsorpM, DesorpM)

}

