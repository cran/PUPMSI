#' @title Langmuir Moisture Sorption Isotherm
#' @description Langmuir Isotherm is a two-parameter model applicable for unimolecular layers with similar sorption sites. Langmuir's isotherm is the most crucial equation among the theoretical models, whose basis are the forces acting between the product surface and the condensed water from the vapor as a monomolecular layer.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the  Langmuir Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples LangmuirMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Andrade, R. D., et al. (2011). Models of sorption isotherms for food: Uses and limitations. Vitae. In Vitae (Vol. 18, Issue 3). Facultad De Qui??mica Farmace??utica, Universidad de Antioquia. http://www.scielo.org.co/scielo.php?script=sci_arttext&pid=S0121-40042011000300012&lng=en&nrm=iso&tlng=en
#' @references Saroyda, J. V., Cruz, et al. (2020) <doi:10.1016/S0001-8686(00)00082> Package "PUPAIM" Type Package Title A Collection of Physical and Chemical Adsorption Isotherm Models Version 0.2.0. <doi:10.1016/S0001-8686(00)00082>
#' @export

LangmuirMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSILangmuir2 <- data.frame(x1, y1, y2)

  start1 <-data.frame(Mo = c(0,1), C = c(-1,1))

  A.eqn1 <- y1 ~ (C*Mo*x1)/((C*x1)+1)
  D.eqn1 <- y2 ~ (C*Mo*x1)/((C*x1)+1)

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSILangmuir2, start = start1,
                 control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings( D.fit1 <- nls2::nls2(D.eqn1, data = MSILangmuir2, start = start1,
                 control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  A.start2 <- list(Mo = stats::coef(A.fit1)[1], C = stats::coef(A.fit1)[2])
  D.start2 <- list(Mo = stats::coef(D.fit1)[1], C = stats::coef(D.fit1)[2])

  A.fit2 <- suppressWarnings(nls2::nls2(A.eqn1, data = MSILangmuir2, start = A.start2,
                 control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  D.fit2 <- suppressWarnings(nls2::nls2(D.eqn1, data = MSILangmuir2, start = D.start2,
                 control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y3 <- stats::predict(A.fit2)
  y4 <- stats::predict(D.fit2)

  Isotherm <- c("Adsorption", "Desorption")

  message("LANGMUIR MOISTURE SORPTION MODEL")

  message("Langmuir Adsorption Isotherm Parameters")
  print(summary(A.fit2))

  message("Langmuir Desorption Isotherm Parameters")
  print(summary(D.fit2))

  message("Predicted Values of Moisture Content from Langmuir Moisture Sorption Model")
  LangmuirPredict <- data.frame(x1, y3, y4)
  names(LangmuirPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(LangmuirPredict, row.names = F, right = F)

  LangmuirMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Langmuir Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit2), stats::AIC(D.fit2))
    MSIBIC <- c(stats::BIC(A.fit2), stats::BIC(D.fit2))
    LangmuirCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(LangmuirCriterion) <- c("Isotherm","AIC", "BIC")
    print(LangmuirCriterion, row.names=F, right = F)

    message("Error Analysis for Langmuir Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1 - y3)^2)/(length(y1)-2)), (sqrt(sum(y2 - y4)^2)/(length(y2)-2)))
    LangmuirError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(LangmuirError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(LangmuirError, row.names=F, right = F)
  }

  LangmuirMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Langmuir Sorption Model")
    Const1 <- c(stats::coef(A.fit2)[1], stats::coef(D.fit2)[1])
    Const2 <- c(stats::coef(A.fit2)[2], stats::coef(D.fit2)[2])
    LangmuirConstant <- data.frame(Isotherm, Const1, Const2)
    names(LangmuirConstant) <- c("Isotherm","Mo", "C")
    print(LangmuirConstant, row.names = F, right = F)
  }

  LangmuirMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    LangmuirPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Langmuir Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "Langmuir (Adsorption)", color = "Langmuir (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "Langmuir (Desorption)", color = "Langmuir (Desorption)"),
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

    suppressWarnings(print(LangmuirPlot))

  }

  LangmuirMSIStats(WaterAct, AdsorpM, DesorpM)
  LangmuirMSIConstant(WaterAct, AdsorpM, DesorpM)
  LangmuirMSIPlot(WaterAct, AdsorpM, DesorpM)

}

