#' @title Peleg Moisture Sorption Isotherm
#' @description Peleg model is an empirical 4-parameter isotherm which describes sigmoidal and non-sigmoidal behavior of isotherm plots.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Peleg Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples PelegMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Abu-Ghannam, N., & McKenna, B. (1997) <doi:10.1016/S0260-8774(97)00034-4> The application of Peleg's equation to model water absorption during the soaking of red kidney beans (Phaseolus vulgaris L.). Journal of Food Engineering, 32(4), 391-401.
#' @references Peleg, M. (1993) <doi:10.1111/j.1745-4530.1993.tb00160.x> Assessment of a semi-empirical four parameter general model for sigmoid moisture sorption isotherms. Journal of Food Process Engineering, 16(1), 21-37.
#' @export

PelegMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSIPeleg1 <- data.frame(x1, y1, y2)

  A.eqn1 <- y1 ~ (k1*(x1^n1))+(k2*(x1^n2))
  D.eqn1 <- y2 ~ (k1*(x1^n1))+(k2*(x1^n2))

  start1 <- data.frame(k1 = c(0.01,10), k2 = c(0.01,10), n1 = c(0.01,1), n2 = c(0.01,1))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSIPeleg1, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSIPeleg1, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  A.start1 <- list(k1 = stats::coef(A.fit1)[1], k2 = stats::coef(A.fit1)[2], n1 = stats::coef(A.fit1)[3], n2 = stats::coef(A.fit1)[4])
  D.start1 <- list(k1 = stats::coef(D.fit1)[1], k2 = stats::coef(D.fit1)[2], n1 = stats::coef(D.fit1)[3], n2 = stats::coef(D.fit1)[4])

  suppressWarnings(A.fit2 <- nls2::nls2(A.eqn1, data = MSIPeleg1, start = A.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit2 <- nls2::nls2(D.eqn1, data = MSIPeleg1, start = D.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))


  y3 <- predict(A.fit2)
  y4 <- predict(D.fit2)

  Isotherm <- c("Adsorption", "Desorption")

  message("PELEG MOISTURE SORPTION MODEL")

  message("Peleg Adsorption Isotherm Parameters")
  print(summary(A.fit2))

  message("Peleg Desorption Isotherm Parameters")
  print(summary(D.fit2))

  message("Predicted Values of Moisture Content from Peleg Moisture Sorption Model")
  PelegPredict <- data.frame(x1, y3, y4)
  names(PelegPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(PelegPredict, row.names = F, right = F)

  PelegMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Peleg Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit2), stats::AIC(D.fit2))
    MSIBIC <- c(stats::BIC(A.fit2), stats::BIC(D.fit2))
    PelegCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(PelegCriterion) <- c("Isotherm","AIC", "BIC")
    print(PelegCriterion, row.names=F, right = F)

    message("Error Analysis for Peleg Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1, y3)^2)/(length(y1)-2)), (sqrt(sum(y2, y4)^2)/(length(y2)-2)))
    PelegError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(PelegError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(PelegError, row.names=F, right = F)
  }

  PelegMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Peleg Sorption Model")
    Const1 <- c(stats::coef(A.fit2)[1], stats::coef(D.fit2)[1])
    Const2 <- c(stats::coef(A.fit2)[2], stats::coef(D.fit2)[2])
    Const3 <- c(stats::coef(A.fit2)[3], stats::coef(D.fit2)[3])
    Const4 <- c(stats::coef(A.fit2)[4], stats::coef(D.fit2)[4])
    PelegConstant <- data.frame(Isotherm, Const1, Const2, Const3, Const4)
    names(PelegConstant) <- c("Isotherm","k1", "k2", "n1", "n2")
    print(PelegConstant, row.names = F, right = F)
  }

  PelegMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    PelegPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Peleg Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "Peleg (Adsorption)", color = "Peleg (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "Peleg (Desorption)", color = "Peleg (Desorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_point(ggplot2::aes(x= x1, y=y1, shape  = "Adsorption"), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x= x1, y=y2, shape  = "Desorption"), size = 2) +
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

    suppressWarnings(print(PelegPlot))

  }

  PelegMSIStats(WaterAct, AdsorpM, DesorpM)
  PelegMSIConstant(WaterAct, AdsorpM, DesorpM)
  PelegMSIPlot(WaterAct, AdsorpM, DesorpM)
}

