#' @title Kuhn Moisture Sorption Isotherm
#' @description Kuhn Isotherm is a two-parameter model which contains many defining characteristics wherein each surface site has a different adsorption potential, as well as cluster formations on each site due to increase in partial pressure.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Kuhn Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples KuhnMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Bi, Y., et al. (2018). The prediction of moisture adsorption isotherm for commercial sodium bicarbonate powder. International Journal of Scientific & Engineering Research, 9(3).
#' @references Kuhn, I. (1967) <doi:10.1016/0021-9797(67)90202-0> A generalized potential theory of adsorption. I. The derivation of a general equation for adsorption isotherms. Journal of Colloid And Interface Science, 23(4), 563-571.
#' @export

KuhnMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  x2 <- 1/log(WaterAct)
  y1 <- log(AdsorpM)
  y2 <- log(DesorpM)
  y3 <- AdsorpM
  y4 <- DesorpM

  MSIKuhn1 <- data.frame(x2, y1, y2)
  MSIKuhn2 <- data.frame(x1, y3, y4)

  A.linear <- stats::lm(y1 ~ x2, data = MSIKuhn1)
  D.linear <- stats::lm(y2 ~ x2, data = MSIKuhn1)

  A.start1 <- list(A = stats::coef(A.linear)[2], B = stats::coef(A.linear)[1])
  D.start1 <- list(A = stats::coef(D.linear)[2], B = stats::coef(D.linear)[1])

  A.eqn1 <- y3 ~ (A/log(x1)) + B
  D.eqn1 <- y4 ~ (A/log(x1)) + B

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSIKuhn2, start = A.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSIKuhn2, start = D.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y5 <- predict(A.fit1)
  y6 <- predict(D.fit1)

  Isotherm <- c("Adsorption", "Desorption")

  message("KUHN MOISTURE SORPTION MODEL")

  message("Kuhn Adsorption Isotherm Parameters")
  print(summary(A.fit1))

  message("Kuhn Desorption Isotherm Parameters")
  print(summary(D.fit1))

  message("Predicted Values of Moisture Content from Kuhn Moisture Sorption Model")
  KuhnPredict <- data.frame(x1, y5, y6)
  names(KuhnPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(KuhnPredict, row.names = F, right = F)

  KuhnMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Kuhn Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit1), stats::AIC(D.fit1))
    MSIBIC <- c(stats::BIC(A.fit1), stats::BIC(D.fit1))
    KuhnCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(KuhnCriterion) <- c("Isotherm","AIC", "BIC")
    print(KuhnCriterion, row.names=F, right = F)

    message("Error Analysis for Kuhn Sorption Model")
    MSRMSE <- c(Metrics::rmse(y3, y5), Metrics::rmse(y4, y6))
    MSMAE <- c(Metrics::mae(y3, y5), Metrics::mae(y4, y6))
    MSMSE <- c(Metrics::mse(y3, y5), Metrics::mse(y4, y6))
    MSRAE <- c(Metrics::rae(y3, y5), Metrics::rae(y4, y6))
    MSSEE <- c((sqrt(sum(y3, y5)^2)/(length(y3)-2)), (sqrt(sum(y4, y6)^2)/(length(y4)-2)))
    KuhnError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(KuhnError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(KuhnError, row.names=F, right = F)
  }

  KuhnMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Kuhn Sorption Model")
    Const1 <- c(stats::coef(A.fit1)[1], stats::coef(D.fit1)[1])
    Const2 <- c(stats::coef(A.fit1)[2], stats::coef(D.fit1)[2])
    KuhnConstant <- data.frame(Isotherm, Const1, Const2)
    names(KuhnConstant) <- c("Isotherm","A", "B")
    print(KuhnConstant, row.names = F, right = F)
  }

  KuhnMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y3,y4,y5,y6)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    KuhnPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Kuhn Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y5, linetype = "Kuhn (Adsorption)", color = "Kuhn (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y6, linetype = "Kuhn (Desorption)", color = "Kuhn (Desorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_point(ggplot2::aes(x= x1, y=y3, shape  = "Adsorption"), size = 2) +
      ggplot2::geom_point(ggplot2::aes(x=x1, y=y4, shape  = "Desorption"), size = 2) +
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

    suppressWarnings(print(KuhnPlot))

  }

  KuhnMSIStats(WaterAct, AdsorpM, DesorpM)
  KuhnMSIConstant(WaterAct, AdsorpM, DesorpM)
  KuhnMSIPlot(WaterAct, AdsorpM, DesorpM)

}

