#' @title Lewicki-2-Parameter Moisture Sorption Isotherm
#' @description Lewicki-2-Parameter MSI is a two-parameter sorption model that was developed based on Raoult's law, which assumes that water is present either as free water or as water of hydration.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for Lewicki-2-Parameter model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples Lewicki2MSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Lewicki, P. P. (2000) <doi:10.1016/S0260-8774(99)00130-2> Raoult's law based food water sorption isotherm. Journal of Food Engineering, 43(1), 31-40.
#' @export

Lewicki2MSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  x2 <- -log(WaterAct)
  y1 <- log(AdsorpM)
  y2 <- log(DesorpM)
  y3 <- AdsorpM
  y4 <- DesorpM

  MSILewicki21 <- data.frame(x2, y1, y2)
  MSILewicki22 <- data.frame(x1, y3, y4)

  A.linear <- stats::lm(y1 ~ x2, data = MSILewicki21)
  D.linear <- stats::lm(y2 ~ x2, data = MSILewicki21)

  A.start1 <- list(A = exp(stats::coef(A.linear)[1]), B = stats::coef(A.linear)[2] + 1)
  D.start1 <- list(A = exp(stats::coef(D.linear)[1]), B = stats::coef(D.linear)[2] + 1)

  A.eqn1 <- y3 ~ A*(((1/x1)-1)^(B-1))
  D.eqn1 <- y4 ~ A*(((1/x1)-1)^(B-1))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSILewicki22, start = A.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSILewicki22, start = D.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y5 <- predict(A.fit1)
  y6 <- predict(D.fit1)

  Isotherm <- c("Adsorption", "Desorption")

  message("LEWICKI 2-PARAMETER MOISTURE SORPTION MODEL")

  message("Lewicki 2-Parameter Adsorption Isotherm Parameters")
  print(summary(A.fit1))

  message("Lewicki 2-Parameter Desorption Isotherm Parameters")
  print(summary(D.fit1))

  message("Predicted Values of Moisture Content from Lewicki 2-Parameter Moisture Sorption Model")
  Lewicki2Predict <- data.frame(x1, y5, y6)
  names(Lewicki2Predict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(Lewicki2Predict, row.names = F, right = F)

  Lewicki2MSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Lewicki 2-Parameter Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit1), stats::AIC(D.fit1))
    MSIBIC <- c(stats::BIC(A.fit1), stats::BIC(D.fit1))
    Lewicki2Criterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(Lewicki2Criterion) <- c("Isotherm","AIC", "BIC")
    print(Lewicki2Criterion, row.names=F, right = F)

    message("Error Analysis for Lewicki 2-Parameter Sorption Model")
    MSRMSE <- c(Metrics::rmse(y3, y5), Metrics::rmse(y4, y6))
    MSMAE <- c(Metrics::mae(y3, y5), Metrics::mae(y4, y6))
    MSMSE <- c(Metrics::mse(y3, y5), Metrics::mse(y4, y6))
    MSRAE <- c(Metrics::rae(y3, y5), Metrics::rae(y4, y6))
    MSSEE <- c((sqrt(sum(y3, y5)^2)/(length(y3)-2)), (sqrt(sum(y4, y6)^2)/(length(y4)-2)))
    Lewicki2Error <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(Lewicki2Error) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(Lewicki2Error, row.names=F, right = F)
  }

  Lewicki2MSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Lewicki 2-Parameter Sorption Model")
    Const1 <- c(stats::coef(A.fit1)[1], stats::coef(D.fit1)[1])
    Const2 <- c(stats::coef(A.fit1)[2], stats::coef(D.fit1)[2])
    Lewicki2Constant <- data.frame(Isotherm, Const1, Const2)
    names(Lewicki2Constant) <- c("Isotherm","A", "B")
    print(Lewicki2Constant, row.names = F, right = F)
  }

  Lewicki2MSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y3,y4,y5,y6)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    Lewicki2Plot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Lewicki 2-Parameter Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y5, linetype = "Lewicki2 (Adsorption)", color = "Lewicki2 (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y6, linetype = "Lewicki2 (Desorption)", color = "Lewicki2 (Desorption)"),
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

    suppressWarnings(print(Lewicki2Plot))

  }

  Lewicki2MSIStats(WaterAct, AdsorpM, DesorpM)
  Lewicki2MSIConstant(WaterAct, AdsorpM, DesorpM)
  Lewicki2MSIPlot(WaterAct, AdsorpM, DesorpM)

}
