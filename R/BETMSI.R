#' @title Brunauer-Emmett-Teller(BET) Moisture Sorption Isotherm
#' @description Brunauer-Emmett-Teller(BET) is a two-parameter isotherm model used for the optimum moisture content determination for drying and storage stability of foods, and in the food's surface area estimation.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Brunauer-Emmett-Teller(BET) Moisture Sorption Isotherm model.
#' @examples \donttest{WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' BETMSI(WaterAct, AdsorpM, DesorpM)}
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Aviara, N. A., et al. (2016). Effect of Temperature and Moisture Sorption Hysteresis on Monolayer Moisture Content of Selected Crops Determined Using BET and GAB Models. 37Th Annual Conference and Annual General Meeting-"Minna 2016," October.
#' @references Staudt, P. B., et al. (2013) <doi:10.1016/j.jfoodeng.2012.07.016> A new method for predicting sorption isotherms at different temperatures using the BET model. Journal of Food Engineering, 114(1), 139-145.
#' @export

BETMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSIBET1 <- data.frame(x1, y1, y2)

  start1 <- data.frame(Mo = c(0, 1), C = c(-200,200))

  A.eqn1 <- y1 ~ (Mo*C*x1)/((1-x1)*(1+ (C*x1) - x1))
  D.eqn1 <- y2 ~ (Mo*C*x1)/((1-x1)*(1+ (C*x1) - x1))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSIBET1, start = start1,
              control = stats::nls.control( maxiter = 100, warnOnly = TRUE),
              algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSIBET1, start = start1,
              control = stats::nls.control( maxiter = 100, warnOnly = TRUE),
              algorithm = "port"))

  A.start1 <- list(Mo = stats::coef(A.fit1)[1], C = stats::coef(A.fit1)[2])
  D.start1 <- list(Mo = stats::coef(D.fit1)[1], C = stats::coef(D.fit1)[2])

  suppressWarnings(A.fit2 <- nls2::nls2(A.eqn1, data = MSIBET1, start = A.start1,
              control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit2 <- nls2::nls2(D.eqn1, data = MSIBET1, start = D.start1,
              control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y3 <- stats::predict(A.fit2)
  y4 <- stats::predict(D.fit2)

  Isotherm <- c("Adsorption", "Desorption")

  message("BET MOISTURE SORPTION MODEL")

  message("BET Adsorption Isotherm Parameters")
  print(summary(A.fit2))

  message("BET Desorption Isotherm Parameters")
  print(summary(D.fit2))

  message("Predicted Values of Moisture Content from BET Moisture Sorption Model")
  BETPredict <- data.frame(x1, y3, y4)
  names(BETPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(BETPredict, row.names = F, right = F)

  BETMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in BET Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit2), stats::AIC(D.fit2))
    MSIBIC <- c(stats::BIC(A.fit2), stats::BIC(D.fit2))
    BETCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(BETCriterion) <- c("Isotherm","AIC", "BIC")
    print(BETCriterion, row.names=F, right = F)

    message("Error Analysis for BET Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1 - y3)^2)/(length(y1)-2)), (sqrt(sum(y2 - y4)^2)/(length(y2)-2)))
    BETError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(BETError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(BETError, row.names=F, right = F)
  }

  BETMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of BET Sorption Model")
    Const1 <- c(stats::coef(A.fit2)[1], stats::coef(D.fit2)[1])
    Const2 <- c(stats::coef(A.fit2)[2], stats::coef(D.fit2)[2])
    BETConstant <- data.frame(Isotherm, Const1, Const2)
    names(BETConstant) <- c("Isotherm","Mo", "C")
    print(BETConstant, row.names = F, right = F)
  }

  BETMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    BETPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "BET Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
           x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "BET (Adsorption)", color = "BET (Adsorption)"),
                 method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "BET (Desorption)", color = "BET (Desorption)"),
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

    suppressWarnings(print(BETPlot))

  }

  BETMSIStats(WaterAct, AdsorpM, DesorpM)
  BETMSIConstant(WaterAct, AdsorpM, DesorpM)
  BETMSIPlot(WaterAct, AdsorpM, DesorpM)
}

