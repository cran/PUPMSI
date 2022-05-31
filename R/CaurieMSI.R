#' @title Caurie Moisture Sorption Isotherm
#' @description Caurie model is a two-parameter isotherm created for calculation of water sorption data of dehydrated foods over a wide range of water activity.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Caurie Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples CaurieMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Caurie, M. (1970) <doi:10.1111/j.1365-2621.1970.tb01571.x> A new model equation for predicting safe storage moisture levels for optimum stability of dehydrated foods. International Journal of Food Science & Technology, 5(3), 301-307.
#' @references Caurie, M. (2007) <doi:10.1111/j.1365-2621.2006.01203.x> Hysteresis phenomenon in foods. International Journal of Food Science and Technology, 42(1), 45-49.
#' @references Caurie, M. (2011) <doi:10.1007/978-90-481-3585-1_71> Hysteresis in foods. In Encyclopedia of Earth Sciences Series: Vol. Part 4 (p. 384). Springer Netherlands.
#' @export

CaurieMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- log(AdsorpM)
  y2 <- log(DesorpM)
  y3 <- AdsorpM
  y4 <- DesorpM

  MSICaurie1 <- data.frame(x1, y1, y2)
  MSICaurie2 <- data.frame(x1, y3, y4)

  A.linear <- stats::lm(y1 ~ x1, data = MSICaurie1)
  D.linear <- stats::lm(y2 ~ x1, data = MSICaurie1)

  A.start1 <- list(A = stats::coef(A.linear)[1], B = stats::coef(A.linear)[2])
  D.start1 <- list(A = stats::coef(D.linear)[1], B = stats::coef(D.linear)[2])

  A.eqn1 <- y3 ~ exp(A + B*x1)
  D.eqn1 <- y4 ~ exp(A + B*x1)

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSICaurie2, start = A.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSICaurie2, start = D.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y5 <- stats::predict(A.fit1)
  y6 <- stats::predict(D.fit1)

  Isotherm <- c("Adsorption", "Desorption")

  message("CAURIE MOISTURE SORPTION MODEL")

  message("Caurie Adsorption Isotherm Parameters")
  print(summary(A.fit1))

  message("Caurie Desorption Isotherm Parameters")
  print(summary(D.fit1))

  message("Predicted Values of Moisture Content from Caurie Moisture Sorption Model")
  CauriePredict <- data.frame(x1, y5, y6)
  names(CauriePredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(CauriePredict, row.names = F, right = F)

  CaurieMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Caurie Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit1), stats::AIC(D.fit1))
    MSIBIC <- c(stats::BIC(A.fit1), stats::BIC(D.fit1))
    CaurieCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(CaurieCriterion) <- c("Isotherm","AIC", "BIC")
    print(CaurieCriterion, row.names=F, right = F)

    message("Error Analysis for Caurie Sorption Model")
    MSRMSE <- c(Metrics::rmse(y3, y5), Metrics::rmse(y4, y6))
    MSMAE <- c(Metrics::mae(y3, y5), Metrics::mae(y4, y6))
    MSMSE <- c(Metrics::mse(y3, y5), Metrics::mse(y4, y6))
    MSRAE <- c(Metrics::rae(y3, y5), Metrics::rae(y4, y6))
    MSSEE <- c((sqrt(sum(y3, y5)^2)/(length(y3)-2)), (sqrt(sum(y4, y6)^2)/(length(y4)-2)))
    CaurieError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(CaurieError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(CaurieError, row.names=F, right = F)
  }

  CaurieMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Caurie Sorption Model")
    Const1 <- c(stats::coef(A.fit1)[1], stats::coef(D.fit1)[1])
    Const2 <- c(stats::coef(A.fit1)[2], stats::coef(D.fit1)[2])
    CaurieConstant <- data.frame(Isotherm, Const1, Const2)
    names(CaurieConstant) <- c("Isotherm","A", "B")
    print(CaurieConstant, row.names = F, right = F)
  }

  CaurieMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y3,y4,y5,y6)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    CauriePlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Caurie Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y5, linetype = "Caurie (Adsorption)", color = "Caurie (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y6, linetype = "Caurie (Desorption)", color = "Caurie (Desorption)"),
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

    suppressWarnings(print(CauriePlot))

  }

  CaurieMSIStats(WaterAct, AdsorpM, DesorpM)
  CaurieMSIConstant(WaterAct, AdsorpM, DesorpM)
  CaurieMSIPlot(WaterAct, AdsorpM, DesorpM)
}

