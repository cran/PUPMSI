#' @title Lewicki-3-Parameter Moisture Sorption Isotherm
#' @description The three-parameter Lewicki model is most suitable for describing the sorption characteristics of raw potato, potato starch, starch-sugar and starch-salt gels within specific temperature and water activity ranges.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for Lewicki-3-Parameter model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples Lewicki3MSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references McMinn, W. A., et al. (2004) <doi:10.1002/jsfa.1866> Assessment of two- and three-parameter Lewicki models for description of sorption phenomena of starch materials. Journal of the Science of Food and Agriculture, 84(13), 1695-1700.
#' @export

Lewicki3MSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  y1 <- AdsorpM
  y2 <- DesorpM

  MSILewicki31 <- data.frame(x1, y1, y2)

  A.eqn1 <- y1 ~ (f/((1-x1)^g))-(f/((1+(x1^h))))
  D.eqn1 <- y2 ~ (f/((1-x1)^g))-(f/((1+(x1^h))))

  start1 <- data.frame(f = c(0,1), g = c(0,1), h = c(0,1))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSILewicki31, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSILewicki31, start = start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  A.start1 <- list(f = stats::coef(A.fit1)[1], g = stats::coef(A.fit1)[2], h = stats::coef(A.fit1)[3])
  D.start1 <- list(f = stats::coef(D.fit1)[1], g = stats::coef(D.fit1)[2], h = stats::coef(D.fit1)[3])

  suppressWarnings(A.fit2 <- nls2::nls2(A.eqn1, data = MSILewicki31, start = A.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit2 <- nls2::nls2(D.eqn1, data = MSILewicki31, start = D.start1,
                                        control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))


  y3 <- predict(A.fit2)
  y4 <- predict(D.fit2)

  Isotherm <- c("Adsorption", "Desorption")

  message("LEWICKI 3-PARAMETER MOISTURE SORPTION MODEL")

  message("Lewicki 3-Parameter Adsorption Isotherm Parameters")
  print(summary(A.fit2))

  message("Lewicki 3-Parameter Desorption Isotherm Parameters")
  print(summary(D.fit2))

  message("Predicted Values of Moisture Content from Lewicki 3-Parameter Moisture Sorption Model")
  Lewicki3Predict <- data.frame(x1, y3, y4)
  names(Lewicki3Predict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(Lewicki3Predict, row.names = F, right = F)


  Lewicki3MSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Lewicki 3-Parameter Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit2), stats::AIC(D.fit2))
    MSIBIC <- c(stats::BIC(A.fit2), stats::BIC(D.fit2))
    Lewicki3Criterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(Lewicki3Criterion) <- c("Isotherm","AIC", "BIC")
    print(Lewicki3Criterion, row.names=F, right = F)

    message("Error Analysis for Lewicki 3-Parameter Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1 - y3)^2)/(length(y1)-2)), (sqrt(sum(y2 - y4)^2)/(length(y2)-2)))
    Lewicki3Error <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(Lewicki3Error) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(Lewicki3Error, row.names=F, right = F)
  }

  Lewicki3MSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Lewicki 3-Parameter Sorption Model")
    Const1 <- c(stats::coef(A.fit2)[1], stats::coef(D.fit2)[1])
    Const2 <- c(stats::coef(A.fit2)[2], stats::coef(D.fit2)[2])
    Const3 <- c(stats::coef(A.fit2)[3], stats::coef(D.fit2)[3])
    Lewicki3Constant <- data.frame(Isotherm, Const1, Const2, Const3)
    names(Lewicki3Constant) <- c("Isotherm","f", "g", "h")
    print(Lewicki3Constant, row.names = F, right = F)
  }

  Lewicki3MSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    Lewicki3Plot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Lewicki 3-Parameter Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "Lewicki3 (Adsorption)", color = "Lewicki3 (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "Lewicki3 (Desorption)", color = "Lewicki3 (Desorption)"),
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

    suppressWarnings(print(Lewicki3Plot))

  }

  Lewicki3MSIStats(WaterAct, AdsorpM, DesorpM)
  Lewicki3MSIConstant(WaterAct, AdsorpM, DesorpM)
  Lewicki3MSIPlot(WaterAct, AdsorpM, DesorpM)
}

