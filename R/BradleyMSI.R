#' @title Bradley Moisture Sorption Isotherm
#' @description Bradley model is a two-parameter isotherm model that measures polar nature of sorptive surfaces.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Bradley Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples BradleyMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Bradley, R. Stevenson (1936) <doi:10.1039/JR9360001467> Polymolecular adsorbed films. Part I. The adsorption of argon on salt crystals at low temperatures, and the determination of surface fields. Journal of the Chemical Society (Resumed), (), 1467-.
#' @export

BradleyMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  x2 <- log(-log(WaterAct))
  y1 <- AdsorpM
  y2 <- DesorpM

  MSIBradley1 <- data.frame(x1, x2, y1, y2)

  A.linear <- stats::lm(y1 ~ x2, data = MSIBradley1)
  D.linear <- stats::lm(y2 ~ x2, data = MSIBradley1)

  A.start1 <- list(k1 = exp(1/stats::coef(A.linear)[2]), k2 = exp(-stats::coef(A.linear)[1]/stats::coef(A.linear)[2]))
  D.start1 <- list(k1 = exp(1/stats::coef(D.linear)[2]), k2 = exp(-stats::coef(D.linear)[1]/stats::coef(D.linear)[2]))

  A.eqn1 <- y1 ~ (((log(-log(x1)))/log(k1)) - (log(k2)/log(k1)))
  D.eqn1 <- y2 ~ (((log(-log(x1)))/log(k1)) - (log(k2)/log(k1)))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSIBradley1, start = A.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "grid-search"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSIBradley1, start = D.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "grid-search"))

  y3 <- stats::predict(A.fit1)
  y4 <- stats::predict(D.fit1)

  Isotherm <- c("Adsorption", "Desorption")

  message("BRADLEY MOISTURE SORPTION MODEL")

  message("Bradley Adsorption Isotherm Parameters")
  print(summary(A.fit1))

  message("Bradley Desorption Isotherm Parameters")
  print(summary(D.fit1))

  message("Predicted Values of Moisture Content from Bradley Moisture Sorption Model")

  BradleyPredict <- data.frame(x1, y3, y4)
  names(BradleyPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(BradleyPredict, row.names = F, right = F)

  BradleyMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Bradley Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit1), stats::AIC(D.fit1))
    MSIBIC <- c(stats::BIC(A.fit1), stats::BIC(D.fit1))
    BradleyCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(BradleyCriterion) <- c("Isotherm","AIC", "BIC")
    print(BradleyCriterion, row.names=F, right = F)

    message("Error Analysis for Bradley Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1 - y3)^2)/(length(y1)-2)), (sqrt(sum(y2 - y4)^2)/(length(y2)-2)))
    BradleyError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(BradleyError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(BradleyError, row.names=F, right = F)
  }

  BradleyMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Bradley Sorption Model")
    Const1 <- c(stats::coef(A.fit1)[1], stats::coef(D.fit1)[1])
    Const2 <- c(stats::coef(A.fit1)[2], stats::coef(D.fit1)[2])
    BradleyConstant <- data.frame(Isotherm, Const1, Const2)
    names(BradleyConstant) <- c("Isotherm","k1", "k2")
    print(BradleyConstant, row.names = F, right = F)
  }

  BradleyMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    BradleyPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Bradley Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y3, linetype = "Bradley (Adsorption)", color = "Bradley (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "Bradley (Desorption)", color = "Bradley (Desorption)"),
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

    suppressWarnings(print(BradleyPlot))

  }

  BradleyMSIStats(WaterAct, AdsorpM, DesorpM)
  BradleyMSIConstant(WaterAct, AdsorpM, DesorpM)
  BradleyMSIPlot(WaterAct, AdsorpM, DesorpM)

}

