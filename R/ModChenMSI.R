#' @title Modified Chen Moisture Sorption Isotherm
#' @description Modified Chen is 2-parameter model related to the drying principle. It is restricted to situations where diffusion is the primary mode of mass transport and is focused on the steady state of the drying equation.
#' @param WaterAct the numerical value of Water Activity, which ranges from 0 to 1.
#' @param AdsorpM the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.
#' @param DesorpM the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.
#' @import nls2
#' @import Metrics
#' @import ggplot2
#' @import stats
#' @return the nonlinear regression, parameters, and graphical visualization for the Modified Chen  Moisture Sorption Isotherm model.
#' @examples WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
#' @examples AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
#' @examples DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
#' @examples ModChenMSI(WaterAct, AdsorpM, DesorpM)
#' @author Benz L. Rivera
#' @author John Carlo F. Panganiban
#' @author Kim M. Villacorte
#' @author Chester C. Deocaris
#' @references Chen, C. (2019) <doi:10.3390/foods8060191> Validation of the Component Model for Prediction of Moisture Sorption Isotherms of Two Herbs and other Products. Foods, 8(6), 191.
#' @references Chen, C. S. (1971) <doi:10.13031/2013.38421> Equilibrium Moisture Curves for Biological Materials. Transactions of the ASAE, 14(5), 0924-0926.
#' @references Chen, C. S. & Clayton, J. T. (1971) <doi:10.13031/2013.38422> The Effect Of Temperature On Sorption Isotherms Of Biological Materials. Transactions of the ASAE, 14(5), 0927-0929.
#' @export

ModChenMSI <- function(WaterAct, AdsorpM, DesorpM)
{
  x1 <- WaterAct
  x2 <- log(-log(WaterAct))
  y1 <- AdsorpM
  y2 <- DesorpM

  MSIModChen1 <- data.frame(x1, x2, y1, y2)

  A.linear <- stats::lm(y1 ~ x2, data = MSIModChen1)
  D.linear <- stats::lm(y2 ~ x2, data = MSIModChen1)

  A.start1 <- list(A = exp(-stats::coef(A.linear)[1]/stats::coef(A.linear)[2]), B = -1/stats::coef(A.linear)[2])
  D.start1 <- list(A = exp(-stats::coef(D.linear)[1]/stats::coef(D.linear)[2]), B = -1/stats::coef(D.linear)[2])

  A.eqn1 <- y1 ~ (1/B)*(log(A/(-log(x1))))
  D.eqn1 <- y2 ~ (1/B)*(log(A/(-log(x1))))

  suppressWarnings(A.fit1 <- nls2::nls2(A.eqn1, data = MSIModChen1, start = A.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))
  suppressWarnings(D.fit1 <- nls2::nls2(D.eqn1, data = MSIModChen1, start = D.start1,
                                  control = stats::nls.control( maxiter = 100, warnOnly = TRUE), algorithm = "port"))

  y3 <- stats::predict(A.fit1)
  y4 <- stats::predict(D.fit1)

  Isotherm <- c("Adsorption", "Desorption")

  message("MODIFIED CHEN MOISTURE SORPTION MODEL")

  message("Modified Chen Adsorption Isotherm Parameters")
  print(summary(A.fit1))

  message("Modified Chen Desorption Isotherm Parameters")
  print(summary(D.fit1))

  message("Predicted Values of Moisture Content from Modified Chen Moisture Sorption Model")
  ModChenPredict <- data.frame(x1, y3, y4)
  names(ModChenPredict) <- c("Water Activity","Adsorption Curve", "Desorption Curve")
  print(ModChenPredict, row.names = F, right = F)

  ModChenMSIStats <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Fitness of Data in Modified Chen Sorption Model")
    MSIAIC <- c(stats::AIC(A.fit1), stats::AIC(D.fit1))
    MSIBIC <- c(stats::BIC(A.fit1), stats::BIC(D.fit1))
    ModChenCriterion <- data.frame(Isotherm, MSIAIC, MSIBIC)
    names(ModChenCriterion) <- c("Isotherm","AIC", "BIC")
    print(ModChenCriterion, row.names=F, right = F)

    message("Error Analysis for Modified Chen Sorption Model")
    MSRMSE <- c(Metrics::rmse(y1, y3), Metrics::rmse(y2, y4))
    MSMAE <- c(Metrics::mae(y1, y3), Metrics::mae(y2, y4))
    MSMSE <- c(Metrics::mse(y1, y3), Metrics::mse(y2, y4))
    MSRAE <- c(Metrics::rae(y1, y3), Metrics::rae(y2, y4))
    MSSEE <- c((sqrt(sum(y1, y3)^2)/(length(y1)-2)), (sqrt(sum(y2, y4)^2)/(length(y2)-2)))
    ModChenError <- data.frame(Isotherm, MSRMSE, MSMAE, MSMSE, MSRAE, MSSEE)
    names(ModChenError) <- c("Isotherm","RMSE","MAE","MSE","RAE","SEE")
    print(ModChenError, row.names=F, right = F)
  }

  ModChenMSIConstant <- function(WaterAct, AdsorpM, DesorpM)
  {
    message("Constants of Modified Chen Sorption Model")
    Const1 <- c(stats::coef(A.fit1)[1], stats::coef(D.fit1)[1])
    Const2 <- c(stats::coef(A.fit1)[2], stats::coef(D.fit1)[2])
    ModChenConstant <- data.frame(Isotherm, Const1, Const2)
    names(ModChenConstant) <- c("Isotherm","A", "B")
    print(ModChenConstant, row.names = F, right = F)
  }

  ModChenMSIPlot <- function(WaterAct, AdsorpM, DesorpM)
  {
    lim1 <- rbind(y1,y2,y3,y4)
    maxlim1 <- max(lim1)+0.01
    minlim1 <- min(lim1)-0.01

    ModChenPlot <- ggplot2::ggplot() +
      ggplot2::labs(title = "Modified Chen Moisture Sorption Model",subtitle = "Combined Sorption Isotherm",
                    x = "Water Activity", y = "Moisture Content (decimal)") +
      ggplot2::geom_smooth(ggplot2::aes(x= x1, y=y1, linetype = "ModChen (Adsorption)", color = "ModChen (Adsorption)"),
                           method = stats::loess, formula = y ~ x, size=1, se =F) +
      ggplot2::geom_smooth(ggplot2::aes(x=x1, y=y4, linetype = "ModChen (Desorption)", color = "ModChen (Desorption)"),
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

    suppressWarnings(print(ModChenPlot))

  }

  ModChenMSIStats(WaterAct, AdsorpM, DesorpM)
  ModChenMSIConstant(WaterAct, AdsorpM, DesorpM)
  ModChenMSIPlot(WaterAct, AdsorpM, DesorpM)

}

