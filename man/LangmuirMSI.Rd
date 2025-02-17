% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LangmuirMSI.R
\name{LangmuirMSI}
\alias{LangmuirMSI}
\title{Langmuir Moisture Sorption Isotherm}
\usage{
LangmuirMSI(WaterAct, AdsorpM, DesorpM)
}
\arguments{
\item{WaterAct}{the numerical value of Water Activity, which ranges from 0 to 1.}

\item{AdsorpM}{the numerical value of the Moisture content of the Adsorption curve, which ranges from 0 to 1.}

\item{DesorpM}{the numerical value of the Moisture content of the Desorption curve, which ranges from 0 to 1.}
}
\value{
the nonlinear regression, parameters, and graphical visualization for the  Langmuir Moisture Sorption Isotherm model.
}
\description{
Langmuir Isotherm is a two-parameter model applicable for unimolecular layers with similar sorption sites. Langmuir's isotherm is the most crucial equation among the theoretical models, whose basis are the forces acting between the product surface and the condensed water from the vapor as a monomolecular layer.
}
\examples{
WaterAct <- c(0.1145,0.2274,0.3265,0.4291,0.6342,0.7385,0.8274,0.9573)
AdsorpM <- c(0.0234, 0.0366, 0.0496, 0.0648, 0.0887, 0.1096, 0.1343, 0.1938)
DesorpM <- c(0.0459, 0.0637, 0.0794, 0.0884, 0.1158, 0.1298,0.1500, 0.1938)
LangmuirMSI(WaterAct, AdsorpM, DesorpM)
}
\references{
Andrade, R. D., et al. (2011). Models of sorption isotherms for food: Uses and limitations. Vitae. In Vitae (Vol. 18, Issue 3). Facultad De Qui??mica Farmace??utica, Universidad de Antioquia. http://www.scielo.org.co/scielo.php?script=sci_arttext&pid=S0121-40042011000300012&lng=en&nrm=iso&tlng=en

Saroyda, J. V., Cruz, et al. (2020) <doi:10.1016/S0001-8686(00)00082> Package "PUPAIM" Type Package Title A Collection of Physical and Chemical Adsorption Isotherm Models Version 0.2.0. <doi:10.1016/S0001-8686(00)00082>
}
\author{
Benz L. Rivera

John Carlo F. Panganiban

Kim M. Villacorte

Chester C. Deocaris
}
