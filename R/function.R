#' Generates \code{\link{model.class}} objects.
#'
#' \code{model} is a function that generates a calibration model and the associated likelihood.
#'
#' @details The different statistical models are: \itemize{\item{Model1:
#' \deqn{for i in [1,...,n]  Yexp_i=f(x_i,\Theta)+\epsilon(x_i)}}
#' \item{Model2:
#' \deqn{for i in [1,...,n]  Yexp_i=F(x_i,\Theta)+\epsilon(x_i)}}
#' \item{Model3:
#' \deqn{for i in [1,...,n]  Yexp_i=f(x_i,\Theta)+\delta(x_i)+\epsilon(x_i)}}
#' \item{Model4:
#' \deqn{for i in [1,...,n]  Yexp_i=F(x_i,\Theta)+\delta(x_i)+\epsilon(x_i)}}
#' }
#' where \eqn{for i in [1,\dots,n] \epsilon(x_i)~N(0,\sigma^2)}, \eqn{F(.,.)~PG(m_1(.,.),c_1{(.,.),(.,.)})}
#'  and \eqn{\delta(.)~PG(m_2(.),c_2(.,.))}.
#' There is four kind of models in calibration. They are properly defined in [1].
#'
#'
#' @param code the computational code (function of X and theta)
#' @param X the matrix of the forced variables
#' @param Yexp the vector of the experiments
#' @param model string of the model chosen ("model1","model2","model3","model4")
#' by default "model1" is choosen. See details for precisions.
#' @param opt.emul is a list containing characteristics ahout emulation option. \itemize{
#' \item{\strong{p}}{ the number of parameter in the model (defaul value 1)}
#' \item{\strong{n.emul}}{ the number of points for contituing the Design Of Experiments (DOE) (default value 100)}
#' \item{\strong{type}}{ type of the chosen kernel (value by default "matern5_2") from \code{\link{km}} function}
#' \item{\strong{binf}{ the lower bound of the parameter vector (default value 0)}}
#' \item{\strong{bsup}{ the upper bound of the parameter vector (default value 1)}}
#' \item{\strong{DOE}{ design of experiments for the surrogate (default value NULL). If NULL the DOE is automatically
#' generated as a maximin LHS}}
#' }
#' @return \code{model} returns a \code{model.class} object. This class contains two main methods:
#' \itemize{
#' \item{$plot(\eqn{\Theta},\eqn{\sigma^2}, points=FALSE)}{ this metod generates the plot for a new
#' \eqn{\Theta}, \eqn{\sigma^2} and a new set of data. The option points allows to vizualize the points from
#' the Design Of Experiments (DOE) used for establishing the surrogate.}
#' \item{$summury()}{ this method presents the main information about the model.}
#' }
#' @author T. Tabouy
#' @references [1] Tabouy et al., SBM with missing data
#' @seealso \code{\link{model.class}},
#' @examples
#' @export
nomADonner <- function(sample, vBlocks, sampling, family, directed)
{




}



