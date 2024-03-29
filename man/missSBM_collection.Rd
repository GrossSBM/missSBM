% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6Class-missSBM_collection.R
\name{missSBM_collection}
\alias{missSBM_collection}
\title{An R6 class to represent a collection of SBM fits with missing data}
\description{
The function \code{\link[=estimateMissSBM]{estimateMissSBM()}} fits a collection of SBM with missing data for
a varying number of block. These models with class \code{\link{missSBM_fit}}  are stored in an instance
of an object with class \code{\link{missSBM_collection}}, described here.

Fields are accessed via active binding and cannot be changed by the user.

This class comes with a set of R6 methods, some of them being useful for the user and exported
as S3 methods. See the documentation for \code{\link[=show]{show()}} and \code{\link[=print]{print()}}
}
\examples{
## Uncomment to set parallel computing with future
## future::plan("multicore", workers = 2)

## Sample 75\% of dyads in  French political Blogosphere's network data
adjacencyMatrix <- missSBM::frenchblog2007 \%>\%
  igraph::delete.vertices(1:100) \%>\%
  igraph::as_adj () \%>\%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
collection <- estimateMissSBM(adjacencyMatrix, 1:5, sampling = "dyad")
class(collection)

}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{models}}{a list of models}

\item{\code{ICL}}{the vector of Integrated Classification Criterion (ICL) associated to the models in the collection (the smaller, the better)}

\item{\code{bestModel}}{the best model according to the ICL}

\item{\code{vBlocks}}{a vector with the number of blocks}

\item{\code{optimizationStatus}}{a data.frame summarizing the optimization process for all models}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-missSBM_collection-new}{\code{missSBM_collection$new()}}
\item \href{#method-missSBM_collection-estimate}{\code{missSBM_collection$estimate()}}
\item \href{#method-missSBM_collection-explore}{\code{missSBM_collection$explore()}}
\item \href{#method-missSBM_collection-plot}{\code{missSBM_collection$plot()}}
\item \href{#method-missSBM_collection-show}{\code{missSBM_collection$show()}}
\item \href{#method-missSBM_collection-print}{\code{missSBM_collection$print()}}
\item \href{#method-missSBM_collection-clone}{\code{missSBM_collection$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-new"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-new}{}}}
\subsection{Method \code{new()}}{
constructor for networkSampling
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$new(partlyObservedNet, sampling, clusterInit, control)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{partlyObservedNet}}{An object with class \code{\link{partlyObservedNetwork}}.}

\item{\code{sampling}}{The sampling design for the modelling of missing data: MAR designs ("dyad", "node") and MNAR designs ("double-standard", "block-dyad", "block-node" ,"degree")}

\item{\code{clusterInit}}{Initial clustering: a list of vectors, each with size \code{ncol(adjacencyMatrix)}.}

\item{\code{control}}{a list of parameters controlling advanced features. Only 'trace' and 'useCov' are relevant here. See \code{\link[=estimateMissSBM]{estimateMissSBM()}} for details.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-estimate"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-estimate}{}}}
\subsection{Method \code{estimate()}}{
method to launch the estimation of the collection of models
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$estimate(control)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{control}}{a list of parameters controlling the variational EM algorithm. See details of function \code{\link[=estimateMissSBM]{estimateMissSBM()}}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-explore"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-explore}{}}}
\subsection{Method \code{explore()}}{
method for performing exploration of the ICL
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$explore(control)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{control}}{a list of parameters controlling the exploration, similar to those found in the regular function \code{\link[=estimateMissSBM]{estimateMissSBM()}}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-plot"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-plot}{}}}
\subsection{Method \code{plot()}}{
plot method for missSBM_collection
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$plot(type = c("icl", "elbo", "monitoring"))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{the type specifies the field to plot, either "icl", "elbo" or "monitoring". Default is "icl"}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-show"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-show}{}}}
\subsection{Method \code{show()}}{
show method for missSBM_collection
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$show()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-print"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-print}{}}}
\subsection{Method \code{print()}}{
User friendly print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-missSBM_collection-clone"></a>}}
\if{latex}{\out{\hypertarget{method-missSBM_collection-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{missSBM_collection$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
