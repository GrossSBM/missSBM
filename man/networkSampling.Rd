% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6Class-networkSampling.R
\name{networkSampling}
\alias{networkSampling}
\title{Definition of R6 Class 'networkSampling'}
\description{
Definition of R6 Class 'networkSampling'

Definition of R6 Class 'networkSampling'
}
\details{
this virtual class is the mother of all subtypes of networkSampling (either sampler or fit)
It is used to define a sampling model for a network.
It has a rSampling method which takes an adjacency matrix as an input and send back an object with class partlyObservedNetwork.
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{type}}{a character for the type of sampling}

\item{\code{parameters}}{the vector of parameters associated with the sampling at play}

\item{\code{df}}{the number of entries in the vector of parameters}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-networkSampling-new}{\code{networkSampling$new()}}
\item \href{#method-networkSampling-show}{\code{networkSampling$show()}}
\item \href{#method-networkSampling-print}{\code{networkSampling$print()}}
\item \href{#method-networkSampling-clone}{\code{networkSampling$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSampling-new"></a>}}
\if{latex}{\out{\hypertarget{method-networkSampling-new}{}}}
\subsection{Method \code{new()}}{
constructor for networkSampling
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSampling$new(type = NA, parameters = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{character for the type of sampling. must be in ("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree")}

\item{\code{parameters}}{the vector of parameters associated to the sampling at play}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSampling-show"></a>}}
\if{latex}{\out{\hypertarget{method-networkSampling-show}{}}}
\subsection{Method \code{show()}}{
show method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSampling$show(
  type = paste0(private$name, "-model for network sampling\\n")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{character used to specify the type of sampling}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSampling-print"></a>}}
\if{latex}{\out{\hypertarget{method-networkSampling-print}{}}}
\subsection{Method \code{print()}}{
User friendly print method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSampling$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSampling-clone"></a>}}
\if{latex}{\out{\hypertarget{method-networkSampling-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSampling$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
