% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6Class-networkSampling_fit.R
\name{doubleStandardSampling_fit}
\alias{doubleStandardSampling_fit}
\title{Class for fitting a double-standard sampling}
\description{
Class for fitting a double-standard sampling

Class for fitting a double-standard sampling
}
\section{Super classes}{
\code{\link[missSBM:networkSampling]{missSBM::networkSampling}} -> \code{\link[missSBM:networkSamplingDyads_fit]{missSBM::networkSamplingDyads_fit}} -> \code{doubleStandardSampling_fit}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{vExpec}}{variational expectation of the sampling}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-doubleStandardSampling_fit-new}{\code{doubleStandardSampling_fit$new()}}
\item \href{#method-doubleStandardSampling_fit-update_parameters}{\code{doubleStandardSampling_fit$update_parameters()}}
\item \href{#method-doubleStandardSampling_fit-update_imputation}{\code{doubleStandardSampling_fit$update_imputation()}}
\item \href{#method-doubleStandardSampling_fit-clone}{\code{doubleStandardSampling_fit$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="missSBM" data-topic="networkSampling" data-id="print"><a href='../../missSBM/html/networkSampling.html#method-networkSampling-print'><code>missSBM::networkSampling$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="missSBM" data-topic="networkSamplingDyads_fit" data-id="show"><a href='../../missSBM/html/networkSamplingDyads_fit.html#method-networkSamplingDyads_fit-show'><code>missSBM::networkSamplingDyads_fit$show()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-doubleStandardSampling_fit-new"></a>}}
\if{latex}{\out{\hypertarget{method-doubleStandardSampling_fit-new}{}}}
\subsection{Method \code{new()}}{
constructor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{doubleStandardSampling_fit$new(partlyObservedNetwork, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{partlyObservedNetwork}}{a object with class partlyObservedNetwork representing the observed data with possibly missing entries}

\item{\code{...}}{used for compatibility}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-doubleStandardSampling_fit-update_parameters"></a>}}
\if{latex}{\out{\hypertarget{method-doubleStandardSampling_fit-update_parameters}{}}}
\subsection{Method \code{update_parameters()}}{
a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{doubleStandardSampling_fit$update_parameters(nu, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nu}}{an adjacency matrix with imputed values (only)}

\item{\code{...}}{use for compatibility}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-doubleStandardSampling_fit-update_imputation"></a>}}
\if{latex}{\out{\hypertarget{method-doubleStandardSampling_fit-update_imputation}{}}}
\subsection{Method \code{update_imputation()}}{
a method to update the imputation of the missing entries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{doubleStandardSampling_fit$update_imputation(nu)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nu}}{the matrix of (uncorrected) imputation for missing entries}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-doubleStandardSampling_fit-clone"></a>}}
\if{latex}{\out{\hypertarget{method-doubleStandardSampling_fit-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{doubleStandardSampling_fit$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
