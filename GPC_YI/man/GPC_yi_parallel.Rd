\name{GPCYI_yi_parallel}
\alias{GPCYI_yi_parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GPC for youden index
}
\description{
scale the posterior for youden index
}
\usage{
GPCYI_yi_parallel(nn, data, theta_boot, data_boot, alpha, M_samp, B_resamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nn}{
	Number of data points.
}
  \item{data}{
	Data.
}
  \item{theta_boot}{
	B bootstrap estimates of the parameter.
}
  \item{data_boot}{
	B bootstrap resampled data sets.
}
  \item{alpha}{
	Nominal coverage probability is given by 100*(1-alpha)%.
}
  \item{M_samp}{
	Number of posterior samples.
}
  \item{B_resamp}{
	Number of bootstrap resamples.
}
}
\value{
 
}

\examples{
\dontrun{
 
}
}

