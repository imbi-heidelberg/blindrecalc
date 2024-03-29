% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Student.R
\name{simulation}
\alias{simulation}
\title{Simulate Rejection Probability and Sample Size for Student's t-Test}
\usage{
simulation(
  design,
  n1,
  nuisance,
  recalculation = TRUE,
  delta_true,
  iters = 1000,
  seed = NULL,
  allocation = c("approximate", "exact"),
  ...
)
}
\arguments{
\item{design}{Object of class \code{Student} created by \code{setupStudent}.}

\item{n1}{Either the sample size of the first stage (if
\code{recalculation = TRUE} or the toal sample size (if
\code{recalculation = FALSE}).}

\item{nuisance}{Value of the nuisance parameter. For the
Student's t-test this is the variance.}

\item{recalculation}{Should the sample size be recalculated after n1
n1 patients are recruited?}

\item{delta_true}{effect measure under which the rejection probabilities are
computed}

\item{iters}{Number of simulation iterations.}

\item{seed}{Random seed for simulation.}

\item{allocation}{Whether the allocation ratio should be preserved
exactly (\code{exact}) or approximately (\code{approximate}).}

\item{...}{Further optional arguments.}
}
\value{
Simulated rejection probabilities and sample sizes for
   each nuisance parameter.
}
\description{
This function simulates the probability that a test defined by
\code{\link{setupStudent}} rejects the null hypothesis.
Note that here the nuisance parameter \code{nuisance} is the variance
of the outcome variable sigma^2.
}
\details{
The implementation follows the algorithm in Lu (2019):
Distribution of the two-sample t-test statistic following blinded
sample size re-estimation.
Pharmaceutical Statistics 15: 208-215.
Since Lu (2019) assumes negative non-inferiority margins, the non-inferiority
margin of \code{design} is multiplied with -1 internally.
}
\examples{
d <- setupStudent(alpha = .025, beta = .2, r = 1, delta = 3.5, delta_NI = 0,
                  alternative = "greater", n_max = 156)
simulation(d, n1 = 20, nuisance = 5.5, recalculation = TRUE, delta_true = 3.5)

}
