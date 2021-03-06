% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contact_outbreak.R
\name{contact_outbreak}
\alias{contact_outbreak}
\title{A contact-based outbreak simulator}
\usage{
contact_outbreak(lambda, R, r_SI, r_contact, max_cases = 100)
}
\arguments{
\item{lambda}{The average number of contacts per infectious case, taken to be
the rate of a Poisson distribution.}

\item{R}{The average reproduction number, i.e. number of secondary cases per
infected individual. Note that this number cannot exceed the average
number of contacts \code{lambda}.}

\item{r_SI}{A function generating random delays from the serial interval
distribution, i.e. the delay between primary and secondary symptoms. See
details for more information on how to specify this distribution.}

\item{r_contact}{A function generating random delays from onset of symptom
to non-infectious contacts.}

\item{max_cases}{The maximum number of new cases the simulation can be run
for.}
}
\description{
Under development, please do not use without contacting the author.
}
\details{
We recommend using the package \code{distcrete} for generating the serial
    interval distribution. See: \url{https://github.com/reconhub/distcrete}.
}
\examples{
if (require(distcrete)) {
## get distributions
 SI <- distcrete("gamma", 1L, w = 0, 10, 0.65)
 d_contacts <- distcrete("exp", 1L, w=0, 0.05)

## simulate outbreak and contacts
 set.seed(1)
 x <- contact_outbreak(3, 1.8, SI$r, d_contacts$r)
 x
 plot(x, group = "case_def")

 if (require(incidence)) {
   plot(incidence(x$linelist$onset, 7))
 }
}


}
\author{
Thibaut Jombart \email{thibautjombart@gmail.com}
}

