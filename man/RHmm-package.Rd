\name{RHmm-package}
\alias{RHmm-package}
\alias{RHmm}
\alias{RHmm}
\docType{package}
\title{
Hidden Markov Models in R
}
\description{
Discrete, univariate or multivariate gaussian, mixture of
        univariate or multivariate gaussian HMM functions for
        simulation and estimation. M. Thrun added the package into the dbt, because it is no longer availible on CRAN.
				
RHmm package version 2.1.1 - August, 10th, 2016
MT: added into dbt, changed namespace, description, Sunspot dataset added, RelDiffWithYesterday.R added, RHmm-package.RD added

RHmm package version 2.1.0 - August, 6th, 2014: added from

}
\details{
\tabular{ll}{
Package: \tab RHmm\cr
Type: \tab Package\cr
Version: \tab 2.1.1\cr
Date: \tab 2016-08-10\cr
License: \tab GPL\cr
}


Thanks to:
- Mait Raag who found a bug in may Multivariate Normal Density function (version <= 1.03).
- Tarmo Leinonen who found some problems in the forwardbackward procedure
- Matteo Re who found bugs in Viterbi path procedure (for discrete distributions)
- Vince Melfi who found a bug in the ForwardBackward procedure.
- Simon Malinowski who found a bug in LLH computation.



}
\author{
Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> and Sebastian Bauer <mail@sebastianbauer.info>

Maintainer: Michael Thrun <mthrun@informatik.uni-marburg.de>
}
\references{
\url{http://r-forge.r-project.org/projects/rhmm/}
}
\examples{
FestgesetzteClustAnz=3;
data(SunSpots)
str(SunSpots)

RelDiffProzent=RelDiffWithYesterday(SunSpots$SmothedMonSunspot,na.rm=TRUE)*100

# Anpassen des HMM mit Baum-Welch
HMMmodell <- HMMFit(obs=RelDiffProzent, nStates=FestgesetzteClustAnz,
control=list(verbose=1, init= 'KMEANS'),# K-means initialisierung
                     asymptCov=TRUE)
############################################################
# analyse des Modells
HMMmodell
summary(HMMmodell)

HMMmeans = round(HMMmodell$HMM$distribution$mean,1)       
print(paste0('HMMmeans ',HMMmeans ))
HMMsdev  = round(sqrt(HMMmodell$HMM$distribution$var),1)
print(paste0('HMMsdev  ',HMMsdev ))

Uebergangsmatrix = HMMmodell$HMM$transMat
UebergangsmatrixProzent= round(Uebergangsmatrix*100)
UebergangsmatrixProzent
HMMmeans
########################################################################
# viterbi Alg zur vorhersage der Zustaende
VitPath <- viterbi(HMMmodell, RelDiffProzent)
HMMcls <- VitPath$states # die vorhergesagten Klassen des HMM

}
\keyword{RHmm}
\keyword{Hidden Markov Models}
\keyword{HMM}
