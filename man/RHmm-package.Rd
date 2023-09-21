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
        simulation and estimation.

RHmm package version 2.1.2 

MCT&LB: updated C++ standard and fixed warnings there, changed namespace, added example and Sunspot dataset, RelDiffWithYesterday.R added, RHmm-package.RD added, updated to CRAN guidelines
			
RHmm package version 2.1.1 - August, 10th, 2016

RHmm package version 2.1.0 - August, 6th, 2014: added from

}
\details{

Thanks to:
- Mait Raag who found a bug in may Multivariate Normal Density function (version <= 1.03).
- Tarmo Leinonen who found some problems in the forwardbackward procedure
- Matteo Re who found bugs in Viterbi path procedure (for discrete distributions)
- Vince Melfi who found a bug in the ForwardBackward procedure.
- Simon Malinowski who found a bug in LLH computation.



}
\author{
Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> and Sebastian Bauer <mail@sebastianbauer.info>

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
# analysis of the model
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
# viterbi alg to predict the states
VitPath <- viterbi(HMMmodell, RelDiffProzent)
HMMcls <- VitPath$states # the predicted classes of the HMM

}
\keyword{RHmm}
\keyword{Hidden Markov Models}
\keyword{HMM}
