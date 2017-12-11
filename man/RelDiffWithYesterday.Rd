\name{RelDiffWithYesterday}
\alias{RelDiffWithYesterday}
\title{ Transformation  von Daten  in Relative Differenzen }
\usage{
RelDiffWithYesterday(Data, minData=0.01,na.rm=F)
}
\arguments{
\item{Data}{Data[1:n,1:d] tageweise in den Spalten.}
\item{minData}{Data = max(Data,Data), default MinPreis = 0.01, verhindert onderflow im denominator
}
\item{na.rm}{
if TRUE, NaN of minDATA set to zero
}
}
\value{
RelDiff[1:n,1:d]
}
\description{
 Transformation  von Daten  in Relative Differenzen bzgl Gestern Spaltenweise, Annahme: jede Zeile entspricht einem Zeitpunkt, 
 jeder Zeitpunkt ist gleiochweit vom naechsten Zeitpunkt entfernt, z.B. taeglicher zeitverlauf
 }
\author{
Michael Thrun
}

