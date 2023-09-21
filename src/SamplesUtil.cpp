/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: SamplesUtil.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

void flatSamples(cDVector* theInVect, uint theNSample, uint theDimObs, uint theNObsAllSamples, cDVector& theOutVect)
{
// Pour un sample Y[i][j]== Y[i + j*DimObs]
// Pour N samples Y[n][i][j] == Y[n][i + j *DimObs] ? remplacer par un seul sample de taille correcte
        for (uint n = 0 ; n < theNSample ; n++)
                for (uint j = 0 ; j < theDimObs ; j++)
                {       uint myNObs = theInVect[n].mSize / theDimObs ;
                        for (uint s = 0 ; s < myNObs ; s++)
                        {       uint myOldInd = s + j*myNObs ;
                                uint myNewInd = n*theNObsAllSamples + myOldInd ;
                                theOutVect[myNewInd] = theInVect[n][myOldInd] ;
                        }
                }
}

void listSamples(cDVector& theInVect, uint theNSample, uint theDimObs, uint* theNObsSample, cDVector* theOutVect)
{
uint myBegInd = 0 ;
        for (uint n = 0 ; n < theNSample ; n++)
        {       for (uint j = 0 ; j < theDimObs ; j++)
                {       for (uint s = 0 ; s < theNObsSample[n] ; s++)
                        {       uint myNewInd = s + j * theDimObs  ;
                                uint myOldInd = myBegInd + myNewInd ;
                                theOutVect[n][myNewInd] = theInVect[myOldInd] ;
                        }
                }
                myBegInd += theNObsSample[n] ;
        }
}
