/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cMultivariateNormal.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***   
 *** Class for multivariate normal distribution                                                      
 **************************************************************/

#include "StdAfxRHmm.h"

/*
 * Constructor for MultivariateNormal (distribution)
 * @param theNClass Number of classes
 * @param theDimObs Dimension of observations (for mutltivariate distributions)
 */
cMultivariateNormal::cMultivariateNormal(uint theNClass, uint theDimObs)
{       MESS_CREAT("cMultivariateNormal")
        mvNClass = theNClass ;
        if ( (mvNClass > 0) && (theDimObs > 0) )
        {       mMean = new cDVector[mvNClass] ;
                mCov = new cDMatrix[mvNClass] ;
                
                for (uint      i = 0 ; i < mvNClass ; i++)
                {       mMean[i].ReAlloc(theDimObs) ;
                        mCov[i].ReAlloc(theDimObs, theDimObs) ;
                }
        }
        else
        {       mMean = NULL ;
                mCov = NULL ;
        }
}

/*
 * Destructor for MultivariateNormal (distribution)
 * @param theNClass Number of classes
 * @param theDimObs Dimension of observations (for mutltivariate distributions)
 */
cMultivariateNormal::~cMultivariateNormal()
{       MESS_DESTR("cMultivariateNormal")
        if (mvNClass > 0)
        {       
        uint i ;
                for (i = 0 ; i < mvNClass ; i++)
                {       mMean[i].Delete() ;
                        mCov[i].Delete() ;
                }
                delete [] mMean ;
                delete [] mCov ;
                mMean = NULL ;
                mCov = NULL ;
                mvNClass = 0 ;
        }
}

void cMultivariateNormal::ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba) 
{
uint i, k ;
uint myDimObs = mMean[0].mSize ;

cDMatrix myInvCov=cDMatrix(myDimObs, myDimObs) ;

double myDet ;
        for (i = 0 ; i < mvNClass ; i++)
        {       SymetricInverseAndDet(mCov[i], myDet, myInvCov) ;

                for (k = 0 ; k < theNSample ; k++)
                        MultivariateNormalDensity(theY[k], mMean[i], myInvCov, myDet, theCondProba[k][i]) ;
        }

}

/*
 * Computation of the derivative of the multivariate normal density (Hessian matrix of density of
 * multivariate normal distributed vector Y), as well as the computation of the gradient.
 * theGrad (Gradient) and the Hessian matrix (theHess) are filled in the function.
 * @param theY Y vector
 * @param theGrad Gradient vector of density, filled in the function
 * @param theHess Hessian matrix of density, filled in the function
 */
void cMultivariateNormal::ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess)
{
uint myDimObs = GetDimObs() ;
uint myT = theY.GetSize()/myDimObs ;
uint myNCovParam = myDimObs * (myDimObs + 1)/2 ;
uint myNParam = myDimObs + myNCovParam ;
cDVector* myGrad = new cDVector[myT] ;
cDMatrix* myHess = new cDMatrix[myT] ;

        for (uint t = 0 ; t < myT ; t++)
        {       myGrad[t].ReAlloc(myNParam) ;
                myHess[t].ReAlloc(myNParam, myNParam) ;
        }

        for (uint j = 0 ; j < mvNClass ; j++)
        {
        cDMatrix myInvCov(myDimObs, myDimObs) ;
        double myDeterminant ;
                ArmadilloInvAndDet(mCov[j], myInvCov, myDeterminant) ;
                MultivariateNormalDensityDeriv(theY, mMean[j], mCov[j], myInvCov, myDeterminant, myGrad, myHess) ;
        uint k = (mvNClass - 1)*(mvNClass + 1) + j*myNParam ;
                for (uint t = 0 ; t < myT ; t++)
                {       theGrad[j][t] = 0.0 ;
                        theHess[j][t] = 0.0 ;
                        for (uint p = 0 ; p < myNParam ; p++)
                        {       theGrad[j][t][p+k] = myGrad[t][p] ;
                                for (uint q = p ; q < myNParam ; q++)
                                        theHess[j][t][p+k][q+k] = theHess[j][t][q+k][p+k] = myHess[t][p][q] ;
                        }
                }
        }

        for (uint t = 0 ; t < myT ; t++)
        {       myGrad[t].Delete() ;
                myHess[t].Delete() ;
        }
        delete [] myGrad ;
        delete [] myHess ;
}

void cMultivariateNormal::ComputeCov(cDMatrix& theCov)
{
}

cDVector cMultivariateNormal::GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd)
{
        return theNumDistrParam ;
}

void cMultivariateNormal::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba)
{       
        for (uint i = 0 ; i < mvNClass ; i++)
        {       uint   n,
                                                t       ;
                double myDenominateur = 0.0 ;
                for (n = 0 ; n < theInParam.mNSample ; n++)
                {       uint myT = theInParam.mY[n].mSize / theInParam.mDimObs ;
                        for (t = 0 ; t < myT ; t++)
                                myDenominateur += theBaumWelch.mGamma[n][i][t] ;
                }
                uint   j,
                                                k       ;
                mMean[i] = 0.0 ;
                mCov[i] = 0.0 ;
                for (n = 0 ; n < theInParam.mNSample ; n++)
                {       
                uint myT = theInParam.mY[n].mSize / theInParam.mDimObs ;
                        for (t = 0 ; t < myT ; t++)
                        {       for (j = 0 ; j < theInParam.mDimObs ; j++)
                                {       mMean[i][j] += theBaumWelch.mGamma[n][i][t] * theInParam.mY[n][t+myT*j] ;
                                        for (k = j ; k < theInParam.mDimObs ; k++)
                                                mCov[i][j][k] += theBaumWelch.mGamma[n][i][t] * theInParam.mY[n][t+myT*j] * theInParam.mY[n][t+myT*k] ;
                                }
                        }
                }
                mMean[i] /= myDenominateur ;
                mCov[i] /= myDenominateur ;
                for (j = 0 ; j < theInParam.mDimObs ; j++)
                        for (k = j ; k < theInParam.mDimObs ; k++)
                        {       mCov[i][j][k] -= mMean[i][j]*mMean[i][k] ; 
                                mCov[i][k][j] = mCov[i][j][k] ;
                        }
        }
}

void cMultivariateNormal::InitParameters(cBaumWelchInParam &theInParam)
{
#ifndef _RDLL_
        if (theInParam.mInitType == eKMeans)
        {       
        uint myT = 0 ;
                for (uint i = 0 ; i < theInParam.mNSample ; i++)
                {       myT +=  theInParam.mY[i].mSize ;
                }
        
        cDVector myY(myT)      ;
        myT /= theInParam.mDimObs ;
                flatSamples(theInParam.mY ,theInParam.mNSample, theInParam.mDimObs, myT, myY) ;

        int* mySeq = new int[myT] ;
        KMeans(myY, mvNClass, theInParam.mDimObs, mySeq) ;
        cMultivariateNormal myLoi = cMultivariateNormal(mvNClass, theInParam.mDimObs) ;
        uint* myNbObs = new uint[mvNClass] ;
                for (uint k = 0 ; k < mvNClass ; k++)
                {       for (uint i = 0 ; i < theInParam.mDimObs ; i++)
                        {       myLoi.mMean[k][i] = 0.0 ;
                                for(uint j= 0 ; j < theInParam.mDimObs ;j++)
                                        myLoi.mCov[k][i][j] = 0.0 ;
                        }
                        myNbObs[k] = 0 ;
                } 
        
                for (uint t = 0 ; t < myT ; t++)
                {       uint k = mySeq[t] ;
                        for (uint i = 0 ; i < theInParam.mDimObs ; i++)
                        {       double myObsPlus = (double)(myNbObs[k]+1) ;
                                myLoi.mMean[k][i] = ((double)myNbObs[k]*myLoi.mMean[k][i] + myY[t + myT*i])/myObsPlus ;
                                for (uint j = 0 ; j < theInParam.mDimObs ; j++)
                                        myLoi.mCov[k][i][j] = ((double)myNbObs[k]*myLoi.mCov[k][i][j] + myY[t + myT*i]*myY[t + myT*j])/myObsPlus ;
                                myNbObs[k]++ ;
                        }
                }

                for (uint k = 0 ; k < mvNClass ; k++)
                        for(uint i = 0 ; i < theInParam.mDimObs ; i++)
                                for(uint j = i ; j < theInParam.mDimObs ; j++)
                                {       myLoi.mCov[k][i][j] -= myLoi.mMean[k][i] * myLoi.mMean[k][j] ;
                                        myLoi.mCov[k][j][i] = myLoi.mCov[k][i][j] ;
                                }
                CopyDistr(&myLoi) ;
                delete mySeq ;
                delete myNbObs ;
//              myY.Delete() ;
                return ;
        }
#else
        GetRNGstate();
#endif //_RDLL_

cDVector       myMoy(theInParam.mDimObs),
                        myVar(theInParam.mDimObs),
                        mySigma(theInParam.mDimObs) ;           

double mys = 0.0       ;

        for (uint n = 0 ; n < theInParam.mNSample ; n++)
        {       uint myT = theInParam.mY[n].mSize/theInParam.mDimObs  ;
                cDVector& myY = theInParam.mY[n] ;
                for (uint t = 0 ; t < myT ; t++)
                {       for(uint i = 0 ; i < theInParam.mDimObs ; i++)
                        {       myMoy[i] = (mys*myMoy[i] + myY[t+i*myT])/(mys+1.0) ;
                                myVar[i] = (mys*myVar[i] + myY[t+i*myT]*myY[t+i*myT])/(mys+1.0) ;
                        }
                        mys += 1.0 ;
                }
        }

        for (uint i = 0 ; i < theInParam.mDimObs ; i++)
        {       myVar[i] -= myMoy[i] * myMoy[i] ;
                mySigma[i] = sqrt(myVar[i]) ;
        }

        for (uint i = 0 ; i < mvNClass ; i++)
                mCov[i] = 0.0 ;

        for (uint i = 0 ; i < mvNClass ; i++)
        {       for (uint j = 0 ; j < theInParam.mDimObs ; j++)
                {       mMean[i][j] = -2*mySigma[j] + myMoy[j] + 2*mySigma[j] * unif_rand() ;
                        mCov[i][j][j] = 0.5*myVar[j] + 3*myVar[j] * unif_rand() ;
                }
        }

/*      myMoy.Delete() ;
        myVar.Delete() ;
        mySigma.Delete() ;      
*/

#ifdef _RDLL_
        PutRNGstate() ;
#endif //_RDLL_
}

void cMultivariateNormal::CopyDistr(cDistribution *theSrc)
{
cMultivariateNormal* mySrc = dynamic_cast<cMultivariateNormal *>(theSrc) ;
        if (mySrc)
        {       mvNClass = mySrc->mvNClass ;
        for (uint i= 0 ; i < mvNClass ; i++)
        {       mMean[i] = mySrc->mMean[i] ;
                        mCov[i] = mySrc->mCov[i] ;
        }
        }
        else
                cOTError("Wrong distribution in cMultivariateNormal") ;
}

cMultivariateNormal::cMultivariateNormal(cDistribution& theSrc)
{
        CopyDistr(&theSrc) ;
}

void cMultivariateNormal::Print()
{
uint myDimObs = mMean[0].mSize ;
        Rprintf("Parameters\n") ;
        for (uint i = 0 ; i < mvNClass ; i++)
        {       Rprintf("State %d\n\tm[%d]:\tCov[%d]:\n", i, i, i) ;
                for (uint j = 0 ; j < myDimObs ; j++)
                {       Rprintf("\t%lf", mMean[i][j]) ;
                        for (uint k = 0 ; k < myDimObs ; k++)
                                Rprintf("\t%lf",mCov[i][j][k]) ;
                        Rprintf("\n") ;
                }
        }
}

uint cMultivariateNormal::GetDimObs(void)
{
        if (mvNClass > 0)
                return mCov[0].mNCol ;
        else
                return 0 ;
}

void cMultivariateNormal::GetParam(uint theDeb, cDVector& theParam)
{
uint myDimObs = GetDimObs() ;
uint k = theDeb ;
        for (uint n = 0 ; n < mvNClass ; n++)
        {       for (uint p = 0 ; p < myDimObs ; p++)
                        theParam[k++] = mMean[n][p] ;
                for (uint p = 0 ; p < myDimObs ; p++)
                        for (uint q = 0 ; q < myDimObs ; q++)
                                theParam[k++] = mCov[n][p][q] ;
        }
}

void cMultivariateNormal::SetParam(uint theDeb, cDVector& theParam)
{
uint myDimObs = GetDimObs() ;
uint k = theDeb ;
        for (uint n = 0 ; n < mvNClass ; n++)
        {       for (uint p = 0 ; p < myDimObs ; p++)
                        mMean[n][p] = theParam[k++]  ;
                for (uint p = 0 ; p < myDimObs ; p++)
                        for (uint q = 0 ; q < myDimObs ; q++)
                                mCov[n][p][q] = theParam[k++] ;
        }
}
