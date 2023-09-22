/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cMixtMultivariateNormal.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

/*
 * Calculates the mixture mutlivariate normal density.
 * The parameter theDens will be filled with the density.
 * @param theY Y vector
 * @param theNMixt Number of mixtures
 * @param theInvCov Inverse covariance matrix
 * @param theDet Determinant of the covariance matrix
 * @param thep Weight vector for mixtures
 * @param theDens Density vector to be filled
 */
static void MixtMultivariateNormalDensity(cDVector& theY, uint theNMixt, cDVector* theMean, cDMatrix* theInvCov, cDVector& theDet, cDVector& thep, double* theDens)
{
uint myT = theY.mSize / theMean[0].mSize ;
double* myDens = new double[myT] ;
        
        for (uint t = 0 ; t < myT ; t++)
                theDens[t] = 0.0 ;
        for (uint j = 0 ; j < theNMixt ; j++)
        {       MultivariateNormalDensity(theY, theMean[j], theInvCov[j], theDet[j], myDens) ;
                for (uint t = 0 ; t < myT ; t++)
                        theDens[t] += thep[j] * myDens[t] ;                     
        }
        for (uint t = 0 ; t < myT ; t++)
                theDens[t] = MAX(theDens[t], 1e-30) ;

        delete[] myDens ;
}

/*
 * Constructor for MixtureMultivariateNormal (distribution)
 * @param theNClass Number of classes
 * @param theNMixt Number of mixtures
 * @param theDimObs Dimension of observations (for mutltivariate distributions)
 */
cMixtMultivariateNormal::cMixtMultivariateNormal(uint theNClass, uint theNMixt, uint theDimObs)
{       MESS_CREAT("cMixtMultivariateNormal")
        mvNClass = theNClass ;
        mvNMixt = theNMixt ;
        mvDimObs = theDimObs ;
        if ( (theNClass > 0) && (theNMixt > 0) && (theDimObs > 0) )
        {       mMean = new cDVector*[theNClass] ;
                mCov = new cDMatrix*[theNClass] ;
                mp = new cDVector[theNClass] ;
                for (uint i = 0 ; i < mvNClass ; i++)
                {       mMean[i] = new cDVector[theNMixt] ;
                        mCov[i] = new cDMatrix[theNMixt] ;
                        mp[i].ReAlloc(theNMixt) ;
                        for (uint j = 0 ; j < theNMixt ; j++)
                        {       mMean[i][j].ReAlloc(theDimObs) ;
                                mCov[i][j].ReAlloc(theDimObs, theDimObs) ;
                        }
                }
        }
        else
        {       mMean = NULL ;
                mp = NULL ;
                mCov = NULL ;
                mvNClass = mvNMixt = mvDimObs = 0 ;
        }
}

/*
 * Deconstructor for MixtureMultivariateNormal (distribution)
 */
cMixtMultivariateNormal::~cMixtMultivariateNormal()
{       MESS_DESTR("cMixtMultivariateNormal")
        for (uint i = 0 ; i < mvNClass ; i++)
        {       for (uint j = 0 ; j < mvNMixt ; j++)
                {       mMean[i][j].Delete() ;
                        mCov[i][j].Delete() ;
                }
                mp[i].Delete() ;
        }
        delete [] mMean ;
        delete [] mCov ;
        delete [] mp ;
        mMean = NULL ;
        mCov = NULL ;
        mp = NULL ;
        mvNClass = mvNMixt = mvDimObs = 0 ;
}

void cMixtMultivariateNormal::ComputeCondProba(cDVector* theY, uint theNSample, cDMatrix* theCondProba)
{
cDMatrix* myInvCov = new cDMatrix[mvNMixt] ;
cDVector myDet = cDVector(mvNMixt) ;

        for (uint i = 0 ; i < mvNMixt ; i++)
                myInvCov[i].ReAlloc(mvDimObs, mvDimObs) ;


        for (uint i = 0 ; i < mvNClass ; i++)
        {       for (uint j = 0 ; j < mvNMixt ; j++)
                        {
                        SymetricInverseAndDet(mCov[i][j], myDet[j], myInvCov[j]) ;
                        }
                for (uint n = 0 ; n < theNSample ; n++)
                        MixtMultivariateNormalDensity(theY[n], mvNMixt, mMean[i], myInvCov, myDet, mp[i], theCondProba[n][i]) ;
        }
        for (uint i = 0 ; i < mvNMixt ; i++)
                myInvCov[i].Delete() ;
        delete [] myInvCov ;
}

/*
 * Computation of the derivative of the mixture multivariate normal density (Hessian matrix of density of
 * mixture multivariate normal distributed vector Y), as well as the computation of the gradient.
 * theGrad (Gradient) and the Hessian matrix (theHess) are filled in the function.
 * @param theY Y vector
 * @param theGrad Gradient vector of density, filled in the function
 * @param theHess Hessian matrix of density, filled in the function
 */
void cMixtMultivariateNormal::ComputeDerivative(cDVector& theY, cDVector** theGrad, cDMatrix** theHess)
{
uint myT = theY.GetSize()/mvDimObs ;
uint myNCovParam = mvDimObs * (mvDimObs + 1)/2 ;
uint myNNormParam = mvDimObs + myNCovParam ;
uint myNParam = mvNMixt*myNNormParam + mvNMixt - 1 ;
cDVector* myGradNorm = new cDVector[myT] ;
cDMatrix* myHessNorm = new cDMatrix[myT] ;
cDVector myDens(myT), myLastDens(myT);

        for (uint t = 0 ; t < myT ; t++)
        {       myGradNorm[t].ReAlloc(myNNormParam) ;
                myHessNorm[t].ReAlloc(myNNormParam, myNNormParam) ;
        }

        for (uint j = 0 ; j < mvNClass ; j++)
        {
                for (uint t = 0 ; t < myT ; t++)
                {       theGrad[j][t] = 0.0 ;
                        theHess[j][t] = 0.0 ;
                }

        cDMatrix myInvCov(mvDimObs, mvDimObs) ;
        double myDeterminant ;
        uint k = (mvNClass - 1)*(mvNClass + 1) + j*myNParam ;
                
                ArmadilloInvAndDet(mCov[j][mvNMixt-1], myInvCov, myDeterminant) ;
                MultivariateNormalDensity(theY, mMean[j][mvNMixt-1], myInvCov, myDeterminant, myLastDens) ;     
                for (uint l = 0 ; l < mvNMixt ; l++)
                {       ArmadilloInvAndDet(mCov[j][l], myInvCov, myDeterminant) ;
                        MultivariateNormalDensity(theY, mMean[j][l], myInvCov, myDeterminant, myDens) ;
                        MultivariateNormalDensityDeriv(theY, mMean[j][l], mCov[j][l], myInvCov, myDeterminant, myDens, myGradNorm, myHessNorm) ;
                        for (uint t = 0 ; t < myT ; t++)
                        {       SetSubVector(mp[j][l]*myGradNorm[t], k, theGrad[j][t]) ;
                        cDMatrix myAuxMat = mp[j][l]*myHessNorm[t] ;    
                                SetSubMatrix(myAuxMat, k, k, theHess[j][t]) ;
                        
                                if (l < mvNMixt - 1)
                                {       theGrad[j][t][k+myNNormParam] = myDens[t] - myLastDens[t] ;
                                        for (uint p = 0 ; p < myNNormParam ; p++)
                                                theHess[j][t][k+myNNormParam][p+k] = theHess[j][t][p+k][k+myNNormParam] = myGradNorm[t][p] ;
                                }                                       
                        }
                        k += myNNormParam ;
                        if (l < mvNMixt - 1)
                                k++ ;
                }
        }

        for (uint t = 0 ; t < myT ; t++)
        {       myGradNorm[t].Delete() ;
                myHessNorm[t].Delete() ;
        }
        delete [] myGradNorm ;
        delete [] myHessNorm ;

}

void cMixtMultivariateNormal::ComputeCov(cDMatrix& theCov)
{
uint myBegIndex = (mvNClass - 1) * (mvNClass + 1) ;
uint myNNormParam = mvDimObs + mvDimObs*(mvDimObs + 1)/2 ;
uint myNFreeMixt = (myNNormParam + 1)*mvNMixt - 1 ;
uint mySizeCour = theCov.GetNCols() ;
cDVector myU(mySizeCour, 0.0) ;

        for (uint n = 0 ; n < mvNClass ; n++)
        {
                for (uint i = myBegIndex + myNNormParam ; i < myBegIndex + myNFreeMixt ; i+=myNNormParam+1)
                        myU[i] = -1.0 ;
                theCov = AddOneVariable(theCov, myU) ;
                mySizeCour++ ;
                myU.ReAlloc(mySizeCour, 0.0) ;
                myBegIndex += myNFreeMixt ;
        }

}

cDVector cMixtMultivariateNormal::GetDistrNumParam(const cDVector& theNumDistrParam, uint& theNextInd)
{
uint myNNormParam = mvDimObs + mvDimObs*(mvDimObs + 1)/2 ;
uint myNFreeParam =(myNNormParam + 1)*mvNMixt - 1 ;
cDVector myNumDistrParam ;
cDVector myNumMixt(myNFreeParam) ;
uint myIndCour = 0 ;
        for (uint j = 0 ; j < mvNClass ; j++)
        {       GetSubVector(theNumDistrParam, myIndCour, myNFreeParam, myNumMixt) ;
                myNumDistrParam = cat(myNumDistrParam, myNumMixt) ;
                myNumDistrParam = cat(myNumDistrParam, (double)theNextInd) ;
                theNextInd++ ;
                myIndCour += myNFreeParam ;
        }
        return myNumDistrParam ;
}

void cMixtMultivariateNormal::UpdateParameters(cInParam& theInParam, cBaumWelch& theBaumWelch, cDMatrix* theCondProba)
{       
cDMatrix*      myInvCov = new cDMatrix[mvNMixt];
double*         myDet = new double[mvNMixt] ;

        for (uint i = 0 ; i < mvNMixt ; i++)
                myInvCov[i].ReAlloc(mvDimObs, mvDimObs) ;

        for (uint i = 0 ; i < mvNClass ; i++)
        {       
        double mySumGammai = 0.0 ;
        uint   n,
                                        t       ;
                
                for (n = 0 ; n < theInParam.mNSample ; n++)
                {       uint myT = theInParam.mY[n].mSize / mvDimObs ;
                        for (t = 0 ; t < myT  ; t++)
                                mySumGammai += theBaumWelch.mGamma[n][i][t] ;
                }

                for (uint j = 0 ; j < mvNMixt ; j++)
                        SymetricInverseAndDet(mCov[i][j], myDet[j], myInvCov[j]) ;
                
        cDVector myMoy = cDVector(mvDimObs) ;
        cDMatrix myCov = cDMatrix(mvDimObs, mvDimObs) ;
                for (uint l = 0 ; l < mvNMixt  ; l++)
                {                               
                        myMoy = 0.0 ;
                        myCov = 0.0 ;
                        double myGammail ;
                        double mySumGammail = 0.0 ;
                        for (n = 0 ; n < theInParam.mNSample ; n++)
                        {       
                        uint myT = theInParam.mY[n].mSize / mvDimObs ;
                        double* myDens = new double[myT] ; 
                                MultivariateNormalDensity(theInParam.mY[n], mMean[i][l], myInvCov[l], myDet[l], myDens) ;
                                for (t = 0 ; t < myT ; t++)
                                {       myGammail = theBaumWelch.mGamma[n][i][t] * mp[i][l] * myDens[t] / theCondProba[n][i][t] ;
                                        mySumGammail += myGammail ; 
                                        for (uint k = 0 ; k < mvDimObs ; k++)
                                        {       myMoy[k] += myGammail * theInParam.mY[n][t+k*myT] ;
                                                for (uint j = k ; j < mvDimObs ; j++)
                                                        myCov[k][j] += myGammail * theInParam.mY[n][t+k*myT] * theInParam.mY[n][t+j*myT] ;
                                        }
                                }
                                delete [] myDens ;
                        }
                        mp[i][l] = mySumGammail / mySumGammai ;
                        mMean[i][l] = myMoy/mySumGammail ;
                        for (int m = 0 ; m < (int)mvDimObs-1 ; m++)
                                for (int p = m+1 ; l < mvDimObs ; p++)
                                        myCov[p][m] = myCov[m][p] ;
                        mCov[i][l] = myCov/mySumGammail ;
                        mCov[i][l] -= mMean[i][l] * Transpose(mMean[i][l]) ;
                }
        }
}               

void cMixtMultivariateNormal::InitParameters(cBaumWelchInParam &theInParam)
{
#ifdef _RDLL_
        GetRNGstate();
#endif //_RDLL_

cDVector       myMoy(mvDimObs),
                        myVar(mvDimObs),
                        myStd(mvDimObs) ;

double mys = 0.0 ;
        for (uint n = 0 ; n < theInParam.mNSample ; n++)
        {
        uint myT = theInParam.mY[n].mSize / mvDimObs ;
                for (uint t = 0 ; t < myT  ; t++)
                {       for (uint i = 0 ; i < mvDimObs ; i++)
                        {       myMoy[i] = (mys*myMoy[i] + theInParam.mY[n][t+i*myT])/(mys+1) ;
                                myVar[i] = (mys*myVar[i] + theInParam.mY[n][t+i*myT]*theInParam.mY[n][t+i*myT])/(mys+1) ;
                        }
                        mys++ ;
                }
        }

        for (uint i = 0 ; i < mvDimObs ; i++)
        {       myVar[i] -= myMoy[i]*myMoy[i] ;
                myStd[i] = sqrt(myVar[i]) ;
        }

        for (uint i = 0 ; i < mvNClass ; i++)
        {       
        double mySomme = 0.0 ;
        uint l ;
                for (l = 0 ; l < mvNMixt ; l++)
                {
                                /* FIXME:
                                 * Zeros returns a global reference, but mCov[i][l] is a instance, so this
                                 * is fine, but weird of course.
                                 */
                                mCov[i][l] = Zeros(mCov[i][l].mNRow,mCov[i][l].mNCol);

                                for (uint k = 0 ; k < mvDimObs ; k++)
                        {       mMean[i][l][k] =  -2*myStd[k] + myMoy[k] + 2*myStd[k] * unif_rand() ;
                                mCov[i][l][k][k] = 0.5*myVar[k] + 3*myVar[k] * unif_rand() ;     
                        }
                        mp[i][l] = unif_rand() ;        
                        mySomme += mp[i][l] ;
                }
                for (l = 0 ; l < mvNMixt ; l++)
                        mp[i][l] /= mySomme ;
        }

                        
#ifndef _RDLL_
        PutRNGstate() ;
#endif //_RDLL_
}

void cMixtMultivariateNormal::CopyDistr(cDistribution* theSrc)
{
cMixtMultivariateNormal* mySrc = dynamic_cast<cMixtMultivariateNormal *>(theSrc) ;
        if (mySrc)
        {       mvNClass = mySrc->mvNClass ;
                mvDimObs = mySrc->mvDimObs ;
                mvNMixt = mySrc->mvNMixt ;
                for (uint i = 0 ; i < mvNClass ; i++)
                {       for (uint l = 0 ; l < mvNMixt ; l++)
                        {       mMean[i][l] = mySrc->mMean[i][l] ;
                                mCov[i][l] = mySrc->mCov[i][l] ;
                        }
                        mp[i] = mySrc->mp[i] ;
                }
        }
        else
                cOTError("Wrong distribution in cMixtMultivariateNormal") ;
}

cMixtMultivariateNormal::cMixtMultivariateNormal(cDistribution& theSrc)
{
        CopyDistr(&theSrc) ;
}

void cMixtMultivariateNormal::Print()
{
        Rprintf("Parameters\n") ;
        for (uint i = 0 ; i < mvNClass ; i++)
        {       Rprintf("State %d\n", i) ;
                for (uint j = 0 ; j < mvNMixt ; j++)
                {       Rprintf("p[%d]=%lf\nEsp[%d]\t\tMatCov[%d]\n", j, mp[i][j], j, j) ;
                        for (uint k = 0 ; k < mvDimObs ; k++)
                        {       Rprintf("%lf\t", mMean[i][j][k]) ;
                                for (uint l = 0 ; l < mvDimObs ; l++)
                                        Rprintf("\t%lf", mCov[i][j][k][l]) ;
                                Rprintf("\n") ;
                        }
                }
                Rprintf("\n") ;
        }
}

void cMixtMultivariateNormal::GetParam(uint theDeb, cDVector& theParam)
{
uint k = theDeb ;
        for (uint n = 0 ; n < mvNClass ; n++)
        {       for (uint p = 0 ; p < mvNMixt ; p++)
                {       for (uint m = 0 ; m < mvDimObs ; m++)
                                theParam[k++] = mMean[n][p][m] ;
                        for (uint m = 0 ; m < mvDimObs ; m++)
                                for (uint l = m ; l < mvDimObs ; l++)
                                        theParam[k++] = mCov[n][p][m][l] ;
                        if (p < mvNMixt-1)
                                theParam[k++] = mp[n][p] ;
                }
        }
}

void cMixtMultivariateNormal::SetParam(uint theDeb, cDVector& theParam)
{
uint k = theDeb ;
        for (uint n = 0 ; n < mvNClass ; n++)
        {       mp[n][mvNMixt-1] = 1.0 ;
                for (uint p = 0 ; p < mvNMixt ; p++)
                {       for (uint i = 0 ; i < mvDimObs ; i++)
                                mMean[n][p][i] = theParam[k++] ;
                        for (uint i = 0 ; i < mvDimObs ; i++)
                                for (uint j = i ; j < mvDimObs ; j++)
                                        mCov[n][p][i][j] = mCov[n][p][j][i] = theParam[k++] ;
                        if (p < mvNMixt-1)
                        {       mp[n][p] = theParam[k++] ;
                                mp[n][mvNMixt-1] -= mp[n][p] ;
                        }
                }
        }
}

