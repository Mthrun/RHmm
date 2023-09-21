/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: MultivariateNormalUtil.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 *** 
 *** Util methods for calculation of multivariate normal density                                                      
 **************************************************************/

#include "StdAfxRHmm.h"

void SymetricInverseAndDet(cDMatrix& theMat, double& theDet, cDMatrix& theInvMat)
{
  ArmadilloInvAndDet(theMat, theInvMat, theDet) ;
}

/*
 * Calculates the mutlivariate normal density, given an vector X, expectation
 * vector mu, the inverse covariance matrix and the determinant of the covariance matrix.
 * The parameter theDens will be filled with the density, where theDens is an double array
 * @param thex X vector
 * @param theMu Expectation vector mu
 * @param theInvCov Inverse covariance matrix
 * @param theDet Determinant of the covariance matrix
 */
void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, double*  theDens)
{
uint i, j, t ;
double  myAux, myRapport ;

uint myDimObs = theMu.mSize ;   
        myRapport = pow(SQRT_TWO_PI, (int)myDimObs)*sqrt(theDet) ;

uint myT = thex.mSize / myDimObs ;

        for ( t = 0 ; t < myT ; t++)
        {       myAux = 0.0 ;
                for (i = 0 ; i < myDimObs ; i++)
                        for (j = 0 ; j < myDimObs ; j++)
                                myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
                theDens[t] = exp(-0.5*myAux)/myRapport ;
        }
}

/*
 * Calculates the mutlivariate normal density, given an vector X, expectation
 * vector mu, the inverse covariance matrix and the determinant of the covariance matrix.
 * The parameter theDens will be filled with the density, where theDens is an cDVector object
 * @param thex X vector
 * @param theMu Expectation vector mu
 * @param theInvCov Inverse covariance matrix
 * @param theDet Determinant of the covariance matrix
 */
void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, cDVector&  theDens)
{
uint i, j, t ;
double  myAux, myRapport ;

uint myDimObs = theMu.mSize ;   
        myRapport = pow(SQRT_TWO_PI, (int)myDimObs)*sqrt(theDet) ;

uint myT = thex.mSize / myDimObs ;

        for ( t = 0 ; t < myT ; t++)
        {       myAux = 0.0 ;
                for (i = 0 ; i < myDimObs ; i++)
                        for (j = 0 ; j < myDimObs ; j++)
                                myAux += (thex[t+i*myT]-theMu[i]) * theInvCov[i][j] * (thex[t+j*myT]-theMu[j]) ;
                theDens[t] = exp(-0.5*myAux)/myRapport ;
        }
}

void InvCovMatDeriv(cDMatrix& theInvCov, cDMatrix* theGrad, cDMatrix** theHess)
{
uint mySize = theInvCov.GetNCols() ;
cDMatrix myCovPrime = Zeros(mySize, mySize) ;
cDMatrix myCovSeconde = Zeros(mySize, mySize) ;
//uint myNParam = mySize*(mySize+1)/2 ;
uint k = 0 ;
        for (uint i = 0 ; i < mySize ; i++)
                for (uint j = i ; j < mySize ; j++)
                {       myCovPrime[i][j] = myCovPrime[j][i] = 1.0 ;
                cDMatrix myAuxMat = myCovPrime * theInvCov ;
                        theGrad[k] = -1.0 * theInvCov * myAuxMat ;
                uint l = 0 ;
                        for (uint p = 0 ; p < mySize ; p++)
                                for (uint q = p ; q < mySize ; q++)
                                {       myCovSeconde[p][q] = myCovSeconde[q][p] = 1.0 ;
                                        theHess[k][l] = theHess[l][k] = -1 * theInvCov * myCovSeconde * theGrad[k] - theGrad[k] * myCovSeconde * theInvCov ;
                                        l++ ;
                                        myCovSeconde[p][q] = myCovSeconde[q][p] = 0.0 ;
                                }
                        k++ ;
                        myCovPrime[i][j] = myCovPrime[j][i] = 0.0 ;
                }
}

void SymDetDeriv(cDMatrix& theMat, cDVector& theGrad, cDMatrix& theHess)
{
// Gradient
cDMatrix myAuxMat = theMat ;
uint myNCol = theMat.GetNCols() ;
uint k = 0 ;
        for (uint i = 0 ; i < myNCol ; i++)
        {       for (uint j = i ; j < myNCol ; j++)
                {       myAuxMat[i][j] = myAuxMat[j][i] = 0 ;
                double myG0 = ArmadilloDet(myAuxMat) ;
                        myAuxMat[i][j] = myAuxMat[j][i] = 1.0 ;
                double myG1 = ArmadilloDet(myAuxMat) ;
                        myAuxMat[i][j] = myAuxMat[j][i] = -1.0 ;
                double myGm1 = ArmadilloDet(myAuxMat) ;
                double myA = (myG1 + myGm1)/2.0 - myG0 ;
                double myB = (myG1 - myGm1)/2.0 ;
                        theGrad[k++] = 2.0*myA*theMat[i][j] + myB ;
                        myAuxMat[i][j] = myAuxMat[j][i] = theMat[i][j] ;
                }
        }
// Hessian
cDMatrix myInvZ = Zeros(9,9) ;  
      myInvZ[0][0] = 0.25 ;
      myInvZ[0][1] = -0.5 ;
      myInvZ[0][2] = 0.25 ;
      myInvZ[0][3] = -0.5 ;
      myInvZ[0][4] = 1.0 ;
      myInvZ[0][5] = -0.5 ;
      myInvZ[0][6] = 0.25 ;
      myInvZ[0][7] = -0.5 ;
      myInvZ[0][8] = 0.25 ;
      myInvZ[1][0] = -0.25 ;
      myInvZ[1][2] = 0.25 ;
      myInvZ[1][3] = 0.5 ;
      myInvZ[1][5] = -0.5 ;
      myInvZ[1][6] = -0.25 ;
      myInvZ[1][8] = 0.25 ;
      myInvZ[2][0] = -0.25 ;
      myInvZ[2][1] = 0.5 ;
      myInvZ[2][2] = -0.25 ;
      myInvZ[2][6] = 0.25 ;
      myInvZ[2][7] = -0.5 ;
      myInvZ[2][8] = 0.25 ;
      myInvZ[3][1] = 0.5 ;
      myInvZ[3][4] = -1.0 ;
      myInvZ[3][7] = 0.5 ;
      myInvZ[4][0] = 0.25 ;
      myInvZ[4][2] = -0.25 ;
      myInvZ[4][6] = -0.25 ;
      myInvZ[4][8] = 0.25 ;
      myInvZ[5][3] = 0.5 ;
      myInvZ[5][4] = -1.0 ;
      myInvZ[5][5] = 0.5 ;
      myInvZ[6][1] = -0.5 ;
      myInvZ[6][7] = 0.5 ;
      myInvZ[7][3] = -0.5 ;
      myInvZ[7][5] = 0.5 ;
      myInvZ[8][4] = 1.0 ;

cDMatrix myInd(9,2) ;
        k = 0 ;
                for (int i =-1 ; i < 2 ; i++)
                        for (int j = -1 ; j < 2 ; j++)
                        {       myInd[k][0] = (double)i ;
                                myInd[k][1] = (double)j ;
                                k++ ;
                        }

        k = 0 ;
        myAuxMat = theMat ;
cDVector myG(9) ;
        for (uint i = 0 ; i < myNCol ; i++)
        {       for (uint j = i ; j < myNCol ; j++)
                {       
                uint l = 0 ;
                double myx = theMat[i][j] ;
                        for (uint p = 0 ; p < myNCol ; p++)
                        {       for (uint q = p ; q < myNCol ; q++)
                                {       
                                double myy = theMat[p][q] ;     
                                        for (uint r = 0 ; r < 9 ; r++)
                                        {       myAuxMat[i][j]= myAuxMat[j][i] = myInd[r][0] ;
                                                myAuxMat[p][q]=myAuxMat[q][p] = myInd[r][1] ;
                                                myG[r] = ArmadilloDet(myAuxMat) ;
                                                myAuxMat[i][j] = myAuxMat[j][i] = myx ;
                                                myAuxMat[p][q] = myAuxMat[q][p] = myy ;
                                        }
                                cDVector myCoef = myInvZ * myG ;
                                        if ( (i == p) && (j == q) )
                                                theHess[k][l] = theHess[l][k] = 2.0*(6.0*myCoef[0]*myx*myx+3.0*(myCoef[1]+myCoef[2])*myx+myCoef[3]+myCoef[4]+myCoef[5]) ;
                                        else
                                                theHess[k][l] = theHess[l][k] = 4.0*myCoef[0]*myx*myy+2.0*(myCoef[1]*myx+myCoef[2]*myy)+myCoef[4] ;
                                        l++ ;
                                } // for q
                        } // for p
                        k++ ;
                } // for j
        } // for i
}

/*
 * Computation of the derivative of the multivariate normal density (Hessian matrix of density of
 * multivariate normal distributed vector X), as well as the computation of the gradient.
 * theGrad (Gradient) and the Hessian matrix (theHess) are filled in the function.
 * This function will calculate the multivariate density before calculating the derivative and is a wrapper for
 * MultivariateNormalDensityDeriv function with the density parameter.
 * @param thex X vector
 * @param theMu Expectation value vector mu
 * @param theCov The covariance matrix
 * @param theInvCov The inverse covariance matrix
 * @param theDet Determinant of the covariance matrix
 * @param theGrad Gradient vector of density, filled in the function
 * @param theHess Hessian matrix of density, filled in the function
 */
void MultivariateNormalDensityDeriv(cDVector& thex, cDVector& theMu, cDMatrix& theCov, cDMatrix& theInvCov, double theDet, cDVector* theGrad, cDMatrix* theHess)
{
uint myDimObs = theMu.mSize ;   
uint myT = thex.mSize / myDimObs ;      

cDVector myDens(myT) ;
        MultivariateNormalDensity(thex, theMu, theInvCov, theDet, myDens) ;
        MultivariateNormalDensityDeriv(thex, theMu, theCov, theInvCov, theDet, myDens, theGrad, theHess) ;
}


/*
 * Computation of the derivative of the multivariate normal density (Hessian matrix of density of
 * multivariate normal distributed vector X), as well as the computation of the gradient.
 * theGrad (Gradient) and the Hessian matrix (theHess) are filled in the function.
 * This function will calculate the multivariate density before calculating the derivative.
 * @param thex X vector
 * @param theMu Expectation value vector mu
 * @param theCov The covariance matrix
 * @param theInvCov The inverse covariance matrix
 * @param theDet Determinant of the covariance matrix
 * @param theDensity The density vector
 * @param theGrad Gradient vector of density, filled in the function
 * @param theHess Hessian matrix of density, filled in the function
 */
void MultivariateNormalDensityDeriv(cDVector& thex, cDVector& theMu, cDMatrix& theCov, cDMatrix& theInvCov, double theDet, cDVector& theDensity, cDVector* theGrad, cDMatrix* theHess)
{
uint myDimObs = theMu.mSize ;   
uint myNCovParam = myDimObs*(myDimObs + 1)/2 ;
uint myT = theDensity.GetSize() ;
        
cDVector myGradDet(myNCovParam) ;
cDMatrix myHessDet(myNCovParam, myNCovParam) ;
        SymDetDeriv(theCov, myGradDet, myHessDet) ;
cDMatrix* myGradInvCov = new cDMatrix[myNCovParam] ;
cDMatrix** myHessInvCov = new cDMatrix*[myNCovParam] ;
        for (uint k = 0 ; k < myNCovParam ; k++)
        {       myGradInvCov[k].ReAlloc(myNCovParam, myNCovParam) ;
                myHessInvCov[k] = new cDMatrix[myNCovParam] ;
                for (uint l = 0 ; l < myNCovParam ; l++)
                        myHessInvCov[k][l].ReAlloc(myNCovParam, myNCovParam) ;
        }
        InvCovMatDeriv(theInvCov, myGradInvCov, myHessInvCov) ;

//
        for (uint t = 0 ; t < myT ; t++)
        {
        // Deriv?e par rapport ? Mu
        cDVector myx(myDimObs) ;
                for (uint i = 0 ; i < myDimObs ; i++)
                        myx[i] = (thex[t+i*myT]-theMu[i]) ;
        
        cDVector myGradMu = theDensity[t]*theInvCov*myx ;
                SetSubVector(myGradMu, 0, theGrad[t]) ;
        
        // D?riv?e par rapport ? Cov
        cDVector myGradLambda(myNCovParam) ;
                for (uint k = 0 ; k < myNCovParam ; k++)
                        myGradLambda[k]  = AsDouble(Transpose(myx) * myGradInvCov[k] * myx) /-2.0 ;
        
                myGradLambda -= 0.5*myGradDet/theDet ;
                myGradLambda *= theDensity[t] ;
                SetSubVector(myGradLambda, myDimObs, theGrad[t]) ;

        // d?riv?e seconde par rapport ? Mu
        cDMatrix myHessMu2 = theInvCov*myx*Transpose(myx)*theInvCov ;
                myHessMu2 -= theInvCov ;
                myHessMu2 *= theDensity[t] ;
                SetSubMatrix(myHessMu2, 0, 0, theHess[t]) ;

   // d?riv?e seconde par rapport ? Mu et ? Cov
        cDMatrix myHessMuLambda(myNCovParam, myDimObs) ;
        cDVector myAuxVect1 = theInvCov * myx ;
                for (uint k = 0 ; k < myNCovParam ; k++)
                {
                cDVector myAuxVect2 = myGradInvCov[k] * myx * theDensity[t] ;
                        for (uint i = 0 ; i < myDimObs ; i++)
                                myHessMuLambda[k][i] = myAuxVect2[i] + myAuxVect1[i]*myGradLambda[k] ;
                }
                SetSubMatrix(myHessMuLambda, myDimObs, 0, theHess[t]) ;
        cDMatrix myHessMuLambdaPrime = Transpose(myHessMuLambda) ;
                SetSubMatrix(myHessMuLambdaPrime, 0, myDimObs, theHess[t]) ;

        // d?riv?e seconde par rapport ? Cov
        cDMatrix myHessLambda2 = -0.5*theDensity[t]*myHessDet/theDet ;
                myHessLambda2 += 0.5*theDensity[t]*myGradDet*Transpose(myGradDet)/(theDet*theDet) ;
                if (theDensity[t] != 0.0)
                        myHessLambda2 += myGradLambda * Transpose(myGradLambda)/theDensity[t] ;
                for (uint k = 0 ; k < myNCovParam ; k++)
                        for (uint l = k ; l < myNCovParam ; l++)
                        {
                        double myDoubleAux =  AsDouble(Transpose(myx)*myHessInvCov[k][l] * myx) * theDensity[t] ;
                                myHessLambda2[k][l] -= 0.5*myDoubleAux ;
                                if (l != k)
                                        myHessLambda2[l][k] -= 0.5*myDoubleAux ;
                        }
                SetSubMatrix(myHessLambda2, myDimObs, myDimObs, theHess[t]) ;
        }

        for (uint i = 0 ; i < myNCovParam ; i++)
        {       myGradInvCov[i].Delete() ;
                for (uint j = 0 ; j < myNCovParam ; j++)
                        myHessInvCov[i][j].Delete() ;
                delete [] myHessInvCov[i] ;
        }

        delete [] myGradInvCov ;
        delete [] myHessInvCov ;
}
