/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cDerivative.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#include "StdAfxRHmm.h"

/*
 * Constructor for cDerivate
 * @param theNSample Number of samples
 * @param theT Index set
 * @param theNClass Number of classes
 * @param theNFreeParam Number of free parameters
 */
cDerivative::cDerivative(uint theNSample, uint* theT, uint theNClass, uint theNFreeParam)
{
        mvNFreeParam = theNFreeParam ;
        mvNClass = theNClass ;
        mvNSample = theNSample ;
        mvT = new uint[theNSample] ;

        mPsi = new cDVector**[theNSample] ;
        mOmega = new cDMatrix**[theNSample] ;
        mScore = new cDVector[theNSample] ;
        mInformation = new cDMatrix[theNSample] ;

        for (uint n = 0 ; n < theNSample ; n++)
    {   mPsi[n] = new cDVector*[mvNClass] ;
                mOmega[n] = new cDMatrix*[mvNClass] ;
                mScore[n].ReAlloc(mvNFreeParam) ;
                mvT[n] = theT[n] ;
                mInformation[n].ReAlloc(mvNFreeParam, mvNFreeParam) ;
                for (uint i = 0 ; i < mvNClass ; i++)
                {       mPsi[n][i] = new cDVector[theT[n]] ;
                        mOmega[n][i] = new cDMatrix[theT[n]] ;
                        for (uint t = 0 ; t < mvT[n] ; t++)
                        {       mPsi[n][i][t].ReAlloc(mvNFreeParam, 0.0) ;
                                mOmega[n][i][t].ReAlloc(mvNFreeParam, mvNFreeParam) ;
                        }
                }
        }
        MESS_CREAT("cDerivative") 
}

/*
 * Constructor for cDerivate
 * @param theInParam cInParam object containing parameters for HMM
 * @param theNFreeParam Number of free parameters
 */
cDerivative::cDerivative(const cInParam &theInParam, uint theNFreeParam)
{
        mvNFreeParam = theNFreeParam ;
        mvNClass = theInParam.mNClass ;
        mvNSample = theInParam.mNSample ;
        mvT = new uint[mvNSample] ;

        mPsi = new cDVector**[mvNSample] ;
        mOmega = new cDMatrix**[mvNSample] ;
        mScore = new cDVector[mvNSample] ;
        mInformation = new cDMatrix[mvNSample] ;

        for (uint n = 0 ; n < mvNSample ; n++)
    {   
        uint myT = theInParam.mY[n].GetSize() / theInParam.mDimObs ;
                mPsi[n] = new cDVector*[mvNClass] ;
                mOmega[n] = new cDMatrix*[mvNClass] ;
                mScore[n].ReAlloc(mvNFreeParam) ;
                mvT[n] = myT ;
                mInformation[n].ReAlloc(mvNFreeParam, mvNFreeParam) ;
                for (uint j = 0 ; j < mvNClass ; j++)
                {       mPsi[n][j] = new cDVector[myT] ;
                        mOmega[n][j] = new cDMatrix[myT] ;      
                        for (uint t = 0 ; t < myT ; t++)
                        {       mPsi[n][j][t].ReAlloc(mvNFreeParam) ;
                                mOmega[n][j][t].ReAlloc(mvNFreeParam, mvNFreeParam) ;
                        }
                }
        }
        MESS_CREAT("cDerivative") 
}

cDerivative::~cDerivative()
{
        for (uint n = 0 ; n < mvNSample ; n++)
        {       for (uint j = 0 ; j < mvNClass ; j++)
                {       for (uint t = 0 ; t < mvT[n] ; t++)
                        {       mPsi[n][j][t].Delete() ;
                                mOmega[n][j][t].Delete() ;
                        }
                        delete [] mPsi[n][j] ;
                        delete [] mOmega[n][j] ;
                }
                delete [] mPsi[n] ;
                delete [] mOmega[n] ;
                mScore[n].Delete() ;
                mInformation[n].Delete() ;
        }
        delete [] mPsi ;
        delete [] mOmega ;
        delete [] mScore ;
        delete [] mInformation ;
        MESS_DESTR("cDerivative") 
}

/*
 * Computation of the derivative of the mixture multivariate normal density (Hessian matrix of density of
 * mixture multivariate normal distributed vector Y), as well as the computation of the gradient.
 * theGrad (Gradient) and the Hessian matrix (theHess) are filled in the function.
 * @param theY Y vector
 * @param theGrad Gradient vector of density, filled in the function
 * @param theHess Hessian matrix of density, filled in the function
 */
void cDerivative::ComputeDerivative(cHmm& theHmm, cInParam& theInParam)
{
cDMatrix* myCondProba = new cDMatrix[theInParam.mNSample] ; 
                
        for (uint n = 0 ; n < mvNSample ; n++)
        {          
        uint mySize = theInParam.mY[n].mSize/theInParam.mDimObs ;
                myCondProba[n].ReAlloc(theInParam.mNClass, mySize) ;
        }
        theHmm.mDistrParam->ComputeCondProba(theInParam.mY, mvNSample, myCondProba) ;

cDMatrix* myLambda = new cDMatrix[theInParam.mNSample] ;
cDVector* mySumLambda = new cDVector[theInParam.mNSample] ;
        
        // CALCUL DES LAMBDAs
        for (uint n = 0 ; n < mvNSample ; n++)
        {
                mySumLambda[n].ReAlloc(mvT[n]) ;
                myLambda[n].ReAlloc(mvNClass, mvT[n]) ;
        
        // t = 0 
                mySumLambda[n][0] = 0.0 ;
                for (uint j = 0 ; j < mvNClass ; j++)
                {
                        myLambda[n][j][0] = theHmm.mInitProba[j] * myCondProba[n][j][0] ;
                        mySumLambda[n][0] += myLambda[n][j][0] ;
                }

        // t > 0
        uint myT = theInParam.mY[n].GetSize() / theInParam.mDimObs ;
                for (uint t = 1 ; t < myT ; t++)
                {       mySumLambda[n][t] = 0.0 ;
                        for (uint j = 0 ; j < mvNClass ; j++)
                        {       
                        double myAux = 0.0 ;
                                for (uint i = 0 ; i < mvNClass ; i++)
                                        myAux += myLambda[n][i][t-1]*theHmm.mTransMatVector[t][i][j] ;
                                myLambda[n][j][t] = myAux * myCondProba[n][j][t]/mySumLambda[n][t-1] ;
                                mySumLambda[n][t] += myLambda[n][j][t] ;
                        }
                }
        }

cDVector* myGradInitProb = new cDVector[mvNClass] ;
cDVector** myGradTransMat = new cDVector*[mvNClass] ;
cDVector** myGradCondProba = new cDVector*[mvNClass] ;
cDMatrix** myHessCondProba = new cDMatrix*[mvNClass] ;

        for (uint n = 0 ; n < mvNClass ; n++)
        {       myGradInitProb[n].ReAlloc(mvNFreeParam, 0.0L) ;
                myGradTransMat[n] = new cDVector[mvNClass] ;
                for (uint i = 0 ; i < mvNClass ; i++)
                        myGradTransMat[n][i].ReAlloc(mvNFreeParam, 0.0L) ;
        }

        /* D?riv?es probabilit?s initiales */
uint myNFreeClass = mvNClass - 1 ;
        for (uint s = 0 ; s < myNFreeClass ; s++)
        {       myGradInitProb[s][s] = 1.0L ;
                myGradInitProb[myNFreeClass][s] = -1.0L ;
        }
        
        /* D?riv?es matrice de transition */
uint myBegIndex = myNFreeClass  ;
        for (uint i = 0 ; i < mvNClass ; i++)
        {       for (uint j = 0 ; j < myNFreeClass  ; j++)
                {       myGradTransMat[i][j][j+myBegIndex] = 1.0L ;
                        myGradTransMat[i][myNFreeClass][j+myBegIndex] = -1.0L ;
                }
                myBegIndex += myNFreeClass ;
        }

        for (uint n = 0 ; n < mvNSample ; n++) // Boucle sur le nombre d'?chantillon
        {
                for (uint j = 0 ; j < mvNClass ; j++)
                {       myGradCondProba[j] = new cDVector[mvT[n]] ;
                        myHessCondProba[j] = new cDMatrix[mvT[n]] ;
                        for (uint t = 0 ; t < mvT[n] ; t++)
                        {       myGradCondProba[j][t].ReAlloc(mvNFreeParam) ;
                                myHessCondProba[j][t].ReAlloc(mvNFreeParam, mvNFreeParam) ;
                        }
                }
                theHmm.mDistrParam->ComputeDerivative(theInParam.mY[n], myGradCondProba, myHessCondProba) ;
        /* t = 0 */

                for (uint j = 0 ; j < mvNClass ; j++)
                {
                        mPsi[n][j][0] = theHmm.mInitProba[j] * myGradCondProba[j][0] + myCondProba[n][j][0] * myGradInitProb[j] ;
                        mOmega[n][j][0] = theHmm.mInitProba[j] * myHessCondProba[j][0] + myGradCondProba[j][0] * Transpose(myGradInitProb[j]) 
                                                        + myGradInitProb[j] * Transpose(myGradCondProba[j][0]) ;                
                }

        /* t > 0 */
        uint myT = mvT[n] ;

                for (uint t = 1 ; t < myT ; t++)
                {       for (uint j = 0 ; j < mvNClass ; j++)
                        {       mPsi[n][j][t] = 0.0 ;
                                mOmega[n][j][t] = 0.0 ;
                                for ( uint i = 0 ; i < mvNClass ; i++)
                                {
                                cDVector myVect1 = myCondProba[n][j][t] * theHmm.mTransMatVector[t][i][j] * mPsi[n][i][t-1] ;
                                cDVector myVect2 = myLambda[n][i][t-1] * theHmm.mTransMatVector[t][i][j] * myGradCondProba[j][t]  ;
                                cDVector myVect3 = myLambda[n][i][t-1] * myCondProba[n][j][t] * myGradTransMat[i][j] ;
                                        
                                        mPsi[n][j][t] += myVect1 + myVect2 + myVect3 ;

                                cDMatrix myMat1 = myCondProba[n][j][t] * theHmm.mTransMatVector[0][i][j] * mOmega[n][i][t-1] ;
                                cDMatrix myMat2 = theHmm.mTransMatVector[t][i][j] * (mPsi[n][i][t-1] * Transpose(myGradCondProba[j][t]) + myGradCondProba[j][t] * Transpose(mPsi[n][i][t-1])) ;
                                cDMatrix myMat3 = myCondProba[n][j][t] * (mPsi[n][i][t-1] * Transpose(myGradTransMat[i][j]) + myGradTransMat[i][j] * Transpose(mPsi[n][i][t-1]))  ;
                                cDMatrix myMat4 = myLambda[n][i][t-1] * (myGradCondProba[j][t] * Transpose(myGradTransMat[i][j]) + myGradTransMat[i][j] * Transpose(myGradCondProba[j][t])) ;
                                cDMatrix myMat5 =  myLambda[n][i][t-1] * theHmm.mTransMatVector[t][i][j] * myHessCondProba[j][t] ;
                                mOmega[n][j][t] += myMat1 + myMat2 + myMat3 + myMat4 + myMat5 ;
                                }
                                mPsi[n][j][t] /= mySumLambda[n][t-1] ;
                                mOmega[n][j][t] /= mySumLambda[n][t-1] ;
                        }
                } /* for t */
                mScore[n] = 0.0 ;
                mInformation[n] = 0.0 ;
                for ( uint j = 0 ; j < mvNClass ; j++)
                {       mScore[n] += mPsi[n][j][myT-1] ;
                        mInformation[n] -= mOmega[n][j][myT-1] ;
                }
                mScore[n] /= mySumLambda[n][myT-1] ;
                mInformation[n] /= mySumLambda[n][myT-1] ;
                mInformation[n] += mScore[n] * Transpose(mScore[n]) ;

                for (uint p = 0 ; p < mvNClass ; p++)
                {       for (uint t = 0 ; t < mvT[n] ; t++)
                        {       myGradCondProba[p][t].Delete() ;
                                myHessCondProba[p][t].Delete() ;
                        }
                        delete [] myGradCondProba[p] ;
                        delete [] myHessCondProba[p] ;
                }

        } /* for n */

        for (uint n = 0 ; n < mvNClass ; n++)
        {       myGradInitProb[n].Delete() ;
                for (uint j = 0 ; j < mvNClass ; j++)
                        myGradTransMat[n][j].Delete() ;
                delete [] myGradTransMat[n] ;
        }

        for (uint n = 0 ; n < mvNSample ; n++)
        {       myLambda[n].Delete() ;
                mySumLambda[n].Delete() ;
        }

        delete [] myGradCondProba ;
        delete [] myHessCondProba ;
        delete [] myGradInitProb ;
        delete [] myGradTransMat ;
        delete [] myLambda ;
        delete [] mySumLambda ;
}

/*
 * Computation of the Score vector and the (Fisher-)Information-matrix
 * @param theScore The score vector to be filled
 * @param theInformation The (Fisher-)Information-matrix to be filled
 */
void cDerivative::ComputeScoreAndInformation(cDVector& theScore, cDMatrix& theInformation)
{
        theScore = 0.0 ;
        theInformation = 0.0 ;
uint myT = 0 ;
        for (uint n = 0 ; n < mvNSample ; n++)
        {       myT += mvT[n] ;
                theScore += mvT[n]*mScore[n] ;
                theInformation +=  mvT[n]*mInformation[n] ;
        }
        theScore /= myT ;
        theInformation /= myT ;
}

void cDerivative::ComputeCov(cHmm& theHmm, cDMatrix& theCov) 
{
uint myNParam = theHmm.GetNParam() ;
//uint* myPlace = new uint[myNParam] ;
cDVector myScore(mvNFreeParam) ;
cDMatrix myInformation(mvNFreeParam, mvNFreeParam) ;
        ComputeScoreAndInformation(myScore, myInformation) ;

        theCov = Inv(myInformation) ;

uint myNFreeClass = mvNClass - 1 ;
uint /*myFreeIndex, myIndexCour, */mySizeCour ;
cDVector myU(mvNFreeParam) ;    

// initProb
        myU = 0.0 ;
        for (uint i = 0 ; i < myNFreeClass ; i++)
                myU[i] = -1.0 ;
        theCov = AddOneVariable(theCov, myU) ;
        mySizeCour = mvNFreeParam + 1 ;
// transMat
uint myBeg = 0 ;        
        for (uint n = 0 ; n < mvNClass ; n++)
        {
                myU.ReAlloc(mySizeCour, 0.0) ;
                myBeg += myNFreeClass ;
                for (uint i = myBeg ; i < myBeg+myNFreeClass ; i++)
                        myU[i] = -1.0 ;
                theCov = AddOneVariable(theCov, myU) ;
                mySizeCour++ ;
        }
// Distribution
        
        theHmm.mDistrParam->ComputeCov(theCov) ;

        // Sorting 
//uint myNumParam = theHmm.GetNParam() ;
cDVector myParamNum(mvNFreeParam) ;
cDVector myNFreeClassVect(myNFreeClass) ;
uint myIndCour ;
uint myNextInd = mvNFreeParam ;
cDVector myResNum ;

        for (uint j = 0 ; j < mvNFreeParam ; j++)
                myParamNum[j] = j ;
        // initProba
        myIndCour = 0 ;
        GetSubVector(myParamNum, myIndCour, myNFreeClass, myNFreeClassVect) ;
        myResNum = cat(myNFreeClassVect, (double)myNextInd) ;
        myNextInd++ ;
        // transMat
        for (uint i = 0 ; i < mvNClass ; i++)
        {       myIndCour += myNFreeClass ;
                GetSubVector(myParamNum, myIndCour, myNFreeClass, myNFreeClassVect) ;
                myResNum = cat(myResNum, myNFreeClassVect) ;
                myResNum = cat(myResNum, (double)myNextInd) ;
                myNextInd++ ;
        }
        // distribution
        myIndCour += myNFreeClass ;
cDVector myParamNumDistr ;
        GetSubVector(myParamNum, myIndCour, mvNFreeParam-myIndCour, myParamNumDistr) ;
cDVector myParamNumDistrAll = theHmm.mDistrParam->GetDistrNumParam(myParamNumDistr, myNextInd) ;
        myResNum = cat(myResNum, myParamNumDistrAll) ;

cDMatrix myCov = theCov ;
        for (uint i = 0 ; i < myNParam ; i++)
                for (uint j = 0 ; j < myNParam ; j++)
                        theCov[i][j] = myCov[(int)myResNum[i]][(int)myResNum[j]] ;
}
