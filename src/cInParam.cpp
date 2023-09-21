/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: cInParam.cpp 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***      
 *** Class for the parameters of the HMM                                                   
 **************************************************************/

#include "StdAfxRHmm.h"

/*
 * Constructor for InParam, containing the parameters of the HMM
 * @param theNSample Number of samples
 * @param theNSample Dimension of observations
 * @param theY Random vector Y
 * @param distrDefinitionEnum Type of distribution
 * @param theNClass Number of classes
 * @param theNMixt Number of mixtures (for mixture distributions)
 * @param theNProba Number of probabilities
 */
cInParam::cInParam(uint theNSample, uint theDimObs, cDVector* theY, distrDefinitionEnum theDistrType, uint theNClass, uint theNMixt, uint theNProba)
{       MESS_CREAT("cInParam")
        mDistrType = theDistrType ;
        mNClass = theNClass ;
        mNMixt = theNMixt ;
        mNProba = theNProba ;
        mNSample = theNSample ;
        mDimObs = theDimObs ;
        if (mNSample > 0)
        {       mY = new cDVector[mNSample] ;
                for (uint i = 0 ; i < mNSample ; i++)
                        mY[i] = theY[i] ;
        }
        else
                mY = (cDVector *)NULL ;
}

/*
 * Destructor for InParam
 */
cInParam::~cInParam()
{       MESS_DESTR("cInParam")
        if (mNSample != 0)
        {       for (uint i = 0 ; i < mNSample ; i++)
                        mY[i].Delete() ;
                delete [] mY ;
                mNSample = 0 ;
        }
}

cInParam &cInParam::operator =(const cInParam &theSrc)
{       
        mDistrType = theSrc.mDistrType ;                
        mNClass = theSrc.mNClass ;
        if (mNSample > 0)
        {       for (uint i = 0 ; i < mNSample ; i++)
                        mY[i].Delete() ;
                delete mY ;
        }
        mNSample = theSrc.mNSample ;
        mY = new cDVector[mNSample] ;
        
        mDimObs = theSrc.mDimObs ;
        mNProba = theSrc.mNProba ;
        mNMixt = theSrc.mNMixt ;
        
        for (uint i = 0 ; i < mNSample ; i++)
                mY[i] = theSrc.mY[i] ;
        
        return(*this) ;
}

void cInParam::Print(void)
{
        Rprintf("NbSample = %d\n", mNSample) ;
        for (uint n = 0 ; n < mNSample ; n++)
                Rprintf("mT[%d]=%d\n", n, mY[n].mSize) ;
}
