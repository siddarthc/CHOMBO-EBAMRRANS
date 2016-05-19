#ifndef _EBGRADDIVFILTERF_F_H_
#define _EBGRADDIVFILTERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEBGDFFACEDIVINCR 
#define GUARDEBGDFFACEDIVINCR 
// Prototype for Fortran procedure ebgdffacedivincr ...
//
void FORTRAN_NAME( EBGDFFACEDIVINCR ,ebgdffacedivincr )(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA(gradvel)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(divdir) );

#define FORT_EBGDFFACEDIVINCR FORTRAN_NAME( inlineEBGDFFACEDIVINCR, inlineEBGDFFACEDIVINCR)
#define FORTNT_EBGDFFACEDIVINCR FORTRAN_NAME( EBGDFFACEDIVINCR, ebgdffacedivincr)

inline void FORTRAN_NAME(inlineEBGDFFACEDIVINCR, inlineEBGDFFACEDIVINCR)(
      CHFp_FRA1(divvel)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA(gradvel)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(divdir) )
{
 CH_TIMELEAF("FORT_EBGDFFACEDIVINCR");
 FORTRAN_NAME( EBGDFFACEDIVINCR ,ebgdffacedivincr )(
      CHFt_FRA1(divvel)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_FRA(gradvel)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(divdir) );
}
#endif  // GUARDEBGDFFACEDIVINCR 

#ifndef GUARDEBGDFGRADVEL 
#define GUARDEBGDFGRADVEL 
// Prototype for Fortran procedure ebgdfgradvel ...
//
void FORTRAN_NAME( EBGDFGRADVEL ,ebgdfgradvel )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_CONST_INT(loworderoneside) );

#define FORT_EBGDFGRADVEL FORTRAN_NAME( inlineEBGDFGRADVEL, inlineEBGDFGRADVEL)
#define FORTNT_EBGDFGRADVEL FORTRAN_NAME( EBGDFGRADVEL, ebgdfgradvel)

inline void FORTRAN_NAME(inlineEBGDFGRADVEL, inlineEBGDFGRADVEL)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(vel)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(divdir)
      ,CHFp_CONST_INT(loworderoneside) )
{
 CH_TIMELEAF("FORT_EBGDFGRADVEL");
 FORTRAN_NAME( EBGDFGRADVEL ,ebgdfgradvel )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(vel)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(divdir)
      ,CHFt_CONST_INT(loworderoneside) );
}
#endif  // GUARDEBGDFGRADVEL 

#ifndef GUARDEBGDFCELLGRAD 
#define GUARDEBGDFCELLGRAD 
// Prototype for Fortran procedure ebgdfcellgrad ...
//
void FORTRAN_NAME( EBGDFCELLGRAD ,ebgdfcellgrad )(
      CHFp_FRA1(graddiv)
      ,CHFp_CONST_FRA1(div)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(imultbylambda) );

#define FORT_EBGDFCELLGRAD FORTRAN_NAME( inlineEBGDFCELLGRAD, inlineEBGDFCELLGRAD)
#define FORTNT_EBGDFCELLGRAD FORTRAN_NAME( EBGDFCELLGRAD, ebgdfcellgrad)

inline void FORTRAN_NAME(inlineEBGDFCELLGRAD, inlineEBGDFCELLGRAD)(
      CHFp_FRA1(graddiv)
      ,CHFp_CONST_FRA1(div)
      ,CHFp_CONST_FRA1(lambda)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(imultbylambda) )
{
 CH_TIMELEAF("FORT_EBGDFCELLGRAD");
 FORTRAN_NAME( EBGDFCELLGRAD ,ebgdfcellgrad )(
      CHFt_FRA1(graddiv)
      ,CHFt_CONST_FRA1(div)
      ,CHFt_CONST_FRA1(lambda)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(imultbylambda) );
}
#endif  // GUARDEBGDFCELLGRAD 

}

#endif
