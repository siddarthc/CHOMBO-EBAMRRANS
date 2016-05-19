#ifndef _INFLOWOUTFLOWF_F_H_
#define _INFLOWOUTFLOWF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETJET2PARAMS 
#define GUARDSETJET2PARAMS 
// Prototype for Fortran procedure setjet2params ...
//
void FORTRAN_NAME( SETJET2PARAMS ,setjet2params )(
      CHFp_CONST_INT(aujet2flowDir)
      ,CHFp_CONST_REAL(aujet2maxVel)
      ,CHFp_CONST_REAL(aujet2radius)
      ,CHFp_CONST_REAL(aujet2viscosity)
      ,CHFp_CONST_VR(aujet2tubeorig)
      ,CHFp_CONST_VR(aujet2tubeend) );

#define FORT_SETJET2PARAMS FORTRAN_NAME( inlineSETJET2PARAMS, inlineSETJET2PARAMS)
#define FORTNT_SETJET2PARAMS FORTRAN_NAME( SETJET2PARAMS, setjet2params)

inline void FORTRAN_NAME(inlineSETJET2PARAMS, inlineSETJET2PARAMS)(
      CHFp_CONST_INT(aujet2flowDir)
      ,CHFp_CONST_REAL(aujet2maxVel)
      ,CHFp_CONST_REAL(aujet2radius)
      ,CHFp_CONST_REAL(aujet2viscosity)
      ,CHFp_CONST_VR(aujet2tubeorig)
      ,CHFp_CONST_VR(aujet2tubeend) )
{
 CH_TIMELEAF("FORT_SETJET2PARAMS");
 FORTRAN_NAME( SETJET2PARAMS ,setjet2params )(
      CHFt_CONST_INT(aujet2flowDir)
      ,CHFt_CONST_REAL(aujet2maxVel)
      ,CHFt_CONST_REAL(aujet2radius)
      ,CHFt_CONST_REAL(aujet2viscosity)
      ,CHFt_CONST_VR(aujet2tubeorig)
      ,CHFt_CONST_VR(aujet2tubeend) );
}
#endif  // GUARDSETJET2PARAMS 

#ifndef GUARDJET2VELINIT 
#define GUARDJET2VELINIT 
// Prototype for Fortran procedure jet2velinit ...
//
void FORTRAN_NAME( JET2VELINIT ,jet2velinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_JET2VELINIT FORTRAN_NAME( inlineJET2VELINIT, inlineJET2VELINIT)
#define FORTNT_JET2VELINIT FORTRAN_NAME( JET2VELINIT, jet2velinit)

inline void FORTRAN_NAME(inlineJET2VELINIT, inlineJET2VELINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_JET2VELINIT");
 FORTRAN_NAME( JET2VELINIT ,jet2velinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDJET2VELINIT 

#ifndef GUARDPOINTJET2VELINIT 
#define GUARDPOINTJET2VELINIT 
// Prototype for Fortran procedure pointjet2velinit ...
//
void FORTRAN_NAME( POINTJET2VELINIT ,pointjet2velinit )(
      CHFp_REALVECT(velocity)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTJET2VELINIT FORTRAN_NAME( inlinePOINTJET2VELINIT, inlinePOINTJET2VELINIT)
#define FORTNT_POINTJET2VELINIT FORTRAN_NAME( POINTJET2VELINIT, pointjet2velinit)

inline void FORTRAN_NAME(inlinePOINTJET2VELINIT, inlinePOINTJET2VELINIT)(
      CHFp_REALVECT(velocity)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTJET2VELINIT");
 FORTRAN_NAME( POINTJET2VELINIT ,pointjet2velinit )(
      CHFt_REALVECT(velocity)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTJET2VELINIT 

#ifndef GUARDJET2GRADPINIT 
#define GUARDJET2GRADPINIT 
// Prototype for Fortran procedure jet2gradPinit ...
//
void FORTRAN_NAME( JET2GRADPINIT ,jet2gradpinit )(
      CHFp_FRA(gP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_JET2GRADPINIT FORTRAN_NAME( inlineJET2GRADPINIT, inlineJET2GRADPINIT)
#define FORTNT_JET2GRADPINIT FORTRAN_NAME( JET2GRADPINIT, jet2gradpinit)

inline void FORTRAN_NAME(inlineJET2GRADPINIT, inlineJET2GRADPINIT)(
      CHFp_FRA(gP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_JET2GRADPINIT");
 FORTRAN_NAME( JET2GRADPINIT ,jet2gradpinit )(
      CHFt_FRA(gP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDJET2GRADPINIT 

#ifndef GUARDPOINTJET2GRADPINIT 
#define GUARDPOINTJET2GRADPINIT 
// Prototype for Fortran procedure pointjet2gradPinit ...
//
void FORTRAN_NAME( POINTJET2GRADPINIT ,pointjet2gradpinit )(
      CHFp_REALVECT(gradP)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTJET2GRADPINIT FORTRAN_NAME( inlinePOINTJET2GRADPINIT, inlinePOINTJET2GRADPINIT)
#define FORTNT_POINTJET2GRADPINIT FORTRAN_NAME( POINTJET2GRADPINIT, pointjet2gradpinit)

inline void FORTRAN_NAME(inlinePOINTJET2GRADPINIT, inlinePOINTJET2GRADPINIT)(
      CHFp_REALVECT(gradP)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTJET2GRADPINIT");
 FORTRAN_NAME( POINTJET2GRADPINIT ,pointjet2gradpinit )(
      CHFt_REALVECT(gradP)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTJET2GRADPINIT 

}

#endif
