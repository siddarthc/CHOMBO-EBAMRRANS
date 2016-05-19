#ifndef _EBSCALARBCF_F_H_
#define _EBSCALARBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETINJECTIONPARAMS 
#define GUARDSETINJECTIONPARAMS 
// Prototype for Fortran procedure setinjectionparams ...
//
void FORTRAN_NAME( SETINJECTIONPARAMS ,setinjectionparams )(
      CHFp_CONST_INT(auinjectDir)
      ,CHFp_CONST_REAL(auscalarInflowVal)
      ,CHFp_CONST_REAL(auinjectRad)
      ,CHFp_CONST_VR(auinjectOrig)
      ,CHFp_CONST_VR(auinjectEnd) );

#define FORT_SETINJECTIONPARAMS FORTRAN_NAME( inlineSETINJECTIONPARAMS, inlineSETINJECTIONPARAMS)
#define FORTNT_SETINJECTIONPARAMS FORTRAN_NAME( SETINJECTIONPARAMS, setinjectionparams)

inline void FORTRAN_NAME(inlineSETINJECTIONPARAMS, inlineSETINJECTIONPARAMS)(
      CHFp_CONST_INT(auinjectDir)
      ,CHFp_CONST_REAL(auscalarInflowVal)
      ,CHFp_CONST_REAL(auinjectRad)
      ,CHFp_CONST_VR(auinjectOrig)
      ,CHFp_CONST_VR(auinjectEnd) )
{
 CH_TIMELEAF("FORT_SETINJECTIONPARAMS");
 FORTRAN_NAME( SETINJECTIONPARAMS ,setinjectionparams )(
      CHFt_CONST_INT(auinjectDir)
      ,CHFt_CONST_REAL(auscalarInflowVal)
      ,CHFt_CONST_REAL(auinjectRad)
      ,CHFt_CONST_VR(auinjectOrig)
      ,CHFt_CONST_VR(auinjectEnd) );
}
#endif  // GUARDSETINJECTIONPARAMS 

#ifndef GUARDSCALARINIT 
#define GUARDSCALARINIT 
// Prototype for Fortran procedure scalarinit ...
//
void FORTRAN_NAME( SCALARINIT ,scalarinit )(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SCALARINIT FORTRAN_NAME( inlineSCALARINIT, inlineSCALARINIT)
#define FORTNT_SCALARINIT FORTRAN_NAME( SCALARINIT, scalarinit)

inline void FORTRAN_NAME(inlineSCALARINIT, inlineSCALARINIT)(
      CHFp_FRA(u)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SCALARINIT");
 FORTRAN_NAME( SCALARINIT ,scalarinit )(
      CHFt_FRA(u)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSCALARINIT 

#ifndef GUARDPOINTSCALARINIT 
#define GUARDPOINTSCALARINIT 
// Prototype for Fortran procedure pointscalarinit ...
//
void FORTRAN_NAME( POINTSCALARINIT ,pointscalarinit )(
      CHFp_REAL(scalarval)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) );

#define FORT_POINTSCALARINIT FORTRAN_NAME( inlinePOINTSCALARINIT, inlinePOINTSCALARINIT)
#define FORTNT_POINTSCALARINIT FORTRAN_NAME( POINTSCALARINIT, pointscalarinit)

inline void FORTRAN_NAME(inlinePOINTSCALARINIT, inlinePOINTSCALARINIT)(
      CHFp_REAL(scalarval)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_POINTSCALARINIT");
 FORTRAN_NAME( POINTSCALARINIT ,pointscalarinit )(
      CHFt_REAL(scalarval)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDPOINTSCALARINIT 

}

#endif
