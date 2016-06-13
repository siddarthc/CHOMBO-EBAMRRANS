#ifndef _EBPATCHRANSMODELF_F_H_
#define _EBPATCHRANSMODELF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPREDADVECTRANS 
#define GUARDPREDADVECTRANS 
// Prototype for Fortran procedure predadvectrans ...
//
void FORTRAN_NAME( PREDADVECTRANS ,predadvectrans )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_FRA(dq)
      ,CHFp_CONST_FRA(velcc)
      ,CHFp_FRA(qlo)
      ,CHFp_FRA(qhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nEqn) );

#define FORT_PREDADVECTRANS FORTRAN_NAME( inlinePREDADVECTRANS, inlinePREDADVECTRANS)
#define FORTNT_PREDADVECTRANS FORTRAN_NAME( PREDADVECTRANS, predadvectrans)

inline void FORTRAN_NAME(inlinePREDADVECTRANS, inlinePREDADVECTRANS)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA(q)
      ,CHFp_CONST_FRA(dq)
      ,CHFp_CONST_FRA(velcc)
      ,CHFp_FRA(qlo)
      ,CHFp_FRA(qhi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nEqn) )
{
 CH_TIMELEAF("FORT_PREDADVECTRANS");
 FORTRAN_NAME( PREDADVECTRANS ,predadvectrans )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA(q)
      ,CHFt_CONST_FRA(dq)
      ,CHFt_CONST_FRA(velcc)
      ,CHFt_FRA(qlo)
      ,CHFt_FRA(qhi)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(nEqn) );
}
#endif  // GUARDPREDADVECTRANS 

#ifndef GUARDPOINTPREDADVECTRANS 
#define GUARDPOINTPREDADVECTRANS 
// Prototype for Fortran procedure pointpredadvectrans ...
//
void FORTRAN_NAME( POINTPREDADVECTRANS ,pointpredadvectrans )(
      CHFp_VR(primit)
      ,CHFp_VR(primlo)
      ,CHFp_VR(primhi)
      ,CHFp_VR(pslope)
      ,CHFp_REALVECT(velcc)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nEqn) );

#define FORT_POINTPREDADVECTRANS FORTRAN_NAME( inlinePOINTPREDADVECTRANS, inlinePOINTPREDADVECTRANS)
#define FORTNT_POINTPREDADVECTRANS FORTRAN_NAME( POINTPREDADVECTRANS, pointpredadvectrans)

inline void FORTRAN_NAME(inlinePOINTPREDADVECTRANS, inlinePOINTPREDADVECTRANS)(
      CHFp_VR(primit)
      ,CHFp_VR(primlo)
      ,CHFp_VR(primhi)
      ,CHFp_VR(pslope)
      ,CHFp_REALVECT(velcc)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nEqn) )
{
 CH_TIMELEAF("FORT_POINTPREDADVECTRANS");
 FORTRAN_NAME( POINTPREDADVECTRANS ,pointpredadvectrans )(
      CHFt_VR(primit)
      ,CHFt_VR(primlo)
      ,CHFt_VR(primhi)
      ,CHFt_VR(pslope)
      ,CHFt_REALVECT(velcc)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(nEqn) );
}
#endif  // GUARDPOINTPREDADVECTRANS 

}

#endif
