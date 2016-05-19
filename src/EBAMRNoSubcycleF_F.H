#ifndef _EBAMRNOSUBCYCLEF_F_H_
#define _EBAMRNOSUBCYCLEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMULTVELINTOFLUX 
#define GUARDMULTVELINTOFLUX 
// Prototype for Fortran procedure MULTVELINTOFLUX ...
//
void FORTRAN_NAME( MULTVELINTOFLUX ,multvelintoflux )(
      CHFp_FRA(flux)
      ,CHFp_FRA(vel)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_MULTVELINTOFLUX FORTRAN_NAME( inlineMULTVELINTOFLUX, inlineMULTVELINTOFLUX)
#define FORTNT_MULTVELINTOFLUX FORTRAN_NAME( MULTVELINTOFLUX, multvelintoflux)

inline void FORTRAN_NAME(inlineMULTVELINTOFLUX, inlineMULTVELINTOFLUX)(
      CHFp_FRA(flux)
      ,CHFp_FRA(vel)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_MULTVELINTOFLUX");
 FORTRAN_NAME( MULTVELINTOFLUX ,multvelintoflux )(
      CHFt_FRA(flux)
      ,CHFt_FRA(vel)
      ,CHFt_BOX(gridBox)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDMULTVELINTOFLUX 

#ifndef GUARDCOMPUTEVORT 
#define GUARDCOMPUTEVORT 
// Prototype for Fortran procedure COMPUTEVORT ...
//
void FORTRAN_NAME( COMPUTEVORT ,computevort )(
      CHFp_FRA1(vort)
      ,CHFp_FRA(vel)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir) );

#define FORT_COMPUTEVORT FORTRAN_NAME( inlineCOMPUTEVORT, inlineCOMPUTEVORT)
#define FORTNT_COMPUTEVORT FORTRAN_NAME( COMPUTEVORT, computevort)

inline void FORTRAN_NAME(inlineCOMPUTEVORT, inlineCOMPUTEVORT)(
      CHFp_FRA1(vort)
      ,CHFp_FRA(vel)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_COMPUTEVORT");
 FORTRAN_NAME( COMPUTEVORT ,computevort )(
      CHFt_FRA1(vort)
      ,CHFt_FRA(vel)
      ,CHFt_BOX(gridBox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDCOMPUTEVORT 

#ifndef GUARDKINETICENERGY 
#define GUARDKINETICENERGY 
// Prototype for Fortran procedure KINETICENERGY ...
//
void FORTRAN_NAME( KINETICENERGY ,kineticenergy )(
      CHFp_FRA1(energy)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_BOX(gridBox) );

#define FORT_KINETICENERGY FORTRAN_NAME( inlineKINETICENERGY, inlineKINETICENERGY)
#define FORTNT_KINETICENERGY FORTRAN_NAME( KINETICENERGY, kineticenergy)

inline void FORTRAN_NAME(inlineKINETICENERGY, inlineKINETICENERGY)(
      CHFp_FRA1(energy)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_BOX(gridBox) )
{
 CH_TIMELEAF("FORT_KINETICENERGY");
 FORTRAN_NAME( KINETICENERGY ,kineticenergy )(
      CHFt_FRA1(energy)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_BOX(gridBox) );
}
#endif  // GUARDKINETICENERGY 

#ifndef GUARDMAGVECT 
#define GUARDMAGVECT 
// Prototype for Fortran procedure MAGVECT ...
//
void FORTRAN_NAME( MAGVECT ,magvect )(
      CHFp_FRA1(magVector)
      ,CHFp_FRA(vect)
      ,CHFp_BOX(gridBox) );

#define FORT_MAGVECT FORTRAN_NAME( inlineMAGVECT, inlineMAGVECT)
#define FORTNT_MAGVECT FORTRAN_NAME( MAGVECT, magvect)

inline void FORTRAN_NAME(inlineMAGVECT, inlineMAGVECT)(
      CHFp_FRA1(magVector)
      ,CHFp_FRA(vect)
      ,CHFp_BOX(gridBox) )
{
 CH_TIMELEAF("FORT_MAGVECT");
 FORTRAN_NAME( MAGVECT ,magvect )(
      CHFt_FRA1(magVector)
      ,CHFt_FRA(vect)
      ,CHFt_BOX(gridBox) );
}
#endif  // GUARDMAGVECT 

#ifndef GUARDCENTERED_DERIV 
#define GUARDCENTERED_DERIV 
// Prototype for Fortran procedure CENTERED_DERIV ...
//
void FORTRAN_NAME( CENTERED_DERIV ,centered_deriv )(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_CENTERED_DERIV FORTRAN_NAME( inlineCENTERED_DERIV, inlineCENTERED_DERIV)
#define FORTNT_CENTERED_DERIV FORTRAN_NAME( CENTERED_DERIV, centered_deriv)

inline void FORTRAN_NAME(inlineCENTERED_DERIV, inlineCENTERED_DERIV)(
      CHFp_FRA1(deriv)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(grid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_CENTERED_DERIV");
 FORTRAN_NAME( CENTERED_DERIV ,centered_deriv )(
      CHFt_FRA1(deriv)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(grid)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDCENTERED_DERIV 

}

#endif
