#ifndef _COEFF_H_
#define _COEFF_H_

#include <AMReX_REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_SET_CC_COEF        set_cc_coef_2d
#define FORT_COEF_TO_EDGES      coef_to_edges_2d
#elif (BL_SPACEDIM == 3)
#define FORT_SET_CC_COEF        set_cc_coef_3d
#define FORT_COEF_TO_EDGES      coef_to_edges_3d
#endif

#else

#if (BL_SPACEDIM == 2)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_SET_CC_COEF        SET_CC_COEF_2d
#define FORT_COEF_TO_EDGES      COEF_TO_EDGES_2d
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_SET_CC_COEF        set_cc_coef_2d
#define FORT_COEF_TO_EDGES      coef_to_edges_2d
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_SET_CC_COEF        set_cc_coef_2d_
#define FORT_COEF_TO_EDGES      coef_to_edges_2d_
#endif
#endif

#if (BL_SPACEDIM == 3)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_SET_CC_COEF        SET_CC_COEF_3d
#define FORT_COEF_TO_EDGES      COEF_TO_EDGES_3d
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_SET_CC_COEF        set_cc_coef_3d
#define FORT_COEF_TO_EDGES      coef_to_edges_3d
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_SET_CC_COEF        set_cc_coef_3d_
#define FORT_COEF_TO_EDGES      coef_to_edges_3d_
#endif
#endif


#include <AMReX_ArrayLim.H>

extern "C"
{
    void FORT_SET_CC_COEF (
        amrex::Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
        const int* lo, const int* hi, const amrex::Real* dx,
        const amrex::Real* prob_lo, const amrex::Real* prob_hi);

    void FORT_COEF_TO_EDGES (int * dir, 
        amrex::Real* edge     , ARLIM_P(edge_lo), ARLIM_P(edge_hi),
        amrex::Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
        const int* lo, const int* hi);

}
#endif

#endif /*_COEF_F_H_*/
