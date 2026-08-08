/* ------------------------------------------------------------------------ */
/* Copyright (c) 2017 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence licensee.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP_Baseband Library                                               */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2017 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
  NatureDSP_Baseband library. Direct Matrix Inversion
    Direct Matrix Inversion
    C code optimized for BBE32
  IntegrIT, 2006-2017
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* NatureDSP_Baseband library API. */
#include "NatureDSP_Baseband_matinv.h"
/* Common utility declarations. */
#include "cmatinvn_common.h"

/*-------------------------------------------------------------------------
Direct Matrix Inversion For Complex Matrices

Description: perform in-place inversion of 2x2 and 4x4 complex matrices. 
2x2 matrices are inverted by Cramer's rule. For 3x3, 4x4 matrices we employ 
the blockwise inversion algorithm encompassed with a suboptimal row/column
permutation that gains better conditioning of the block structure.

Data format and order options:
  Suffix   Data Order                 Data Format   
    n        Block     16-bit signed fixed-point
    nf       Block     IEEE-754 Std single precision floating-point
    s        Stream    16-bit signed fixed-point
    sf       Stream    IEEE-754 Std single precision floating-point

Notes:
1. In general, accuracy of a matrix inversion algorithm implementation is a
   function of input matrix condition number. Thus it is user's responsibility
   to qualify the reliability of numeric results. Refer to NatureDSP Baseband 
   Library Reference for details. 
2. For blockwise inversion of 4x4 fixed-point matrices, it is reasonable to
   limit the dynamic range of input data by 11..13 significant bits. This
   measure reduces the possibility of an overflow at internal computations.

Parameters:
Temporary:
pScr      Scratch memory area. To determine the scratch area size required by
          a function <fun>, use the respective helper function 
          <fun>_getScratchSize()
Input:
L         Number of matrices
qA        Number of fractional bits for fixed-point input/output data
Input/Output:
A[L][SA]  Sequence of L NxN complex input/result matrices. SA is the number
          of data elements occupied by a single NxN matrix in a block 
          (stream) ordered sequence, see function specifications.

Restrictions:
pScr,A    Aligned on 32-byte boundary

Specification of a particular function may impose additional restrictions.
-------------------------------------------------------------------------*/

#if !(HAVE_MATINV2X2)
DISCARD_FUN(void,cmatinv2x2n,(void * restrict pScr, complex_fract16 * restrict A, int qA, int L))
DISCARD_FUN(size_t,cmatinv2x2n_getScratchSize,(int L))
#else

#define STEP  (BBE_SIMD_WIDTH/2)

void cmatinv2x2n(void * restrict pScr, complex_fract16 * restrict A, int qA, int L)
{
    int l;
    xb_vecNx16 a0, a1, a2, a3, t;
    xb_vecNx16 det, det2, recip_det,mrecip_det;
    
    xb_vecNx40 DET, DET2, Z, _1q27; 
    vsaN vsa_sh, ve_det, _16;    

    const xb_vecNx16 * restrict A_rd =   (const xb_vecNx16*)A;
          xb_vecNx16 * restrict A_wr =   (      xb_vecNx16*)A;
          xb_vecNx16 * restrict pDet =   (      xb_vecNx16*)pScr;

    NASSERT_ALIGN32(A);
    NASSERT_ALIGN32(pScr);

    if (L <= 0) return;

    A_rd = (const xb_vecNx16*)A;
    A_wr = (      xb_vecNx16*)A;
    _16 = BBE_MOVVSA32(16);

    for (l=0; l<L/STEP; l++)
    {
        // Load SIMD_WIDTH/2 2x2  matrices
        BBE_LVNX16_IP( a0, A_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a1, A_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a2, A_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a3, A_rd, 2*BBE_SIMD_WIDTH );

        // Convert data to streaming format 
        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2 );

        // Overwrite original input data with streaming data to avoid conversion
        // in the second loop.
        BBE_SVNX16_XP( a0, A_wr, 2*BBE_SIMD_WIDTH );
        BBE_SVNX16_XP( a1, A_wr, 2*BBE_SIMD_WIDTH );
        BBE_SVNX16_XP( a2, A_wr, 2*BBE_SIMD_WIDTH );
        BBE_SVNX16_XP( a3, A_wr, 2*BBE_SIMD_WIDTH );

        // compute SIMD_WIDTH/2 complex determinants
        DET=BBE_MULNX16C(a0,a3);// Q(2a)
        BBE_MULSNX16C(DET,a1,a2);//sub Q(2a)
        //sh=52-2*qA-e_det;
        ve_det=BBE_NSANX40C(DET);
        vsa_sh =BBE_SUBSAVSN(52-2*qA,ve_det);
        ve_det =BBE_SUBSAVSN(24,ve_det);
            
        DET=BBE_RNDADJNX40(DET, ve_det); 
        det=BBE_PACKVNX40(DET,ve_det);// ->Q(2a-e_det)
        DET2= BBE_MULRNX16J (det, det, _16);// -> Q(2a-e_det) - 16
        det2=BBE_PACKVNX40(DET2,_16);
        /* divide half of vector (even elements only) */
        _1q27=BBE_MOVWA32(0x8000000); 
        BBE_DIVNX32S_5STEP0_0(_1q27,det2);   
        BBE_DIVNX16S_4STEP_0(det2);         
        BBE_DIVNX16S_4STEP_0(det2);         
        recip_det=BBE_DIVNX16S_3STEPN_0(det2); 
        recip_det = BBE_SHFLNX16I( recip_det, BBE_SHFLI_DUPLICATE_1_EVEN );
        recip_det = BBE_MULNX16PACKQ(det, recip_det); 

        BBE_SVNX16_IP(recip_det , pDet, 2*BBE_SIMD_WIDTH);
        t=BBE_MOVVVS(vsa_sh);
        BBE_SVNX16_IP(t , pDet, 2*BBE_SIMD_WIDTH);
    }

    pDet = (xb_vecNx16*)pScr;
    A_rd = (xb_vecNx16*)A;
    A_wr = (xb_vecNx16*)A;

    for (l=0; l<L/STEP; l++)
    {
        // Load SIMD_WIDTH/2 2x2 matrices in streaming format
        BBE_LVNX16_IP( a0, A_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a1, A_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a2, A_rd, 2*BBE_SIMD_WIDTH );
        BBE_LVNX16_IP( a3, A_rd, 2*BBE_SIMD_WIDTH );

        BBE_LVNX16_IP(recip_det , pDet, 2*BBE_SIMD_WIDTH);
        BBE_LVNX16_IP(t         , pDet, 2*BBE_SIMD_WIDTH); vsa_sh=BBE_MOVVSV(t,0);

        mrecip_det= BBE_NEGSNX16(recip_det);
        Z = BBE_MULRNX16J(a0, recip_det, vsa_sh ); a0 = BBE_PACKVNX40( Z, vsa_sh );
        Z = BBE_MULRNX16J(a1,mrecip_det, vsa_sh ); a1 = BBE_PACKVNX40( Z, vsa_sh );    
        Z = BBE_MULRNX16J(a2,mrecip_det, vsa_sh ); a2 = BBE_PACKVNX40( Z, vsa_sh );        
        Z = BBE_MULRNX16J(a3, recip_det, vsa_sh ); a3 = BBE_PACKVNX40( Z, vsa_sh );

        //converts matrices representing in the streaming order to block (packed) order
        BBE_DSELNX16I( a0, a1, a0, a1, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a2, a3, a2, a3, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a0, a2, a0, a2, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a1, a3, a1, a3, BBE_DSELI_INTERLEAVE_2 );

        //Save SIMD_WIDTH/2 inverted 2x2 matrices.
        BBE_SVNX16_IP(a3 , A_wr, 2*BBE_SIMD_WIDTH);
        BBE_SVNX16_IP(a1 , A_wr, 2*BBE_SIMD_WIDTH); 
        BBE_SVNX16_IP(a2 , A_wr, 2*BBE_SIMD_WIDTH); 
        BBE_SVNX16_IP(a0 , A_wr, 2*BBE_SIMD_WIDTH); 
    }

    L-=(L/STEP)*STEP;

    // last 4, 8 or 12 matrices
    if (L>0)
    {
        valign rd_va, wr_va;
        int rd_nb, wr_nb;

        _1q27=BBE_MOVWA32(0x8000000); 

        // last 1x, 2x or 3xSIMD_WIDTH/8 matrices
        rd_nb = L*4*4;
        rd_va = BBE_LAVNX16_PP( A_rd );
        BBE_LAVNX16_XP( a0, rd_va, A_rd, rd_nb ); rd_nb -= 2*BBE_SIMD_WIDTH;
        BBE_LAVNX16_XP( a1, rd_va, A_rd, rd_nb ); rd_nb -= 2*BBE_SIMD_WIDTH;
        BBE_LAVNX16_XP( a2, rd_va, A_rd, rd_nb );
        a3=0;

        // Convert data to streaming format 
        BBE_DSELNX16I( a1, a0, a1, a0, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a2, a3, a2, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a3, a1, a3, a1, BBE_DSELI_DEINTERLEAVE_2 );
        BBE_DSELNX16I( a2, a0, a2, a0, BBE_DSELI_DEINTERLEAVE_2 );
        // compute SIMD_WIDTH/2 complex determinants
        DET=BBE_MULNX16C(a0,a3);// Q(2a)
        BBE_MULSNX16C(DET,a1,a2);//sub Q(2a)
        //sh=52-2*qA-e_det;
        ve_det=BBE_NSANX40C(DET);
        vsa_sh =BBE_SUBSAVSN(52-2*qA,ve_det);
        ve_det =BBE_SUBSAVSN(24,ve_det);
        
        DET=BBE_RNDADJNX40(DET, ve_det); 
        det=BBE_PACKVNX40(DET,ve_det);// ->Q(2a-e_det)
        DET2= BBE_MULRNX16J (det, det, _16);// -> Q(2a-e_det) - 16
        det2=BBE_PACKVNX40(DET2,_16);
        /* divide half of vector (even elements only) */
        BBE_DIVNX32S_5STEP0_0(_1q27,det2);   
        BBE_DIVNX16S_4STEP_0(det2);         
        BBE_DIVNX16S_4STEP_0(det2);         
        recip_det=BBE_DIVNX16S_3STEPN_0(det2); 
        recip_det = BBE_SHFLNX16I( recip_det, BBE_SHFLI_DUPLICATE_1_EVEN );
        recip_det = BBE_MULNX16PACKQ(det, recip_det); 
        mrecip_det= BBE_NEGSNX16(recip_det);
        Z = BBE_MULRNX16J(a0, recip_det, vsa_sh ); a0 = BBE_PACKVNX40( Z, vsa_sh );
        Z = BBE_MULRNX16J(a1,mrecip_det, vsa_sh ); a1 = BBE_PACKVNX40( Z, vsa_sh );    
        Z = BBE_MULRNX16J(a2,mrecip_det, vsa_sh ); a2 = BBE_PACKVNX40( Z, vsa_sh );        
        Z = BBE_MULRNX16J(a3, recip_det, vsa_sh ); a3 = BBE_PACKVNX40( Z, vsa_sh );

        //converts matrices representing in the streaming order to block (packed) order
        BBE_DSELNX16I( a0, a1, a0, a1, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a2, a3, a2, a3, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a0, a2, a0, a2, BBE_DSELI_INTERLEAVE_2 );
        BBE_DSELNX16I( a1, a3, a1, a3, BBE_DSELI_INTERLEAVE_2 );

        //Save 1x, 2x or 3xSIMD_WIDTH/8 inverted 2x2 matrices.
        wr_nb = L*4*4;
        wr_va = BBE_ZALIGN();
        BBE_SAVNX16_XP( a3, wr_va, A_wr, wr_nb ); wr_nb -= 2*BBE_SIMD_WIDTH;
        BBE_SAVNX16_XP( a1, wr_va, A_wr, wr_nb ); wr_nb -= 2*BBE_SIMD_WIDTH;
        BBE_SAVNX16_XP( a2, wr_va, A_wr, wr_nb );
        BBE_SAVNX16POS_FP( wr_va, A_wr );
    }

} /* cmatinv2x2n() */

/* Return the scratch area size, in bytes. */
size_t cmatinv2x2n_getScratchSize ( int L )
{
    return ( L>0 ? 2*(L>>(LOG2_BBE_SIMD_WIDTH-1))*2*BBE_SIMD_WIDTH : 0 );
}
#endif
