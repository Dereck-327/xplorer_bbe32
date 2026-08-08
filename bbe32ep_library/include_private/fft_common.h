#ifndef _FFT_COMMON_H_
#define _FFT_COMMON_H_
#include "NatureDSP_types.h"
#include "common.h"
#if HAVE_FFT
#include "fft_macro.h"
#include "fft_tw.h"

#ifdef IS_INV_FFT
#define R2_DFT5xIN_5_unroll_all iR2_DFT5xIN_5_unroll_all
#define R3_DFT5xI3xIv iR3_DFT5xI3xIv
#define R2_tDFT4xI4_U2 iR2_tDFT4xI4_U2
#define R1_tDFT4_L64_16 iR1_tDFT4_L64_16
#define R3_I3xDFT3xIv iR3_I3xDFT3xIv
#define R2_DFT3xIN_3_compress65 iR2_DFT3xIN_3_compress65
#define R2_DFT5xI4_expand56 iR2_DFT5xI4_expand56
#define R2_DFT9xIN_9 iR2_DFT9xIN_9
#define R2_DFT6xI4_N8 iR2_DFT6xI4_N8
#define R1_DFT4_L64_16_mr iR1_DFT4_L64_16_mr
#define R2_DFT4xIN_4_rm iR2_DFT4xIN_4_rm
#define R2_DFT4xI4_rr iR2_DFT4xI4_rr
#define R1_DFT4_L64_16_norm iR1_DFT4_L64_16_norm
#define R2_DFT2xI4_N8 iR2_DFT2xI4_N8
#define R2_DFT5xIN_5_compress iR2_DFT5xIN_5_compress
#define R2_DFT3xIN_3_compress iR2_DFT3xIN_3_compress
#define R2_DFT3xI4_expand iR2_DFT3xI4_expand
#define R1_DFT4_L64_16_N4 iR1_DFT4_L64_16_N4
#define R2_DFT3xI4_expand iR2_DFT3xI4_expand
#define R3_DFT3xIv        iR3_DFT3xIv  
#define R1_DFT3xIN_compress iR1_DFT3xIN_compress
#define R3_DFT2xIv iR3_DFT2xIv
#define R1_DFT5xIN_compress iR1_DFT5xIN_compress
#define R3_DFT5xIv iR3_DFT5xIv
#define R1_DFT2xIN_compress iR1_DFT2xIN_compress
#define R1_DFT4_L64_16 iR1_DFT4_L64_16
#define R2_DFT4xI4 iR2_DFT4xI4
#define R2_DFT4xIN_4 iR2_DFT4xIN_4
#define R2_DFT8xIN_8 iR2_DFT8xIN_8
#define R3_DFT4xIv iR3_DFT4xIv
#define R2_DFT5xIN_5 iR2_DFT5xIN_5
#define R2_DFT2xI4 iR2_DFT2xI4
#define R1_DFT4_L64_16_N16 iR1_DFT4_L64_16_N16
#define R2_DFT4xI4_N16 iR2_DFT4xI4_N16
#define R2_DFT4xI4_wSEL iR2_DFT4xI4_wSEL
#define R3_DFT4xIv_2U iR3_DFT4xIv_2U
#define R1_DFT4_L64_16_fs iR1_DFT4_L64_16_fs
#define R2_DFT4xI4_U2 iR2_DFT4xI4_U2
#define R2_DFT6xIN_6 iR2_DFT6xIN_6
#define R2_DFT6xI4 iR2_DFT6xI4
#define R2_DFT3xIN_3 iR2_DFT3xIN_3
#define R2_DFT4xI4_N4 iR2_DFT4xI4_N4
#endif 





inline_ const int16_t* R3_DFT2xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       )ATTRIBUTE_ALWAYS_INLINE; 

inline_ void  R1_DFT2xIN_compress( 
                                       int16_t *x, 
                                       int16_t *y, 
                                       const int N )ATTRIBUTE_ALWAYS_INLINE;

//inline_ int R1_DFT4_L64_16_N4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE; 

inline_ int R1_DFT4_L64_16_N16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE; 
inline_ int R1_DFT4_L64_16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE;
inline_ int R2_DFT4xI4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE; 
inline_ int R2_DFT4xI4_N16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE;
inline_ int R2_DFT4xIN_4(int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE;
inline_ int R2_DFT8xIN_8(int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE; 
inline_ int R2_DFT2xI4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)ATTRIBUTE_ALWAYS_INLINE; 



inline_ const int16_t* R3_DFT2xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       )
{

    const int N1 = 2;
    const int stride = 2*N/N1; /*in int16_t */
    int i, j; 
//    unsigned int tmp;

    int num_bfls = N/N1/v;

    xb_vecN_2xcq15 y0, y1, x0, x1; 
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
 //   xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 tw1; 

    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    xb_vecNx16 zero = BBE_ZERONX16();  

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 

    BBE_MOVSAV(zero);    
    BBE_MOVSBV(zero);    

    for(i=0; i<num_bfls; i++)
    {
        xb_vecNx16 tmp0; 
        BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
        tmp0 = BBE_REPNX16C(tmp0, 0); 
        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 

        for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
        {

            BBE_LVN_2XCQ15_IP(x0, px0, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(x1, px1, sizeof(*px1)); 
            /* y2 = x0 + x1 + A + B; A==0 && B==0 */		
            y0 = BBE_FFTADD4SABN_2XCQ15(x0, x1, 0, 0);				
            /* y2 = x0 - x1 + A - B; A==0 && B==0 */				
            y1 = BBE_FFTADD4SABN_2XCQ15(x0, x1, 2, 0);				
									\
            MUL(y1, tw1);		

            BBE_SVRN_2XCQ15_X (y1, py, 1*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

        }
        py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
    }    
    return (int16_t*)tw; 
}





inline_ void  R1_DFT2xIN_compress( 
                                       int16_t *x, 
                                       int16_t *y, 
                                       const int N )
{
    const int v = 16;       // vector size
    const int N1 = 2;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_out = 2*N*3/4/N1; /*in int16_t */
    
    int i; 

    int   num_bfls = N/N1/v;
//    int   num_frac_bfls = N/N1 - v * num_bfls;
    xb_vecN_2xcq15 y0, y1, x0, x1; 

    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);


    xb_vecN_2xcq15 *py0 = (xb_vecN_2xcq15 *)(y);
    xb_vecN_2xcq15 *py1 = (xb_vecN_2xcq15 *)(y + 1*stride_out);


    xb_vecNx16 zero = BBE_ZERONX16();  
    valign a0 =  BBE_ZALIGN();
    valign a1 =  BBE_ZALIGN();

    BBE_MOVSAV(zero);    
    BBE_MOVSBV(zero); 


    ASSERT((N/N1 - v * num_bfls)==0); 

    for(i=0; i<num_bfls; i++)
    {
    								
        BBE_LVN_2XCQ15_IP(x0, px0, sizeof(*px1));
        BBE_LVN_2XCQ15_IP(x1, px1, sizeof(*px1)); 
        /* y2 = x0 + x1 + A + B; A==0 && B==0 */		
        y0 = BBE_FFTADD4SABN_2XCQ15(x0, x1, 0, 0);				
        /* y2 = x0 - x1 + A - B; A==0 && B==0 */				
        y1 = BBE_FFTADD4SABN_2XCQ15(x0, x1, 2, 0);			

        BBE_SAVRN_2XCQ15_XP(y0, a0, py0, 12*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y1, a1, py1, 12*2*sizeof(int16_t)); 

    }

    BBE_SAN_2XCQ15POS_FC(a0, py0);	
    BBE_SAN_2XCQ15POS_FC(a1, py1);
}

/*
    N must be multiple of 16
*/

inline_ int R1_DFT4_L64_16_N16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{

    int stride = 2*N/4; 
    int count = N/4/(BBE_SIMD_WIDTH/2) + (N%(4*BBE_SIMD_WIDTH/2) > 0 );                                        
    int i;        
    VT * p_tw = (VT *)(tw);  

    xb_vecNx16 * px0 = (xb_vecNx16 *)(x);
    xb_vecNx16 * px1 = (xb_vecNx16 *)(x+stride*1);
    xb_vecNx16 * px2 = (xb_vecNx16 *)(x+stride*2);
    xb_vecNx16 * px3 = (xb_vecNx16 *)(x+stride*3);

    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;
    valign v1, v3; 
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  


    int scaling = bexp;
    vsaN shift; 
  
    bexp = (bexp>3)? 3: bexp;
    scaling = bexp-scaling;

    shift = BBE_MOVVSA32(scaling); 
    // Number of int16 which must be filled by zeros
    // for avoid influence unknown data to FFT scaling
    int num_to_fill = count * BBE_SIMD_WIDTH + 2*N/4*3; 

    ASSERT( N%(BBE_SIMD_WIDTH) == 0 ); 
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    //Filling input buffer by zeros beyond actual data
    if( num_to_fill > 0 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 0);
    }
    if( num_to_fill > 4 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 8);
    }
    if( num_to_fill > 8 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 16);
    }
    if( num_to_fill > 12 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 24);
    }    
    __Pragma("ymemory(p_dst)"); 

    v1 = BBE_LA_PP(px1); 
    v3 = BBE_LA_PP(px3); 
    BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
    BBE_LANX16_IP(_t1, v1, px1);                                    
    BBE_LVNX16_XP(_t2, px2, BBE_SIMD_WIDTH*2);                                    
    BBE_LANX16_IP(_t3, v3, px3);                                
                                                                    
    _t0 = BBE_SRANX16(_t0, shift);                                   
    _t1 = BBE_SRANX16(_t1, shift);                                   
    _t2 = BBE_SRANX16(_t2, shift);                                   
    _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                    
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_MOVSAV(_t2);                                                  
    BBE_MOVSBV(_t3);                                                  
          
    BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                  
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2); //INTLV(t0, t2);                                                
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2); //INTLV(t1, t3);         
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);  
    uu0 = BBE_MOVUVR(_t0);                             
    uu1 = BBE_MOVUVR(_t2);                             
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*2*BBE_SIMD_WIDTH, 0);     
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -2*BBE_SIMD_WIDTH, 0);     

    for (i=0; i<count-1; i++) 
    {   
        BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LVNX16_XP(_t2, px2, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t3, v3, px3);   

        _t0 = BBE_SRANX16(_t0, shift);                                   
        _t1 = BBE_SRANX16(_t1, shift);                                   
        _t2 = BBE_SRANX16(_t2, shift);                                   
        _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                        
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                      
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           
                                                                      
        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*2*BBE_SIMD_WIDTH, 0);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,    -2*BBE_SIMD_WIDTH, 0);
    }

    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                              
    BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*2*BBE_SIMD_WIDTH);                      
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*2*BBE_SIMD_WIDTH);

    return scaling; 
}

inline_ int R1_DFT4_L64_16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    int stride = N/4*2*sizeof(int16_t); 
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  

    

    int scaling = bexp;
    vsaN shift; 
  
    bexp = (bexp>3)? 3: bexp;
    scaling = bexp-scaling;

    shift = BBE_MOVVSA32(scaling); 


    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    
    __Pragma("ymemory(p_dst)"); 

    BBE_LVNX16_XP(_t0, p_src, stride);                                    
    BBE_LVNX16_XP(_t1, p_src, stride);                                    
    BBE_LVNX16_XP(_t2, p_src, stride);                                    
    BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                    
    _t0 = BBE_SRANX16(_t0, shift);                                   
    _t1 = BBE_SRANX16(_t1, shift);                                   
    _t2 = BBE_SRANX16(_t2, shift);                                   
    _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                    
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_MOVSAV(_t2);                                                  
    BBE_MOVSBV(_t3);                                                  
          
    BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                  
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2); //INTLV(t0, t2);                                                
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2); //INTLV(t1, t3);         
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);  
    uu0 = BBE_MOVUVR(_t0);                             
    uu1 = BBE_MOVUVR(_t2);                             
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*2*BBE_SIMD_WIDTH, 0);     
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -2*BBE_SIMD_WIDTH, 0);     

    for (i=0; i<count-1; i++) 
    {   
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t0, p_src, stride);    
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t1, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t2, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                        
        _t0 = BBE_SRANX16(_t0, shift);                                   
        _t1 = BBE_SRANX16(_t1, shift);                                   
        _t2 = BBE_SRANX16(_t2, shift);                                   
        _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                        
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                      
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           
                                                                      
        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*2*BBE_SIMD_WIDTH, 0);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,    -2*BBE_SIMD_WIDTH, 0);
    }

    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                              
    BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*2*BBE_SIMD_WIDTH);                      
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*2*BBE_SIMD_WIDTH);

    return scaling; 
}



inline_ int R2_DFT4xI4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/4*2*sizeof(int16_t); 

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
    valign align =  BBE_LAN_2XCQ15_PP(p_tw);

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

   {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                            
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);

        uu0 = BBE_MOVUVR(_t0);                             
        uu1 = BBE_MOVUVR(_t2);                             
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  2*BBE_SIMD_WIDTH, 1); 

    }

    for (i=0; i<count-1; i++) 
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }

    BBE_SALIGNVRNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH);        
    return scaling;
}

inline_ int R2_DFT4xI4_N16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = 2*N/4; 

    int count = N/4/(BBE_SIMD_WIDTH/2) + (N%(4*BBE_SIMD_WIDTH/2) > 0 );                                       
    int i;        
    VT * p_tw = (VT *)(tw);  

    xb_vecNx16 * px0 = (xb_vecNx16 *)(x);
    xb_vecNx16 * px1 = (xb_vecNx16 *)(x+stride*1);
    xb_vecNx16 * px2 = (xb_vecNx16 *)(x+stride*2);
    xb_vecNx16 * px3 = (xb_vecNx16 *)(x+stride*3);

    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
   // valign align =  BBE_LAN_2XCQ15_PP(p_tw);
    valign v1, v3;
    ASSERT( N%(BBE_SIMD_WIDTH) == 0 ); 

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

   {   
        v1 = BBE_LA_PP(px1); 
        v3 = BBE_LA_PP(px3); 

        BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LVNX16_XP(_t2, px2, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t3, v3, px3);                    
         
        LD_IX4_TW_UNPACK_3(0, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);

        uu0 = BBE_MOVUVR(_t0);                             
        uu1 = BBE_MOVUVR(_t2);                             
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  2*BBE_SIMD_WIDTH, 1); 

    }

    for (i=0; i<count-1; i++) 
    {   
        BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LVNX16_XP(_t2, px2, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t3, v3, px3);                           
         
        LD_IX4_TW_UNPACK_3(0, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }

    BBE_SALIGNVRNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH);        
    return scaling;
}

inline_ int R2_DFT4xIN_4(int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/4*2*sizeof(int16_t); 
    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
                                      
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
                                               
    VT t0, t1, t2, t3, t4, t5;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
     // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 0, bexp, scaling);
    __Pragma("no_reorder");

    for (i=0; i<count; i++) 
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVNX16_XP(_t2, p_src, stride);                                    
        BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
                
        DFT4(t0, t1, t2, t3, t4, t5, A, B, 0);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SVRNX16_XP(_t0, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t1, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t2, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t3, p_dst, -3*stride + 2*BBE_SIMD_WIDTH); 
    }
    return scaling;
}

inline_ int R2_DFT8xIN_8(int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/8*2*sizeof(int16_t); 
    int count = N/8/(BBE_SIMD_WIDTH/2);                                        
    int i;        
                                         
    xb_vecNx16 * p_src0 = (xb_vecNx16 *)(x);   
    xb_vecNx16 * p_src1 = (xb_vecNx16 *)(x+stride/2); 
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
                                               
    VT t0, t1, t2, t3, t4, t5, t6, t7, e4, e8;    
    VT x0, x1; 
    
    xb_vecNx16 _t0, _t1, _t2, _t3, _t4, _t5, _t6, _t7; 
    xb_vecNx16 a0, a1, b0, b1;  

    __Pragma("ymemory(p_src0)"); 
    __Pragma("ymemory(p_src1)");

    // RANGE_END(bexp);
    RANGE_BEGIN(8, -1, 0, bexp, scaling); 
    __Pragma("no_reorder");

    RX8_TW(-1, e4, e8);   

    for (i=0; i<count; i++) 
    {   
        BBE_LVNX16_XP(a0, p_src0, 2*stride);                                    
        BBE_LVNX16_XP(b0, p_src1, 2*stride);                                    
        BBE_LVNX16_XP(a1, p_src0, 2*stride);                                    
        BBE_LVNX16_XP(b1, p_src1, 2*stride);  
        BBE_LVA_XP(p_src0, 2*stride);                       
        BBE_LVC_XP(p_src1, 2*stride);                       
        BBE_LVB_XP(p_src0, -6*stride + 2*BBE_SIMD_WIDTH);   
        BBE_LVD_XP(p_src1, -6*stride + 2*BBE_SIMD_WIDTH);   
                                                
        x0 = BBE_MOVN_2XCQ15_FROMNX16(a0);                               
        x1 = BBE_MOVN_2XCQ15_FROMNX16(a1);  
                
        DFT4(t0, t1, t2, t3, x0, x1, A, B, 1);
                                                          
        x0 = BBE_MOVN_2XCQ15_FROMNX16(b0);                               
        x1 = BBE_MOVN_2XCQ15_FROMNX16(b1);  

        DFT4(t4, t5, t6, t7, x0, x1, C, D, 1);

        MUL(t5, e8);                
        MUL(t6, e4);                
        MUL_CONJ(t7, e8);           
                                    
        BFLY(t0, t4);               
        BFLY(t1, t5);               
        BFLY(t2, t6);               
        BFLY(t3, t7);               

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);        
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                                                      
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   
        _t4 = BBE_MOVNX16_FROMN_2XCQ15(t4);                           
        _t5 = BBE_MOVNX16_FROMN_2XCQ15(t5);        
        _t6 = BBE_MOVNX16_FROMN_2XCQ15(t6);                                                      
        _t7 = BBE_MOVNX16_FROMN_2XCQ15(t7);   

        BBE_SVRNX16_XP(_t0, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t1, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t2, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t7, p_dst, stride);     

        BBE_SVRNX16_XP(_t4, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t5, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t6, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t3, p_dst, -7*stride + 2*BBE_SIMD_WIDTH); 
    }
    return scaling;
}

inline_ ATTRIBUTE_ALWAYS_INLINE int R3_DFT3xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v, /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       int bexp
                                       )
{

    const int N1 = 3;
    const int stride = 2*N/N1; /*in int16_t */
    int i, j; 
//    unsigned int tmp;

    int num_bfls = N/N1/v;
//    xb_vecN_2xcq15 x0, x1, x2, t0, t1, t2; 
    xb_vecN_2xcq15 y0, y1, y2; 
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 tw1, tw2, r3_tw; 

    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    valign vtw = BBE_LAVNX16_PP(tw); 

    // int bexp; 
    int scaling = 0;
    // RANGE_END(bexp);
    RANGE_BEGIN(N1, -1, 1, bexp, scaling);

    BBE_MOVSBV(0); 
    BBE_MOVSDV(0); 

    __Pragma("no_reorder");

    RX3_TW(-1, r3_tw);	

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 
    if(v==(BBE_SIMD_WIDTH/2))
    {
        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1; 
            BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));

            tmp0 = BBE_REPNX16C(tmp0, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 0); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 

            BBE_LVA_IP(px0, 0); 
            BBE_LVC_IP(px0,  sizeof(*px0));

            BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2)); 

            DFT3(y0, y1, y2, r3_tw, 0);													
            MUL(y1, tw1);						
            MUL(y2, tw2);

            BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));
            BBE_SVRN_2XCQ15_IP(y1, py, sizeof(*py));
            BBE_SVRN_2XCQ15_IP(y2, py, sizeof(*py));
        }    
    }
    else if(v==(BBE_SIMD_WIDTH))
    {
        

        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1;
#if 0
            BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));
            tmp0 = BBE_REPNX16C(tmp0, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 0); 
#else
            BBE_LAVNX16_XP(tmp1, vtw, tw, 2*4); 
            tmp0 = BBE_REPNX16C(tmp1, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 1); 
#endif

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            {

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride); 
                                
                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride); 

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_XP(y0, py, sizeof(*py)+(N1-1)*2*v*sizeof(int16_t));
            }     
        }
    }
    else if(v==(BBE_SIMD_WIDTH/2*3))
    {
        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1; 
            BBE_LAVNX16_XP(tmp1, vtw, tw, 2*4); 
            tmp0 = BBE_REPNX16C(tmp1, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 1);

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 

            
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride); 

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);  

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride); 

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_XP(y0, py, sizeof(*py)+(N1-1)*2*v*sizeof(int16_t));

            }
         
        }
    }
    else if (v==(BBE_SIMD_WIDTH*3))
    {
        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1; 

            BBE_LAVNX16_XP(tmp1, vtw, tw, 2*4); 
            tmp0 = BBE_REPNX16C(tmp1, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 1); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 

            
            {            
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));
            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_XP(y0, py, sizeof(*py)+(N1-1)*2*v*sizeof(int16_t));
            }
         
        }
    }
    else if(v == BBE_SIMD_WIDTH*4)
    {
        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1; 

            BBE_LAVNX16_XP(tmp1, vtw, tw, 2*4); 
            tmp0 = BBE_REPNX16C(tmp1, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 1); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 

            
            {            
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));
            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            
                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));

            }
            {
            

                BBE_LVN_2XCQ15_XP(y2, px2, -2*stride);
                BBE_LVN_2XCQ15_XP(y1, px2, -2*stride);
                BBE_LVC_IP(px2,  0);
                BBE_LVA_XP(px2, 2*BBE_SIMD_WIDTH + 4*stride);

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_XP(y0, py, sizeof(*py)+(N1-1)*2*v*sizeof(int16_t));
            }
         
        }
    }
    else
    {
        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1; 
            BBE_LAVNX16_XP(tmp1, vtw, tw, 2*4); 
            tmp0 = BBE_REPNX16C(tmp1, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 1);

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 

            for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {            
                BBE_LVA_IP(px0, 0); 
                BBE_LVC_IP(px0,  sizeof(*px0));

                BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1));
                BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2)); 

                DFT3(y0, y1, y2, r3_tw, 0);													
                MUL(y1, tw1);						
                MUL(y2, tw2);

                BBE_SVRN_2XCQ15_X(y1, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(y2, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(y0, py, sizeof(*py));
            }
            py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
        }
    }
    return scaling; 
}



inline_ int R2_DFT2xI4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/2*2*sizeof(int16_t); 

    // RANGE_END(bexp);
    RANGE_BEGIN(2, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    int count = N/2/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    xb_vecNx16 * p_tw = (xb_vecNx16 *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0;                                                  
  //  VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
   // xb_vecNx16 _t0, _t1, _t2, _t3; 
    xb_vecNx16 x0, x1, y0, y1, tw1; 
   
    valign align =  BBE_LANX16_PP(p_tw);

    BBE_MOVSAV(0);                                                  
    BBE_MOVSBV(0); 
   {   
        NASSERT_ALIGN32(p_src); 
        BBE_LVNX16_XP(x0, p_src, stride);
        NASSERT_ALIGN32(p_src); 
        BBE_LVNX16_XP(x1, p_src, -stride + 2*BBE_SIMD_WIDTH);                           
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 
                                                                            
        uu0 = BBE_MOVUVR(y0);                                                          
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
    }

    for (i=0; i<count-1; i++) 
    {   

        BBE_LVNX16_XP(x0, p_src, stride);                                    
        BBE_LVNX16_XP(x1, p_src, -stride + 2*BBE_SIMD_WIDTH);                             
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
    }
    BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);

    return scaling;
}
inline_ ATTRIBUTE_ALWAYS_INLINE int R3_DFT4xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v, /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       int bexp
                                       )
{

    const int N1 = 4;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_bytes = 2*stride; 
    int i, j; 


    int num_bfls = N/N1/v;
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    const xb_vecNx16 *p_src = (xb_vecNx16*)x;
    
    xb_vecN_2xcq15 tw1, tw2, tw3; 
    xb_vecNx16 _t0, _t1;
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5 ;
    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    xb_vecNx16 tmp0, tmp1, tmp2; 

    // int bexp; 
    int scaling = 0;
    valign u = BBE_LAVNX16_PP(tw); 

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0);
    if(v == (BBE_SIMD_WIDTH/2) )
    {
        for(i=0; i<num_bfls; i++)
        {
           
            BBE_LAVNX16_XP(tmp2, u, tw, 3*4); 
            tmp0 = BBE_REPNX16C(tmp2, 0); 
            tmp1 = BBE_REPNX16C(tmp2, 1); 
            tmp2 = BBE_REPNX16C(tmp2, 2); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);

            // for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
            

                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                 BBE_SVRN_2XCQ15_IP(t0, py, 2*BBE_SIMD_WIDTH);
                 BBE_SVRN_2XCQ15_XP(t1, py, 2*BBE_SIMD_WIDTH); 
                 BBE_SVRN_2XCQ15_XP(t2, py, 2*BBE_SIMD_WIDTH); 
                 BBE_SVRN_2XCQ15_XP(t3, py, 2*BBE_SIMD_WIDTH);                  
            }

        }
    }//if(v == (BBE_SIMD_WIDTH/2) )
    else if(v == (BBE_SIMD_WIDTH) )
    {
        

        for(i=0; i<num_bfls; i++)
        {           
            BBE_LAVNX16_XP(tmp2, u, tw, 3*4); 
            tmp0 = BBE_REPNX16C(tmp2, 0); 
            tmp1 = BBE_REPNX16C(tmp2, 1); 
            tmp2 = BBE_REPNX16C(tmp2, 2); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);
 
        //    for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 =  BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_XP(t0, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t2, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t3, py, 2*BBE_SIMD_WIDTH-3*2*v*sizeof(int16_t));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP( p_src, stride_bytes);                                    
                BBE_LVD_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH); 

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_XP(t0, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t2, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t3, py, 2*BBE_SIMD_WIDTH);

            }
//            py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
        }

    }//if(v == (BBE_SIMD_WIDTH) )
    else if(v == 3*BBE_SIMD_WIDTH/2)
    {
        for(i=0; i<num_bfls; i++)
        {
           
            BBE_LAVNX16_XP(tmp2, u, tw, 3*4); 
            tmp0 = BBE_REPNX16C(tmp2, 0); 
            tmp1 = BBE_REPNX16C(tmp2, 1); 
            tmp2 = BBE_REPNX16C(tmp2, 2); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);

           // for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP(  p_src, stride_bytes);                                    
                BBE_LVD_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP(  p_src, stride_bytes);                                    
                BBE_LVB_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP( p_src, stride_bytes);                                    
                BBE_LVD_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t0, py, (N1-1)*2*v*sizeof(int16_t)+2*BBE_SIMD_WIDTH);

            }
        }
    }
    else if(v == 2*BBE_SIMD_WIDTH)
    {
        for(i=0; i<num_bfls; i++)
        {
           
            BBE_LAVNX16_XP(tmp2, u, tw, 3*4); 
            tmp0 = BBE_REPNX16C(tmp2, 0); 
            tmp1 = BBE_REPNX16C(tmp2, 1); 
            tmp2 = BBE_REPNX16C(tmp2, 2); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);

           // for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP(  p_src, stride_bytes);                                    
                BBE_LVB_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP(  p_src, stride_bytes);                                    
                BBE_LVB_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t0, py, (N1-1)*2*v*sizeof(int16_t)+2*BBE_SIMD_WIDTH);

            }

        }
    }
    else if(v == 4*BBE_SIMD_WIDTH)
    {
        for(i=0; i<num_bfls; i++)
        {
           
            BBE_LAVNX16_XP(tmp2, u, tw, 3*4); 
            tmp0 = BBE_REPNX16C(tmp2, 0); 
            tmp1 = BBE_REPNX16C(tmp2, 1); 
            tmp2 = BBE_REPNX16C(tmp2, 2); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);

           // for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP(  p_src, stride_bytes);                                    
                BBE_LVD_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP(  p_src, stride_bytes);                                    
                BBE_LVD_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP(  p_src, stride_bytes);                                    
                BBE_LVD_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP(  p_src, stride_bytes);                                    
                BBE_LVB_XP(  p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  


                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            {
                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVC_XP( p_src, stride_bytes);                                    
                BBE_LVD_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_XP(t0, py, (N1-1)*2*v*sizeof(int16_t)+2*BBE_SIMD_WIDTH);

            }

        }
    }
    else
    {
        for(i=0; i<num_bfls; i++)
        {
            BBE_LAVNX16_XP(tmp2, u, tw, 3*4); 
            tmp0 = BBE_REPNX16C(tmp2, 0); 
            tmp1 = BBE_REPNX16C(tmp2, 1); 
            tmp2 = BBE_REPNX16C(tmp2, 2); 


            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);

            for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
            

                BBE_LVNX16_XP(_t0, p_src, stride_bytes);                                    
                BBE_LVNX16_XP(_t1, p_src, stride_bytes);                                    
                BBE_LVA_XP( p_src, stride_bytes);                                    
                BBE_LVB_XP( p_src, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

                t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
                t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

                __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

                 BBE_SVRN_2XCQ15_X(t1, py, 1*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t2, py, 2*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_X(t3, py, 3*2*v*sizeof(int16_t)); 
                 BBE_SVRN_2XCQ15_IP(t0, py, sizeof(*py));

            }
            py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
        }

    }
    return scaling; 
}
inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT5xIN_5( int16_t *x   /*input*/, 
                          int16_t *y   /*output*/, 
                          const int N, int bexp)
{
    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    const int v = N/N1; 
    int j; 

    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15  r5_tw1, r5_tw2, r5_tw3; 
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;

    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    int scaling = 0;

    RANGE_BEGIN(N1, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 

    for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
    {
        BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
        BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
        BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
        BBE_LVB_IP(px3, sizeof(*px1)); 
        BBE_LVA_IP(px4, sizeof(*px1)); 
        t0 = BBE_FFTSRAN_2XCQ15(t0);

        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        

        BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(t5));							
        BBE_MOVSDV(BBE_MOVNX16_FROMN_2XCQ15(t6));							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								       
        BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

    }
  
    return scaling; 
}   //R2_DFT5xIN_5

inline_ ATTRIBUTE_ALWAYS_INLINE int R3_DFT5xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v, /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       int bexp
                                       )
{

    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    int i, j; 
//    unsigned int tmp;

    int num_bfls = N/N1/v;
//    xb_vecN_2xcq15 x0, x1, x2, t0, t1, t2; 
//    xb_vecN_2xcq15 y0, y1, y2; 
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15 tw1, tw2, tw3, tw4, r5_tw1, r5_tw2, r5_tw3; 
     xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;

    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    valign vtw = BBE_LAVNX16_PP(tw); 

    // int bexp; 
    int scaling = 0;
    // RANGE_END(bexp);
    RANGE_BEGIN(N1, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 
    if (v==16)
    {
        for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0 ,tmp1, tmp2, tmp3; 
#if 0
                
            BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp2, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp3, tw, 2*sizeof(int16_t));

            tmp0 = BBE_REPNX16C(tmp0, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 0); 
            tmp2 = BBE_REPNX16C(tmp2, 0); 
            tmp3 = BBE_REPNX16C(tmp3, 0); 
#else
            BBE_LAVNX16_XP(tmp0, vtw, tw, 4*4); 
            tmp1 = BBE_REPNX16C(tmp0, 1); 
            tmp2 = BBE_REPNX16C(tmp0, 2); 
            tmp3 = BBE_REPNX16C(tmp0, 3); 
            tmp0 = BBE_REPNX16C(tmp0, 0); 

#endif
            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
            tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 

           // for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
                BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
                BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
                BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
                BBE_LVB_IP(px3, sizeof(*px1)); 
                BBE_LVA_IP(px4, sizeof(*px1)); 


                t0 = BBE_FFTSRAN_2XCQ15(t0);

                /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
                DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
                								
                t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
                {  
                    xb_vecN_2xcq15 tmpc;
                    tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                    t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
                }                                    
                MUL(t4, r5_tw1);							
                MUL(t5, r5_tw2);							
                MUL(t6, r5_tw3);							
                
                tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
                tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

                BBE_MOVSCV(tmp0); //  MV ## C (t5);							
                BBE_MOVSDV(tmp1); //  MV ## D (t6);							
            									
                t3 = t4;								
                								
                DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
                								
                MUL(t4, tw1);								
                MUL(t5, tw2);								
                MUL(t7, tw3);								
                MUL(t6, tw4);								
 
                BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

            }
            
            {
                BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
                BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
                BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
                BBE_LVB_IP(px3, sizeof(*px1)); 
                BBE_LVA_IP(px4, sizeof(*px1)); 


                t0 = BBE_FFTSRAN_2XCQ15(t0);

                /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
                DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
                								
                t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
                {  
                    xb_vecN_2xcq15 tmpc;
                    tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                    t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
                }                                    
                MUL(t4, r5_tw1);							
                MUL(t5, r5_tw2);							
                MUL(t6, r5_tw3);							
                
                tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
                tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

                BBE_MOVSCV(tmp0); //  MV ## C (t5);							
                BBE_MOVSDV(tmp1); //  MV ## D (t6);							
            									
                t3 = t4;								
                								
                DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
                								
                MUL(t4, tw1);								
                MUL(t5, tw2);								
                MUL(t7, tw3);								
                MUL(t6, tw4);								
 
                BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

            }
            py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
        } 
    } 
    else //if (v==16)
    { 
      for(i=0; i<num_bfls; i++)
        {
            xb_vecNx16 tmp0, tmp1, tmp2, tmp3; 
            BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp2, tw, 2*sizeof(int16_t));
            BBE_LPNX16_IP(tmp3, tw, 2*sizeof(int16_t));

            tmp0 = BBE_REPNX16C(tmp0, 0); 
            tmp1 = BBE_REPNX16C(tmp1, 0); 
            tmp2 = BBE_REPNX16C(tmp2, 0); 
            tmp3 = BBE_REPNX16C(tmp3, 0); 

            tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
            tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
            tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
            tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 

            for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
            {
                BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
                BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
                BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
                BBE_LVB_IP(px3, sizeof(*px1)); 
                BBE_LVA_IP(px4, sizeof(*px1)); 

    //888888888888888888888888888888888888888888888888888888888888888
                /*
                LD_INC(t0, p_src, is);						
                LD_INC(t1, p_src, is);						
                LD_INC(t2, p_src, is);						
                LD_ ## B (p_src, is);						
                LD_ ## A (p_src, -4*is + ivs + iinc);				
                */		

                t0 = BBE_FFTSRAN_2XCQ15(t0);

                /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
                DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
                								
                t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
                {  
                    xb_vecN_2xcq15 tmpc;
                    tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                    t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
                }                                    
                MUL(t4, r5_tw1);							
                MUL(t5, r5_tw2);							
                MUL(t6, r5_tw3);							
                
                tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
                tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

                BBE_MOVSCV(tmp0); //  MV ## C (t5);							
                BBE_MOVSDV(tmp1); //  MV ## D (t6);							
            									
                t3 = t4;								
                								
                DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
                								
                MUL(t4, tw1);								
                MUL(t5, tw2);								
                MUL(t7, tw3);								
                MUL(t6, tw4);								
                /*							
                RST_INC(t1, p_dst, os);						
                RST_INC(t4, p_dst, os);						
                RST_INC(t5, p_dst, os);						
                RST_INC(t7, p_dst, os);						
                RST_INC(t6, p_dst, -4*os + ovs + oinc);				*/
                //888888888888888888888888888888888888888888888888888888888888888
                /*            MUL(y1, tw1);						
                MUL(y2, tw2);
                MUL(y3, tw3);						
                MUL(y4, tw4);*/

                BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
                BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

            }
            py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
        } 
    }
    return scaling; 
}


 

inline_ int R1_DFT4_L64_16_fs(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
ATTRIBUTE_ALWAYS_INLINE; 
inline_ int R2_DFT4xI4_wSEL(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
ATTRIBUTE_ALWAYS_INLINE; 

/* Fast scaling (scaling  without additional left shift) */
inline_ int R1_DFT4_L64_16_fs(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    int stride = N/4*2*sizeof(int16_t); 
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * px0 = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  
    int scaling = 0;  
    bexp = (bexp>4)? 4: bexp;
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    
    __Pragma("ymemory(p_dst)"); 

    BBE_LVNX16_XP(_t0, px0, stride);                                    
    BBE_LVNX16_XP(_t1, px0, stride);  

    BBE_LVA_XP(px0, stride); 
    BBE_LVB_XP(px0, -3*stride + 2*BBE_SIMD_WIDTH);   
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
     
    BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                  
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2); //INTLV(t0, t2);                                                
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2); //INTLV(t1, t3);         
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);  
    uu0 = BBE_MOVUVR(_t0);                             
    uu1 = BBE_MOVUVR(_t2);                             
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*2*BBE_SIMD_WIDTH, 0);     
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -2*BBE_SIMD_WIDTH, 0);     

    for (i=0; i<count-1; i++) 
    {   
        BBE_LVNX16_XP(_t0, px0, stride);    
        BBE_LVNX16_XP(_t1, px0, stride);

        BBE_LVA_XP(px0, stride); 
        BBE_LVB_XP(px0, -3*stride + 2*BBE_SIMD_WIDTH);   
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                      
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           
                                                                      
        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*2*BBE_SIMD_WIDTH, 0);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,    -2*BBE_SIMD_WIDTH, 0);
    }

    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                              
    BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*2*BBE_SIMD_WIDTH);                      
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*2*BBE_SIMD_WIDTH);

    return scaling; 
}

#define UNPK_TW(tw1, tw2, tw3, packed_tw, offset) \
{                                                               \
    xb_vecNx16 tmp0 = BBE_REPNX16C(packed_tw, 0+offset);        \
    xb_vecNx16 tmp1 = BBE_REPNX16C(packed_tw, 1+offset);        \
    xb_vecNx16 tmp2 = BBE_REPNX16C(packed_tw, 2+offset);        \
    tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0);                       \
    tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1);                       \
    tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2);                       \
}



/*
    Outer loop unrolled twice 
*/
inline_ ATTRIBUTE_ALWAYS_INLINE int R3_DFT4xIv_2U(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v, /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       int bexp
                                       )
{

    const int N1 = 4;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_bytes = 2*stride; 
    int j;
    int i;

    int num_bfls = N/N1/v;
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    const xb_vecNx16 *px0 = (xb_vecNx16*)(x);
    const xb_vecNx16 *px2 = (xb_vecNx16*)(x+1*2*v);
    
    xb_vecN_2xcq15 tw1, tw2, tw3; 
    xb_vecNx16 _t0, _t1;
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5 ;

    xb_vecN_2xcq15 *py0 = (xb_vecN_2xcq15 *)(y);
    xb_vecN_2xcq15 *py2 = (xb_vecN_2xcq15 *)(y+1*N1*2*v);

    // int bexp; 
    int scaling = 0;
    valign u = BBE_LAVNX16_PP(tw); 

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0);


    for(i=0; i<num_bfls/2; i++)
    {
        xb_vecNx16 tw50; 
        xb_vecN_2xcq15 tw11, tw21, tw31; 
        BBE_LAVNX16_XP(tw50,  u, tw, 6*4);
        
        UNPK_TW(tw1,  tw2, tw3, tw50,   0); 
        UNPK_TW(tw11, tw21, tw31, tw50, 3);

        for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
        {  
           
            BBE_LVNX16_XP(_t0, px0, stride_bytes);                                    
            BBE_LVNX16_XP(_t1, px0, stride_bytes);                                    
            BBE_LVA_XP( px0, stride_bytes);                                    
            BBE_LVB_XP( px0, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

             BBE_SVRN_2XCQ15_XP(t0, py0, 2*v*sizeof(int16_t));
             BBE_SVRN_2XCQ15_XP(t1, py0, 2*v*sizeof(int16_t)); 
             BBE_SVRN_2XCQ15_XP(t2, py0, 2*v*sizeof(int16_t)); 
             BBE_SVRN_2XCQ15_XP(t3, py0, 2*BBE_SIMD_WIDTH - 3*2*v*sizeof(int16_t)); 
            
            BBE_LVNX16_XP(_t0, px2, stride_bytes);                                    
            BBE_LVNX16_XP(_t1, px2, stride_bytes);                                    
            BBE_LVC_XP( px2, stride_bytes);                                    
            BBE_LVD_XP( px2, -3*stride_bytes + 2*BBE_SIMD_WIDTH);

            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);    

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw11, tw21, tw31);   

             BBE_SVRN_2XCQ15_XP(t0, py2, 2*v*sizeof(int16_t));
             BBE_SVRN_2XCQ15_XP(t1, py2, 2*v*sizeof(int16_t)); 
             BBE_SVRN_2XCQ15_XP(t2, py2, 2*v*sizeof(int16_t)); 
             BBE_SVRN_2XCQ15_XP(t3, py2, 2*BBE_SIMD_WIDTH - 3*2*v*sizeof(int16_t)); 
        }

        py0 += (2*N1-1)*2*v*sizeof(int16_t)/sizeof(*py0);
        py2 += (2*N1-1)*2*v*sizeof(int16_t)/sizeof(*py0);
        px0 += 2*v*sizeof(int16_t)/sizeof(*py0);
        px2 += 2*v*sizeof(int16_t)/sizeof(*py0);
    }

    return scaling; 
}  
 
inline_ int R2_DFT4xI4_wSEL(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/4*2*sizeof(int16_t); 

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
    valign align =  BBE_LAN_2XCQ15_PP(p_tw);

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");


    for (i=0; i<count; i++) 
    {   
        xb_vecNx16 y0, y1, y2, y3;
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        y0 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_LO_HALVES);		
        y1 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_LO_HALVES);	
        y2 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_HI_HALVES);		
        y3 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_HI_HALVES);	

        BBE_SVRNX16_XP(y0, p_dst, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y1, p_dst, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y2, p_dst, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y3, p_dst, 2*BBE_SIMD_WIDTH); 
    }
      
    return scaling;
}



inline_ int R2_DFT4xI4_U2(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
ATTRIBUTE_ALWAYS_INLINE; 

inline_ int R2_DFT4xI4_U2(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/4*2*sizeof(int16_t); 

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
    valign align =  BBE_LAN_2XCQ15_PP(p_tw);

    ASSERT((count&1)==0); 

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

   {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                            
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);

        uu0 = BBE_MOVUVR(_t0);                             
        uu1 = BBE_MOVUVR(_t2);                             
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  2*BBE_SIMD_WIDTH, 1); 

    }

    for (i=0; i < (count-2)/2; i++) 
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           

        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);

        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVC_XP( p_src, stride);                                    
        BBE_LVD_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);

    }
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }
    BBE_SALIGNVRNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH);        
    return scaling;
}

inline_ int R2_DFT6xIN_6(int16_t *x, int16_t *y,  int N, int bexp)
ATTRIBUTE_ALWAYS_INLINE; 

inline_ int R2_DFT6xIN_6(int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/6*2*sizeof(int16_t); 
    int count = N/6/(BBE_SIMD_WIDTH/2);                                        
    int i;        
                                      
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
       
    VT r3_tw; 
    xb_vecNx16 _t0, _t1, _t2, _t3, _t4, _t5; 
    // RANGE_END(bexp);
    RANGE_BEGIN(6, -1, 1, bexp, scaling);
    __Pragma("no_reorder");
    BBE_MOVSBV(0);	
    BBE_MOVSDV(0);	
    RX3_TW(-1, r3_tw);  
    for (i=0; i<count; i++) 
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVNX16_XP(_t2, p_src, stride);      
        BBE_LVNX16_XP(_t3, p_src, stride);                                    
        BBE_LVNX16_XP(_t4, p_src, stride);      
        BBE_LVNX16_XP(_t5, p_src, -5*stride + 2*BBE_SIMD_WIDTH);

        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw);

        BBE_SVNX16_XP(_t0, p_dst, stride);                                    
        BBE_SVNX16_XP(_t1, p_dst, stride);                                    
        BBE_SVNX16_XP(_t2, p_dst, stride);     
        BBE_SVNX16_XP(_t3, p_dst, stride);                                    
        BBE_SVNX16_XP(_t4, p_dst, stride);                                    
        BBE_SVNX16_XP(_t5, p_dst, -5*stride + 2*BBE_SIMD_WIDTH); 
    }
    return scaling;
}




inline_ int R2_DFT6xI4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
ATTRIBUTE_ALWAYS_INLINE; 

inline_ int R2_DFT6xI4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/6*2*sizeof(int16_t); 

    int count = N/6/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    xb_vecNx16 * p_tw = (xb_vecNx16 *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1, uu2;                                                  
    xb_vecNx16 tw1, tw2, tw3, tw4, tw5;                         
    xb_vecNx16 _t0, _t1, _t2, _t3, _t4, _t5; 
    //xb_vecNx16 _tw1, _tw2, _tw3, _tw4, _tw5;   
    VT r3_tw; 

    valign align =  BBE_LANX16_PP(p_tw);

    BBE_MOVSBV(0);	
    BBE_MOVSDV(0);	
    RX3_TW(-1, r3_tw);  
    // RANGE_END(bexp);
    RANGE_BEGIN(6, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

   {   
        
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVNX16_XP(_t2, p_src, stride);      
        BBE_LVNX16_XP(_t3, p_src, stride);                                    
        BBE_LVNX16_XP(_t4, p_src, stride);      
        BBE_LVNX16_XP(_t5, p_src, -5*stride + 2*BBE_SIMD_WIDTH);                          
         
        __LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 
        __LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5); 

        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw); 

        _MUL(_t1,  tw1);	
        _MUL(_t2,  tw2);	
        _MUL(_t3,  tw3);	
        _MUL(_t4,  tw4);	
        _MUL(_t5,  tw5);


        uu0 = BBE_MOVUVR(_t0);                             
        uu1 = BBE_MOVUVR(_t2);  
        uu2 = BBE_MOVUVR(_t4); 

        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  2*BBE_SIMD_WIDTH, 1); 
        BBE_SVINTLARNX16_XP(_t5, uu2, p_dst,  2*BBE_SIMD_WIDTH, 1); 

    }
 //   __Pragma("no_reorder");
    for (i=0; i<count-1; i++) 
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVNX16_XP(_t2, p_src, stride);      
        BBE_LVNX16_XP(_t3, p_src, stride);                                    
        BBE_LVNX16_XP(_t4, p_src, stride);      
        BBE_LVNX16_XP(_t5, p_src, -5*stride + 2*BBE_SIMD_WIDTH);                          
         
        __LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 
        __LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5); 

        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw); 
                                                                            
        _MUL(_t1, tw1);	
        _MUL(_t2, tw2);	
        _MUL(_t3, tw3);	
        _MUL(_t4, tw4);	
        _MUL(_t5, tw5);	

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);  
        BBE_SALIGNVRNX16_XP(_t4, uu2, p_dst, 2*BBE_SIMD_WIDTH);  
        
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
        BBE_SVINTLARNX16_XP(_t5, uu2, p_dst, 2*BBE_SIMD_WIDTH, 1);

    }

    BBE_SALIGNVRNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH); 
    BBE_SALIGNVRNX16_XP(_t5, uu2, p_dst, 2*BBE_SIMD_WIDTH);

    return scaling;
}


inline_ int R2_DFT3xIN_3(int16_t *x, int16_t *y,  int N, int bexp)
ATTRIBUTE_ALWAYS_INLINE; 

inline_ int R2_DFT3xIN_3(int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/3*2*sizeof(int16_t); 
    int count = N/3/(BBE_SIMD_WIDTH/2);                                        
    int i;        
                                      
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
                                               
//    VT t0, t1, t2, t3, t4, t5;                         
    xb_vecNx16 _t0, _t1, _t2; 
    xb_vecN_2xcq15 y0, y1, y2, r3_tw;

    // RANGE_END(bexp);
    RANGE_BEGIN(3, -1, 1, bexp, scaling);

    RX3_TW(-1, r3_tw);	

    BBE_MOVSBV(0); 
    BBE_MOVSDV(0); 
    __Pragma("no_reorder");
    for (i=0; i<count; i++) 
    {
#if 0
        BBE_LVA_IP(p_src, 0); 
        BBE_LVC_XP(p_src,  stride);

       // BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVNX16_XP(_t2, p_src, -2*stride + 2*BBE_SIMD_WIDTH);                           
/*
        BBE_MOVSAV(_t0);                                                  
        BBE_MOVSCV(_t0);
*/
#else
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVNX16_XP(_t2, p_src, -2*stride + 2*BBE_SIMD_WIDTH);                           

        BBE_MOVSAV(_t0);                                                  
        BBE_MOVSCV(_t0);
#endif
        y1 = BBE_MOVN_2XCQ15_FROMNX16(_t1); 
        y2 = BBE_MOVN_2XCQ15_FROMNX16(_t2); 

        DFT3(y0, y1, y2, r3_tw, 0);	              
        
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(y0);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(y1);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(y2);                           

        BBE_SVNX16_XP(_t0, p_dst, stride);                                    
        BBE_SVNX16_XP(_t1, p_dst, stride);                                    
        BBE_SVNX16_XP(_t2, p_dst, -2*stride + 2*BBE_SIMD_WIDTH); 
    }
    return scaling;
}
#if 1
inline_ ATTRIBUTE_ALWAYS_INLINE int R1_DFT4_L64_16_N4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{

    int stride = 2*N/4; 
    int count = N/4/(BBE_SIMD_WIDTH/2) + (N%(4*BBE_SIMD_WIDTH/2) > 0 );                                        
    int i;        
    VT * p_tw = (VT *)(tw);  

    xb_vecNx16 * px0 = (xb_vecNx16 *)(x);
    xb_vecNx16 * px1 = (xb_vecNx16 *)(x+stride*1);
    xb_vecNx16 * px2 = (xb_vecNx16 *)(x+stride*2);
    xb_vecNx16 * px3 = (xb_vecNx16 *)(x+stride*3);

    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0;
    valign v1, v2, v3; 
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  

    // Number of int16 which must be filled by zeros
    // for avoid influence unknown data to FFT scaling
    int num_to_fill = count * BBE_SIMD_WIDTH + 2*N/4*3 - 2*N; 


    int scaling = bexp;
    vsaN shift; 
  
    bexp = (bexp>3)? 3: bexp;
    scaling = bexp-scaling;
    shift = BBE_MOVVSA32(scaling); 
    
    //Filling input buffer by zeros beyond actual data

    if( num_to_fill > 0 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 0);
    }
    if( num_to_fill > 4 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 8);
    }
    if( num_to_fill > 8 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 16);
    }
    if( num_to_fill > 12 )
    {
         BBE_SV4X16_I(0, (void*)(x+2*N), 24);
    }


    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    
    __Pragma("ymemory(p_dst)"); 

    v1 = BBE_LA_PP(px1); 
    v2 = BBE_LA_PP(px2); 
    v3 = BBE_LA_PP(px3); 


    if(count<5)
    {
     {   
            BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
            BBE_LANX16_IP(_t1, v1, px1);                                    
            BBE_LANX16_IP(_t2, v2, px2);                                  
            BBE_LANX16_IP(_t3, v3, px3);   

            _t0 = BBE_SRANX16(_t0, shift);                                   
            _t1 = BBE_SRANX16(_t1, shift);                                   
            _t2 = BBE_SRANX16(_t2, shift);                                   
            _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                            
            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            BBE_MOVSAV(_t2);                                                  
            BBE_MOVSBV(_t3);                                                  
                  
            BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                          
            BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
            BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

            _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
            _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
            _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
            _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

     
            uu0 = BBE_MOVUVR(_t0);                  
            BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);                
            BBE_SALIGNVRNX16_XP(_t2, uu0, p_dst,    2*BBE_SIMD_WIDTH);            
            BBE_SVINTLARNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);
        }   

        if(count>1)
        {   
            BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
            BBE_LANX16_IP(_t1, v1, px1);                                    
            BBE_LANX16_IP(_t2, v2, px2);                                  
            BBE_LANX16_IP(_t3, v3, px3);   

            _t0 = BBE_SRANX16(_t0, shift);                                   
            _t1 = BBE_SRANX16(_t1, shift);                                   
            _t2 = BBE_SRANX16(_t2, shift);                                   
            _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                            
            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            BBE_MOVSAV(_t2);                                                  
            BBE_MOVSBV(_t3);                                                  
                  
            BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                          
            BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
            BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

            _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
            _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
            _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
            _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

     
            BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*BBE_SIMD_WIDTH);                   
            BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);                
            BBE_SALIGNVRNX16_XP(_t2, uu0, p_dst,    2*BBE_SIMD_WIDTH);            
            BBE_SVINTLARNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);
        }
        if(count>2)
        {   
            BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
            BBE_LANX16_IP(_t1, v1, px1);                                    
            BBE_LANX16_IP(_t2, v2, px2);                                  
            BBE_LANX16_IP(_t3, v3, px3);   

            _t0 = BBE_SRANX16(_t0, shift);                                   
            _t1 = BBE_SRANX16(_t1, shift);                                   
            _t2 = BBE_SRANX16(_t2, shift);                                   
            _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                            
            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            BBE_MOVSAV(_t2);                                                  
            BBE_MOVSBV(_t3);                                                  
                  
            BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                          
            BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
            BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

            _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
            _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
            _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
            _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

     
            BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*BBE_SIMD_WIDTH);                   
            BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);                
            BBE_SALIGNVRNX16_XP(_t2, uu0, p_dst,    2*BBE_SIMD_WIDTH);            
            BBE_SVINTLARNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);
        }
        if(count>3)
        {   
            BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
            BBE_LANX16_IP(_t1, v1, px1);                                    
            BBE_LANX16_IP(_t2, v2, px2);                                  
            BBE_LANX16_IP(_t3, v3, px3);   

            _t0 = BBE_SRANX16(_t0, shift);                                   
            _t1 = BBE_SRANX16(_t1, shift);                                   
            _t2 = BBE_SRANX16(_t2, shift);                                   
            _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                            
            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            BBE_MOVSAV(_t2);                                                  
            BBE_MOVSBV(_t3);                                                  
                  
            BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                          
            BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
            BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

            _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
            _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
            _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
            _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

     
            BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*BBE_SIMD_WIDTH);                   
            BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);                
            BBE_SALIGNVRNX16_XP(_t2, uu0, p_dst,    2*BBE_SIMD_WIDTH);            
            BBE_SVINTLARNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);
        }
    }
    else //if(count==3)
    {
        {   
            BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
            BBE_LANX16_IP(_t1, v1, px1);                                    
            BBE_LANX16_IP(_t2, v2, px2);                                  
            BBE_LANX16_IP(_t3, v3, px3);   

            _t0 = BBE_SRANX16(_t0, shift);                                   
            _t1 = BBE_SRANX16(_t1, shift);                                   
            _t2 = BBE_SRANX16(_t2, shift);                                   
            _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                            
            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            BBE_MOVSAV(_t2);                                                  
            BBE_MOVSBV(_t3);                                                  
                  
            BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                          
            BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
            BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

            _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
            _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
            _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
            _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

     
            uu0 = BBE_MOVUVR(_t0);                  
            BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);                
            BBE_SALIGNVRNX16_XP(_t2, uu0, p_dst,    2*BBE_SIMD_WIDTH);            
            BBE_SVINTLARNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);
        }   

        for (i=0; i<count-1; i++) 
        {   
            BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
            BBE_LANX16_IP(_t1, v1, px1);                                    
            BBE_LANX16_IP(_t2, v2, px2);                                  
            BBE_LANX16_IP(_t3, v3, px3);   

            _t0 = BBE_SRANX16(_t0, shift);                                   
            _t1 = BBE_SRANX16(_t1, shift);                                   
            _t2 = BBE_SRANX16(_t2, shift);                                   
            _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                            
            t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
            t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

            BBE_MOVSAV(_t2);                                                  
            BBE_MOVSBV(_t3);                                                  
                  
            BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
            BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

            __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                          
            BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
            BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

            _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
            _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
            _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
            _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

     
            BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*BBE_SIMD_WIDTH);                   
            BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);                
            BBE_SALIGNVRNX16_XP(_t2, uu0, p_dst,    2*BBE_SIMD_WIDTH);            
            BBE_SVINTLARNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH, 0);
        }
    } // if(count==3).. else
        BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst,    2*BBE_SIMD_WIDTH); 
    

    return scaling; 
}

#else
inline_ ATTRIBUTE_ALWAYS_INLINE int R1_DFT4_L64_16_N4(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
     

    int stride = 2*N/4; 

    int count = N/4/(BBE_SIMD_WIDTH/2) + (N%(4*BBE_SIMD_WIDTH/2) > 0 );                                       
    int i;        
    VT * p_tw = (VT *)(tw);  

    xb_vecNx16 * px0 = (xb_vecNx16 *)(x);
    xb_vecNx16 * px1 = (xb_vecNx16 *)(x+stride*1);
    xb_vecNx16 * px2 = (xb_vecNx16 *)(x+stride*2);
    xb_vecNx16 * px3 = (xb_vecNx16 *)(x+stride*3);

    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
                                                      
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
   // valign align =  BBE_LAN_2XCQ15_PP(p_tw);
    valign v1, v2, v3;
    ASSERT( N%(4) == 0 ); 

    int scaling = bexp;
    vsaN shift; 
  
    bexp = (bexp>3)? 3: bexp;
    scaling = bexp-scaling;
    shift = BBE_MOVVSA32(scaling); 


    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");


    v1 = BBE_LA_PP(px1); 
    v2 = BBE_LA_PP(px2);
    v3 = BBE_LA_PP(px3); 

    for (i=0; i<count; i++) 
    {   
        BBE_LVNX16_XP(_t0, px0, BBE_SIMD_WIDTH*2);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LANX16_IP(_t2, v2, px2);                                    
        BBE_LANX16_IP(_t3, v3, px3);                         

        _t0 = BBE_SRANX16(_t0, shift);                                   
        _t1 = BBE_SRANX16(_t1, shift);                                   
        _t2 = BBE_SRANX16(_t2, shift);                                   
        _t3 = BBE_SRANX16(_t3, shift);   

        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        BBE_DSELN_2XCQ15I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_2);          
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_4);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

        BBE_SVRNX16_XP(_t0, p_dst,    2*BBE_SIMD_WIDTH);                   
        BBE_SVRNX16_XP(_t2, p_dst,    2*BBE_SIMD_WIDTH);            
        BBE_SVRNX16_XP(_t1, p_dst,    2*BBE_SIMD_WIDTH);                
        BBE_SVRNX16_XP(_t3, p_dst,    2*BBE_SIMD_WIDTH);
    }
    
    return scaling;
}
#endif

inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT3xI4_expand(const int16_t *ptw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/3*2;  /*in int16_t */ 

    int count = (N+3*BBE_SIMD_WIDTH/2-1)/(3*BBE_SIMD_WIDTH/2);                                        
    int i;        
    const xb_vecN_2xcq15  *tw = (xb_vecN_2xcq15 *)ptw;
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 tw1, tw2, r3_tw; 
    valign v1 = BBE_LAN_2XCQ15_PP(px1);
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
//    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
    xb_vecN_2xcq15 y0, y1, y2; 
   // valign align =  BBE_LAN_2XCQ15_PP(tw);
    xb_vecN_2xcq15 __tw; 

    BBE_MOVSBV(0); 
    BBE_MOVSDV(0);    


    // RANGE_END(bexp);
    RANGE_BEGIN(3, -1, 1, bexp, scaling);
    RX3_TW(-1, r3_tw);	

    __Pragma("no_reorder");
        
    _t3 = 0; 
   {   
        

        BBE_LVN_2XCQ15_XP(__tw, tw, 8*4);                
        tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));

        BBE_LAN_2XCQ15_IP(y1, v1, px1);
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2)); 

        DFT3(y0, y1, y2, r3_tw, 0);													
        MUL(y1, tw1);						
        MUL(y2, tw2);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(y0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(y2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(y1);                           
        

        uu0 = BBE_MOVUVR(_t0);                             
        uu1 = BBE_MOVUVR(_t2);                             
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  2*BBE_SIMD_WIDTH, 1); 
    }
#if 0
    for (i=0; i<count-1; i++) 
    {   
        
        LDA_IX4_TW_UNPACK_2(align, tw, tw1, tw2); 

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));

        BBE_LAN_2XCQ15_IP(y1, v1, px1);
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2)); 

        DFT3(y0, y1, y2, r3_tw, 0);													
        MUL(y1, tw1);						
        MUL(y2, tw2);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(y0);    
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(y1);                                   
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(y2);                           

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }
#else
    __Pragma("ymemory(px0)");
    __Pragma("ymemory(px1)");
    __Pragma("ymemory(px2)");

    for (i=0; i<(count-1)/2; i++) 
    {   


        tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);
        tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_3);

        ///////
        BBE_LVA_IP(px0, 0); 
        BBE_LVC_XP(px0,  sizeof(*px0));

        BBE_LAN_2XCQ15_IP(y1, v1, px1);
        BBE_LVN_2XCQ15_XP(y2, px2, sizeof(*px2)); 

        DFT3(y0, y1, y2, r3_tw, 0);													
        MUL(y1, tw1);						
        MUL(y2, tw2);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(y0);    
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(y1);                                   
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(y2);                           

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
        ///////
        BBE_LVN_2XCQ15_XP(__tw, tw, 8*4);                
        tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_1);

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_XP(px0,  sizeof(*px0));

        BBE_LAN_2XCQ15_IP(y1, v1, px1);
        BBE_LVN_2XCQ15_XP(y2, px2, sizeof(*px2)); 

        DFT3(y0, y1, y2, r3_tw, 0);													
        MUL(y1, tw1);						
        MUL(y2, tw2);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(y0);    
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(y1);                                   
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(y2);                           

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }

    if((count-1)%2)
    {   
        
        tw1 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_2);
        tw2 = BBE_SHFLN_2XCQ15I(__tw,  BBE_SHFLI_REP_2X4_OFFSET_3);

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));

        BBE_LAN_2XCQ15_IP(y1, v1, px1);
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2)); 

        DFT3(y0, y1, y2, r3_tw, 0);													
        MUL(y1, tw1);						
        MUL(y2, tw2);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(y0);    
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(y1);                                   
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(y2);                           

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }


#endif
    BBE_SALIGNVRNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH);        
    return scaling;
}

inline_ int ATTRIBUTE_ALWAYS_INLINE R2_DFT3xIN_3_compress( 
                                       int16_t *x, 
                                       int16_t *y, 
                                       const int N, int bexp )
{
    const int v = 8;       // vector size
    const int N1 = 3;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_out = 2*N*3/4/N1; /*in int16_t */
    int scaling = 0; 
    int i; 

    int num_bfls = N/N1/v;
   // int   num_frac_bfls = N/N1 - v * num_bfls;
    xb_vecN_2xcq15 y0, y1, y2; 

    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);

    xb_vecN_2xcq15 *py0 = (xb_vecN_2xcq15 *)(y);
    xb_vecN_2xcq15 *py1 = (xb_vecN_2xcq15 *)(y + 1*stride_out);
    xb_vecN_2xcq15 *py2 = (xb_vecN_2xcq15 *)(y + 2*stride_out);

    xb_vecN_2xcq15 r3_tw;
 
    valign a0 =  BBE_ZALIGN();
    valign a1 =  BBE_ZALIGN();
    valign a2 =  BBE_ZALIGN();
  
    RX3_TW(-1, r3_tw);

    BBE_MOVSBV(0); 
    BBE_MOVSDV(0); 

    // RANGE_END(bexp);
    RANGE_BEGIN(3, -1, 0, bexp, scaling);

    ASSERT((N/N1 - v * num_bfls)==0); 
    ASSERT(num_bfls%2==0); 

    for(i=0; i<num_bfls/2; i++)
    {
    								
        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));  
        BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1) );
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2) ); 

        DFT3(y0, y1, y2, r3_tw, 0);	

        BBE_SAVRN_2XCQ15_XP(y0, a0, py0, 8*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y1, a1, py1, 8*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y2, a2, py2, 8*2*sizeof(int16_t)); 

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));  
        BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1) );
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2) ); 

        DFT3(y0, y1, y2, r3_tw, 0);	

        BBE_SAVRN_2XCQ15_XP(y0, a0, py0, 4*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y1, a1, py1, 4*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y2, a2, py2, 4*2*sizeof(int16_t)); 
    }

    BBE_SAN_2XCQ15POS_FC(a0, py0);	
    BBE_SAN_2XCQ15POS_FC(a1, py1);
    BBE_SAN_2XCQ15POS_FC(a2, py2);

    return scaling; 
}



inline_ int  ATTRIBUTE_ALWAYS_INLINE  R2_DFT5xIN_5_compress( 
                                       int16_t *x, 
                                       int16_t *y, 
                                       const int N, int bexp )
{
    const int v = 8;       // vector size
    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_out = 2*N*3/4/N1; /*in int16_t */
    
    int i; 

    int   num_bfls = N/N1/v;
    int scaling = 0 ;
   
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15 r5_tw1, r5_tw2, r5_tw3; 
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;

    xb_vecN_2xcq15 *py0 = (xb_vecN_2xcq15 *)(y);
    xb_vecN_2xcq15 *py1 = (xb_vecN_2xcq15 *)(y + 1*stride_out);
    xb_vecN_2xcq15 *py2 = (xb_vecN_2xcq15 *)(y + 2*stride_out);
    xb_vecN_2xcq15 *py3 = (xb_vecN_2xcq15 *)(y + 3*stride_out);
    xb_vecN_2xcq15 *py4 = (xb_vecN_2xcq15 *)(y + 4*stride_out);

    valign a0 =  BBE_ZALIGN();
    valign a1 =  BBE_ZALIGN();
    valign a2 =  BBE_ZALIGN();
    valign a3 =  BBE_ZALIGN();
    valign a4 =  BBE_ZALIGN();

    // RANGE_END(bexp);
    RANGE_BEGIN(N1, -1, 1, bexp, scaling);

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);



    ASSERT((N/N1 - v * num_bfls)==0); 
    ASSERT(num_bfls%2 == 0); 

    for(i=0; i<num_bfls/2; i++)
    {
            xb_vecNx16 tmp0, tmp1; 
            BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
            BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
            BBE_LVB_IP(px3, sizeof(*px1)); 
            BBE_LVA_IP(px4, sizeof(*px1)); 

            t0 = BBE_FFTSRAN_2XCQ15( t0);						
            								
            /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
            DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
            								
            t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
            {  
                xb_vecN_2xcq15 tmpc;
                tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
            }                                    
            MUL(t4, r5_tw1);							
            MUL(t5, r5_tw2);							
            MUL(t6, r5_tw3);							
            
            tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
            tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

            BBE_MOVSCV(tmp0); //  MV ## C (t5);							
            BBE_MOVSDV(tmp1); //  MV ## D (t6);							
        									
            t3 = t4;								
            								
            DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);	
	

            BBE_SAVRN_2XCQ15_XP(t1, a0, py0, 8*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t4, a1, py1, 8*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t5, a2, py2, 8*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t7, a3, py3, 8*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t6, a4, py4, 8*2*sizeof(int16_t));

            // ============== Store only 4 complex samples from earch butterfly's output ==============================
            BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
            BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
            BBE_LVB_IP(px3, sizeof(*px1)); 
            BBE_LVA_IP(px4, sizeof(*px1)); 

            t0 = BBE_FFTSRAN_2XCQ15( t0);						
            								
            /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
            DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
            								
            t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
            {  
                xb_vecN_2xcq15 tmpc;
                tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
            }                                    
            MUL(t4, r5_tw1);							
            MUL(t5, r5_tw2);							
            MUL(t6, r5_tw3);							
            
            tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
            tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

            BBE_MOVSCV(tmp0); //  MV ## C (t5);							
            BBE_MOVSDV(tmp1); //  MV ## D (t6);							
        									
            t3 = t4;								
            								
            DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);	
	
            BBE_SAVRN_2XCQ15_XP(t1, a0, py0, 4*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t4, a1, py1, 4*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t5, a2, py2, 4*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t7, a3, py3, 4*2*sizeof(int16_t)); 
            BBE_SAVRN_2XCQ15_XP(t6, a4, py4, 4*2*sizeof(int16_t));
    }

    BBE_SAN_2XCQ15POS_FC(a0, py0);	
    BBE_SAN_2XCQ15POS_FC(a1, py1);
    BBE_SAN_2XCQ15POS_FC(a2, py2);
    BBE_SAN_2XCQ15POS_FC(a3, py3);
    BBE_SAN_2XCQ15POS_FC(a4, py4);

    return scaling; 
}


inline_ int ATTRIBUTE_ALWAYS_INLINE R2_DFT2xI4_N8(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/2*2*sizeof(int16_t); 

    // RANGE_END(bexp);
    RANGE_BEGIN(2, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    int count = (N+BBE_SIMD_WIDTH-1)/BBE_SIMD_WIDTH;                                        
    int i;        
    xb_vecNx16 * p_tw = (xb_vecNx16 *)(tw);                                        
    xb_vecNx16 * px0 = (xb_vecNx16 *)(x);
    xb_vecNx16 * px1 = (xb_vecNx16 *)(x+stride/2); 
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);          
    valign v1 = BBE_LANX16_PP(px1);
    valign uu0;                                                  

    xb_vecNx16 x0, x1, y0, y1, tw1; 
   
    valign align =  BBE_LANX16_PP(p_tw);

    BBE_MOVSAV(0);                                                  
    BBE_MOVSBV(0); 


   {   

        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                           
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 
                                                                            
        uu0 = BBE_MOVUVR(y0);                                                          
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
    }
#if 0
    for (i=0; i<count-1; i++) 
    {   

        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
    }
#else
    for (i=0; i<(count-1)/4; i++) 
    {   
        xb_vecNx16 tmp;
        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        BBE_LAVNX16_XP(tmp, align, p_tw, 8*4);                
        tw1 = BBE_SHFLNX16I(tmp,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
//1
        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        tw1 = BBE_SHFLNX16I(tmp,  BBE_SHFLI_REP_2X4_OFFSET_1);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
//2
              
        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        tw1 = BBE_SHFLNX16I(tmp,  BBE_SHFLI_REP_2X4_OFFSET_2);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 
    
        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
// 3
        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        tw1 = BBE_SHFLNX16I(tmp,  BBE_SHFLI_REP_2X4_OFFSET_3);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                

    }
    if((count-1)%4 > 0)
    {   

        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
    }
    if((count-1)%4 > 1)
    {   

        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
    }
    if((count-1)%4 > 2)
    {   

        BBE_LVNX16_XP(x0, px0, 2*BBE_SIMD_WIDTH);
        BBE_LANX16_IP(x1, v1, px1);                                                      
         
        BBE_LAVNX16_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLNX16I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);        
	
        y0 = BBE_FFTADD4SABNX16( x0, x1, 0, 0 );
        y1 = BBE_FFTADD4SABNX16( x0, x1, 2, 0 );
                                                               
        _MUL(y1, tw1); 

        BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SVINTLARNX16_XP(y1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
    }

#endif
    BBE_SALIGNVRNX16_XP(y0,  uu0, p_dst, 2*BBE_SIMD_WIDTH);

    return scaling;
}


inline_ ATTRIBUTE_ALWAYS_INLINE int R1_DFT4_L64_16_norm(const int16_t *tw, int16_t *x, int16_t *y,  int N)
{

    int scaling = 0;
    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    int stride = N/4*2*sizeof(int16_t); 
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  

    RANGE_BEGIN(4, -1, 1, 0, scaling);
    
    __Pragma("ymemory(p_dst)"); 

    BBE_LVNX16_XP(_t0, p_src, stride);                                    
    BBE_LVNX16_XP(_t1, p_src, stride);                                    
    BBE_LVNX16_XP(_t2, p_src, stride);                                    
    BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                                                                                                                      
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_MOVSAV(_t2);                                                  
    BBE_MOVSBV(_t3);                                                  
          
    BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
    BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                  
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2); //INTLV(t0, t2);                                                
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2); //INTLV(t1, t3);         
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);  
    uu0 = BBE_MOVUVR(_t0);                             
    uu1 = BBE_MOVUVR(_t2);                             
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*2*BBE_SIMD_WIDTH, 0);     
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -2*BBE_SIMD_WIDTH, 0);     

    for (i=0; i<count-1; i++) 
    {   
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t0, p_src, stride);    
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t1, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t2, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                                                                                                                                 
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                      
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           
                                                                      
        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*2*BBE_SIMD_WIDTH, 0);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,    -2*BBE_SIMD_WIDTH, 0);
    }

    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                              
    BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*2*BBE_SIMD_WIDTH);                      
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*2*BBE_SIMD_WIDTH);

    return scaling; 
}


/* First stage radix 4, 
   suffix _mr mean : memory is input, registers is output  */

inline_ ATTRIBUTE_ALWAYS_INLINE int R1_DFT4_L64_16_mr(const int16_t *tw, int16_t *x, xb_vecNx16 *y,  int N, int bexp)
{

    int count = N/4/(BBE_SIMD_WIDTH/2);                                             
    int stride = N/4*2*sizeof(int16_t); 
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                                                                       
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  

    

    int scaling = bexp;
    vsaN shift; 
  
    bexp = (bexp>3)? 3: bexp;
    scaling = bexp-scaling;

    shift = BBE_MOVVSA32(scaling); 


    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    

    ASSERT(count==2); 
     if(count==2)
    {   
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t0, p_src, stride);    
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t1, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t2, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                        
        _t0 = BBE_SRANX16(_t0, shift);                                   
        _t1 = BBE_SRANX16(_t1, shift);                                   
        _t2 = BBE_SRANX16(_t2, shift);                                   
        _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                        
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

        BBE_RANGENX16(_t0); 
        BBE_RANGENX16(_t1); 
        BBE_RANGENX16(_t2); 
        BBE_RANGENX16(_t3);


     
        BBE_DSELN_2XCQ15I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_2);          
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_4);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           


        y[0]=_t0; //BBE_SVRNX16_XP(_t0, p_dst,    2*BBE_SIMD_WIDTH);                   
        y[1]=_t2; //BBE_SVRNX16_XP(_t2, p_dst,    2*BBE_SIMD_WIDTH);            
        y[2]=_t1; //BBE_SVRNX16_XP(_t1, p_dst,    2*BBE_SIMD_WIDTH);                
        y[3]=_t3; //BBE_SVRNX16_XP(_t3, p_dst,    2*BBE_SIMD_WIDTH);
    }
    {   
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t0, p_src, stride);    
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t1, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t2, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                        
        _t0 = BBE_SRANX16(_t0, shift);                                   
        _t1 = BBE_SRANX16(_t1, shift);                                   
        _t2 = BBE_SRANX16(_t2, shift);                                   
        _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                        
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        BBE_LVN_2XCQ15_IP(tw1, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw2, p_tw, 2*BBE_SIMD_WIDTH); 
        BBE_LVN_2XCQ15_IP(tw3, p_tw, 2*BBE_SIMD_WIDTH); 

        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
 
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

        BBE_RANGENX16(_t0); 
        BBE_RANGENX16(_t1); 
        BBE_RANGENX16(_t2); 
        BBE_RANGENX16(_t3);

        BBE_DSELN_2XCQ15I(t1, t0, t1, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t2, t3, t2, BBE_DSELI_INTERLEAVE_2);          
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_4);
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_4);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           

        y[4]=_t0; //BBE_SVRNX16_XP(_t0, p_dst,    2*BBE_SIMD_WIDTH);                   
        y[5]=_t2; //BBE_SVRNX16_XP(_t2, p_dst,    2*BBE_SIMD_WIDTH);            
        y[6]=_t1; //BBE_SVRNX16_XP(_t1, p_dst,    2*BBE_SIMD_WIDTH);                
        y[7]=_t3; //BBE_SVRNX16_XP(_t3, p_dst,    2*BBE_SIMD_WIDTH);
    }

    return scaling; 
}

inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT4xI4_rr(const int16_t *tw, xb_vecNx16 *x, xb_vecNx16 *y,  int N, int bexp)
{
    int scaling = 0;

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    VT * p_tw = (VT *)(tw);                                        
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
    valign align =  BBE_LAN_2XCQ15_PP(p_tw);

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);


    ASSERT(count==2); //for (i=0; i<count; i++) 
    if(count==2)
    {   
        xb_vecNx16 y0, y1, y2, y3; 

        _t0 = x[0];                            
        _t1 = x[2];                            
        BBE_MOVSAV(x[4]);               
        BBE_MOVSBV(x[6]);                            
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_RANGENX16(_t0); 
        BBE_RANGENX16(_t1); 
        BBE_RANGENX16(_t2); 
        BBE_RANGENX16(_t3);

        y0 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_LO_HALVES);
        y1 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_LO_HALVES);
        y2 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_HI_HALVES);
        y3 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_HI_HALVES);
              
        y[0] = y0; 
        y[1] = y1; 
        y[2] = y2; 
        y[3] = y3; 
    }
    {   
        xb_vecNx16 y0, y1, y2, y3; 
        _t0 = x[1];                           
        _t1 = x[3];                           
        BBE_MOVSAV(x[5]);              
        BBE_MOVSBV(x[7]);                           
         
        LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_RANGENX16(_t0); 
        BBE_RANGENX16(_t1); 
        BBE_RANGENX16(_t2); 
        BBE_RANGENX16(_t3);

        y0 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_LO_HALVES);
        y1 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_LO_HALVES);
        y2 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_HI_HALVES);
        y3 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_HI_HALVES);
              
        y[0+4] = y0;  
        y[1+4] = y1;  
        y[2+4] = y2;  
        y[3+4] = y3;  
    }

    return scaling;
}

inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT4xIN_4_rm(xb_vecNx16 *x, int16_t *y,  int N, int bexp)
{

    int scaling = 0;
    int stride = N/4*2*sizeof(int16_t); 



    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
                                               
    VT t0, t1, t2, t3, t4, t5;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 0, bexp, scaling);


    ASSERT(count==2); 

  //  for (i=0; i<count; i++) 
     if(count==2)
    {   
        _t0 = x[0];       
        _t1 = x[2];       
        BBE_MOVSAV(x[4]); 
        BBE_MOVSBV(x[6]); 
         
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);                                                                                                     
                
        DFT4(t0, t1, t2, t3, t4, t5, A, B, 0);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SVRNX16_XP(_t0, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t1, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t2, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t3, p_dst, -3*stride + 2*BBE_SIMD_WIDTH); 
    }
    {   
        _t0 = x[1];                    
        _t1 = x[3];                    
        BBE_MOVSAV(x[5]);              
        BBE_MOVSBV(x[7]);                           
         
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                
        DFT4(t0, t1, t2, t3, t4, t5, A, B, 0);

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SVRNX16_XP(_t0, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t1, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t2, p_dst, stride);                                    
        BBE_SVRNX16_XP(_t3, p_dst, -3*stride + 2*BBE_SIMD_WIDTH); 
    }
    return scaling;
}
#if 1
inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT6xI4_N8(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    int scaling = 0;
    int stride = N/6*2*sizeof(int16_t); 

    int count = (N+6*BBE_SIMD_WIDTH/2-1)/(6*BBE_SIMD_WIDTH/2);                                        
    int i;        
    const xb_vecNx16 * p_tw = (xb_vecNx16 *)(tw);                                        
    const xb_vecNx16 * px1 = (xb_vecNx16 *)(x + 1*stride/2); 

    const xb_vecNx16 * px3 = (xb_vecNx16 *)(x + 3*stride/2); 
    const xb_vecNx16 * px4 = (xb_vecNx16 *)(x + 4*stride/2); 
    const xb_vecNx16 * px5 = (xb_vecNx16 *)(x + 5*stride/2); 

    NASSERT_ALIGN32(px4); 

    valign v1 = BBE_LA_PP(px1); 
    valign v3 = BBE_LA_PP(px3); 
    valign v5 = BBE_LA_PP(px5); 

    xb_vecNx16 * py0  = (xb_vecNx16 *)(y);
    xb_vecNx16 * py1  = (xb_vecNx16 *)(y+3*BBE_SIMD_WIDTH);

    valign uu0;                                                  
    xb_vecNx16 tw1, tw2, tw3, tw4, tw5;                         
    xb_vecNx16 _t0, _t1, _t2, _t3, _t4, _t5; 

    VT r3_tw; 

    BBE_MOVSBV(0);	
    BBE_MOVSDV(0);	
    RX3_TW(-1, r3_tw);  

    RANGE_BEGIN(6, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    ASSERT(N%8==0); 
    __Pragma("ymemory(px1)");
    __Pragma("ymemory(px3)");
    __Pragma("ymemory(px4)");
    __Pragma("ymemory(px5)");
    __Pragma("ymemory(p_tw)");  
   {   
        xb_vecNx16 y0; 

        BBE_LANX16_IP(_t5, v5, px5);
        BBE_LANX16_IP(_t3, v3, px3);                                    
        BBE_LANX16_IP(_t1, v1, px1); 

        BBE_LVNX16_XP(_t4, px4, -2*stride);   
        BBE_LVNX16_XP(_t2, px4, -2*stride);   
        BBE_LVNX16_XP(_t0, px4,  4*stride + 2*BBE_SIMD_WIDTH);                     
         
        BBE_LVNX16_IP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw1 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      
        tw2 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_1);      
        tw3 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_2);      
        tw4 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_3);

        BBE_LVNX16_IP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw5 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      


        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw); 
                                                                            
        _MUL(_t1, tw1);	
        _MUL(_t2, tw2);	
        _MUL(_t3, tw3);	
        _MUL(_t4, tw4);	
        _MUL(_t5, tw5);	

        uu0 = BBE_MOVUVR(_t0);                  
        BBE_SVINTLARNX16_XP(_t1, uu0, py0,    2*BBE_SIMD_WIDTH, 1);                
        BBE_SALIGNVRNX16_XP(_t2, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t3, uu0, py0,    2*BBE_SIMD_WIDTH, 1);
        BBE_SALIGNVRNX16_XP(_t4, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t5, uu0, py0,    4*2*BBE_SIMD_WIDTH, 1);

      }

   
    for (i=0; i<count-1; i++) 
    {   
        xb_vecNx16 y0; 
        BBE_LANX16_IP(_t5, v5, px5);
        BBE_LANX16_IP(_t3, v3, px3);                                    
        BBE_LANX16_IP(_t1, v1, px1); 

        BBE_LVNX16_XP(_t4, px4, -2*stride);   
        BBE_LVNX16_XP(_t2, px4, -2*stride);   
        BBE_LVNX16_XP(_t0, px4,  4*stride + 2*BBE_SIMD_WIDTH);                                    
 
        BBE_LVNX16_XP(y0, p_tw, 2*BBE_SIMD_WIDTH);

        tw1 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      
        tw2 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_1);      
        tw3 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_2);      
        tw4 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_3);

        BBE_LVNX16_XP(tw5, p_tw, 2*BBE_SIMD_WIDTH);         
        tw5 = BBE_SHFLNX16I(tw5,  BBE_SHFLI_REP_2X4_OFFSET_0);      
         
        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw);  
                                                                            
        _MUL(_t1, tw1);	
        _MUL(_t2, tw2);	
        _MUL(_t3, tw3);	
        _MUL(_t4, tw4);	
        _MUL(_t5, tw5);	

        BBE_SALIGNVRNX16_XP(_t0, uu0, py1,    4*2*BBE_SIMD_WIDTH);                  
        BBE_SVINTLARNX16_XP(_t1, uu0, py0,    2*BBE_SIMD_WIDTH, 1);                
        BBE_SALIGNVRNX16_XP(_t2, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t3, uu0, py0,    2*BBE_SIMD_WIDTH, 1);
        BBE_SALIGNVRNX16_XP(_t4, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t5, uu0, py0,    4*2*BBE_SIMD_WIDTH, 1);

      }
      BBE_SALIGNVRNX16_XP(_t5, uu0, py1,    2*BBE_SIMD_WIDTH); 

    return scaling;
}  

#else
inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT6xI4_N8(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/6*2*sizeof(int16_t); 

    int count = (N+6*BBE_SIMD_WIDTH/2-1)/(6*BBE_SIMD_WIDTH/2);                                        
    int i;        
    const xb_vecNx16 * p_tw = (xb_vecNx16 *)(tw);                                        
   // xb_vecNx16 * p_src = (xb_vecNx16 *)(x);  

    const xb_vecNx16 * px0 = (xb_vecNx16 *)(x + 0*stride/2); 
    const xb_vecNx16 * px1 = (xb_vecNx16 *)(x + 1*stride/2); 
    const xb_vecNx16 * px2 = (xb_vecNx16 *)(x + 2*stride/2); 
    const xb_vecNx16 * px3 = (xb_vecNx16 *)(x + 3*stride/2); 
    const xb_vecNx16 * px4 = (xb_vecNx16 *)(x + 4*stride/2); 
    const xb_vecNx16 * px5 = (xb_vecNx16 *)(x + 5*stride/2); 

    NASSERT_ALIGN32(px0); 
  //  NASSERT_ALIGN32(px1); 
    NASSERT_ALIGN32(px2); 
    //NASSERT_ALIGN32(px3); 
    NASSERT_ALIGN32(px4); 
    //NASSERT_ALIGN32(px5); 


    valign v1 = BBE_LA_PP(px1); 
    valign v3 = BBE_LA_PP(px3); 
    valign v5 = BBE_LA_PP(px5); 

  //  xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);  
    xb_vecNx16 * py0  = (xb_vecNx16 *)(y);
    xb_vecNx16 * py1  = (xb_vecNx16 *)(y+3*BBE_SIMD_WIDTH);

   // valign uu0;                                                  
    xb_vecNx16 tw1, tw2, tw3, tw4, tw5;                         
    xb_vecNx16 _t0, _t1, _t2, _t3, _t4, _t5; 
    
    //xb_vecNx16 _tw1, _tw2, _tw3, _tw4, _tw5;   
    VT r3_tw; 

   // valign align =  BBE_LANX16_PP(p_tw);

    BBE_MOVSBV(0);	
    BBE_MOVSDV(0);	
    RX3_TW(-1, r3_tw);  
    // RANGE_END(bexp);
    RANGE_BEGIN(6, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    ASSERT(N%8==0); 
#if 0
  
        
    __Pragma("ymemory(px0)");
    __Pragma("ymemory(px1)");
    __Pragma("ymemory(px2)");
    __Pragma("ymemory(px3)");
    __Pragma("ymemory(px4)");
    __Pragma("ymemory(px5)");
    __Pragma("ymemory(p_tw)");  
   {   
        xb_vecNx16 y0; 

        BBE_LVNX16_IP(_t0, px0, 2*BBE_SIMD_WIDTH);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LVNX16_IP(_t2, px2, 2*BBE_SIMD_WIDTH);      
        BBE_LANX16_IP(_t3, v3, px3);                                    
        BBE_LVNX16_IP(_t4, px4, 2*BBE_SIMD_WIDTH);      
        BBE_LANX16_IP(_t5, v5, px5);                        
         
        BBE_LVNX16_IP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw1 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      
        tw2 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_1);      
        tw3 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_2);      
        tw4 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_3);

        BBE_LVNX16_IP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw5 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      


       //__LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 
       // __LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5); 

        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw); 
                                                                            
        _MUL(_t1, tw1);	
        _MUL(_t2, tw2);	
        _MUL(_t3, tw3);	
        _MUL(_t4, tw4);	
        _MUL(_t5, tw5);	

        uu0 = BBE_MOVUVR(_t0);                  
        BBE_SVINTLARNX16_XP(_t1, uu0, py0,    2*BBE_SIMD_WIDTH, 1);                
        BBE_SALIGNVRNX16_XP(_t2, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t3, uu0, py0,    2*BBE_SIMD_WIDTH, 1);
        BBE_SALIGNVRNX16_XP(_t4, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t5, uu0, py0,    4*2*BBE_SIMD_WIDTH, 1);

      }

    for (i=0; i<count-1; i++) 
    {   
        xb_vecNx16 y0; 

        BBE_LVNX16_XP(_t0, px0, 2*BBE_SIMD_WIDTH);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LVNX16_XP(_t2, px2, 2*BBE_SIMD_WIDTH);      
        BBE_LANX16_IP(_t3, v3, px3);                                    
        BBE_LVNX16_XP(_t4, px4, 2*BBE_SIMD_WIDTH);      
        BBE_LANX16_IP(_t5, v5, px5);                        
         
        BBE_LVNX16_XP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw1 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      
        tw2 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_1);      
        tw3 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_2);      
        tw4 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_3);

        BBE_LVNX16_XP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw5 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      


       //__LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 
       // __LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5); 

        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw);  
                                                                            
        _MUL(_t1, tw1);	
        _MUL(_t2, tw2);	
        _MUL(_t3, tw3);	
        _MUL(_t4, tw4);	
        _MUL(_t5, tw5);	

        BBE_SALIGNVRNX16_XP(_t0, uu0, py1,    4*2*BBE_SIMD_WIDTH);                  
        BBE_SVINTLARNX16_XP(_t1, uu0, py0,    2*BBE_SIMD_WIDTH, 1);                
        BBE_SALIGNVRNX16_XP(_t2, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t3, uu0, py0,    2*BBE_SIMD_WIDTH, 1);
        BBE_SALIGNVRNX16_XP(_t4, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t5, uu0, py0,    4*2*BBE_SIMD_WIDTH, 1);

      }
      BBE_SALIGNVRNX16_XP(_t5, uu0, py1,    2*BBE_SIMD_WIDTH); 
#else
    __Pragma("ymemory(px0)");
    __Pragma("ymemory(px1)");
    __Pragma("ymemory(px2)");
    __Pragma("ymemory(px3)");
    __Pragma("ymemory(px4)");
    __Pragma("ymemory(px5)");
    __Pragma("ymemory(p_tw)");  

    for (i=0; i<count; i++) 
    {   
        xb_vecNx16 y0, y1, y2; 

        BBE_LVNX16_IP(_t0, px0, 2*BBE_SIMD_WIDTH);                                    
        BBE_LANX16_IP(_t1, v1, px1);                                    
        BBE_LVNX16_IP(_t2, px2, 2*BBE_SIMD_WIDTH);      
        BBE_LANX16_IP(_t3, v3, px3);                                    
        BBE_LVNX16_IP(_t4, px4, 2*BBE_SIMD_WIDTH);      
        BBE_LANX16_IP(_t5, v5, px5);                        
         
        BBE_LVNX16_XP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw1 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      
        tw2 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_1);      
        tw3 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_2);      
        tw4 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_3);

        BBE_LVNX16_XP(y0, p_tw, 2*BBE_SIMD_WIDTH);
        tw5 = BBE_SHFLNX16I(y0,  BBE_SHFLI_REP_2X4_OFFSET_0);      


       //__LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 
       // __LDA_IX4_TW_UNPACK_2(align, p_tw, tw4, tw5); 

        __DFT6(_t0, _t1, _t2, _t3, _t4, _t5, r3_tw); 
                                                                            
        _MUL(_t1, tw1);	
        _MUL(_t2, tw2);	
        _MUL(_t3, tw3);	
        _MUL(_t4, tw4);	
        _MUL(_t5, tw5);	

        y0 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_LO_HALVES);		
        y1 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_LO_HALVES);
        y2 = BBE_SELNX16I(_t5, _t4, BBE_SELI_EXTRACT_LO_HALVES);

        BBE_SVRNX16_XP(y0, py0, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y1, py0, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y2, py0, 4*2*BBE_SIMD_WIDTH); 

        y0 = BBE_SELNX16I(_t1, _t0, BBE_SELI_EXTRACT_HI_HALVES);		
        y1 = BBE_SELNX16I(_t3, _t2, BBE_SELI_EXTRACT_HI_HALVES);
        y2 = BBE_SELNX16I(_t5, _t4, BBE_SELI_EXTRACT_HI_HALVES);

        BBE_SVRNX16_XP(y0, py1, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y1, py1, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y2, py1, 4*2*BBE_SIMD_WIDTH);

      }


#endif
    return scaling;
} // R2_DFT6xI4_N8 Old variant
#endif

inline_  ATTRIBUTE_ALWAYS_INLINE  int R2_DFT9xIN_9(    int16_t *x,/*input*/ 
                               int16_t *y/*output*/, 
                               const int N, 
                               int bexp)
{
    

    const int v = 8;       // vector size
    const int N1 = 9;
    
    const int stride_bytes = 2*2*N/N1;
    int i; 

    const int   num_bfls = (N+N1*v-1)/(N1*v);

    int shift;

    VT*  px0 = ( VT*)(x + 8*stride_bytes/2);
    VT*  py0 = ( VT*)(y );
 
    VT r3_tw, 
        t11, t12, t22; // twiddles for radix 9 bfl 
    RX3_TW(-1, r3_tw); 
    RX9_TW(t11, t12, t22); 

    
    NASSERT_ALIGN32(x); 
    NASSERT_ALIGN32(y); 
    BBE_MOVSBV(0);	
    BBE_MOVSDV(0);
    BBE_MOVSCV(0);	
    BBE_MOVSAV(0);
    //RANGE_BEGIN 
    shift = 4-bexp; 
    shift = (shift<0)? 0: shift;

    {
        int mode; 
        if(shift>3)
            mode = 16+4+shift-1;
        else
            mode = 16+shift;

       
        BBE_FFTWMODE(mode); 
    }
   // __Pragma("ymemory(px0)");
    __Pragma("ymemory(py0)");
    __Pragma("no_reorder");  
    for (i = 0; i < num_bfls; i++) 
    {
         VT   A0, A1, A2, B0, B1, B2, C0, C1, C2; 
           // radix 3 bfl 3/0
        BBE_LVN_2XCQ15_XP(C2, px0, -3*stride_bytes/*8*stride_bytes*/ );
        BBE_LVN_2XCQ15_XP(C1, px0, -3*stride_bytes/*5*stride_bytes*/ );
        BBE_LVC_IP(px0, 0); //C0 = BBE_LVN_2XCQ15_X(px0, 2*stride_bytes);
        BBE_LVA_XP(px0, 5*stride_bytes);

        //BBE_MOVSAV(BBE_MOVNX16_FROMN_2XCQ15(C0)); BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(C0));		
        DFT3(C0, C1, C2, r3_tw, 0);  

        // radix 3 bfl 2/0
        BBE_LVN_2XCQ15_XP(B2, px0, -3*stride_bytes/*7*stride_bytes*/);
        BBE_LVN_2XCQ15_XP(B1, px0, -3*stride_bytes/*4*stride_bytes*/);
        BBE_LVC_IP(px0, 0); //C0 = BBE_LVN_2XCQ15_X(px0, 2*stride_bytes);
        BBE_LVA_XP(px0, 5*stride_bytes/*2*BBE_SIMD_WIDTH - 2*stride_bytes*/);
        //BBE_MOVSAV(BBE_MOVNX16_FROMN_2XCQ15(B0)); BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(B0));		
        DFT3(B0, B1, B2, r3_tw, 0);  

        // radix 3 bfl 1/0
        BBE_LVN_2XCQ15_XP(A2, px0, -3*stride_bytes/*6*stride_bytes*/);
        BBE_LVN_2XCQ15_XP(A1, px0, -3*stride_bytes/*3*stride_bytes*/);
        BBE_LVC_IP(px0, 0); //B0 = BBE_LVN_2XCQ15_X(px0, 1*stride_bytes);
        BBE_LVA_XP(px0, 8*stride_bytes+BBE_SIMD_WIDTH*2); 
        //BBE_MOVSAV(BBE_MOVNX16_FROMN_2XCQ15(A0)); BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(A0));		
        DFT3(A0, A1, A2, r3_tw, 0);  

        

        MUL(B1, t11);
        MUL(B2, t12);
        MUL(C1, t12);
        MUL(C2, t22);

        // radix 3 bfl 1/1
        BBE_MOVSAV(BBE_MOVNX16_FROMN_2XCQ15(A0)); BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(A0));		
        DFT3(A0, B0, C0, r3_tw, 1);  
        BBE_SVN_2XCQ15_I(A0, py0, 0);
        BBE_SVN_2XCQ15_X(B0, py0, 3*stride_bytes);
        BBE_SVN_2XCQ15_X(C0, py0, 6*stride_bytes);

  
        // radix 3 bfl 2/1
        BBE_MOVSAV( BBE_MOVNX16_FROMN_2XCQ15(A1)); BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(A1));		
        DFT3(A1, B1, C1, r3_tw, 1);  
        BBE_SVN_2XCQ15_X(A1, py0, 1*stride_bytes);
        BBE_SVN_2XCQ15_X(B1, py0, 4*stride_bytes);
        BBE_SVN_2XCQ15_X(C1, py0, 7*stride_bytes);

        // radix 3 bfl 2/1
        BBE_MOVSAV(BBE_MOVNX16_FROMN_2XCQ15(A2)); BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(A2));		
        DFT3(A2, B2, C2, r3_tw, 1);  //x0  = A2;  x1 = B2; x2 = C2; 
        BBE_SVN_2XCQ15_X(A2, py0, 2*stride_bytes);
        BBE_SVN_2XCQ15_X(B2, py0, 5*stride_bytes);
        BBE_SVN_2XCQ15_X(C2, py0, 8*stride_bytes);

        py0 ++; 
    //    px0 ++; 
    }

    return shift;
} //R2_DFT9xIN_9


#if 1
// Interleaving based on BBE_SALIGN 
inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT5xI4_expand56(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       int bexp
                                       )
{

    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    int i;

    int num_bfls = (N+N1*BBE_SIMD_WIDTH/2-1)/(N1*BBE_SIMD_WIDTH/2);
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    
    //xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    //xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15 tw1, tw2, tw3, tw4, r5_tw1, r5_tw2, r5_tw3; 
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;
    xb_vecNx16 * py0  = (xb_vecNx16 *)(y);
    xb_vecNx16 * py1  = (xb_vecNx16 *)(y+3*BBE_SIMD_WIDTH);

  
    valign v1 = BBE_LAN_2XCQ15_PP(px1), 
           v3 = BBE_LAN_2XCQ15_PP(px3); 
    valign uu0;    
    int scaling = 0;

    RANGE_BEGIN(N1, -1, 1, bexp, scaling);
    __Pragma("no_reorder");
    

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);

    //NASSERT_ALIGN32(px0); 
    //NASSERT_ALIGN32(px1); 
    //NASSERT_ALIGN32(px2); 
   // NASSERT_ALIGN32(px3); 
    NASSERT_ALIGN32(px4); 

    __Pragma("ymemory(py0)"); 
    __Pragma("ymemory(py1)"); 

    {
        VT tmp_; 
        xb_vecNx16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5 = 0 ; 
       // xb_vecNx16 y0, y1, y2; 

        BBE_LVNX16_IP(tmp3, tw, 2*BBE_SIMD_WIDTH); 

        tmp0 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tmp1 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_1);
        tmp2 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_2);
        tmp3 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_3);
        
        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
        tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
        tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 

        BBE_LAN_2XCQ15_IP(t1, v1, px1);            
        BBE_LAN_2XCQ15_IP(tmp_, v3, px3);          
        BBE_MOVSBV(BBE_MOVNX16_FROMN_2XCQ15(tmp_));

        BBE_LVA_XP(px4, -4*stride); 
        BBE_LVN_2XCQ15_XP(t2, px4, -4*stride);
        BBE_LVN_2XCQ15_XP(t0, px4, 8*stride + 2*BBE_SIMD_WIDTH);

        t0 = BBE_FFTSRAN_2XCQ15(t0);
        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        
        tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
        tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

        BBE_MOVSCV(tmp0); //  MV ## C (t5);							
        BBE_MOVSDV(tmp1); //  MV ## D (t6);							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								
        MUL(t4, tw1);								
        MUL(t5, tw2);								
        MUL(t7, tw3);								
        MUL(t6, tw4);	

        tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t1);
        tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t4);
        tmp2 = BBE_MOVNX16_FROMN_2XCQ15(t5);
        tmp3 = BBE_MOVNX16_FROMN_2XCQ15(t7);
        tmp4 = BBE_MOVNX16_FROMN_2XCQ15(t6);
        
         uu0 = BBE_MOVUVR(tmp0);                  
        BBE_SVINTLARNX16_XP(tmp1, uu0, py0,    2*BBE_SIMD_WIDTH, 1);                
        BBE_SALIGNVRNX16_XP(tmp2, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(tmp3, uu0, py0,    2*BBE_SIMD_WIDTH, 1);
        BBE_SALIGNVRNX16_XP(tmp4, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(tmp5, uu0, py0,    4*2*BBE_SIMD_WIDTH, 1);

      
    } 
    
    for(i=0; i<num_bfls-1; i++)
    {
        // 15 Cycles
        VT tmp_; 
        xb_vecNx16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5 = 0 ; 
       // xb_vecNx16 y0, y1, y2; 

        BBE_LVNX16_IP(tmp3, tw, 2*BBE_SIMD_WIDTH); 

        tmp0 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tmp1 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_1);
        tmp2 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_2);
        tmp3 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_3);
        
        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
        tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
        tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 

        BBE_LAN_2XCQ15_IP(t1, v1, px1);            
        BBE_LAN_2XCQ15_IP(tmp_, v3, px3);          
        BBE_MOVSBV(BBE_MOVNX16_FROMN_2XCQ15(tmp_));

        BBE_LVA_XP(px4, -4*stride); 
        BBE_LVN_2XCQ15_XP(t2, px4, -4*stride);
        BBE_LVN_2XCQ15_XP(t0, px4, 8*stride + 2*BBE_SIMD_WIDTH);
 
        t0 = BBE_FFTSRAN_2XCQ15(t0);
        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        
        tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
        tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

        BBE_MOVSCV(tmp0); //  MV ## C (t5);							
        BBE_MOVSDV(tmp1); //  MV ## D (t6);							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								
        MUL(t4, tw1);								
        MUL(t5, tw2);								
        MUL(t7, tw3);								
        MUL(t6, tw4);	

        tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t1);
        tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t4);
        tmp2 = BBE_MOVNX16_FROMN_2XCQ15(t5);
        tmp3 = BBE_MOVNX16_FROMN_2XCQ15(t7);
        tmp4 = BBE_MOVNX16_FROMN_2XCQ15(t6);
        
        BBE_SALIGNVRNX16_XP(tmp0, uu0, py1,    4*2*BBE_SIMD_WIDTH);                  
        BBE_SVINTLARNX16_XP(tmp1, uu0, py0,    2*BBE_SIMD_WIDTH, 1);                
        BBE_SALIGNVRNX16_XP(tmp2, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(tmp3, uu0, py0,    2*BBE_SIMD_WIDTH, 1);
        BBE_SALIGNVRNX16_XP(tmp4, uu0, py1,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(tmp5, uu0, py0,    4*2*BBE_SIMD_WIDTH, 1);
     
    } 
    BBE_SALIGNVRNX16_XP(0, uu0, py1,    2*BBE_SIMD_WIDTH); 

    return scaling; 
}
#else
// Interleaved based on BBE_SEL
inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT5xI4_expand56(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       int bexp
                                       )
{

    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    int i;



    int num_bfls = (N+N1*BBE_SIMD_WIDTH/2-1)/(N1*BBE_SIMD_WIDTH/2);
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15 tw1, tw2, tw3, tw4, r5_tw1, r5_tw2, r5_tw3; 
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;
    xb_vecNx16 *py = (xb_vecNx16 *)(y);
  
    valign v1 = BBE_LAN_2XCQ15_PP(px1), 
           v3 = BBE_LAN_2XCQ15_PP(px3); 

    int scaling = 0;

    RANGE_BEGIN(N1, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);

    NASSERT_ALIGN32(px0); 
    //NASSERT_ALIGN32(px1); 
    NASSERT_ALIGN32(px2); 
   // NASSERT_ALIGN32(px3); 
    NASSERT_ALIGN32(px4); 
    // 17 Cycles
    for(i=0; i<num_bfls; i++)
    {
        VT tmp_; 
        xb_vecNx16 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5 = 0 ; 
        xb_vecNx16 y0, y1, y2; 

        BBE_LVNX16_IP(tmp3, tw, 2*BBE_SIMD_WIDTH); 

        tmp0 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tmp1 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_1);
        tmp2 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_2);
        tmp3 = BBE_SHFLNX16I(tmp3,  BBE_SHFLI_REP_2X4_OFFSET_3);
        
        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
        tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
        tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 

        BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
        BBE_LAN_2XCQ15_IP(t1, v1, px1);             //BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
        BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));

        BBE_LAN_2XCQ15_IP(tmp_, v3, px3);           //BBE_LVB_IP(px3, sizeof(*px1)); 
        BBE_MOVSBV(BBE_MOVNX16_FROMN_2XCQ15(tmp_)); 
        BBE_LVA_IP(px4, sizeof(*px1)); 

        t0 = BBE_FFTSRAN_2XCQ15(t0);
        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        
        tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
        tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

        BBE_MOVSCV(tmp0); //  MV ## C (t5);							
        BBE_MOVSDV(tmp1); //  MV ## D (t6);							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								
        MUL(t4, tw1);								
        MUL(t5, tw2);								
        MUL(t7, tw3);								
        MUL(t6, tw4);	

        tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t1);
        tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t4);
        tmp2 = BBE_MOVNX16_FROMN_2XCQ15(t5);
        tmp3 = BBE_MOVNX16_FROMN_2XCQ15(t7);
        tmp4 = BBE_MOVNX16_FROMN_2XCQ15(t6);
        
        y0 = BBE_SELNX16I(tmp1, tmp0, BBE_SELI_EXTRACT_LO_HALVES);		
        y1 = BBE_SELNX16I(tmp3, tmp2, BBE_SELI_EXTRACT_LO_HALVES);
        y2 = BBE_SELNX16I(tmp5, tmp4, BBE_SELI_EXTRACT_LO_HALVES);

        BBE_SVRNX16_XP(y0, py, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y1, py, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y2, py, 2*BBE_SIMD_WIDTH); 

        y0 = BBE_SELNX16I(tmp1, tmp0, BBE_SELI_EXTRACT_HI_HALVES);		
        y1 = BBE_SELNX16I(tmp3, tmp2, BBE_SELI_EXTRACT_HI_HALVES);
        y2 = BBE_SELNX16I(tmp5, tmp4, BBE_SELI_EXTRACT_HI_HALVES);

        BBE_SVRNX16_XP(y0, py, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y1, py, 2*BBE_SIMD_WIDTH); 
        BBE_SVRNX16_XP(y2, py, 2*BBE_SIMD_WIDTH);     
    } 

    return scaling; 
}
#endif

inline_ int ATTRIBUTE_ALWAYS_INLINE R2_DFT3xIN_3_compress65( 
                                       int16_t *x, 
                                       int16_t *y, 
                                       const int N, int bexp )
{
    const int v = 8;       // vector size
    const int N1 = 3;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_out = 2*N*5/6/N1; /*in int16_t */
    int scaling = 0; 
    int i; 

    int num_bfls = N/N1/v;
   // int   num_frac_bfls = N/N1 - v * num_bfls;
    xb_vecN_2xcq15 y0, y1, y2; 

    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);

    xb_vecN_2xcq15 *py0 = (xb_vecN_2xcq15 *)(y);
    xb_vecN_2xcq15 *py1 = (xb_vecN_2xcq15 *)(y + 1*stride_out);
    xb_vecN_2xcq15 *py2 = (xb_vecN_2xcq15 *)(y + 2*stride_out);

    xb_vecN_2xcq15 r3_tw;
 
    valign a0 =  BBE_ZALIGN();
    valign a1 =  BBE_ZALIGN();
    valign a2 =  BBE_ZALIGN();
  
    RX3_TW(-1, r3_tw);

    BBE_MOVSBV(0); 
    BBE_MOVSDV(0); 

    // RANGE_END(bexp);
    RANGE_BEGIN(3, -1, 0, bexp, scaling);

    ASSERT((N/N1 - v * num_bfls)==0); 
    ASSERT(num_bfls%3==0); 

    for(i=0; i<num_bfls/3; i++)
    {
    								
        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));  
        BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1) );
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2) ); 

        DFT3(y0, y1, y2, r3_tw, 0);	

        BBE_SAVRN_2XCQ15_XP(y0, a0, py0, 8*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y1, a1, py1, 8*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y2, a2, py2, 8*2*sizeof(int16_t)); 

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));  
        BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1) );
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2) ); 

        DFT3(y0, y1, y2, r3_tw, 0);	

        BBE_SAVRN_2XCQ15_XP(y0, a0, py0, 8*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y1, a1, py1, 8*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y2, a2, py2, 8*2*sizeof(int16_t)); 

        BBE_LVA_IP(px0, 0); 
        BBE_LVC_IP(px0,  sizeof(*px0));  
        BBE_LVN_2XCQ15_IP(y1, px1, sizeof(*px1) );
        BBE_LVN_2XCQ15_IP(y2, px2, sizeof(*px2) ); 

        DFT3(y0, y1, y2, r3_tw, 0);	

        BBE_SAVRN_2XCQ15_XP(y0, a0, py0, 4*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y1, a1, py1, 4*2*sizeof(int16_t)); 
        BBE_SAVRN_2XCQ15_XP(y2, a2, py2, 4*2*sizeof(int16_t)); 
    }

    BBE_SAN_2XCQ15POS_FC(a0, py0);	
    BBE_SAN_2XCQ15POS_FC(a1, py1);
    BBE_SAN_2XCQ15POS_FC(a2, py2);

    return scaling; 
}


inline_ ATTRIBUTE_ALWAYS_INLINE int R3_DFT5xI3xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v, /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       int bexp
                                       )
{

    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    int /*i,*/ j; 
//    unsigned int tmp;

    int num_bfls = N/N1/v;
//    xb_vecN_2xcq15 x0, x1, x2, t0, t1, t2; 
//    xb_vecN_2xcq15 y0, y1, y2; 
    const xb_vecNx16 *tw = (xb_vecNx16*)ptw;
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15 tw1, tw2, tw3, tw4, r5_tw1, r5_tw2, r5_tw3; 
     xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;

    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    
    valign vtw = BBE_LAVNX16_PP(tw); 

    // int bexp; 
    int scaling = 0;
    // RANGE_END(bexp);
    RANGE_BEGIN(N1, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 
    ASSERT(num_bfls==3);(void)num_bfls;

    
    {
        xb_vecNx16 tmp0, tmp1,  tmp3; 
        BBE_LAVNX16_XP(tmp3, vtw, tw, 4*4);
        //  Iteration 1 twidles always equal to ones
        for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
        {
            BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
            BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
            BBE_LVB_IP(px3, sizeof(*px1)); 
            BBE_LVA_IP(px4, sizeof(*px1)); 


            t0 = BBE_FFTSRAN_2XCQ15(t0);

            /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
            DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
            								
            t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
            {  
                xb_vecN_2xcq15 tmpc;
                tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
            } 

            MUL(t4, r5_tw1);							
            MUL(t5, r5_tw2);							
            MUL(t6, r5_tw3);							
          
            tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
            tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

            BBE_MOVSCV(tmp0); //  MV ## C (t5);							
            BBE_MOVSDV(tmp1); //  MV ## D (t6);							
        									
            t3 = t4;								
            								
            DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
/*            								
            MUL(t4, tw1);								
            MUL(t5, tw2);								
            MUL(t7, tw3);								
            MUL(t6, tw4);								
*/
            BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

        }
        py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
    } 
    {
        xb_vecNx16 tmp0, tmp1, tmp2, tmp3; 

        BBE_LAVNX16_XP(tmp3, vtw, tw, 4*4);
        tmp0 = BBE_REPNX16C(tmp3, 0); 
        tmp1 = BBE_REPNX16C(tmp3, 1); 
        tmp2 = BBE_REPNX16C(tmp3, 2); 
        tmp3 = BBE_REPNX16C(tmp3, 3); 

        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
        tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
        tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 

        //  Iteration 2
        for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
        {
            BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
            BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
            BBE_LVB_IP(px3, sizeof(*px1)); 
            BBE_LVA_IP(px4, sizeof(*px1)); 


            t0 = BBE_FFTSRAN_2XCQ15(t0);

            /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
            DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
            								
            t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
            {  
                xb_vecN_2xcq15 tmpc;
                tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
            }                                    
            MUL(t4, r5_tw1);							
            MUL(t5, r5_tw2);							
            MUL(t6, r5_tw3);							
            
            tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
            tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

            BBE_MOVSCV(tmp0); //  MV ## C (t5);							
            BBE_MOVSDV(tmp1); //  MV ## D (t6);							
        									
            t3 = t4;								
            								
            DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
            								
            MUL(t4, tw1);								
            MUL(t5, tw2);								
            MUL(t7, tw3);								
            MUL(t6, tw4);								

            BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

        }
        py += (N1-1)*2*v*sizeof(int16_t)/sizeof(*py); 
    } 
    {
        xb_vecNx16 tmp0, tmp1, tmp2, tmp3; 

        BBE_LAVNX16_XP(tmp3, vtw, tw, 4*4);
        tmp0 = BBE_REPNX16C(tmp3, 0); 
        tmp1 = BBE_REPNX16C(tmp3, 1); 
        tmp2 = BBE_REPNX16C(tmp3, 2); 
        tmp3 = BBE_REPNX16C(tmp3, 3); 

        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1); 
        tw3 = BBE_MOVN_2XCQ15_FROMNX16(tmp2); 
        tw4 = BBE_MOVN_2XCQ15_FROMNX16(tmp3); 
        //  Iteration 3
        for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
        {
            BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
            BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
            BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
            BBE_LVB_IP(px3, sizeof(*px1)); 
            BBE_LVA_IP(px4, sizeof(*px1)); 


            t0 = BBE_FFTSRAN_2XCQ15(t0);

            /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
            DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
            								
            t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
            {  
                xb_vecN_2xcq15 tmpc;
                tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
                t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
            }                                    
            MUL(t4, r5_tw1);							
            MUL(t5, r5_tw2);							
            MUL(t6, r5_tw3);							
            
            tmp0 = BBE_MOVNX16_FROMN_2XCQ15(t5); 
            tmp1 = BBE_MOVNX16_FROMN_2XCQ15(t6); 

            BBE_MOVSCV(tmp0); //  MV ## C (t5);							
            BBE_MOVSDV(tmp1); //  MV ## D (t6);							
        									
            t3 = t4;								
            								
            DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
            								
            MUL(t4, tw1);								
            MUL(t5, tw2);								
            MUL(t7, tw3);								
            MUL(t6, tw4);								

            BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
            BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

        }
    } 

    return scaling; 
} // R3_DFT5xI3xIv

inline_ ATTRIBUTE_ALWAYS_INLINE int R3_I3xDFT3xIv(     const int16_t *ptw,   
                                       int16_t *x   /*input*/, 
                                       int16_t *y   /*output*/, 
                                       const int N, 
                                       const int v, /*vector length must be multiple of BBE_SIMD_WIDTH/2 */
                                       int bexp
                                       )
{

    const int N1 = 3;
    const int stride = 2*N/N1; /*in int16_t */
    const int stride_bytes = 2*stride; 
    int i, j; 
//    unsigned int tmp;

    int num_bfls = N/N1/v;
//    xb_vecN_2xcq15 x0, x1, x2, t0, t1, t2; 
    xb_vecN_2xcq15 y0, y1, y2; 
    const xb_vecNx16 *tw = (xb_vecNx16*)(ptw+4);
    
    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*2*v + 2*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*2*v + 2*stride);
    xb_vecN_2xcq15 r3_tw; 

    xb_vecN_2xcq15 *py0 = (xb_vecN_2xcq15 *)(y);
    xb_vecN_2xcq15 *py1 = (xb_vecN_2xcq15 *)(y + 1*2*v*N1);    
    xb_vecN_2xcq15 *py2 = (xb_vecN_2xcq15 *)(y + 2*2*v*N1);    


    // int bexp; 
    int scaling = 0;
    // RANGE_END(bexp);
    RANGE_BEGIN(N1, -1, 1, bexp, scaling);

    BBE_MOVSBV(0); 
    BBE_MOVSDV(0); 

    __Pragma("no_reorder");

    RX3_TW(-1, r3_tw);	

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 

//    ASSERT(v==48);
    ASSERT(num_bfls==3);

  
    for(i=0; i<num_bfls/3; i++)
    {
        xb_vecNx16 tmp0, tmp1;
        xb_vecN_2xcq15 tw11, tw21, tw12, tw22; 
#if 0
        BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
        BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));

        tmp0 = BBE_REPNX16C(tmp0, 0); 
        tmp1 = BBE_REPNX16C(tmp1, 0); 

        tw1 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw2 = BBE_MOVN_2XCQ15_FROMNX16(tmp1);
#endif
        BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
        BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));

        tmp0 = BBE_REPNX16C(tmp0, 0); 
        tmp1 = BBE_REPNX16C(tmp1, 0); 

        tw11 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw21 = BBE_MOVN_2XCQ15_FROMNX16(tmp1);
      
        BBE_LPNX16_IP(tmp0, tw, 2*sizeof(int16_t));
        BBE_LPNX16_IP(tmp1, tw, 2*sizeof(int16_t));

        tmp0 = BBE_REPNX16C(tmp0, 0); 
        tmp1 = BBE_REPNX16C(tmp1, 0); 

        tw12 = BBE_MOVN_2XCQ15_FROMNX16(tmp0); 
        tw22 = BBE_MOVN_2XCQ15_FROMNX16(tmp1);

        for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
        {
            BBE_LVN_2XCQ15_XP(y2, px0, -stride_bytes); 
            BBE_LVN_2XCQ15_XP(y1, px0, -stride_bytes);
            
            BBE_LVC_IP(px0, 0); 
            BBE_LVA_XP(px0,  2*stride_bytes+2*BBE_SIMD_WIDTH);

            DFT3(y0, y1, y2, r3_tw, 0);													
          //  MUL(y1, tw1);						
          //  MUL(y2, tw2);

            BBE_SVRN_2XCQ15_XP(y0, py0, 2*v*2);
            BBE_SVRN_2XCQ15_XP(y1, py0, 2*v*2); 
            BBE_SVRN_2XCQ15_XP(y2, py0, 2*BBE_SIMD_WIDTH - 4*v*2); 
/////////////////////////////
            BBE_LVN_2XCQ15_XP(y2, px1, -stride_bytes); 
            BBE_LVN_2XCQ15_XP(y1, px1, -stride_bytes);
            
            BBE_LVC_IP(px1, 0); 
            BBE_LVA_XP(px1,  2*stride_bytes+2*BBE_SIMD_WIDTH);

            DFT3(y0, y1, y2, r3_tw, 0);													
            MUL(y1, tw11);						
            MUL(y2, tw21);

            BBE_SVRN_2XCQ15_XP(y0, py1, 2*v*2);
            BBE_SVRN_2XCQ15_XP(y1, py1, 2*v*2); 
            BBE_SVRN_2XCQ15_XP(y2, py1, 2*BBE_SIMD_WIDTH - 4*v*2); 
////////////
            BBE_LVN_2XCQ15_XP(y2, px2, -stride_bytes); 
            BBE_LVN_2XCQ15_XP(y1, px2, -stride_bytes);

            BBE_LVC_IP(px2, 0); 
            BBE_LVA_XP(px2,      2*stride_bytes+2*BBE_SIMD_WIDTH);                        

            DFT3(y0, y1, y2, r3_tw, 0);			 										
            MUL(y1, tw12);						
            MUL(y2, tw22);
   
            BBE_SVRN_2XCQ15_XP(y0, py2, 2*v*2);
            BBE_SVRN_2XCQ15_XP(y1, py2, 2*v*2); 
            BBE_SVRN_2XCQ15_XP(y2, py2, 2*BBE_SIMD_WIDTH - 4*v*2);
        }

        py0 += 2*(N1)*2*v*sizeof(int16_t)/sizeof(*py0); 
        py1 += 2*(N1)*2*v*sizeof(int16_t)/sizeof(*py0); 
        py2 += 2*(N1)*2*v*sizeof(int16_t)/sizeof(*py0); 
    }
    return scaling; 
}


inline_ ATTRIBUTE_ALWAYS_INLINE int R1_tDFT4_L64_16(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    int stride = N/4*2*sizeof(int16_t); 
    xb_vecNx16 * p_tw = (xb_vecNx16 *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3;  
    valign v = BBE_LAVNX16_PP(p_tw); 
    VT r1;//, r2, r3; 
 
    int scaling = bexp;
    vsaN shift; 
  
    bexp = (bexp>3)? 3: bexp;
    scaling = bexp-scaling;

    shift = BBE_MOVVSA32(scaling); 

    BBE_LAVNX16_XP(_t1, v, p_tw, 2*BBE_SIMD_WIDTH); 

    r1 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_LAVNX16_XP(_t3, v, p_tw, 1*4); 
    tw1 = r1; 
    tw2 = r1;
    tw2 = BBE_MULN_2XCQ15PACKQ(tw2, r1);

    tw3 = tw2; 
    tw3 = BBE_MULN_2XCQ15PACKQ(tw3, r1);

    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    
    __Pragma("ymemory(p_dst)"); 

    BBE_LVNX16_XP(_t0, p_src, stride);                                    
    BBE_LVNX16_XP(_t1, p_src, stride);                                    
    BBE_LVNX16_XP(_t2, p_src, stride);                                    
    BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
                                                                    
    _t0 = BBE_SRANX16(_t0, shift);                                   
    _t1 = BBE_SRANX16(_t1, shift);                                   
    _t2 = BBE_SRANX16(_t2, shift);                                   
    _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                    
    t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
    t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

    BBE_MOVSAV(_t2);                                                  
    BBE_MOVSBV(_t3);                                                            

    __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                  
    BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2); //INTLV(t0, t2);                                                
    BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2); //INTLV(t1, t3);         
    _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
    _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
    _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);  
    uu0 = BBE_MOVUVR(_t0);                             
    uu1 = BBE_MOVUVR(_t2);                             
    BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*2*BBE_SIMD_WIDTH, 0);     
    BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  -2*BBE_SIMD_WIDTH, 0);     
  
    for (i=0; i<count-1; i++) 
    {

        BBE_LAVNX16_XP(_t3, v, p_tw,  1*4); 

        _t1 = BBE_REPNX16C(_t3, 0); 
        tw1 = BBE_MOVN_2XCQ15_FROMNX16(_t1);    
        tw1 = BBE_MULN_2XCQ15PACKQ(tw1, r1);
        tw2 = tw1; 
        tw2 = BBE_MULN_2XCQ15PACKQ(tw2, tw2);
        tw3 = tw2; 
        tw3 = BBE_MULN_2XCQ15PACKQ(tw3, tw1);

        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t0, p_src, stride);    
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t1, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t2, p_src, stride);
        NASSERT_ALIGN32(p_src);
        BBE_LVNX16_XP(_t3, p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
 
        _t0 = BBE_SRANX16(_t0, shift);                                   
        _t1 = BBE_SRANX16(_t1, shift);                                   
        _t2 = BBE_SRANX16(_t2, shift);                                   
        _t3 = BBE_SRANX16(_t3, shift);                                   
                                                                        
        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  

        BBE_MOVSAV(_t2);                                                  
        BBE_MOVSBV(_t3);                                                  
              
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                      
        BBE_DSELN_2XCQ15I(t2, t0, t2, t0, BBE_DSELI_INTERLEAVE_2);                                                 
        BBE_DSELN_2XCQ15I(t3, t1, t3, t1, BBE_DSELI_INTERLEAVE_2);          

        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                           
                                                                      
        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst,    2*2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst,    2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,    2*2*BBE_SIMD_WIDTH, 0);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,    -2*BBE_SIMD_WIDTH, 0);
    }

    _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);                              
    BBE_SALIGNVRNX16_XP(_t3, uu0, p_dst, 2*2*BBE_SIMD_WIDTH);                      
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*2*BBE_SIMD_WIDTH);

    return scaling; 
}


inline_ ATTRIBUTE_ALWAYS_INLINE int R2_tDFT4xI4_U2(const int16_t *tw, int16_t *x, int16_t *y,  int N, int bexp)
{
    // int bexp; 
    int scaling = 0;
    int stride = N/4*2*sizeof(int16_t); 

    int count = N/4/(BBE_SIMD_WIDTH/2);                                        
    int i;        
    VT * p_tw = (VT *)(tw);                                        
    xb_vecNx16 * p_src = (xb_vecNx16 *)(x);                                       
    xb_vecNx16 * p_dst = (xb_vecNx16 *)(y);                      
    valign uu0, uu1;                                                  
    VT t0, t1, t2, t3, t4, t5,  tw1, tw2, tw3;                         
    xb_vecNx16 _t0, _t1, _t2, _t3; 
   
    valign align =  BBE_LAN_2XCQ15_PP(p_tw);

    ASSERT((count&1)==0); 

    // RANGE_END(bexp);
    RANGE_BEGIN(4, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

   {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
    //    LDA_IX4_TW_UNPACK_3(align, p_tw, tw1, tw2, tw3); 
        BBE_LAVN_2XCQ15_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLN_2XCQ15I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tw2 = BBE_MULN_2XCQ15PACKQ(tw1, tw1); 
        tw3 = BBE_MULN_2XCQ15PACKQ(tw2, tw1);

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                            
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);

        uu0 = BBE_MOVUVR(_t0);                             
        uu1 = BBE_MOVUVR(_t2);                             
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst,  2*BBE_SIMD_WIDTH, 1);     
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst,  2*BBE_SIMD_WIDTH, 1); 

    }

    for (i=0; i < (count-2)/2; i++) 
    {   
        xb_vecN_2xcq15 tmp; 
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           

        BBE_LAVN_2XCQ15_XP(tmp, align, p_tw, 4*4);                
        tw1 = BBE_SHFLN_2XCQ15I(tmp,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tw2 = BBE_MULN_2XCQ15PACKQ(tw1, tw1); 
        tw3 = BBE_MULN_2XCQ15PACKQ(tw2, tw1);


        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);

        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVC_XP( p_src, stride);                                    
        BBE_LVD_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
             
        tw1 = BBE_SHFLN_2XCQ15I(tmp,  BBE_SHFLI_REP_2X4_OFFSET_1);
        tw2 = BBE_MULN_2XCQ15PACKQ(tw1, tw1); 
        tw3 = BBE_MULN_2XCQ15PACKQ(tw2, tw1);

        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, C, D, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);

    }
    {   
        BBE_LVNX16_XP(_t0, p_src, stride);                                    
        BBE_LVNX16_XP(_t1, p_src, stride);                                    
        BBE_LVA_XP( p_src, stride);                                    
        BBE_LVB_XP( p_src, -3*stride + 2*BBE_SIMD_WIDTH);                           
         
        BBE_LAVN_2XCQ15_XP(tw1, align, p_tw, 2*4);                
        tw1 = BBE_SHFLN_2XCQ15I(tw1,  BBE_SHFLI_REP_2X4_OFFSET_0);
        tw2 = BBE_MULN_2XCQ15PACKQ(tw1, tw1); 
        tw3 = BBE_MULN_2XCQ15PACKQ(tw2, tw1);


        t4 = BBE_MOVN_2XCQ15_FROMNX16(_t0);                               
        t5 = BBE_MOVN_2XCQ15_FROMNX16(_t1);  
                                                        
        __DFT4_SKEW_T(t0, t1, t2, t3, t4, t5, A, B, tw1, tw2, tw3);   
                                                                            
        _t0 = BBE_MOVNX16_FROMN_2XCQ15(t0);                           
        _t2 = BBE_MOVNX16_FROMN_2XCQ15(t2);                           
        _t1 = BBE_MOVNX16_FROMN_2XCQ15(t1);                           
        _t3 = BBE_MOVNX16_FROMN_2XCQ15(t3);   

        BBE_SALIGNVRNX16_XP(_t0, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
        BBE_SALIGNVRNX16_XP(_t2, uu1, p_dst, 2*BBE_SIMD_WIDTH);            
        BBE_SVINTLARNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH, 1);                
        BBE_SVINTLARNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH, 1);
    }
    BBE_SALIGNVRNX16_XP(_t1, uu0, p_dst, 2*BBE_SIMD_WIDTH);                   
    BBE_SALIGNVRNX16_XP(_t3, uu1, p_dst, 2*BBE_SIMD_WIDTH);        
    return scaling;
}



inline_ ATTRIBUTE_ALWAYS_INLINE int R2_DFT5xIN_5_unroll_all( int16_t *x   /*input*/, 
                          int16_t *y   /*output*/, 
                          const int N, int bexp)
{
    const int N1 = 5;
    const int stride = 2*N/N1; /*in int16_t */
    const int v = N/N1; 
    

    xb_vecN_2xcq15 *px0 = (xb_vecN_2xcq15 *)(x);
    xb_vecN_2xcq15 *px1 = (xb_vecN_2xcq15 *)(x + 1*stride);
    xb_vecN_2xcq15 *px2 = (xb_vecN_2xcq15 *)(x + 2*stride);
    xb_vecN_2xcq15 *px3 = (xb_vecN_2xcq15 *)(x + 3*stride);
    xb_vecN_2xcq15 *px4 = (xb_vecN_2xcq15 *)(x + 4*stride);
    xb_vecN_2xcq15  r5_tw1, r5_tw2, r5_tw3; 
    xb_vecN_2xcq15 t0, t1, t2, t3, t4, t5, t6, t7;

    xb_vecN_2xcq15 *py = (xb_vecN_2xcq15 *)(y);
    int scaling = 0;

    RANGE_BEGIN(N1, -1, 1, bexp, scaling);
    __Pragma("no_reorder");

    RX5_TW(-1, r5_tw1, r5_tw2, r5_tw3);

    ASSERT(v%(BBE_SIMD_WIDTH/2) == 0); 

    
    //for(j=0; j<v; j+=(BBE_SIMD_WIDTH/2))
    {
        BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
        BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
        BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
        BBE_LVB_IP(px3, sizeof(*px1)); 
        BBE_LVA_IP(px4, sizeof(*px1)); 
        t0 = BBE_FFTSRAN_2XCQ15(t0);

        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        

        BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(t5));							
        BBE_MOVSDV(BBE_MOVNX16_FROMN_2XCQ15(t6));							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								       
        BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

    }
    if(v>8)
     {
        BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
        BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
        BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
        BBE_LVB_IP(px3, sizeof(*px1)); 
        BBE_LVA_IP(px4, sizeof(*px1)); 
        t0 = BBE_FFTSRAN_2XCQ15(t0);

        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        

        BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(t5));							
        BBE_MOVSDV(BBE_MOVNX16_FROMN_2XCQ15(t6));							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								       
        BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

    }
    if(v>16)
    {
        BBE_LVN_2XCQ15_IP(t0, px0, sizeof(*px1));
        BBE_LVN_2XCQ15_IP(t1, px1, sizeof(*px2)); 
        BBE_LVN_2XCQ15_IP(t2, px2, sizeof(*px1));
        BBE_LVB_IP(px3, sizeof(*px1)); 
        BBE_LVA_IP(px4, sizeof(*px1)); 
        t0 = BBE_FFTSRAN_2XCQ15(t0);

        /* both DFT4 are forward, expecting BBE_MODE[4] set to 1 */		
        DFT4(t3, t4, t5, t6, t1, t2, A, B, 1);				
        								
        t1 = BBE_FFTADDSSRN_2XCQ15(t0, t3);					
        {  
            xb_vecN_2xcq15 tmpc;
            tmpc = BBE_SRAIN_2XCQ15(t3, 2);   
            t2 = BBE_SUBSN_2XCQ15(t0, tmpc);   /* t2 = t0 - (t3 >> 2);*/                            
        }                                    
        MUL(t4, r5_tw1);							
        MUL(t5, r5_tw2);							
        MUL(t6, r5_tw3);							
        

        BBE_MOVSCV(BBE_MOVNX16_FROMN_2XCQ15(t5));							
        BBE_MOVSDV(BBE_MOVNX16_FROMN_2XCQ15(t6));							
    									
        t3 = t4;								
        								
        DFT4(t4, t5, t6, t7, t2, t3, C, D, 0);				
        								       
        BBE_SVRN_2XCQ15_X(t4, py, 1*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t5, py, 2*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t7, py, 3*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_X(t6, py, 4*2*v*sizeof(int16_t)); 
        BBE_SVRN_2XCQ15_IP(t1, py, sizeof(*py));

    }
    return scaling; 
}   // R2_DFT5xIN_5_unroll_all

#endif // HAVE_FFT
#endif //#ifndef _FFT_COMMON_H_
