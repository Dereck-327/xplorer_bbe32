/**
    @file       dsp_ops.h
    @brief      dsp 算子重命名 指向`bbe/`下的集成到项目的算子

    @version    0.1.0
    @date       2026-08-10
    @author     hjk

*/

#ifndef SIMULATE_DSP_OPS_H__
#define SIMULATE_DSP_OPS_H__

#include "prepareData_U16_U16_vs.h"
#include "nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs.h"
#include "zero_padding.h"
#include "copyConvRes_S16_S16_v.h"
#include "calc_peaks_info.h"
#include "peak_idx_comp.h"
#include "calc_peaks_idx_frac.h"

#define DSP_PrepareData			prepareData_U16_U16_vs
#define DSP_NciFourFft          nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs
#define DSP_MulS16U16           mul_S16_U16_vws
#define DSP_ConvS16Vw15         conv_S16_S16_vw_size15
#define DSP_CopyConvRes         copyConvRes_S16_S16_v
#define DSP_CalcPeaksInfo       calc_peaks_info
#define DSP_PeakIdxComp         peak_idx_comp
#define DSP_CalcPeaksIdxFrac    calc_peaks_idx_frac

#define DSP_ZeroPad             zero_padding


#endif /* SIMULATE_DSP_OPS_H__ */
