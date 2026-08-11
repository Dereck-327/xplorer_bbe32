/*
 * nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs.h
 *
 *  Created on: 2026.3.21
 *      Author: qirui
 */

#ifndef NCI_U16_S16_WITHMAXBEXPSCALING_VSM__FOURFFTS_H_
#define NCI_U16_S16_WITHMAXBEXPSCALING_VSM__FOURFFTS_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void nci_U16_S16_withMaxBexpScaling_vsm__fourFFTs(
    const uint16_t *restrict mag2,
	const int16_t *restrict mag2_bexp,
    int16_t *restrict mag2_acc,
	int16_t *restrict mag2_acc_bexp,
    const uint16_t samples);

#endif /* NCI_U16_S16_WITHMAXBEXPSCALING_VSM__FOURFFTS_H_ */
