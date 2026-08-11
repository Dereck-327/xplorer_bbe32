/*
 * prepareData_U16_U16_vs.h
 *
 *  Created on: 2026.3.21
 *      Author: qirui
 */

#ifndef PREPAREDATA_U16_U16_VS_H_
#define PREPAREDATA_U16_U16_VS_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void prepareData_U16_U16_vs(
    const uint16_t *restrict input,
    uint16_t *restrict output,
    const uint16_t samples);

#endif /* PREPAREDATA_U16_U16_VS_H_ */
