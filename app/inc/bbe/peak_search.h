/*
 * peak_search.h
 *
 *  Created on: 2025.11.7
 *      Author: qirui
 */

#ifndef PEAK_SEARCH_H_
#define PEAK_SEARCH_H_

//==========================[ Headers ]=========================
#include "NatureDSP_types.h"
#include "common.h"

//==========================[ Function Prototypes ]=========================
void peak_search(
	int16_t *restrict m,
	uint16_t *restrict idx,
    const int16_t *restrict x,
	const uint16_t N);

#endif /* PEAK_SEARCH_H_ */
