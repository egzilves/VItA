/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VesselFilterByTerminal.h
 *
 *  Created on: 09/02/2024
 *      Author: Eduardo G. Zilves
 */

#ifndef FILTERS_VESSELFILTERBYTERMINAL_H_
#define FILTERS_VESSELFILTERBYTERMINAL_H_

#include "AbstractVesselFilter.h"

class VesselFilterByTerminal: public AbstractVesselFilter {
	int stage;
public:
	VesselFilterByTerminal();
	virtual ~VesselFilterByTerminal();

	vector<SingleVessel *> apply(vector<SingleVessel *> vessels);
};

#endif /* FILTERS_VESSELFILTERBYTERMINAL_H_ */
