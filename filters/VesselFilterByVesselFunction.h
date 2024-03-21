/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2024 Eduardo Zilves & Gonzalo Maso Talou */
/*
 * VesselFilterByVesselFunction.h
 *
 *  Created on: 21/03/2024
 *      Author: Eduardo G. Zilves
 */

#ifndef FILTERS_VESSELFILTERBYVESSELFUNCTION_H_
#define FILTERS_VESSELFILTERBYVESSELFUNCTION_H_

#include "AbstractVesselFilter.h"
#include "../structures/vascularElements/AbstractVascularElement.h"

class VesselFilterByVesselFunction : public AbstractVesselFilter{
	AbstractVascularElement::VESSEL_FUNCTION function;
public:
	VesselFilterByVesselFunction(AbstractVascularElement::VESSEL_FUNCTION function);
	virtual ~VesselFilterByVesselFunction();

	vector<SingleVessel *> apply(vector<SingleVessel *> vessels);
};

#endif /* FILTERS_VESSELFILTERBYBRANCHINGMODE_H_ */
