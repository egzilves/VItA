/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2024 Eduardo Zilves & Gonzalo Maso Talou */
/*
 * VesselFilterByVesselFunction.cpp
 *
 *  Created on: 21/03/2024
 *      Author: Eduardo G. Zilves
 */

#include "VesselFilterByVesselFunction.h"

VesselFilterByVesselFunction::VesselFilterByVesselFunction(AbstractVascularElement::VESSEL_FUNCTION function) : AbstractVesselFilter() {
	this->function = function;
}

VesselFilterByVesselFunction::~VesselFilterByVesselFunction(){
	// TODO Auto-generated destructor stub
}

vector<SingleVessel*> VesselFilterByVesselFunction::apply(vector<SingleVessel*> vessels){
	vector<SingleVessel*> filteredVessels;
	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if( (*it)->vesselFunction == function){
			filteredVessels.push_back(*it);
		}
	}
	return filteredVessels;
}
