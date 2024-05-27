/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2024 Eduardo Zilves & Gonzalo Maso Talou */
/*
 * VesselFilterByTerminal.cpp
 *
 *  Created on: 09/02/2024
 *      Author: Eduardo G. Zilves
 */

#include "VesselFilterByTerminal.h"

VesselFilterByTerminal::VesselFilterByTerminal() : AbstractVesselFilter(){
}

VesselFilterByTerminal::~VesselFilterByTerminal(){
}

vector<SingleVessel*> VesselFilterByTerminal::apply(vector<SingleVessel*> vessels){
	vector<SingleVessel*> filteredVessels;
	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
		if( (*it)->children.size() == 0){
			filteredVessels.push_back(*it);
		}
	}
	return filteredVessels;
}
