/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * ProjectedCostEstimator.cpp
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 * 				Jul 12, 2023
 * 				Eduardo Guerreiro Zilves
 * 
 */

#include "ProjectedCostEstimator.h"
#include "../vascularElements/SingleVessel.h"

ProjectedCostEstimator::ProjectedCostEstimator(double volumeFactor, double proteolyticFactor, double diffusionFactor): AbstractCostEstimator(){
	previousVolume = 0.0;
	this->proteolyticFactor = proteolyticFactor;
	this->diffusionFactor = diffusionFactor;
	this->volumeFactor = volumeFactor;
	distToParent = 1.0;
	parentRadius = 1.0;
	// projectedSquareDistance = 0.0;
}

ProjectedCostEstimator::~ProjectedCostEstimator(){
}

void ProjectedCostEstimator::previousState(AbstractObjectCCOTree* tree, AbstractVascularElement* parent, point iNew, point iTest, double dLim){
	previousVolume = ((SingleVessel *) tree->getRoot())->treeVolume;

	point a = ((SingleVessel *)parent)->xProx;
	point b = ((SingleVessel *)parent)->xDist;
	//	Parent-to-iNew distance
	//	Parent vessel slope
	point m = b - a;
	//	Parameter for closer projection
	double t = (m ^ (iNew - a)) / (m^m);
	//	Confine t into [0,1] interval
	if (t < 0){
		t = 0;
	}
	else if( t > 1.0){
		t = 1.0;
	}
	//	Closest segment between iNew and parent vessel
	point proj = (iNew - a) - m * t;
	distToParent = sqrt(proj ^ proj);

	parentRadius = ((SingleVessel *)parent)->radius;
}

double ProjectedCostEstimator::computeCost(AbstractObjectCCOTree* tree){
	double volCost = volumeFactor * (computeTreeCost(tree->getRoot()) - previousVolume);
	double proteolysisCost = proteolyticFactor * parentRadius; // 500.0
	double stimulusCost = diffusionFactor * (distToParent * distToParent);
	// cout << "Volumetric cost = " << volCost << ", Protease degradation cost = " << proteolysisCost << ", VEGF/FGF difussion cost = " << stimulusCost << endl;
	return volCost + proteolysisCost + stimulusCost ;
}

double ProjectedCostEstimator::computeTreeCost(AbstractVascularElement* root) {
	double currentCost = ((SingleVessel *)root)->getVolume();
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

AbstractCostEstimator* ProjectedCostEstimator::clone(){
	return (new ProjectedCostEstimator(volumeFactor, proteolyticFactor, diffusionFactor));
}

double ProjectedCostEstimator::getVolumeFactor()
{
	return this->volumeFactor;
}

double ProjectedCostEstimator::getProteolyticFactor()
{
	return this->proteolyticFactor;
}

double ProjectedCostEstimator::getDiffusionFactor()
{
	return this->diffusionFactor;
}

void ProjectedCostEstimator::logCostEstimator(FILE *fp) {
	double v_fac, p_fac, d_fac;        
    v_fac = this->getVolumeFactor();
    p_fac = this->getProteolyticFactor();
    d_fac = this->getDiffusionFactor();
    fprintf(fp, "This domain uses ProjectedCostEstimator.\n");
    fprintf(fp, "Volume factor = %f.\n", v_fac);
    fprintf(fp, "Proteolytic factor = %f.\n", p_fac);
    fprintf(fp, "Diffusion factor = %f.\n", d_fac);
}