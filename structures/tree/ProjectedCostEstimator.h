/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * ProjectedCostEstimator.h
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 * 				Jul 12, 2023
 * 				Eduardo Guerreiro Zilves
 * 
 */

#ifndef TREE_PROJECTEDCOSTESTIMATOR_H_
#define TREE_PROJECTEDCOSTESTIMATOR_H_

#include "AbstractCostEstimator.h"

/**
 * Cost estimator that considers projection onto parent vessel direction
 * and minimizes the projection, giving priority to perpendicular vessels.
 */
class ProjectedCostEstimator: public AbstractCostEstimator {
	/**	Volume at the previous step. */
	double previousVolume;
	/**	Vessel wall proteolytic degradation factor. */
	double proteolyticFactor;
	/**	Factor of media diffusion of angiogenic growth factors. */
	double diffusionFactor;
	/**	Weigh for tree volume. */
	double volumeFactor;
	/**	Distance between sprouting source and candidate parent vessel. */
	double distToParent;
	/**	Radius of the candidate parent vessel. */
	double parentRadius;
public:
	/**
	 * Common constructor.
	 */
	ProjectedCostEstimator(double volumeFactor, double proteolyticFactor, double diffusionFactor);
	/**
	 * Common destructor.
	 */
	virtual ~ProjectedCostEstimator();
	/**
	 * Clones the current estimator instance.
	 * @return Cloned instance.
	 */
	AbstractCostEstimator *clone();

	/**
	 * Extracts information of the tree at the previous step.
	 * @param root	Root of the tree at the previous step.
	 * @param parent	Vascular element where the new vessel will be connected.
	 * @param iNew	Distal position of the new vessel.
	 * @param iTest	Proximal position of the new vessel.
	 * @param dLim	Minimum radius distance from the new vessel to the tree.
	 */
	void previousState(AbstractObjectCCOTree *tree, AbstractVascularElement *parent, point iNew, point iTest, double dLim);

	/**
	 * Computes the functional cost of the given tree.
	 * @param tree	Tree at the current step.
	 * @return Cost of the given tree.
	 */
	double computeCost(AbstractObjectCCOTree* tree);
	/**
	 * @return @p volumeFactor
	 */
	double getVolumeFactor();
	/**
	 * @return @p proteolyticFactor
	 */
	double getProteolyticFactor();
	/**
	 * @return @p diffusionFactor
	 */
	double getDiffusionFactor();

	void logCostEstimator(FILE *fp);
private:
	/**
	 * Computes the volume for the tree with root @p root.
	 * @param root	Root of the tree.
	 * @return	Volume of the tree.
	 */
	double computeTreeCost(AbstractVascularElement* root);

};

#endif /* TREE_SPROUTINGVOLUMETRICCOSTESTIMATOR_H_ */
