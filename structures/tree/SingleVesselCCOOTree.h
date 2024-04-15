/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SingleVesselCCOOTree.h
 *
 *	Jan 28, 2021: Important! This class needs cleaning and refactoring.
 *
 *  Created on: Mar 29, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef TREE_SINGLEVESSELCCOOTREE_H_
#define TREE_SINGLEVESSELCCOOTREE_H_

#include <iostream>
#include <string>
#include <vector>
#include<unordered_set>

#include "../../constrains/AbstractConstraintFunction.h"
#include "../vascularElements/SingleVessel.h"
#include "AbstractObjectCCOTree.h"

using namespace std;

/**
 * N-furcation tree with only SingleVessel elements as vascular elements.
 */
class SingleVesselCCOOTree: public AbstractObjectCCOTree {
	string filenameCCO;
	/**	Root radius. */
	double rootRadius;
	/** Convergence tolerance. */
	double variationTolerance;
	/**	Amount of non-common terminals. */
	long long int nCommonTerminals;
	/** Use gamma based on stage. */
	bool isGammaStage;
	friend class PruningCCOOTree;
	friend class BreadthFirstPruning;
	friend class TreeMerger;
	friend class StagedFRROTreeGeneratorLogger;
	// friend class PenetratingVesselTreeGenerator;
public:
	/**
	 * Common tree creator.
	 * @param xi	Root point.
	 * @param rootRadius	Root radius.
	 * @param qi	Flow at the root.
	 * @param gam	Murray law function.
	 * @param epsLim	Sibling vessels ratio function.
	 * @param nu	Viscosity function with respect to the tree level.
	 * @param minAngle	Minimum angle allowed.
	 * @param resistanceVariationTolerance	Convergence tolerance for the iterative viscosity scheme.
	 */
	SingleVesselCCOOTree(point xi, double rootRadius, double qi, AbstractConstraintFunction<double,int> *gam, AbstractConstraintFunction<double,int> *epsLim,
			AbstractConstraintFunction<double,int> *nu, double refPressure, double resistanceVariationTolerance, GeneratorData *instanceData);

//	/**
//	 * Creates a new tree from the .cco file @p filename. The obtained tree has not vtkLine objects of the vessels.
//	 * @param filenameCCO Path to the .cco file.
//	 * @param filenameVTK Path to the .vtk file.
//	 */
//	SingleVesselCCOOTree(string filenameCCO, string filenameVTK, GeneratorData *instanceData);

	/**
	 * Creates a new tree from the .cco file @p filename in VItA format.
	 * @param filenameCCO Path to the .cco file.
	 * @param gam	Murray law function.
	 * @param epsLim	Sibling vessels ratio function.
	 * @param nu	Viscosity function with respect to the tree level.
	 */
	SingleVesselCCOOTree(string filenameCCO, GeneratorData *instanceData, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
			AbstractConstraintFunction<double, int> *nu);
	/**
	 * Creates a new tree from the .cco file @p filename in HeMoLab format.
	 * @param filenameCCO Path to the .cco file.
	 * @param qi	Flow at the root.
	 * @param gam	Murray law function.
	 * @param epsLim	Sibling vessels ratio function.
	 * @param nu	Viscosity function with respect to the tree level.
	 * @param minAngle	Minimum angle allowed.
	 * @param refPressure	Convergence tolerance for the iterative viscosity scheme.
	 * @param viscosityTolerance	Convergence tolerance for the iterative viscosity scheme.
	 */
	SingleVesselCCOOTree(string filenameCCO, GeneratorData* instanceData, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
			AbstractConstraintFunction<double, int> *nu, double refPressure, double viscosityTolerance);
	/**
	 * Common destructor.
	 */
	~SingleVesselCCOOTree();
	/**
	 *	Creates a copy from the tree without the vtk structures (neither from the tree nor from the vessels).
	 * @return	Copy from the tree
	 */
	SingleVesselCCOOTree *clone();
	/**
	 * Returns the closest point in the CCOTree with respect to @p xNew point.
	 * @param xNew	Point from which the minimum distance is computed.
	 * @param xBif	Closest point in the tree to @p xNew.
	 * @param dist	Minimum distance between @p xNew and the tree.
	 */
	void getClosestTreePoint(point xNew, point *xBif, double *dist);

	/**
	 * Return the segments in a close neighborhood of @p xNew. The neighborhood is computed based on the perfusion
	 * volume indicated by @p domain.
	 * @param xNew	Center point of the neighborhood of interest.
	 * @param domain	Domain of the segments.
	 * @param nFound	Amount of segments in the neighborhood.
	 * @return	Array of segments in the neighborhood of @p xNew.
	 */
	vector<AbstractVascularElement *> getCloseSegments(point xNew, AbstractDomain *domain, int *nFound);

	/**
	 * Adds a new vessel to the CCO tree. @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 * @param vesselFunction Vessel function of the added vessel.
	 */
	void addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction);

	void addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, unordered_set<vtkIdType>* partVessels, long long int *termPart, const vector<double> qPart) override;

	/**
	 * Adds a new vessel to the CCO tree. @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 * @param vesselFunction Vessel function of the added vessel.
	 * @param branchingMode Branching mode of the added vessel.
	 */
	void addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, 
					AbstractVascularElement::BRANCHING_MODE branchingMode);
	void addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, 
					AbstractVascularElement::BRANCHING_MODE branchingMode, unordered_set<vtkIdType>* partVessels, long long int *termPart, const vector<double> qPart) override;

	/**
	 * Adds a new vessel to the CCO tree without updating viscosities and radii. 
	 * This is useful for fast generating steps e.g. merging the tree or generating penetrating vessels.
	 * Must call updateMassiveTree afterwards to generate a tree with valid vessels.
	 * @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 * @param vesselFunction Vessel function of the added vessel.
	 * @param branchingMode Branching mode of the added vessel.
	 */
	void addVesselNoUpdate(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, 
					AbstractVascularElement::BRANCHING_MODE branchingMode);

	/**
	 * Adds a new vessel to the CCO tree without updating viscosities and radii. 
	 * This is useful for fast generating steps e.g. merging the tree or generating penetrating vessels.
	 * Must call updateMassiveTree afterwards to generate a tree with valid vessels.
	 * @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 * @param vesselFunction Vessel function of the added vessel.
	 * @param branchingMode Branching mode of the added vessel.
	 * @param addedVesselID return value of the vessel id added.
	 */
	void addVesselNoUpdate(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, 
					AbstractVascularElement::BRANCHING_MODE branchingMode, vtkIdType &addedVesselID);

	/**
	 * Runs the updateTree and updateViscosities functions manually after the addVesselNoUpdate.
	 */
	void updateMassiveTree();
	/**
	 * Runs the updateTree and updateViscosities functions for a subtree given a root vessel.
	 * @param subtreeRoot The root vessel of the subtree.
	 */
	void updateSubtree(SingleVessel* subtreeRoot);
	/**
	 * Runs the updateTree and updateViscosities functions for a subtree given a root vessel.
	 * @param subtreeRoot The root vessel of the subtree. Can be "this->root".
	 * @param tolerance Viscosity tolerance for convergence.
	 */
	void updateSubtree(SingleVessel* subtreeRoot, double tolerance);
	/**
	 * Scale the tree root radius by a factor.
	 * @param scaleFactor Factor to scale the diameters. (1.0 to not scale)
	 */
	void scaleTreeRootRadius(double scaleFactor);
	/**
	 * Scale the entire tree radii by a factor. Do not use with the scaleTreeRoot.
	 * @param scaleFactor Factor to scale the diameters. (1.0 to not scale)
	 */
	void scaleTreeRadius(double scaleFactor);

	void addVesselMergeFast(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, int stage, unordered_map<string, SingleVessel *>* stringToPointer);

	void addVesselMerge(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, int savedStage, unordered_map<string, SingleVessel *>* stringToPointer);
	/** 
	 * Adds a vessel that has already been validated. This function is used by BreadthFirstPrunning.
	 * @param newVessel Vessel that will be added.
	 * @param originalVessel Vessel that has been validaded in a previous tree.
	 * @param copiedTo Mapping such that copiedTo[originalVessel] = copiedVessel;
	*/
	//	FIXME This function probably should be part of other class
	void addValitatedVessel(SingleVessel *newVessel, SingleVessel *originalVessel, unordered_map<SingleVessel *, SingleVessel *>& copiedTo);

	/** 
	 * Adds a vessel that has already been validated, but do not update the tree. This function is used by BreadthFirstPrunning.
	 * @param newVessel Vessel that will be added.
	 * @param originalVessel Vessel that has been validaded in a previous tree.
	 * @param copiedTo Mapping such that copiedTo[originalVessel] = copiedVessel;
	*/
	//	FIXME This function probably should be part of other class
	void addValitatedVesselFast(SingleVessel *newVessel, SingleVessel *originalVessel, unordered_map<SingleVessel *, SingleVessel *>& copiedTo);

	/**
	 * Adds a new pregenerated subtree as a replacement to a terminal vessel, or as is appended to its parent.
	 * TODO: rename this to replaceSubtree and create new addSubtree without terminal.
	 * TODO: adapt this to replace entire subtrees.
	 * The terminal vessel is passed as a SingleVessel(AbstractVascularElement) where the replacing will occur.
	 * TODO: If terminal is null, tree is appended to parent with no scaling.
	 * The subtree is passed as another SingleVesselCCOOTree object with the correct point coordinates mapped *BEFORE* this step.
	 * This ensures no resource wasted with allocation.
	 * The CCO tree diameters are not updated, must use the Update Massive Tree method to validate the vessels.
	 * TODO: if extend this function, use this signature as a wrapper for extended method.
	 * @param subtree The Subtree to be appended to the tree, replacing the terminal vessel.
	 * @param terminalVessel The terminal vessel to be replaced in the operation.
	 * -@param parent The parent vessel to the new subtree. // get parent inside with vessel->getParent()
	 * @param nNewTerms Number of new terminals to add
	 */
	void addSubtree(AbstractObjectCCOTree *newSubtree, AbstractVascularElement *oldTerminalVessel, int nNewTerms);

//	/**
//	 * Adds a new vessel to the CCO tree as continuation of the pre-existent vessel @p parent. @param xDist is the distal nodes of the new
//	 * vessel and @param parent is the proximal attachment parent vessel.
//	 * @param xDist Distal point of the new vessel.
//	 * @param parent Parent to the new vessel.
//	 * @param mode Branching mode of the added vessel.
//	 * @param vesselFunction Vessel function of the added vessel.
//	 */
//	void addVessel(point xDist, AbstractVascularElement *parent, AbstractVascularElement::BRANCHING_MODE mode, AbstractVascularElement::VESSEL_FUNCTION vesselFunction);

	/**
	 * For a given spatial point @p xNew test its connection with @p parent vessel. It must evaluate if the restrictions
	 * of geometry and symmetry are satisfied and also if it do not intersects with other vessel of this tree. It returns
	 * in @p cost the functional variation by the inclusion of such vessel.
	 * @param xNew	Distal point for the new vessel to test.
	 * @param parent	Parent vessel to test the @p xNew point connection.
	 * @param domain	Tree domain.
	 * @param neighbors	Close neighbors used for intersection test.
	 * @param dlim	Not used in the current implementation.
	 * @param xBif	Proximal point for the new vessel that present the lower function cost.
	 * @param cost	Functional cost variation for the best bifurcation position.
	 * @return	If the connection of the tree with xNew is possible. If not @p cost is INFINITY.
	 */
	int testVessel(point xNew, AbstractVascularElement *parent, AbstractDomain *domain, vector<AbstractVascularElement *> neighbors, double dlim, point *xBif, double *cost);

	int testVessel(point xNew, AbstractVascularElement *parent, AbstractDomain *domain, vector<AbstractVascularElement *> neighbors, double dlim, point *xBif, double *cost, unordered_set<vtkIdType>* partVessels, long long int *termPart, const vector<double> qPart) override;

	/**
	 * Prints the current tree node by node.
	 */
	void print();

	/**
	 * Returns the class name identifier.
	 * @return Class name.
	 */
	string getTreeName();

	/**
	 * Returns the root radius.
	 * @return Root radius.
	 */
	double getRootRadius();

	/**
	 * Removes all branches at stage @p stage such that do not branch in the next stages.
	 * @param stage
	 */
	void removeWitheredBranches(int stage);

	/**
	 * Removes a vessel and its subtree from the current tree.
	 * @param vessel Vessel to be removed.
	 */
	void remove(SingleVessel *vessel);

	/**
	 * Returns if the @p vessel at stage @p stage has no connexions with vessels at stage @p stage + 1.
	 * @param vessel Vessel to evaluate if it is withered or not.
	 * @param stage Stage of the vessel @p vessl
	 * @return If it is withered or not.
	 */
	bool isWithered(SingleVessel *vessel);

	string getFilenameCCO();

	void updateAll();

	/**
	 * Set if the tree updates gamma using stage or level.
	 */
	void setIsGammaStage(bool isGammaStage);

	/**
	 * Set the root vessel
	 */
	void setRoot(SingleVessel* newRoot);

	/**
	 * Bypass the vessel function when generating partitioned domains, allows bifurcating from partly-outside vessels
	 */
	bool bypassFunctionIfMidpointInside;

	/**
	 * Bypass the isValidOpeningAngle() function when the parent vessel is rigid
	 * and the planeNormal is {0,0,0}, and openingAngle returns -nan(0x8000...) NaN
	 */
	bool bypassRigidParentOpeningAnglePlaneNormal;

protected:
	/**
	 * Returns a string with the tree atributes to create the .cco file.
	 * @param outFile	File writer for the .cco file.
	 */
	void saveTree(ofstream *outFile);

private:
	/**
	 * Clones the subtree with parent vessel @p levels .
	 * @param root	Root of the tree to clone.
	 * @param segments	Segments of the tree.
	 * @return Cloned subtree.
	 */
	SingleVesselCCOOTree *cloneUpTo(int levels, SingleVessel *parent);
	/**
	 * Clones the subtree with parent vessel @p root recursively.
	 * @param root	Root of the tree to clone.
	 * @param segments	Segments of the tree.
	 * @return Cloned subtree.
	 */
	SingleVessel *cloneTree(SingleVessel *root, unordered_map<long long, AbstractVascularElement *> *segments);
	/**
	 * Returns a partial variation of the cost functional due to the new segment inclusion.
	 * @param xNew	Proximal point of the new vessel.
	 * @param xTest Distal point of the new vessel.
	 * @param parent Parent to the new vessel.
	 * @param dLim Minimum distance from the new vessel to the tree.
	 */
	double evaluate(point xNew, point xTest, SingleVessel *parent, double dLim);
	double evaluate(point xNew, point xTest, SingleVessel *parent, double dLim, unordered_set<vtkIdType>* partVessels, long long int *termPart, const vector<double> qPart);
	/**
	 * Returns a partial variation of the cost functional due to the new segment inclusion. This method is only used for DISTAL_BRANCHING vessels.
	 * @param xNew	Proximal point of the new vessel.
	 * @param parent Parent to the new vessel.
	 * @param dLim Minimum distance from the new vessel to the tree.
	 */
	double evaluate(point xNew, SingleVessel *parent, double dLim);
	double evaluate(point xNew, SingleVessel *parent, double dLim, unordered_set<vtkIdType>* partVessels, long long int *termPart, const vector<double> qPart);
	/**
	 * Updates the tree values for the current topology in only one tree "in order" swept (O(N)).
	 * As the recursion deepens, the level number is computed for each element. As the
	 * recursion is returning, it computes the flow and resistance for the current node and the
	 * radius ratio for its childs.
	 * @param root Root vessel for the tree to update.
	 * @param tree Tree to update.
	 */
	void updateTree(SingleVessel *root, SingleVesselCCOOTree *tree);

	void updateTree(SingleVessel *root, SingleVesselCCOOTree *tree, unordered_set<vtkIdType>* partVessels, long long int *nPartTotal, const vector<double> qPart);
	/**
	 * For a giving pair of beta between sibling of a parent vessel, it analyze the symmetry constrain given by
	 * epsLim function.
	 *
	 * @param beta1	Beta for the 1st sibling.
	 * @param beta2	Beta for the 2nd sibling.
	 * @param nLevel Tree level of the bifurcation.
	 * @return	If the simmetry constrain is satisfied.
	 */
	int isSymmetricallyValid(double beta1, double beta2, int nLevel);
	/**
	 * Checks if the angles of the parent vessel and the new vessel do not violate the minimum angle constraint.
	 * @param xBif	Bifurcation point between the new vessel and the distal part of the parent vessel (iCon).
	 * @param xNew	Distal point of the new vessel.
	 * @param parent	Parent vessel.
	 * @param minAngle 	Minimum angle constraint.
	 * @return	If the angles are higher than the minimum allowed.
	 */
	int areValidAngles(point xBif, point xNew, SingleVessel *parent, double minAngle);
	/**
	 * Checks if the angle between the plane of the parent and sibling vessel and the new vessel satisfies the opening angle constraint.
	 * @param xBif	Bifurcation point between the new vessel and the distal part of the parent vessel (iCon).
	 * @param xNew	Distal point of the new vessel.
	 * @param parent	Parent vessel.
	 * @param minPlaneAngle 	Minimum opening angle constraint.
	 * @return	If the angles satisfies the minimum plane angle.
	 */
	int isValidOpeningAngle(point xBif, point xNew, SingleVessel *parent, double minPlaneAngle);
	/**
	 * It returns if the line @param p1 - @param p2 intersects any vessel of the tree beside parent.
	 * @param p1	Extreme point 1 of the line.
	 * @param p2	Extreme point 2 of the line.
	 * @param parent	Vessel excluded from the checking.
	 * @param boundaryTol	Factor of line contraction at the extremes to avoid false intersections due to contact with anastomose.
	 * @return	If the segment intersects any segment of the tree.
	 */
	int isIntersectingVessels(point p1, point p2, SingleVessel *parent, vector<AbstractVascularElement *> neighbors);
	/**
	 * Determines if the segment xBif-xNew is closer than dLim to any other segment of the tree (without being its parent
	 * vessel). Not used.
	 * @param xBif	Proximal point for the new segment.
	 * @param xNew	Distal point for the new segment.
	 * @param parent	Parent vessel of the new segment.
	 * @param dLim	Perfusion volume for each terminal at the current state of the tree.
	 * @return	If the middle point of the segment is sufficiently distant from the tree.
	 */
	int isOverlapped(point xBif, point xNew, SingleVessel* parent, double dLim);
	/**
	 * Updates the tree beta values for the current vessel diameters.
	 * @param root	Tree root.
	 * @param parentRadius	1.0 if root is actually the tree root, or parent radius of root otherwise.
	 * @param maxBetaVariation	The maximum variation of the beta value due to the update.
	 */
	void updateTreeViscositiesBeta(SingleVessel *root, double *maxBetaVariation);
	/**
	 * Returns the viscosity estimated with the Fahraeus-Lindquist model.
	 * @param radius
	 * @return Viscosity value.
	 */
	double getNuFL(double radius);

	/**
	 * Creates the VTK lines and points associated to a HeMoLab file loaded.
	 */
	void createSegmentVtkLines(AbstractVascularElement *rootVessel);

	double getVariationTolerance();

	/**
	 * Checks if length/radius > 2.
	 */
	bool isValidAspectRatio(SingleVessel *vessel);

	/**
	 * Checks if length ratio is not smaller than 10% of parent, receives parent iPar and continuation iCon, not new vessel
	 */
	bool isValidLengthRatio(SingleVessel *vessel1, SingleVessel *vessel2);

	/*
	* Returns the gamma.
	*/
	double getGamma(SingleVessel* vessel);

	/** 
	 * Get vessel resistance.
	 */
	double getRealViscosity(SingleVessel *vessel);
};

#endif /* TREE_SINGLEVESSELCCOOTREE_H_ */
