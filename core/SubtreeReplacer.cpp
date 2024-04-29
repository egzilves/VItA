/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2024 Eduardo Guerreiro Zilves & Gonzalo Maso Talou */
/*
 * SubtreeReplacer.cpp
 *
 *  Created on: Mar 20, 2024
 *      Author: Eduardo G. Zilves
 */

#include "SubtreeReplacer.h"

#include <cmath>
#include <fstream>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <string>

#include <vtkPolyDataReader.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataNormals.h>

#include "../io/VTKObjectTreeElementalWriter.h"
#include "../io/VTKObjectTreeSplinesNodalWriter.h"

#include "../io/VTKConverter.h"
#include "../io/VTKObjectTreeNodalWriter.h"
#include "../io/VTKObjectTreeSplinesNodalWriter.h"

#include "../structures/CCOCommonStructures.h"
#include "../structures/tree/SingleVesselCCOOTree.h"
#include "../structures/vascularElements/AbstractVascularElement.h"
#include "../utils/MemoryMonitor.h"

#include "../filters/AbstractVesselFilter.h"
#include "../filters/VesselFilterComposite.h"
#include "../filters/VesselFilterByBranchingMode.h"
#include "../filters/VesselFilterByTerminal.h"
#include "../filters/VesselFilterByVesselFunction.h"
#include "../filters/VesselFilterByStage.h"

#include "GeneratorData.h"

#include <omp.h>


#include "../io/task/AbstractSavingTask.h"
#include "../io/task/VisualizationSavingTask.h"
#include "../io/task/CheckpointSavingTask.h"
#include "../stats/StructuredTreeStatsManager.h"
#include "../stats/VesselStructHandler.h"

#include "../structures/tree/FRRVaViOptCCOSTree.h"
#include "../structures/tree/SingleVesselCCOOTree.h"
#include "../structures/tree/VolumetricCostEstimator.h"
#include "../structures/tree/SproutingVolumetricCostEstimator.h"
#include "../structures/tree/AdimSproutingVolumetricCostEstimator.h"
#include "../structures/vascularElements/AbstractVascularElement.h"

#include "../structures/domain/AbstractDomain.h"
#include "../structures/domain/SimpleDomain.h"
#include "../structures/domain/PartiallyVascularizedDomain.h"
#include "../structures/domain/StagedDomain.h"
#include "../structures/domain/TreeProjector.h"

#include "StagedFRROTreeGenerator.h"
#include "FRRSTreeGenerator.h"
#include "GeneratorData.h"

#include "../constrains/AbstractConstraintFunction.h"
#include "../constrains/ConstantConstraintFunction.h"
#include "../constrains/ConstantPiecewiseConstraintFunction.h"


SubtreeReplacer::SubtreeReplacer(
		StagedDomain* domain, string projectionDomainFile, AbstractObjectCCOTree* tree, long long nTerm,
		vector<AbstractConstraintFunction<double,int> *>gam,
		vector<AbstractConstraintFunction<double,int> *>epsLim,
		vector<AbstractConstraintFunction<double,int> *>nu) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	SingleVessel::bifurcationTests = instanceData->nBifurcationTest;
	this->nTerminals = nTerm;
	this->tree = (SingleVesselCCOOTree*) tree;
	//	Stage can be loaded from file
	this->stage = domain->getCurrentStage();
	this->tree->setCurrentStage(domain->getCurrentStage());

	this->gams = gam;
	this->epsLims = epsLim;
	this->nus = nu;

	this->tree->setGam(gam[0]);
	this->tree->setEpsLim(epsLim[0]);
	this->tree->setNu(nu[0]);

	this->dLim = domain->getDLim(1, instanceData->perfusionAreaFactor);

	this->isGeneratingConfFile = 0;
	this->confFilename = "";

	this->dataMonitor = new GeneratorDataMonitor(domain);
	this->monitor = new MemoryMonitor(MemoryMonitor::MEGABYTE);

	this->didAllocateTree = false;

	// this->descendingOffset = 0;
	// this->endpointOffset = 0;

	// this->projectionDomainFile = projectionDomainFile;
	// this->domainFile = projectionDomainFile;

/*
	this->maxIterationsLimit = 1000000;


	this->gen_data_0 = new GeneratorData(16000, 2000, 0.95, 1.0, 1.0, 0.25, 7, 0, false, new VolumetricCostEstimator());
    this->gam_0 = new ConstantConstraintFunction<double, int>(3.0);
    this->eps_lim_1 = new ConstantPiecewiseConstraintFunction<double, int>({0.0, 0.0},{0, 2});
    this->nu = new ConstantConstraintFunction<double, int>(3.6); //cP

	this->toAppendVesselData.clear();
*/

}

SubtreeReplacer::~SubtreeReplacer() {
	delete this->dataMonitor;
	delete this->monitor;
	if (this->didAllocateTree) {
		delete this->tree;
	}
}


AbstractObjectCCOTree *SubtreeReplacer::appendSubtree(long long int saveInterval, string tempDirectory, string subtreeFilename, string filenameData, int max_iterations_count){
	if (!allowThisClass) {
		cout << "FATAL: experimental class, set bool 'allowThisClass' to true to use this" << endl;
		exit(1);
	}
	this->beginTime = time(nullptr);
	generatesConfigurationFile(ios::out);

	int generatedVessels = 0;
	string modelsFolder = "./";
	string outputDir = "./";
	string prefix = "output";



	// pass parameters, longTreeList.cco, shortTreeList.cco, percentages
	vector<vector<string>> populations; // each element is a population, contains a list of cco files
	vector<double> accumulatedPercentages; // the ACCUMULATED distributions for each population

	// Filter vessels by type
	// filter the penetrating vessels, distalbranching, etc.
	// terminal && function=penetrating && mode=distal
	AbstractVascularElement::VESSEL_FUNCTION vesselfunction = AbstractVascularElement::VESSEL_FUNCTION::PERFORATOR; //penetrating
	AbstractVascularElement::BRANCHING_MODE branchingmode = AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING; //penetrating
	AbstractVesselFilter *replacedFilters = new VesselFilterComposite({
		new VesselFilterByTerminal(), 
		new VesselFilterByVesselFunction(vesselfunction), 
		new VesselFilterByBranchingMode(branchingmode)
		});
	vector<SingleVessel *> treeVessels = this->tree->getVessels();
	vector<SingleVessel *> replacedVessels = replacedFilters->apply(treeVessels);

	/// TODO: for each (SingleVessel *) vessel
	int maxIterations = max_iterations_count;
	int itCount = 0;
    GeneratorData *gen_data_0 {new GeneratorData(16000, 2000, 0.95,
        1.0, 1.0, 0.25, 7, 0, false, new VolumetricCostEstimator())};
    AbstractConstraintFunction<double,int> *gam_0 {new ConstantConstraintFunction<double, int>(3.0)};
    AbstractConstraintFunction<double, int> *eps_lim_1 {new ConstantPiecewiseConstraintFunction<double, int>({0.0, 0.0},{0, 2})};
    AbstractConstraintFunction<double,int> *nu {new ConstantConstraintFunction<double, int>(3.6)}; //cP

	// Where to save the data to be read.
	vector<ProxySegment> proxySegments;
	proxySegments.clear();
	// Load the segment ID of parent + coordinates for subtree, to find the correct when appending the tree.
	// This was done via a new method that returns the ID via parameter. Now the function:
	ifstream streamIn;
	streamIn.open(filenameData.c_str(), ios::in);
	if (!streamIn.is_open()) {
		cout << "ERROR: file was not open!" << endl;
		return tree;
	}
	cout << "Loading data from filename " << filenameData << " and saving to memory." << endl;
	vtkIdType parentVesselID;
	double x1, y1, z1, x2, y2, z2;
	while (streamIn >> parentVesselID >> x1 >> y1 >> z1 >> x2 >> y2 >> z2) {
		point pProx;
		pProx.p[0] = x1;
		pProx.p[1] = y1;
		pProx.p[2] = z1;
		point pDist;
		pDist.p[0] = x2;
		pDist.p[1] = y2;
		pDist.p[2] = z2;
		// this->toAppendVesselData[parentVesselID] = vector<point> {pProx, pDist};
		proxySegments.push_back(ProxySegment(parentVesselID, pProx, pDist));
	}
	streamIn.close();
	cout << "Data loaded." << endl;

	cout << "WARNING: testing for 1 subtree only. Filename: " << subtreeFilename << " " << endl;

	/// TODO: sort type of tree
	cout << "WARNING: limiting max iterations to " << maxIterations << endl;
	// for (vector<SingleVessel *>::iterator it = replacedVessels.begin(); it != replacedVessels.end() && itCount<maxIterations; ++it, ++itCount) {
	// I should instead iterate over the list of proxy vessels, and not the subtree. This ensures we don't try to access an invalid vessel.
	for (vector<ProxySegment>::iterator it = proxySegments.begin(); it != proxySegments.end() && itCount<maxIterations; ++it, ++itCount) {
		// SingleVessel* parentVessel = (*it);
		vtkIdType parentSegmentID = (*it).parentID;
		cout << "iteration number " << itCount << " adding subtree for parent segment id " << parentSegmentID << "\n";
		// TODO: get properties, get distal (coordinates), get radius, length
		// point vesselProx = this->toAppendVesselData[parentVessel->vtkSegmentId][0]; // xProx
		point vesselProx = (*it).xProx; // xProx
		// point vesselDist = this->toAppendVesselData[parentVessel->vtkSegmentId][1]; // xDist
		point vesselDist = (*it).xDist; // xDist

		// Instantiate a new subtree
		// string subtreeFilename = subtreeFilename;
		SingleVesselCCOOTree *newSubtree {new SingleVesselCCOOTree(subtreeFilename, gen_data_0, gam_0, eps_lim_1, nu)};
		(*newSubtree).setIsInCm(true);
		vector<SingleVessel *> subtreeVessels = newSubtree->getVessels();


		// Map geometry of subtree
		// map: 0,0,0 -> proximal, map 0,0,h -> distal with h=0.25cm
		// map xyz-translation, map xy-rotation, map z-scale
		// keep xy-scale (or use z-scale)
		// TODO: random z-rotation
		// cout << "tree imported, calculating basing characteristics" << "\n";
		// NOTE: assuming subtree is generated from (0,0,0) to (0,0,h)
		point originSubtree = {0,0,0};
		double heightSubtree = 0.25; // NOTE: assuming h = 2.5mm, and shorter vessels (1.0mm) will be short penetrating
		point terminalSubtree = {0,0,heightSubtree};
		point displacementSubtree = terminalSubtree-originSubtree;
		double lengthSubtree = sqrt(displacementSubtree^displacementSubtree);
		point unitSubtree = displacementSubtree/lengthSubtree;

		point displacement = vesselDist-vesselProx; 
		double length = sqrt(displacement^displacement);
		point unitDirection = displacement/length;

		// cout << "linear mapping the tree..." << "\n";
		// SCALING
		double scaleFactor = length / heightSubtree;
		// scale the tree, for each terminal scale distal/proximal points
		for (vector<SingleVessel *>::iterator itVessel = subtreeVessels.begin(); itVessel != subtreeVessels.end(); ++itVessel) {
			(*itVessel)->xProx.p[0] = (*itVessel)->xProx.p[0]*scaleFactor;
			(*itVessel)->xDist.p[0] = (*itVessel)->xDist.p[0]*scaleFactor;
			(*itVessel)->xProx.p[1] = (*itVessel)->xProx.p[1]*scaleFactor;
			(*itVessel)->xDist.p[1] = (*itVessel)->xDist.p[1]*scaleFactor;
			(*itVessel)->xProx.p[2] = (*itVessel)->xProx.p[2]*scaleFactor;
			(*itVessel)->xDist.p[2] = (*itVessel)->xDist.p[2]*scaleFactor;
			(*itVessel)->length = (*itVessel)->length * scaleFactor;
			(*itVessel)->radius = (*itVessel)->radius * scaleFactor;
		}
		// cout << "scaled, now rotating" << "\n";
		// ROTATION
		// Rodrigues formula.
		point u {unitSubtree};
		point Ru = {unitDirection};
		matrix Identity = {{1,0,0, 0,1,0, 0,0,1}};
		matrix Rotation;
		matrix Krotation;
		double smallValue = 1e-5;
		bool calculateRotation = true;
		double cosineAngle = u^Ru;
		if (abs(cosineAngle-1) < smallValue){
			Rotation = Identity;
			calculateRotation = false;
		}
		if (abs(cosineAngle+1) < smallValue){
			Rotation = Identity*(-1);
			calculateRotation = false;
		}
		if (calculateRotation){
			Krotation = outer(Ru,u) - outer(u,Ru);
			Rotation = Identity + Krotation + (Krotation*Krotation)/(1+cosineAngle);
		}
		// TODO: add random rotation, add matrix to rotate in xy plane, z axis, random angle.
		// Rotate the points for each terminal, multiply R*v for every point v
		for (vector<SingleVessel *>::iterator itVessel = subtreeVessels.begin(); itVessel != subtreeVessels.end(); ++itVessel) {
			(*itVessel)->xProx = Rotation*(*itVessel)->xProx;
			(*itVessel)->xDist = Rotation*(*itVessel)->xDist;
		}
		// cout << "rotated, now translating" << "\n";
		// TRANSLATION
		point translationVector = vesselProx - originSubtree;
		// translate for each point
		for (vector<SingleVessel *>::iterator itVessel = subtreeVessels.begin(); itVessel != subtreeVessels.end(); ++itVessel) {
			(*itVessel)->xProx = (*itVessel)->xProx + translationVector;
			(*itVessel)->xDist = (*itVessel)->xDist + translationVector;
		}
		cout << "subtree ready for append" << "\n";
		// Now the subtree is geometrically located in the correct point. Time to replace the subtree.

		// TODO: make subtree and append
		// Read CCO, map root, and childs recursively
		// map proximal and distal of subtrees, recursively for every child.
		// update radius, update tree
		// int newTerms = 101;
		// cout << "WARNING: hardcode for " << newTerms << " new terms in subtree (ignore)" << "\n";
		// tree->addSubtree(newSubtree, parentVessel, newTerms);
		tree->appendSubtree(newSubtree, newSubtree->getRoot(), nullptr /* parentVessel */, 0 /* level */, parentSegmentID);

		// DO NOT CLEAR ELEMENTS, i reverted the NoAlloc, tree is copyed instead of moved, and deletion should occur normally.
		// newSubtree->clearElements();

		delete newSubtree;
	}
	


	cout << "iterated through all vessels" << endl;


	// yes run it because we dont want to update everything after EVERY vessel.

    // Update tree
	cout << "updating the tree" << endl;
	((SingleVesselCCOOTree*) tree)->updateMassiveTree();
	cout << "tree updated" << endl;


	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());


	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;


}

vector<SingleVessel *> SubtreeReplacer::getFilteredVessels() {
	// Filter vessels by type
	// filter the penetrating vessels, distalbranching, etc.
	// terminal && function=penetrating && mode=distal
	AbstractVascularElement::VESSEL_FUNCTION vesselfunction = AbstractVascularElement::VESSEL_FUNCTION::PERFORATOR; //penetrating
	AbstractVascularElement::BRANCHING_MODE branchingmode = AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING; //penetrating
	int stage = this->stageToAppend;
	AbstractVesselFilter *replacedFilters = new VesselFilterComposite({
		new VesselFilterByTerminal(), 
		new VesselFilterByVesselFunction(vesselfunction), 
		new VesselFilterByBranchingMode(branchingmode),
		new VesselFilterByStage(stage)
		});
	vector<SingleVessel *> treeVessels = this->tree->getVessels();
	vector<SingleVessel *> replacedVessels = replacedFilters->apply(treeVessels);
	delete replacedFilters;
	return replacedVessels;
}

SingleVesselCCOOTree *SubtreeReplacer::buildNewSubtree(string subtreeFilename) {
		// Instantiate a new subtree
		// string subtreeFilename = subtreeFilename;
		SingleVesselCCOOTree *newSubtree {new SingleVesselCCOOTree(subtreeFilename, this->gen_data_0, this->gam_0, this->eps_lim_1, this->nu)};
		return newSubtree;
}

void SubtreeReplacer::mapSubtree(SingleVesselCCOOTree *subtree, point xProx, point xDist) {
	// Map geometry of subtree
	// map: 0,0,0 -> proximal, map 0,0,h -> distal with h=0.25cm
	// map xyz-translation, map xy-rotation, map z-scale
	// keep xy-scale (or use z-scale)
	// TODO: random z-rotation
	// NOTE: assuming subtree is generated from (0,0,0) to (0,0,h)
	cout << "tree imported, calculating basic characteristics" << endl;
	// NOTE: assuming subtree is generated from (0,0,0) to (0,0,h)
	point originSubtree = {0,0,0};
	double heightSubtree = 0.25; // NOTE: assuming h = 2.5mm, and shorter vessels (1.0mm) will be short penetrating
	point terminalSubtree = {0,0,heightSubtree};
	cout << "WARNING: hardcode for " << heightSubtree << " [cm] height of subtree" << endl;
	// find the versors of subtree
	point displacementSubtree = terminalSubtree-originSubtree;
	double lengthSubtree = sqrt(displacementSubtree^displacementSubtree);
	point unitVersorSubtree = displacementSubtree/lengthSubtree;
	// find the versors of the target position
	point displacement = xProx-xDist; 
	double length = sqrt(displacement^displacement);
	point unitVersorDirection = displacement/length;

	cout << "linear mapping of the tree..." << endl;
	vector<SingleVessel *> subtreeVessels = subtree->getVessels();
	// SCALING
	double scaleFactor = length / heightSubtree;
	// scale the tree, for each terminal scale distal/proximal points
	for (vector<SingleVessel *>::iterator itVessel = subtreeVessels.begin(); itVessel != subtreeVessels.end(); ++itVessel) {
		(*itVessel)->xProx.p[0] = (*itVessel)->xProx.p[0]*scaleFactor;
		(*itVessel)->xDist.p[0] = (*itVessel)->xDist.p[0]*scaleFactor;
		(*itVessel)->xProx.p[1] = (*itVessel)->xProx.p[1]*scaleFactor;
		(*itVessel)->xDist.p[1] = (*itVessel)->xDist.p[1]*scaleFactor;
		(*itVessel)->xProx.p[2] = (*itVessel)->xProx.p[2]*scaleFactor;
		(*itVessel)->xDist.p[2] = (*itVessel)->xDist.p[2]*scaleFactor;
		// and other dimensional properties
		(*itVessel)->length = (*itVessel)->length * scaleFactor;
		(*itVessel)->radius = (*itVessel)->radius * scaleFactor;
	}
	cout << "scaled, now rotating" << endl;
	// ROTATION
	// Rodrigues formula.
	point u {unitVersorSubtree};
	point Ru = {unitVersorDirection};
	matrix Identity = {{1,0,0, 0,1,0, 0,0,1}};
	matrix Rotation;
	matrix Krotation;
	double smallValue = 1e-5;
	bool calculateRotation = true;
	double cosineAngle = u^Ru;
	if (abs(cosineAngle-1) < smallValue){
		Rotation = Identity;
		calculateRotation = false;
	}
	if (abs(cosineAngle+1) < smallValue){
		Rotation = Identity*(-1);
		calculateRotation = false;
	}
	if (calculateRotation){
		Krotation = outer(Ru,u) - outer(u,Ru);
		Rotation = Identity + Krotation + (Krotation*Krotation)/(1+cosineAngle);
	}
	// TODO: add random rotation, add matrix to rotate in xy plane, z axis, random angle.
	// Rotate the points for each terminal, multiply R*v for every point v
	for (vector<SingleVessel *>::iterator itVessel = subtreeVessels.begin(); itVessel != subtreeVessels.end(); ++itVessel) {
		(*itVessel)->xProx = Rotation*(*itVessel)->xProx;
		(*itVessel)->xDist = Rotation*(*itVessel)->xDist;
	}

	cout << "rotated, now translating" << endl;
	// TRANSLATION
	point translationVector = xProx - originSubtree;
	// translate for each point
	for (vector<SingleVessel *>::iterator itVessel = subtreeVessels.begin(); itVessel != subtreeVessels.end(); ++itVessel) {
		(*itVessel)->xProx = (*itVessel)->xProx + translationVector;
		(*itVessel)->xDist = (*itVessel)->xDist + translationVector;
	}
	cout << "subtree ready for replacement" << endl;
	// Now the subtree is geometrically located in the correct point. Time to replace the subtree.
	return;
}

AbstractObjectCCOTree *SubtreeReplacer::replaceSegments(long long int saveInterval, string tempDirectory, string subtreeFilename){
	if (!allowThisClass) {
		cout << "FATAL: experimental class, set bool 'allowThisClass' to true to use this" << endl;
		exit(1);
	}

	// pass parameters, longTreeList.cco, shortTreeList.cco, percentages

	vector<SingleVessel *> parentVessels = getFilteredVessels();
	int itCount = 0;

	/// TODO: for each (SingleVessel *) vessel
	/// TODO: sort type of tree
	cout << "WARNING: limiting max iterations to " << maxIterationsLimit << endl;
	for (vector<SingleVessel *>::iterator it = parentVessels.begin(); it != parentVessels.end() && itCount<maxIterationsLimit; ++it, ++itCount) {
		SingleVessel* parentTerminal = (*it);
		vtkIdType parentID = parentTerminal->vtkSegmentId;
		/// TODO: get properties, get distal (coordinates), get radius, length
		point subtreeProx = toAppendVesselData[parentID][0];
		point subtreeDist = toAppendVesselData[parentID][1];

		SingleVesselCCOOTree *newSubtree = buildNewSubtree(subtreeFilename);
		

		mapSubtree(newSubtree, subtreeProx, subtreeDist);
		// Now the subtree is geometrically located in the correct point. Time to replace the subtree.


		/// TODO: make subtree and append
		// Read CCO, map root, and childs recursively
		// map proximal and distal of subtrees, recursively for every child.
		// update radius, update tree
		int newSegments = 101;
		/// TODO: smarter way to pass nTerms, maybe with:
		// int qtyTerms = newSubtree->getNTerms();
		// but for the amount of segments, not terminals;
		cout << "WARNING: hardcode for " << newSegments << " new segments in subtree" << endl;
		tree->addSubtree(newSubtree, parentTerminal, newSegments);

		
		// append tree, use SVCCOOT::addSubtree() method.
		// transfer ownership, use addChild, removeChildren, setParent, 

		// update the tree, terms, VTK_ID, etc.

		delete newSubtree;
	}
	


	cout << "iterated through all vessels" << endl;


    // Update tree
	cout << "updating the tree" << endl;
	((SingleVesselCCOOTree*) tree)->updateMassiveTree();
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());
	saveStatus(nTerminals-1);
	cout << "tree updated" << endl;

	return tree;

}

int SubtreeReplacer::loadData(string filename) {
	// Load the segment ID of parent + coordinates for subtree, to find the correct when appending the tree.
	// This was done via a new method that returns the ID via parameter. Now the function:
	/// I need to save the @param toAppendVesselData to a .txt file or binary...
	ifstream inStream;
	inStream.open(filename, ios::in);
	vtkIdType vesselID;
	double x1, y1, z1, x2, y2, z2;
	if (!inStream.is_open()) {
		cout << "ERROR: file was not open!" << endl;
		return 1;
	}
	while (inStream >> vesselID >> x1 >> y1 >> z1 >> x2 >> y2 >> z2) {
		point pProx;
		pProx.p[0] = x1;
		pProx.p[1] = y1;
		pProx.p[2] = z1;
		point pDist;
		pDist.p[0] = x2;
		pDist.p[1] = y2;
		pDist.p[2] = z2;
		this->toAppendVesselData[vesselID] = vector<point> {pProx, pDist};
	}
	inStream.close();
	return 0;
}

int SubtreeReplacer::isValidSegment(point xNew, int iTry) {

	if (iTry % instanceData->nTerminalTrial == 0) {
		dLim *= instanceData->dLimReductionFactor;
		cout << "DLim reduced." << endl;
	}
	point xBif;
	double dist;
	tree->getClosestTreePoint(xNew, &xBif, &dist);
//	cout << iTry << ": Closest point to xNew=" << xNew << " is " << xBif << " at " << dist << " (distance limit " << dLim << ")" << endl;

	return dist > dLim;
}

StagedDomain * SubtreeReplacer::getDomain() {
	return domain;
}

vector<AbstractObjectCCOTree *> SubtreeReplacer::getTrees() {
	vector<AbstractObjectCCOTree *> trees;
	trees.push_back(tree);
	return trees;
}

void SubtreeReplacer::enableConfigurationFile(
		string filename) {
	this->isGeneratingConfFile = 1;
	this->confFilename = filename;
}

void SubtreeReplacer::generatesConfigurationFile(ios::openmode mode) {

	if (isGeneratingConfFile) {
		confFile.open(confFilename.c_str(), mode);
		confFile.setf(ios::scientific, ios::floatfield);
		confFile.precision(16);

		confFile << "*GlobalParameters" << endl;
		confFile << "N_LEVELS_SCALING_TEST " << instanceData->nLevelTest << endl;
		confFile << "N_TRIAL " << instanceData->nTerminalTrial << endl;
		confFile << "DLIM_REDUCTION_FACTOR " << instanceData->dLimReductionFactor << endl;
		confFile << "PERFUSION_AREA_FACTOR " << instanceData->perfusionAreaFactor << endl;
		confFile << "CLOSE_NEIGHBORHOOD_FACTOR " << instanceData->closeNeighborhoodFactor << endl;
		confFile << "MIDPOINT_DLIM_FACTOR " << instanceData->midPointDlimFactor << endl;
		confFile << "N_BIF_TRIES " << instanceData->nBifurcationTest << endl;
		confFile << endl;

		confFile << "*TreeParameters" << endl;
		confFile << "TREE_CLASS " << this->tree->getTreeName() << endl;
		confFile << "N_TERMINALS " << this->nTerminals << endl;
		confFile << "INPUT_POSITION " << this->tree->getXProx() << endl;
		confFile << "INPUT_FLOW " << this->tree->getQProx() << endl;
		confFile << "PRESSURE_DROP " << this->tree->getDp() << endl;
		confFile << endl;

		confFile << "*Timestamps" << endl;
		beginningTime = chrono::steady_clock::now();
	}
}

void SubtreeReplacer::markTimestampOnConfigurationFile(
		string label) {
	if (isGeneratingConfFile) {
		confFile << (chrono::duration_cast<chrono::microseconds>(
				chrono::steady_clock::now() - beginningTime).count())
				/ 1000000.0 << ": " << label << endl;
	}
}

void SubtreeReplacer::closeConfigurationFile() {

	confFile << endl << "DOMAIN_POINTS_GENERATED " << this->domain->getPointCounter() << endl;
	if (isGeneratingConfFile) {
		confFile.flush();
		confFile.close();
	}
}

SingleVesselCCOOTree*& SubtreeReplacer::getTree() {
	return this->tree;
}

void SubtreeReplacer::setSavingTasks(const vector<AbstractSavingTask*>& savingTasks){
	this->savingTasks = savingTasks;
}

void SubtreeReplacer::saveStatus(long long int terminals){
	tree->setPointCounter(domain->getPointCounter());
	for (std::vector<AbstractSavingTask *>::iterator it = savingTasks.begin(); it != savingTasks.end(); ++it) {
		(*it)->execute(terminals,tree);
	}
//			nodalWriter->write(tempDirectory+ "/step" + to_string(i) + "_view.vtp",tree);
	markTimestampOnConfigurationFile("Generating vessel #" + to_string(terminals));
	markTimestampOnConfigurationFile("Total RAM consumption: " + to_string(monitor->getProcessMemoryConsumption()) + " MB.");
}

void SubtreeReplacer::observableModified(IDomainObservable* observableInstance) {
	cout << "Changing instance parameters from " << endl << instanceData;
	instanceData = ((AbstractDomain *) observableInstance)->getInstanceData();
	cout << "To " << endl << instanceData << endl;
	tree->setInstanceData(instanceData);

	SingleVessel::bifurcationTests = instanceData->nBifurcationTest;
	this->stage = ((StagedDomain *) observableInstance)->getCurrentStage();
	cout << "Changing to stage " << stage << endl;
	this->tree->setCurrentStage(this->stage);
	this->tree->setGam(gams[stage]);
	this->tree->setEpsLim(epsLims[stage]);
	this->tree->setNu(nus[stage]);

	this->dataMonitor->reset();
}

vector<AbstractConstraintFunction<double, int> *>* SubtreeReplacer::getGams()
{
	return &(this->gams);
}

vector<AbstractConstraintFunction<double, int> *>* SubtreeReplacer::getEpsLims()
{
	return &(this->epsLims);
}
vector<AbstractConstraintFunction<double, int> *>* SubtreeReplacer::getNus()
{
	return &(this->nus);
}

double SubtreeReplacer::getDLim() {
	return this->dLim;
}

void SubtreeReplacer::setDLim(double newDLim) {
	this->dLim = newDLim;
}

double SubtreeReplacer::getDLimInitial() {
	return this->dLimInitial;
}

double SubtreeReplacer::getDLimLast() {
	return this->dLimLast;
}

time_t SubtreeReplacer::getBeginTime() {
	return this->beginTime;
}

time_t SubtreeReplacer::getEndTime() {
	return this->endTime;
}


// int PenetratingVesselTreeGenerator::isValidRootSegment(point xNew,
// 		int iTry) {
// 	if (iTry % instanceData->nTerminalTrial == 0) {
// 		dLim *= instanceData->dLimReductionFactor;
// //		cout << "DLim reduced." << endl;
// 	}
// 	if( !(domain->isSegmentInside(tree->getXProx(),xNew)) )
// 		return false;
// 	point dVect = xNew - tree->getXProx();
// 	return sqrt(dVect ^ dVect) > dLim;
// }

// AbstractObjectCCOTree *PenetratingVesselTreeGenerator::resume(long long int saveInterval, string tempDirectory) {

// 	this->beginTime = time(nullptr);
// 	this->dLimInitial = this->dLim;

// //	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
// 	generatesConfigurationFile(ios::out);

// 	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
// 		domain->getRandomPoint();
// 	}

// 	//	Compute current nTerm
// 	long long currentTerminals = tree->getNTerms();
// 	point xNew;

// 	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
// 	//	Be careful nTerminals may differ from the current amount of terminals since vessel-tip conexions are allowed in some cases.
// 	for (long long i = currentTerminals; i < nTerminals; i = tree->getNTerms()) {

// 		dataMonitor->update();

// 		if (i % saveInterval == 0 || i == currentTerminals) {
// 			saveStatus(i);
// 		}

// 		int invalidTerminal = true;
// 		int iTry = 0;
// 		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
// 		while (invalidTerminal) {

// 			do {
// 				xNew = domain->getRandomPoint();
// 			} while (!isValidSegment(xNew, ++iTry));

// 			int nNeighbors;
// 			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
// 			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

// 			double minCost = INFINITY;
// 			point minBif;
// 			AbstractVascularElement *minParent = NULL;
// #pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
// 			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
// 				point xBif;
// 				double cost;
// 				tree->testVessel(xNew, neighborVessels[j], domain,
// 						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
// #pragma omp critical
// 				{
// 					if (cost < minCost) {
// 						minCost = cost;
// 						minBif = xBif;
// 						minParent = neighborVessels[j];
// 					}
// 				}
// 			}
// 			//	end for trees

// 			if (minCost < INFINITY) {
// 				cout << "Added with a cost of " << minCost << " with a total cost of " << ((SingleVessel *) tree->getRoot())->treeVolume << endl;
// 				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction,
// 								(AbstractVascularElement::BRANCHING_MODE) instanceData->branchingMode);
// 				invalidTerminal = false;
// 			}
// 		}
// 		dataMonitor->addDLimValue(dLim,i);
// 		domain->update();
// 		//	TODO Reset data monitor!!
// 		//	Terminal added.
// 	}
// 	tree->computePressure(tree->getRoot());
// 	tree->setPointCounter(domain->getPointCounter());

// 	this->endTime = time(nullptr);
// 	this->dLimLast = this->dLim;

// 	saveStatus(nTerminals-1);
// 	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
// 	markTimestampOnConfigurationFile("Tree successfully generated.");
// 	closeConfigurationFile();

// 	return tree;

// }

// static void originalVesselsRecursive(SingleVessel *root, unordered_set<vtkIdType>* originals, vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints, long long int &partTerm) {
// 	if(!root) {
// 		return;
// 	}
// 	point midPoint = (root->xProx + root->xDist) / 2.;
// 	if (enclosedPoints->IsInsideSurface(midPoint.p)) {
// 		originals->insert(root->vtkSegmentId);
// 		if ((root->children).empty()) {
// 			++partTerm;
// 		}
// 	}
// 	for(auto it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
// 		originalVesselsRecursive(static_cast<SingleVessel *>((*it)), originals, enclosedPoints, partTerm);
// 	}
// }


// static unordered_set<vtkIdType>* originalVessels(SingleVesselCCOOTree *tree, vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints, long long int &partTerm) {
// 	unordered_set<vtkIdType>* originalSV = new unordered_set<vtkIdType>;
// 	originalVesselsRecursive(static_cast<SingleVessel *>(tree->getRoot()), originalSV, enclosedPoints, partTerm);
// 	return originalSV;
// }

// AbstractObjectCCOTree *PenetratingVesselTreeGenerator::resumeSavePointsMidPoint(long long int saveInterval, string tempDirectory, FILE *fp) {
// 	// This function checks for vessels differently, instead of checking if parent is in domain/region/partition, it checks if
// 	// its midpoint is inside the region, and then considers it as a candidate
// 	// Does not discard parent vessel if parameter
// 	/// @param bypassFunctionIfMidpointInside
// 	// is set to true (default: false)

// 	this->beginTime = time(nullptr);
// 	this->dLimInitial = this->dLim;

// 	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = domain->getEnclosedPoints();
// 	/* It holds all the segments such that the midpoint is inside the partition. */
// 	long long int partTerm {0};
// 	unordered_set<vtkIdType> *partVessels = originalVessels(static_cast<SingleVesselCCOOTree *>(this->tree), enclosedPoints, partTerm);
// 	// Assuming that all terminals have the same outflow.
// 	const double qOutsideTerm {this->tree->getQProx() / this->tree->getNTerms()};
// 	const double qPart {partTerm * qOutsideTerm};
// 	const vector<double> qPartitions {qPart, qOutsideTerm};

// //	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
// 	generatesConfigurationFile(ios::out);

// 	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
// 		domain->getRandomPoint();
// 	}

// 	//	Compute current nTerm
// 	long long currentTerminals = tree->getNTerms();
// 	point xNew;

// 	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
// 	//	Be careful nTerminals may differ from the current amount of terminals since vessel-tip conexions are allowed in some cases.
// 	for (long long i = currentTerminals; i < nTerminals; i = tree->getNTerms()) {

// 		dataMonitor->update();

// 		if (i % saveInterval == 0 || i == currentTerminals) {
// 			saveStatus(i);
// 		}

// 		int invalidTerminal = true;
// 		int iTry = 0;
// 		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
// 		while (invalidTerminal) {

// 			do {
// 				xNew = domain->getRandomPoint();
// 			} while (!isValidSegment(xNew, ++iTry));

// 			int nNeighbors;
// 			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
// 			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

// 			double minCost = INFINITY;
// 			point minBif;
// 			AbstractVascularElement *minParent = NULL;
// 			((SingleVesselCCOOTree*) tree)->bypassFunctionIfMidpointInside = this->bypassFunctionIfMidpointInside;
// #pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
// 			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
// 				point xBif;
// 				double cost;
// 				if (partVessels->find(static_cast<SingleVessel *>(neighborVessels[j])->vtkSegmentId) == partVessels->end()) {
// 					continue;
// 				}
// 				tree->testVessel(xNew, neighborVessels[j], domain,
// 						neighborVessels, dLim, &xBif, &cost, partVessels, &partTerm, qPartitions); //	Inf cost stands for invalid solution
// #pragma omp critical
// 				{
// 					if (cost < minCost) {
// 						minCost = cost;
// 						minBif = xBif;
// 						minParent = neighborVessels[j];
// 					}
// 				}
// 			}
// 			//	end for trees

// 			if (minCost < INFINITY) {
// 				cout << "Added with a cost of " << minCost << " with a total cost of " << ((SingleVessel *) tree->getRoot())->treeVolume << endl;
// 				SingleVessel *minParentSV = (SingleVessel *) minParent;
// 				fwrite(&(minBif.p[0]), sizeof(double), 1, fp);
// 				fwrite(&(minBif.p[1]), sizeof(double), 1, fp);
// 				fwrite(&(minBif.p[2]), sizeof(double), 1, fp);
// 				fwrite(&(xNew.p[0]), sizeof(double), 1, fp);
// 				fwrite(&(xNew.p[1]), sizeof(double), 1, fp);
// 				fwrite(&(xNew.p[2]), sizeof(double), 1, fp);
// 				fwrite(&(minParentSV->xProx.p[0]), sizeof(double), 1, fp);
// 				fwrite(&(minParentSV->xProx.p[1]), sizeof(double), 1, fp);
// 				fwrite(&(minParentSV->xProx.p[2]), sizeof(double), 1, fp);
// 				fwrite(&(minParentSV->xDist.p[0]), sizeof(double), 1, fp);
// 				fwrite(&(minParentSV->xDist.p[1]), sizeof(double), 1, fp);
// 				fwrite(&(minParentSV->xDist.p[2]), sizeof(double), 1, fp);
// 				fwrite(&(instanceData->vesselFunction), sizeof(int), 1, fp);
// 				fwrite(&(this->stage), sizeof(int), 1, fp);
// 				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction,
// 								(AbstractVascularElement::BRANCHING_MODE) instanceData->branchingMode, partVessels, &partTerm, qPartitions);
// 				invalidTerminal = false;
// 			}
// 		}
// 		dataMonitor->addDLimValue(dLim,i);
// 		domain->update();
// 		//	TODO Reset data monitor!!
// 		//	Terminal added.
// 	}
// 	tree->computePressure(tree->getRoot());
// 	tree->setPointCounter(domain->getPointCounter());

// 	this->endTime = time(nullptr);
// 	this->dLimLast = this->dLim;

// 	partVessels->clear();
// 	delete partVessels;

// 	saveStatus(nTerminals-1);
// 	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
// 	markTimestampOnConfigurationFile("Tree successfully generated.");
// 	closeConfigurationFile();

// 	return tree;
// }


// void PenetratingVesselTreeGenerator::projectTerminals(vector<SingleVessel *> vessels){
// 	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
// 		point terminal = (*it)->xDist;
// 		point projection,normal;
// 		vtkIdType closeCellId;
// 		int subId;
// 		double distance;
// 		locator->FindClosestPoint(terminal.p,projection.p,closeCellId,subId,distance);

// 		//	Get normal of the projected element - // TESTME
// 		vtkSmartPointer<vtkDataArray> cellNormalsRetrieved = vtkGeometry->GetCellData()->GetNormals();
// 		cellNormalsRetrieved->GetTuple(closeCellId,normal.p);

// 		point displacement = projection-terminal;
// 		//	Displacement versor
// 		displacement = displacement / sqrt(displacement^displacement);
// 		//	Projection + offset - Inner product checks that the offset is applied towards the geometry interior.
// 		if((normal^displacement)<0)
// 			projection = projection + displacement * offset;
// 		else
// 			projection = projection - displacement * offset;
// 		(*it)->xDist = projection;	//	Need to update the VTK segment!
// 	}

// }


// void PenetratingVesselTreeGenerator::projectVessel(vector<SingleVessel *> vessels){
// 	for (std::vector<SingleVessel *>::iterator it = vessels.begin(); it != vessels.end(); ++it) {
// 		point distal = (*it)->xDist;
// 		point proximal = (*it)->xProx;
// 		point projectionProx, projectionDist, normalProx, normalDist,normal;
// 		vtkIdType closeCellIdProx;
// 		vtkIdType closeCellIdDist;
// 		int subId;
// 		double distance;
// 		locator->FindClosestPoint(distal.p,projectionProx.p,closeCellIdProx,subId,distance);
// 		locator->FindClosestPoint(distal.p,projectionDist.p,closeCellIdDist,subId,distance);

// 		//	Get normal of the projected element
// 		vtkSmartPointer<vtkDoubleArray> cellNormalsRetrieved = vtkDoubleArray::SafeDownCast(vtkGeometry->GetCellData()->GetNormals());
// 		cellNormalsRetrieved->GetTuple(closeCellIdProx,normalProx.p);
// 		cellNormalsRetrieved->GetTuple(closeCellIdDist,normalDist.p);
// 		normal = (normalProx + normalDist)/2;

// 		point displacementDist = projectionDist-distal;
// 		point displacementProx = projectionProx-proximal;
// 		point displacement = (displacementDist + displacementProx ) / 2;
// 		//	Displacement versor
// 		displacement = displacement / sqrt(displacement^displacement);
// 		//	Projection + offset - Inner product checks that the offset is applied towards the geometry interior.
// 		if((normal^displacement)<0){
// 			projectionDist = projectionDist + displacement * offset;
// 			projectionProx = projectionProx + displacement * offset;
// 		}
// 		else{
// 			projectionDist = projectionDist - displacement * offset;
// 			projectionProx = projectionProx - displacement * offset;
// 		}
// 		(*it)->xDist = projectionDist;	//	Need to update the VTK segment!
// 		(*it)->xProx = projectionProx;  //	Need to update the VTK segment!
// 	}

// }


// void PenetratingVesselTreeGenerator::setProjectionDomain(string projectionDomainFile){
// 	//	Read all the data from the file
// 	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
// 	reader->SetFileName(projectionDomainFile.c_str());
// 	reader->Update();
// 	vtkGeometry = reader->GetOutput();

// 	//	Generate normals for the geometry
// 	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
// 	normalGenerator->SetInputData(vtkGeometry);
// 	normalGenerator->ComputePointNormalsOff();
// 	normalGenerator->ComputeCellNormalsOn();
// 	normalGenerator->Update();
// 	vtkGeometry = normalGenerator->GetOutput();

// 	//	Create the tree
// 	locator = vtkSmartPointer<vtkCellLocator>::New();
// 	locator->SetDataSet(vtkGeometry);
// 	locator->BuildLocator();
// }

// void PenetratingVesselTreeGenerator::setDescendingOffset(double descendingOffset){
// 	this->descendingOffset = descendingOffset;
// }
// void PenetratingVesselTreeGenerator::setPenetrationOffset(double penetrationOffset){
// 	this->penetrationOffset = penetrationOffset;
// }




// void standardCCOPaperCase(string inputDir, string prefix, string outputDir){

// 	string modelsFolder = "./";
// 	int nTree = 1;
// 	// point *xi = new point[nTree];
// 	// xi[0] = {28.98,-1517.55, -140.01}; // in mm
// 	// double *ri = new double[nTree];
// 	// ri[0] = 0.198;		//	in mm
// 	double *qi = new double[nTree];
// 	qi[0] = 2;				//	in mm^3/s
// 	// int terminal1 = 10;
// 	// int terminal2 = 10;
// 	// int terminal3 = 10;
// 	// int terminal4 = 10;
// 	int terminals = 10; //terminal1 + terminal2 + terminal3 + terminal4;

// 	GeneratorData *instanceData = new GeneratorData(32000, 500, 0.9, 1.0, 2, 0.25, 7,
// 			AbstractVascularElement::VESSEL_FUNCTION::DISTRIBUTION,false,new VolumetricCostEstimator());
// 	SimpleDomain *domain = new SimpleDomain(modelsFolder + "sphereR1.vtk", 10000l, 5, instanceData);

// 	StagedDomain *stagedDomain = new StagedDomain();
// 	stagedDomain->addStage(terminals, domain);

// 	AbstractConstraintFunction<double, int> *gam = new ConstantConstraintFunction<double, int>(3);
// 	AbstractConstraintFunction<double, int> *epsLim = new ConstantPiecewiseConstraintFunction<double, int>( { 0.4, 0.0 }, { 0, 5 });
// 	AbstractConstraintFunction<double, int> *nu = new ConstantConstraintFunction<double, int>(0.036);

// 	SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(modelsFolder + "toProject.cco",
// 			instanceData, qi[0], gam, epsLim, nu, 0.0, 1E-5);

// 	cout << "Connecting arteries to the cortex surface..." << endl;
// 	TreeProjector *projector = new TreeProjector(modelsFolder + "sphereR1.vtk", 0.05);
// 	AbstractVesselFilter *filter = new VesselFilterComposite({new VesselFilterByBranchingMode(AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING)});
// 	vector<SingleVessel *> treeVessels = tree->getVessels();
// 	vector<SingleVessel *> terminalVessels = filter->apply(treeVessels);
// 	projector->projectTerminals(terminalVessels);

// 	cout << "Trying to generate the tree..." << endl;

// 	StagedFRROTreeGenerator * treeGenerator = new StagedFRROTreeGenerator(stagedDomain, tree,
// 			terminals, {gam}, {epsLim}, {nu});

// 	 treeGenerator->setSavingTasks({new VisualizationSavingTask(outputDir,prefix + "_vis_TERM" + to_string(terminals)+"_step"),
// 	 	new CheckpointSavingTask(outputDir,prefix + "_chk_TERM" + to_string(terminals)+"_step")});

// 	treeGenerator->enableConfigurationFile(
// 			outputDir + prefix + "_TERM" + to_string(terminals) + "TRIES" + to_string(instanceData->nTerminalTrial) + "FNEIGHBOR" + to_string(
// 					instanceData->closeNeighborhoodFactor) + "BIF_TR" + to_string(instanceData->nBifurcationTest) + ".cfg");
// 	tree = (SingleVesselCCOOTree *) treeGenerator->resume(20,outputDir);
// 	cout << "Constructor executed." << endl;

// 	cout << "Writing "
// 			<< outputDir + prefix + "_TERM" + to_string(terminals) + "TRIES" + to_string(instanceData->nTerminalTrial) + "FNEIGHBOR" + to_string(
// 					instanceData->closeNeighborhoodFactor) + "BIF_TR" + to_string(instanceData->nBifurcationTest) + ".vtp" << endl;

// 	VTKObjectTreeNodalWriter *nodalWriter = new VTKObjectTreeNodalWriter();
// 	nodalWriter->write(outputDir + prefix + "_nodal_TERM" +  to_string(terminals) + "TRIES" + to_string(instanceData->nTerminalTrial) + "FNEIGHBOR" + to_string(
// 			instanceData->closeNeighborhoodFactor) + "BIF_TR" + to_string(instanceData->nBifurcationTest) + ".vtp",tree);

// 	VTKObjectTreeSplinesNodalWriter *splinesWriter = new VTKObjectTreeSplinesNodalWriter();
// 	splinesWriter->write(outputDir + prefix + "_view_TERM" + to_string(terminals) + "TRIES" + to_string(instanceData->nTerminalTrial) + "FNEIGHBOR" + to_string(
// 			instanceData->closeNeighborhoodFactor) + "BIF_TR" + to_string(instanceData->nBifurcationTest) + ".vtp",tree);
// 	((SingleVesselCCOOTree *)tree)->save(
// 			outputDir + prefix + "_TERM" + to_string(terminals) + "TRIES" + to_string(instanceData->nTerminalTrial) + "FNEIGHBOR" + to_string(
// 					instanceData->closeNeighborhoodFactor) + "BIF_TR" + to_string(instanceData->nBifurcationTest) + ".cco");
// 	//	FIXME Potential bug, the program crash before writing vtk files.

// }
// // standardCCOPaperCase("./","original", "./");





// string coordToString(const point &xProx,const point &xDist) {
//     double coordArray[6] = {xProx.p[0], xProx.p[1], xProx.p[2],
// 		xDist.p[0], xDist.p[1], xDist.p[2]};
//     char *coordCString = (char *) malloc(6 * sizeof(double));
//     memcpy(coordCString, &coordArray, 6 * sizeof(double));
//     string coordString(coordCString, (6 * sizeof(double) / sizeof(char)));
//     free(coordCString);
//     return coordString;
// }


// void TreeMerger::createMapping(SingleVessel *vessel) {
//     if(!vessel) {
//         return;
//     }
//     printf("Added vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
//         vessel->xProx.p[0], vessel->xProx.p[1], vessel->xProx.p[2],
//         vessel->xDist.p[0], vessel->xDist.p[1], vessel->xDist.p[2]);
//     string coordString = vessel->coordToString();
//     pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert = this->stringToPointer->insert(pair<string, SingleVessel *>(coordString, vessel));
//     if(!didInsert.second) {
//         printf("Failed to add element!\n");
//     }
//     printf("stringToPoint.size = %lu\n", this->stringToPointer->size());
//     for (auto it = vessel->children.begin(); it != vessel->children.end(); ++it) {
//         SingleVessel *child = (SingleVessel *) (*it);
//         createMapping(child);
//     }
// }




// void tTreeMerger(SingleVesselCCOOTree *baseTree, vector<string>& derivedTreePoints) {

//         this->tree = baseTree;

//         this->skipFailedMerges = false;

//         this->vesselToMerge = new vector<vector<ReadData> *>;

//         for (auto it = derivedTreePoints.begin(); it != derivedTreePoints.end(); ++it) {

//             FILE *fp = fopen((*it).c_str(), "rb");
//             if (!fp) {
//                 fprintf(stderr, "Failed to open derived tree file!\n");
//                 exit(EXIT_FAILURE);
//             }

//             vector<ReadData> *toMerge = new vector<ReadData>;
//             this->vesselToMerge->push_back(toMerge);

//             ReadData readLine;
//             double tempPoints[12];
//             int tempFunction;
//             int tempStage;
//             fread(&tempPoints, sizeof(double), 12, fp);
//             fread(&tempFunction, sizeof(int), 1, fp);
//             fread(&tempStage, sizeof(int), 1, fp);
//             readLine.xBif.p[0] = tempPoints[0];
//             readLine.xBif.p[1] = tempPoints[1];
//             readLine.xBif.p[2] = tempPoints[2];
//             readLine.xNew.p[0] = tempPoints[3];
//             readLine.xNew.p[1] = tempPoints[4];
//             readLine.xNew.p[2] = tempPoints[5];
//             readLine.xPProx.p[0] = tempPoints[6];
//             readLine.xPProx.p[1] = tempPoints[7];
//             readLine.xPProx.p[2] = tempPoints[8];
//             readLine.xPDist.p[0] = tempPoints[9];
//             readLine.xPDist.p[1] = tempPoints[10];
//             readLine.xPDist.p[2] = tempPoints[11];
//             readLine.function = static_cast<AbstractVascularElement::VESSEL_FUNCTION>(tempFunction);
//             readLine.stage = tempStage;
//             while (!feof(fp)) {
//                 toMerge->push_back(readLine);
//                 fread(&tempPoints, sizeof(double), 12, fp);
//                 fread(&tempFunction, sizeof(int), 1, fp);
//                 fread(&tempStage, sizeof(int), 1, fp);
//                 readLine.xBif.p[0] = tempPoints[0];
//                 readLine.xBif.p[1] = tempPoints[1];
//                 readLine.xBif.p[2] = tempPoints[2];
//                 readLine.xNew.p[0] = tempPoints[3];
//                 readLine.xNew.p[1] = tempPoints[4];
//                 readLine.xNew.p[2] = tempPoints[5];
//                 readLine.xPProx.p[0] = tempPoints[6];
//                 readLine.xPProx.p[1] = tempPoints[7];
//                 readLine.xPProx.p[2] = tempPoints[8];
//                 readLine.xPDist.p[0] = tempPoints[9];
//                 readLine.xPDist.p[1] = tempPoints[10];
//                 readLine.xPDist.p[2] = tempPoints[11];
//                 readLine.function = static_cast<AbstractVascularElement::VESSEL_FUNCTION>(tempFunction);
//                 readLine.stage = tempStage;
//             }
//             fclose(fp);
//         }

//         // for (auto it1 = vesselToMerge->begin(); it1 != vesselToMerge->end(); ++it1) {
//         //     for (auto it2 = (*it1)->begin(); it2 != (*it1)->end(); ++it2) {
//         //         printf("xBif = (%.16e %.16e %.16e)\txNew = (%.16e %.16e %.16e)\txProx = (%.16e %.16e %.16e)\t xDist = (%.16e %.16e %.16e)\n",
//         //             (*it2).xBif.p[0], (*it2).xBif.p[1], (*it2).xBif.p[2],
//         //             (*it2).xNew.p[0], (*it2).xNew.p[1], (*it2).xNew.p[2],
//         //             (*it2).xPProx.p[0], (*it2).xPProx.p[1], (*it2).xPProx.p[2],
//         //             (*it2).xPDist.p[0], (*it2).xPDist.p[1],(*it2).xPDist.p[2]);
//         //     }
//         // }

//         this->stringToPointer = new unordered_map<string, SingleVessel *>();
//         this->createMapping((SingleVessel *) this->tree->getRoot());

//         // printf("The mapping:\n");
//         // for (auto it = this->stringToPointer->begin(); it != this->stringToPointer->end(); ++it) {
//         //     printf("key = %s\n, pointer = %p\n", (*it).first.c_str(), (*it).second);
//         // }
// }

// TreeMerger::~TreeMerger() {
//     delete this->stringToPointer;
//     for (auto it = this->vesselToMerge->begin(); it != this->vesselToMerge->end(); ++it) {
//         delete (*it);
//     }
//     delete this->vesselToMerge;
// }

// SingleVesselCCOOTree* PenetratingVesselTreeGenerator::mergeFast() {
//     for (auto itTree = vesselToMerge->begin(); itTree != vesselToMerge->end(); ++itTree) {
//         for (auto itVessels = (*itTree)->begin(); itVessels != (*itTree)->end(); ++itVessels) {
//             printf("Trying to access parent at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
//                 (*itVessels).xPProx.p[0], (*itVessels).xPProx.p[1], (*itVessels).xPProx.p[2],
//                 (*itVessels).xPDist.p[0], (*itVessels).xPDist.p[1], (*itVessels).xPDist.p[2]);
//             // printf("key_accessed = %s\n", coordToString((*itVessels).xPProx, (*itVessels).xPDist).c_str());
//             try {
//                 // A junção
//                 SingleVessel *parentPointer = this->stringToPointer->at(coordToString((*itVessels).xPProx, (*itVessels).xPDist));
//                 tree->addVesselMergeFast((*itVessels).xBif, (*itVessels).xNew, parentPointer, (*itVessels).function, (*itVessels).stage, this->stringToPointer);
//             }
//             catch(std::out_of_range const&) {
//                 std::cout << "This merge failed!" << std::endl;
//                 if(this->skipFailedMerges) continue;
//                 throw; // default behaviour to throw the exception
//             }
//         }
//     }

//     // Update tree
//     this->tree->updateTree(((SingleVessel *) this->tree->getRoot()), this->tree);

// 	double maxVariation = INFINITY;
// 	while (maxVariation > this->tree->variationTolerance) {
// 			this->tree->updateTreeViscositiesBeta(((SingleVessel *) this->tree->getRoot()), &maxVariation);
// 	}

//     this->tree->computePressure(this->tree->getRoot());

//     return this->tree;
// }




// typedef struct {
//     ReadData *readData;
//     size_t noVessels;
//     size_t cur;
//     bool isDone;
// } ReadDataArray;

// typedef struct {
//     ReadDataArray *vessels;
//     size_t size;
// } TreesArray;

// bool didEnd(const TreesArray &treesArray) {
//     bool didEnd = true;
//     for (size_t i = 0; i < treesArray.size; ++i) {
//         didEnd = didEnd && (treesArray.vessels)[i].isDone;
//     }
//     return didEnd;
// }

// SingleVesselCCOOTree* TreeMerger::merge() {
//     TreesArray treesArray;
//     treesArray.size = this->vesselToMerge->size();
//     treesArray.vessels = new ReadDataArray[treesArray.size];
//     for (size_t i = 0; i < treesArray.size; ++i) {
//         (treesArray.vessels)[i].readData = this->vesselToMerge->at(i)->data();
//         (treesArray.vessels)[i].noVessels = this->vesselToMerge->at(i)->size();
//         (treesArray.vessels)[i].cur = 0;
//         (treesArray.vessels)[i].isDone = false;
//     }

//     size_t i = 0;
//     while(!didEnd(treesArray)) {
//         if (!(treesArray.vessels)[i].isDone) {
//             ReadData curData = ((treesArray.vessels)[i].readData)[(treesArray.vessels)[i].cur];
//             printf("Trying to access parent at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
//                 curData.xPProx.p[0], curData.xPProx.p[1], curData.xPProx.p[2],
//                 curData.xPDist.p[0], curData.xPDist.p[1], curData.xPDist.p[2]);
//             // printf("key_accessed = %s\n", coordToString((*itVessels).xPProx, (*itVessels).xPDist).c_str());
//             SingleVessel *parentPointer = this->stringToPointer->at(coordToString(curData.xPProx, curData.xPDist));
//             tree->addVesselMerge(curData.xBif, curData.xNew, parentPointer, curData.function, curData.stage, this->stringToPointer);
//             if ((treesArray.vessels)[i].cur + 1 == (treesArray.vessels)[i].noVessels) {
//                 (treesArray.vessels)[i].isDone = true;
//             }
//             ++((treesArray.vessels)[i].cur);
//         }
//         i = (i + 1) % treesArray.size;
//     }

//     this->tree->updateTree((SingleVessel *) this->tree->getRoot(), this->tree);
// 	double maxVariation = INFINITY;
// 	while (maxVariation > this->tree->variationTolerance) {
// 		this->tree->updateTreeViscositiesBeta(((SingleVessel *) this->tree->getRoot()), &maxVariation);
//     }
//     this->tree->computePressure(this->tree->getRoot());

//     delete treesArray.vessels;
//     return this->tree;
// }