/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2024 Eduardo Guerreiro Zilves & Gonzalo Maso Talou */
/*
 * PenetratingVesselTreeGenerator.cpp
 *
 *  Created on: Feb 06, 2024
 *      Author: Eduardo G. Zilves
 */

#include "PenetratingVesselTreeGenerator.h"

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


PenetratingVesselTreeGenerator::PenetratingVesselTreeGenerator(
		StagedDomain* domain, string projectionDomainFile, AbstractObjectCCOTree* tree, long long nTerm,
		vector<AbstractConstraintFunction<double,int> *>gam,
		vector<AbstractConstraintFunction<double,int> *>epsLim,
		vector<AbstractConstraintFunction<double,int> *>nu) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	SingleVessel::bifurcationTests = instanceData->nBifurcationTest;
	this->nTerminals = nTerm;
	this->tree = tree;
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

	this->descendingOffset = 0;
	this->endpointOffset = 0;

	this->projectionDomainFile = projectionDomainFile;
	this->domainFile = projectionDomainFile;


	// Penetration parameters
	/// TODO: pass penetrationFactor as parameter
	this->descendingOffset = max<double>(descendingOffset, 1.0E-4);
	this->endpointOffset = max<double>(endpointOffset, 1.0E-4);
	this->maxDistanceToClosestPoint = 0.4; // cm
	this->maxPenetratingVesselLength = 0.25; //cm

	this->penetrationFactor = 1.0;
	this->maxGenerateLimit = 1000000; 

	// get final point of bifurcation per new segment. this technique uses a ray-cast and
	// finds the first intersection to a surface, and considers it as the end of domain.
	this->xRayDisplacement = 100.0;
	// find the endpoint for vessel
	this->intersectionTolerance = 1e-5;

	// set CTOR as terminal and generate, set 0.5 for midpoint and generate again
	this->parametricT = 1.0;

	// function and mode for the addvessel
	this->vesselFunction = (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction;
	this->branchingMode = (AbstractVascularElement::BRANCHING_MODE) instanceData->branchingMode;

	this->penetratingData.clear();

}

PenetratingVesselTreeGenerator::~PenetratingVesselTreeGenerator() {
	delete this->dataMonitor;
	delete this->monitor;
	if (this->didAllocateTree) {
		delete this->tree;
	}
}

inline void PenetratingVesselTreeGenerator::failsafeCheck() {
	if (!allowThisClass) {
		cout << "FATAL: experimental class, set bool 'allowThisClass' to true to use this" << endl;
		exit(1);
	}
	return;
}

vector<SingleVessel *> PenetratingVesselTreeGenerator::filterTerminalVessels() {
	// filter the bifurcable vessels, stage, distalbranching, etc.
	AbstractVesselFilter *terminalFilters = new VesselFilterComposite({new VesselFilterByTerminal()});
	vector<SingleVessel *> treeVessels = this->tree->getVessels();
	vector<SingleVessel *> terminalVessels = terminalFilters->apply(treeVessels);
	return terminalVessels;
}

void PenetratingVesselTreeGenerator::buildNormalLocator() {
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> projectionReader = vtkSmartPointer<vtkPolyDataReader>::New();
	projectionReader->SetFileName(this->projectionDomainFile.c_str());
	projectionReader->Update();
	this->vtkGeometryProjection = projectionReader->GetOutput();
	//	Generate normals for the geometry
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(vtkGeometryProjection);
	normalGenerator->ComputePointNormalsOff();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();
	this->vtkGeometryProjection = normalGenerator->GetOutput();
	//	Create the tree locator
	this->locatorProjection = vtkSmartPointer<vtkCellLocator>::New();
	this->locatorProjection->SetDataSet(vtkGeometryProjection);
	this->locatorProjection->BuildLocator();
	// prepare the normal geometries
	this->cellNormalsRetrieved = this->vtkGeometryProjection->GetCellData()->GetNormals();
}

void PenetratingVesselTreeGenerator::buildIntersectionLocator() {
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> readerIntersect = vtkSmartPointer<vtkPolyDataReader>::New();
	readerIntersect->SetFileName(this->domainFile.c_str());
	readerIntersect->Update();
	this->vtkGeometryIntersect = readerIntersect->GetOutput();
	//	Create the tree locator
	this->locatorIntersect = vtkSmartPointer<vtkCellLocator>::New();
	this->locatorIntersect->SetDataSet(vtkGeometryIntersect);
	this->locatorIntersect->BuildLocator();
}

point PenetratingVesselTreeGenerator::findDistalDescending(point terminal, point& normalT, double& distance2T) {
	// Find the endpoint for the terminal
	point projectionT;
	vtkIdType closeCellIdT;
	int subIdT;
	// double distance2T; // This is returned by reference
	// find closest cell from point of bifurcation
	this->locatorProjection->FindClosestPoint(terminal.p, projectionT.p, closeCellIdT, subIdT, distance2T); // the distance is squared!

	// get normal of the closest cell at closest point and the normal versor
	cellNormalsRetrieved->GetTuple(closeCellIdT, normalT.p);

	// define displacement versor for the offset (+-something%), get projection and add/subtract the displacement, always below surface.
	// this displacement is NOT the normal, it is in the direction from terminal to projection
	point displacementT = projectionT - terminal;
	displacementT = displacementT / sqrt(displacementT^displacementT);

	// add if same direction, subtract if going outwards
	// however, start from terminal if point is inside, to avoid zig-zag pattern
	// it checks if point is inside domain and then change how the first step is added, to avoid weird zig-zag patterns.
	bool terminalIsInside = false;
	if ((normalT^displacementT)<0) {
		projectionT = projectionT + displacementT * descendingOffset;
	} else { // here the point is inside the domain
		// projectionT = projectionT - displacementT * descendingOffset;
		projectionT = terminal - displacementT * descendingOffset;
		terminalIsInside = true;
	}

	return projectionT;
}

bool PenetratingVesselTreeGenerator::isDescendingValid(double distance2T) {
	// check if the distance is under a maximum value, to verify if there is gray matter below it
	// if not, then abort this vessel/point from bifurcating
	bool outputCheck = true;
	if (distance2T > pow(this->maxDistanceToClosestPoint, 2)) {
		cout << "WARNING: Terminal aborted, closest point beyond maximum distance." << endl;
		outputCheck = false;
	}
	return outputCheck;
}

point PenetratingVesselTreeGenerator::findDistalPenetrating(point terminal, point normal, int& foundRaycast) {
	/// NOTE: @param terminal is the xNew1/distal of the descending, the xBif/proximal of the penetrating.
	// get final point of bifurcation per new segment. this technique uses a ray-cast and
	// finds the first intersection to a surface, and considers it as the end of domain.
	point xRayDisplacementT = (normal*xRayDisplacement)*(-1);
	point xRaycastT = terminal + xRayDisplacementT;

	// find the endpoint for vessel
	double tParamT = 0.0;
	point hitpointT;
	point pcoordsT;
	vtkIdType endCellIdT;
	int endSubIdT = 0;

	foundRaycast = this->locatorIntersect->IntersectWithLine(terminal.p, xRaycastT.p, intersectionTolerance,
		tParamT, hitpointT.p, pcoordsT.p, endSubIdT, endCellIdT);
	
	point directionT = hitpointT - terminal;
	directionT = directionT / sqrt(directionT^directionT);

	// Check if segment is too long, and clamp it to a maximum value
	// point is always from inside the surface, we subtract to make the length shorter.
	point endpointT;
	endpointT = hitpointT - directionT * endpointOffset;

	point penetratingT = endpointT - terminal;
	double lengthT = sqrt(penetratingT^penetratingT);
	lengthT = adjustMaxPenetratingLength(lengthT);

	point xNew2 = terminal + (penetratingT/sqrt(penetratingT^penetratingT))*lengthT*penetrationFactor;

	return xNew2;
}

double PenetratingVesselTreeGenerator::adjustMaxPenetratingLength(double length) {
	if (length > maxPenetratingVesselLength){
		cout << "NOTE: Terminal shortened, length beyond maximum distance." << "\n";
		length = maxPenetratingVesselLength;
	}
	return length;
}

bool PenetratingVesselTreeGenerator::isPenetratingValid(int foundRaycast) {
	bool outputCheck = true;
	if (!foundRaycast){
		cout << "WARNING: no intersection found! Terminal aborted." << endl;
		outputCheck = false;
	}
	return outputCheck;
}

bool PenetratingVesselTreeGenerator::isPenetratingInside(point xProx, point xDist) {
	bool isInside = true;
	// check if the displaced new point/segment is indeed inside domain.
	isInside = domain->isSegmentInside(xProx, xDist);
	if (!isInside) {
		cout << "WARNING: Terminal penetrating step 2 outside domain. Aborted this vessel. Is the domain thick enough?" << "\n";
		isInside = false;
	}
	return isInside;
}

void PenetratingVesselTreeGenerator::setParentMode(SingleVessel* parent, AbstractVascularElement::BRANCHING_MODE mode) {
	parent->branchingMode = mode; 
}

AbstractObjectCCOTree *PenetratingVesselTreeGenerator::generateData(
		long long int saveInterval, string tempDirectory){
	failsafeCheck();

	buildNormalLocator();
	buildIntersectionLocator();

	// get list of bifurcable vessels, the parent vessels
	vector<SingleVessel *> vesselsList = filterTerminalVessels();
	long long int vesselcount = 0;
	// parametricT = 1.0;
	// cout << "WARNING: Hardcode of parametricT to " << parametricT << " . " << endl;
	if (parametricT > 1.0) {
		cout << "ERROR: invalid parametricT " << parametricT << " ." << endl;
		return tree;
	}

	cout << "INFO: Limiting to " << maxGenerateLimit << " vessels." << endl;
	cout << "iterating all " << vesselsList.size() << " terminal segments, bifurcating from terminals and midpoint" << "\n";
	for (vector<SingleVessel *>::iterator it = vesselsList.begin(); it != vesselsList.end() && vesselcount<maxGenerateLimit; ++it, ++vesselcount) {
		dataMonitor->update();
		cout << "Adding: Vessel count " << vesselcount << "; vtksegmentid " << (*it)->vtkSegmentId << endl;
		SingleVessel* parent = (*it);
		
		if (parametricT == 1.0) {
			parent->AbstractVascularElement::BRANCHING_MODE::DISTAL_BRANCHING;
		} else if (parametricT < 1.0) {
			parent->AbstractVascularElement::BRANCHING_MODE::RIGID_PARENT;
		}
		// get the bifurcation point from the terminal
		// the bifurcation points are the terminal and midpoint points
		// Generic point terminal is a combination of terminal (t=1.0) and midpoint (t=0.5)
		point xBifT = parent->xProx + (parent->xDist - parent->xProx)*parametricT;
		
		// define the new terminal point as the projection point.
		// we do not do (*it)->xDist = projectionT; because it changes the tree
		// we need to addVessel(...) instead.

		double length2;
		point normal;
		// get intersection point from bifurcation to surface along versor
		point xNew1 = findDistalDescending(xBifT, normal, length2);

		int foundRaycast;
		point xNew2 = findDistalPenetrating(xNew1, normal, foundRaycast);
		
		bool generateFromTerminal = (isDescendingValid(length2) && 
									 isPenetratingValid(foundRaycast) && 
									 isPenetratingInside(xNew1, xNew2));

		// generate for terminal, generate from midpoint, save the results in a variable:
		// if generate from terminal: save the following:
		// parent (ID), xbif, xnew1, xnew2.
		if (generateFromTerminal) {
			cout << "Saving coordinates to local variable" << endl;
			penetratingData[parametricT][parent->vtkSegmentId] = {xBifT, xNew1, xNew2};
		}
	}
	cout << "iterated through all " << vesselsList.size() << " vessels" << endl;
	return tree;
}

/* 
int temporaryfunction(){

	// INSIDE LOOP
		/// TODO: move this to another function, 2 functions.
		if (generateFromTerminal) {
			cout << "add descending 1st step" << endl;
			// Move to descend();
			((SingleVesselCCOOTree*) tree)->addVesselNoUpdate(xBifT, xNew1, parent, this->vesselFunction, this->branchingMode);

			// Move to penetrate();
			// The following assumes parent vessel was a terminal and has a single child at the distal tip
			cout << "add descending 2nd step" << endl;
			SingleVessel* firstStepVesselT = (SingleVessel *) parent->getChildren()[0];
			((SingleVesselCCOOTree*) tree)->addVesselNoUpdate(xNew1, xNew2, firstStepVesselT, this->vesselFunction, this->branchingMode);
			SingleVessel* secondStepVesselT = (SingleVessel *) firstStepVesselT->getChildren()[0];
			// cout << "added!." << endl;
		}

	// OUTSIDE LOOP
/// TODO: move this to another function
    // Update tree
	cout << "updating the tree" << endl;
	((SingleVesselCCOOTree*) tree)->updateMassiveTree();
	cout << "tree updated" << endl;

	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());
	saveStatus(nTerminals-1);

}
 */

AbstractObjectCCOTree *PenetratingVesselTreeGenerator::descend(
		long long int saveInterval, string tempDirectory){
	//
	//
	return tree;
}

AbstractObjectCCOTree *PenetratingVesselTreeGenerator::penetrate(
		long long int saveInterval, string tempDirectory){
	//
	//
	return tree;
}

int PenetratingVesselTreeGenerator::saveData(/*file type*/) {
	/// TODO: write function
	return 0;
}

int PenetratingVesselTreeGenerator::isValidSegment(point xNew, int iTry) {

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

StagedDomain * PenetratingVesselTreeGenerator::getDomain() {
	return domain;
}

vector<AbstractObjectCCOTree *> PenetratingVesselTreeGenerator::getTrees() {
	vector<AbstractObjectCCOTree *> trees;
	trees.push_back(tree);
	return trees;
}

void PenetratingVesselTreeGenerator::enableConfigurationFile(
		string filename) {
	this->isGeneratingConfFile = 1;
	this->confFilename = filename;
}

void PenetratingVesselTreeGenerator::generatesConfigurationFile(ios::openmode mode) {

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

void PenetratingVesselTreeGenerator::markTimestampOnConfigurationFile(
		string label) {
	if (isGeneratingConfFile) {
		confFile << (chrono::duration_cast<chrono::microseconds>(
				chrono::steady_clock::now() - beginningTime).count())
				/ 1000000.0 << ": " << label << endl;
	}
}

void PenetratingVesselTreeGenerator::closeConfigurationFile() {

	confFile << endl << "DOMAIN_POINTS_GENERATED " << this->domain->getPointCounter() << endl;
	if (isGeneratingConfFile) {
		confFile.flush();
		confFile.close();
	}
}

AbstractObjectCCOTree*& PenetratingVesselTreeGenerator::getTree() {
	return tree;
}

void PenetratingVesselTreeGenerator::setSavingTasks(const vector<AbstractSavingTask*>& savingTasks){
	this->savingTasks = savingTasks;
}

void PenetratingVesselTreeGenerator::saveStatus(long long int terminals){
	tree->setPointCounter(domain->getPointCounter());
	for (std::vector<AbstractSavingTask *>::iterator it = savingTasks.begin(); it != savingTasks.end(); ++it) {
		(*it)->execute(terminals,tree);
	}
//			nodalWriter->write(tempDirectory+ "/step" + to_string(i) + "_view.vtp",tree);
	markTimestampOnConfigurationFile("Generating vessel #" + to_string(terminals));
	markTimestampOnConfigurationFile("Total RAM consumption: " + to_string(monitor->getProcessMemoryConsumption()) + " MB.");
}

void PenetratingVesselTreeGenerator::observableModified(IDomainObservable* observableInstance) {
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

vector<AbstractConstraintFunction<double, int> *>* PenetratingVesselTreeGenerator::getGams()
{
	return &(this->gams);
}

vector<AbstractConstraintFunction<double, int> *>* PenetratingVesselTreeGenerator::getEpsLims()
{
	return &(this->epsLims);
}
vector<AbstractConstraintFunction<double, int> *>* PenetratingVesselTreeGenerator::getNus()
{
	return &(this->nus);
}

double PenetratingVesselTreeGenerator::getDLim() {
	return this->dLim;
}

void PenetratingVesselTreeGenerator::setDLim(double newDLim) {
	this->dLim = newDLim;
}

double PenetratingVesselTreeGenerator::getDLimInitial() {
	return this->dLimInitial;
}

double PenetratingVesselTreeGenerator::getDLimLast() {
	return this->dLimLast;
}

time_t PenetratingVesselTreeGenerator::getBeginTime() {
	return this->beginTime;
}

time_t PenetratingVesselTreeGenerator::getEndTime() {
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