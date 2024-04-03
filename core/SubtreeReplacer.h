/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2024 Eduardo Guerreiro Zilves & Gonzalo Maso Talou */
/*
 * SubtreeReplacer.h
 *
 *  Created on: Mar 20, 2024
 *      Author: Eduardo G. Zilves
 */

#ifndef SUBTREEREPLACER_H_
#define SUBTREEREPLACER_H_

#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <unordered_map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellLocator.h>

#include "../constrains/AbstractConstraintFunction.h"
#include "../io/task/AbstractSavingTask.h"
#include "../structures/domain/IDomainObserver.h"
#include "../structures/domain/StagedDomain.h"
#include "../structures/tree/AbstractObjectCCOTree.h"
#include "../structures/tree/SingleVesselCCOOTree.h"
#include "../utils/MemoryMonitor.h"
#include "../structures/vascularElements/SingleVessel.h"
#include "../structures/CCOCommonStructures.h"

#include "GeneratorDataMonitor.h"
#include "StagedFRROTreeGenerator.h"



using namespace std;


// struct ReadData {
//     point xBif, xNew, xPProx, xPDist;
//     AbstractVascularElement::VESSEL_FUNCTION function;
//     int stage;
// };

/**
 * Generator for a projected perfusion from existing tree with many stages.
 * Generates a segment from terminals into the projectiondomain, and continues inside the domain
 * with the same normal versor.
 */
class SubtreeReplacer : public IDomainObserver  {

	/** Time of the beggining of the tree generation process*/
	time_t beginTime;
	/** Time of the end of the tree generation process*/
	time_t endTime;
	/** Wrapper with parameters associated to a tree generation process. */
	GeneratorData *instanceData;
	/**	Monitor of the @p instanceData . */
	GeneratorDataMonitor *dataMonitor;
	/**	Monitors the memory usage. */
	MemoryMonitor *monitor;
	/**	Action executed at each save interval during generate and resume methods */
	vector<AbstractSavingTask *> savingTasks;

	/** Amount of terminals in the trees.*/
	long long int nTerminals;
	/** Perfusion volume.*/
	double dLim;
	/** Initial dLim value*/
	double dLimInitial;
	/** Final dLim value*/
	double dLimLast;
	/** Perfusion domain. */
	StagedDomain *domain;

	/** Generated trees. */
	// AbstractObjectCCOTree *tree;
	SingleVesselCCOOTree *tree;

	vector<AbstractConstraintFunction<double,int> *> gams;
	vector<AbstractConstraintFunction<double,int> *> epsLims;
	vector<AbstractConstraintFunction<double,int> *> nus;

	/**	Current stage of generation.*/
	int stage;
	/**	If the current generation saves the configuration used.*/
	int isGeneratingConfFile;
	/**	File name for the configuration file.*/
	string confFilename;

	bool didAllocateTree;
	/**
	 * Bypass vessel function when bifurcating, to allow bifurcation from vessels with midpoint in partitioned domain.
	 * Used in the resume midpoint function.
	 */
	// bool bypassFunctionIfMidpointInside;

	// /** Domain. */
	// string domainFile;
	// /** Projection domain. */
	// string projectionDomainFile;
	// /** vtkPolydata description of the domain. */
	// vtkSmartPointer<vtkPolyData> vtkGeometryIntersect;
	// /** vtkPolydata description of the projection domain. */
	// vtkSmartPointer<vtkPolyData> vtkGeometryProjection;
	// /** Cell locator responsible to determine if a segment is inside the domain. */
	// vtkSmartPointer<vtkCellLocator> locatorIntersect;
	// /** Cell locator responsible to determine if a segment is inside the domain, check closest cell to project. */
	// vtkSmartPointer<vtkCellLocator> locatorProjection;
	// /**	Descending offset for projected points */
	// double descendingOffset;
	// /**	Penetration offset for projected points */
	// double endpointOffset;
	// /** Maximum distance between terminal and closest point/cell on projection surface to allow generating a vessel. */
	// double maxDistanceToClosestPoint;
	// /** Maximum segment length inside domain. */
	// double maxPenetratingVesselLength;

    // SingleVesselCCOOTree *treee;
    // vector<vector<ReadData>*> *vesselToMerge;
    // unordered_map<string, SingleVessel *> *stringToPointer;

public:
	/**
	 * Constructor to resume from a pre-existent tree.
	 * @param domain	Perfusion domain.
	 * @param tree		Pre-existent tree.
	 * @param nTerm		Total number of tree terminals.
	 * @param gam		Murray's power law function.
	 * @param epsLim	Symmetry constraint function.
	 * @param nu		Viscosity function.
	 */
	SubtreeReplacer(StagedDomain *domain, string projectionDomainFile, AbstractObjectCCOTree *tree, long long int nTerm, 
		vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, 
		vector<AbstractConstraintFunction<double,int> *>nu);
	/**
	 * Common destructor.
	 */
	~SubtreeReplacer();

	/**
	 * Replace the segment with a subtree given the probability of picking
	 */
	AbstractObjectCCOTree *replaceSegments(long long int saveInterval, string tempDirectory);
	
	// // /**
	// //  * Resumes the tree generation.
	// //  * @param saveInterval Number of iterations performed between saved steps.
	// //  * @param tempDirectory Directory where intermediate solutions are saved.
	// //  * @return	Perfusion tree.
	// //  */
	// // AbstractObjectCCOTree *resume(long long int saveInterval, string tempDirectory);
	// // /**
	// //  * Resumes the tree generation process and saves the optimal xNew and xBif in @param fp.
	// //  */
	// // AbstractObjectCCOTree *resumeSavePointsMidPoint(long long int saveInterval, string tempDirectory, FILE *fp);
	// /**
	//  * Generates the tree penetration into domain.
	//  * @param saveInterval Number of iterations performed between saved steps.
	//  * @param tempDirectory Directory where intermediate solutions are saved.
	//  * @return	Perfusion tree.
	//  */
	// AbstractObjectCCOTree *generatePenetrating(long long int saveInterval, string tempDirectory);
	
	/**
	 * Returns the perfusion domain.
	 * @return Perfusion domain.
	 */
	StagedDomain* getDomain();
	/**
	 * Returns the generated tree.
	 * @return Generated tree.
	 */
	vector<AbstractObjectCCOTree *> getTrees();
	/**
	 * Enables the configuration file generation capabilities.
	 * @param filename	File name where the generator stores the configuration data.
	 */
	void enableConfigurationFile(string filename);
	/**
	 * Updates internal data from this instance after the domain has been modified.
	 * @param observableInstance Domain instance modified.
	 */
	void observableModified(IDomainObservable * observableInstance);
	/**
	 * Return the currently generated @p tree
	 * @return Generated tree.
	 */
	SingleVesselCCOOTree*& getTree();
	/**
	 * Saves the current generation status using the @p savingTasks and the memory monitor.
	 * @param terminals	Current iteration number.
	 */
	void saveStatus(long long int terminals);
	/**
	 * Sets a set of saving tasks to be produced each time the generation process saves.
	 * @param savingTasks	Set of tasks to be performed each time that the tree is saved.
	 */
	void setSavingTasks(const vector<AbstractSavingTask*>& savingTasks);
	/**
	 * Returns a pointer to the vector of gams.
	 */
	vector<AbstractConstraintFunction<double, int> *>* getGams();
	/**
	 * Returns a pointer to the vector of eps_lims.
	 */
	vector<AbstractConstraintFunction<double, int> *>* getEpsLims();
	/**
	 * Returns a pointer to the vector of nus.
	 */
	vector<AbstractConstraintFunction<double, int> *>* getNus();
	/**
	 * Returns the current DLim value.
	 */
	double getDLim();
	/**
	 * Set the DLim value. Use this only to resume process from a previous generation.
	 */
	void setDLim(double newDLim);
	
	time_t getBeginTime();

	time_t getEndTime();

	double getDLimInitial();

	double getDLimLast();
	
	/**
	 * Projects the terminals of the vessels list to the closest point within the domain.
	 * @param vessels
	 */
	void projectTerminals(vector<SingleVessel *> vessels);
	// /**
	//  * Projects the vessels in the vessels list to the closest cell within the domain.
	//  * @param vessels
	//  */
	// void projectVessel(vector<SingleVessel *> vessels);
	/**
	 * Sets the projection surface associated to the domain
	 * @param projectionDomainFile projection domain file, VTK file which contains the surface of the projection domain.
	 */
	void setProjectionDomain(string projectionDomainFile);
	/**
	 * Sets the descending offset of the projector into the domain. All projected points will be pushed inside.
	 * @param descendingOffset projection offset.
	 */
	void setDescendingOffset(double descendingOffset);
	/**
	 * Sets the penetration offset the projector into the domain. All projected points will be pushed inside.
	 * @param penetrationOffset projection offset.
	 */
	void setPenetrationOffset(double penetrationOffset);


    // void TreeMerger(SingleVesselCCOOTree *baseTree, vector<string>& derivedTreePoints);
    // void /*~*/dTreeMerger();
    // SingleVesselCCOOTree *mergeFast();


protected:
	/**	Configuration file stream. */
	ofstream confFile;
	/** Initial timestamp of the generation process. */
	chrono::steady_clock::time_point beginningTime;
	/**
	 * Returns if the length of the segment defined by xProx and xNew is higher than dLim (distance criterion).
	 * @param xNew	Proposed distal point.
	 * @param iTry	Number of trial.
	 * @return If the root segment is valid.
	 */
	// int isValidRootSegment(point xNew, int iTry);
	/**
	 * Returns if the closest point to the whole try is greater than dLim (distance criterion for a new segment).
	 * @param xNew Proposed distal point.
	 * @param iTry Number of trial.
	 * @return If the segment is valid.
	 */
	int isValidSegment(point xNew, int iTry);
	/**
	 * Generates the configuration file for the current tree generation.
	 * @param mode Is the openmode used for the generated file (ios::out for generation, ios::app for resume).
	 */
	void generatesConfigurationFile(ios::openmode mode);
	/**
	 * Saves the timestamp for event @p label.
	 */
	void markTimestampOnConfigurationFile(string label);
	/**
	 * Closes the configuration file for the current tree generation.
	 */
	void closeConfigurationFile();


private:
    // void createMapping(SingleVessel *vessel);

};
#endif