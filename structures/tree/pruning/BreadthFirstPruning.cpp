#include"BreadthFirstPruning.h"

#include<queue>

#include"../SingleVesselCCOOTree.h"
#include"../../vascularElements/SingleVessel.h"
#include"AbstractPruningRule.h"
#include"../AbstractCostEstimator.h"

using namespace std;

BreadthFirstPruning::BreadthFirstPruning(const SingleVesselCCOOTree *ogTree) {
    this->ogTree = ogTree;
    this->cost = ogTree->instanceData->costEstimator->clone();
    this->genData = new GeneratorData(ogTree->instanceData->nLevelTest, ogTree->instanceData->nTerminalTrial,
        ogTree->instanceData->dLimReductionFactor, ogTree->instanceData->perfusionAreaFactor, ogTree->instanceData->closeNeighborhoodFactor,
        ogTree->instanceData->midPointDlimFactor, ogTree->instanceData->nBifurcationTest, ogTree->instanceData->vesselFunction,
        ogTree->instanceData->resetsDLim, this->cost);
    this->gam = ogTree->gam->clone();
    this->epsLim = ogTree->epsLim->clone();
    this->nu = ogTree->nu->clone();
    this->prunedTree = new SingleVesselCCOOTree(this->ogTree->xPerf, this->ogTree->rootRadius, this->ogTree->qProx
        ,this->gam, this->epsLim, this->nu, this->ogTree->refPressure, this->ogTree->variationTolerance, this->genData);
}

BreadthFirstPruning::~BreadthFirstPruning() {
    delete this->cost;
    delete this->genData;
    delete this->gam;
    delete this->epsLim;
    delete this->nu;
    delete this->prunedTree;
}

bool isVesselToBePruned(SingleVessel *vessel, vector<AbstractPruningRule *>& rules) {
    for (auto it = rules.begin(); it != rules.end(); ++it) {
        bool flag = (*it)->needsPruning(vessel);
        if (flag) {
            return true;
        }
    }
    return false;
}

queue<SingleVessel *>* BreadthFirstPruning::vesselsToPreserve(vector<AbstractPruningRule *>& rules) {
    queue<SingleVessel *> *toPreserve = new queue<SingleVessel *>;
    queue<SingleVessel *> *searchable = new queue<SingleVessel *>;
    SingleVessel *root = static_cast<SingleVessel *>(this->ogTree->root);
    searchable->push(root);
    while (!searchable->empty()) {
        SingleVessel *curVessel = searchable->front();
        searchable->pop();
        if(!isVesselToBePruned(curVessel, rules)) {
            toPreserve->push(curVessel);
        }
        vector<AbstractVascularElement *> children = curVessel->getChildren();
        for (auto it = children.begin(); it != children.end(); ++it) {
            searchable->push((SingleVessel *) (*it));
        }
    }
    delete searchable;
    return toPreserve;
}

void BreadthFirstPruning::pruneTree(vector<AbstractPruningRule *>& rules) {
    unordered_map<SingleVessel*, SingleVessel*> copiedTo;
    // This represents the root parent
    copiedTo.insert(pair<SingleVessel*, SingleVessel*>(static_cast<SingleVessel *>(nullptr), static_cast<SingleVessel *>(nullptr)));
    queue<SingleVessel *>* toPreserve = this->vesselsToPreserve(rules);
    while(!toPreserve->empty()) {
        SingleVessel *vesselToCopy = toPreserve->front();
        toPreserve->pop();
        SingleVessel *vesselToAdd = new SingleVessel();
        copiedTo.insert(pair<SingleVessel*, SingleVessel*>(static_cast<SingleVessel *>(vesselToCopy), static_cast<SingleVessel *>(vesselToAdd)));
        this->prunedTree->addValitatedVessel(vesselToAdd, vesselToCopy, copiedTo);
    }
    delete toPreserve;
}

void BreadthFirstPruning::pruneTreeFast(vector<AbstractPruningRule *>& rules) {
    unordered_map<SingleVessel*, SingleVessel*> copiedTo;
    // This represents the root parent
    copiedTo.insert(pair<SingleVessel*, SingleVessel*>(static_cast<SingleVessel *>(nullptr), static_cast<SingleVessel *>(nullptr)));
    queue<SingleVessel *>* toPreserve = this->vesselsToPreserve(rules);
    while(!toPreserve->empty()) {
        SingleVessel *vesselToCopy = toPreserve->front();
        toPreserve->pop();
        SingleVessel *vesselToAdd = new SingleVessel();
        copiedTo.insert(pair<SingleVessel*, SingleVessel*>(static_cast<SingleVessel *>(vesselToCopy), static_cast<SingleVessel *>(vesselToAdd)));
        this->prunedTree->addValitatedVesselFast(vesselToAdd, vesselToCopy, copiedTo);
    }

   	//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
	this->prunedTree->updateTree(static_cast<SingleVessel *>(this->prunedTree->root), this->prunedTree);

    //	Update resistance, pressure and betas
	double maxVariation = INFINITY;
	while (maxVariation > this->prunedTree->variationTolerance) {
	    this->prunedTree->updateTreeViscositiesBeta(static_cast<SingleVessel *>(this->prunedTree->root), &maxVariation);
	}    

    delete toPreserve;
}

SingleVesselCCOOTree* BreadthFirstPruning::getPrunedTree() {
    return this->prunedTree;
}