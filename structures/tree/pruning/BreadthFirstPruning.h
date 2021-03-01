#ifndef BREADTH_FIRST_PRUNING_
#define BREADH_FIRT_PRUNING_H_

#include"../SingleVesselCCOOTree.h"
#include"AbstractPruningRule.h"
#include"../AbstractCostEstimator.h" 
#include<queue>

//	TODO The class is not following the abstract prunnig class.
//	TODO Two implementations in one class, Split this class into two.
class BreadthFirstPruning {
private:
    const SingleVesselCCOOTree *ogTree;
    SingleVesselCCOOTree *prunedTree;
    GeneratorData *genData;
    AbstractConstraintFunction<double, int> *gam;
    AbstractConstraintFunction<double, int> *epsLim;
    AbstractConstraintFunction<double, int> *nu;
    AbstractCostEstimator *cost;

    queue<SingleVessel *>* vesselsToPreserve(vector<AbstractPruningRule *>& rules);
public:
    BreadthFirstPruning(const SingleVesselCCOOTree *ogTree);
    ~BreadthFirstPruning();
    SingleVesselCCOOTree* getPrunedTree();
    void pruneTree(vector<AbstractPruningRule *>& rules);
    void pruneTreeFast(vector<AbstractPruningRule *>& rules);
};

#endif
