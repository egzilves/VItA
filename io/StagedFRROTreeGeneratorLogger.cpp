// Standard library dependencies
#include<cstdio>
#include<typeinfo>

// VItA dependencies
#include"../core/StagedFRROTreeGenerator.h"
#include"../core/GeneratorData.h"
#include"../structures/tree/AbstractCostEstimator.h"
#include"../structures/tree/VolumetricCostEstimator.h"
#include"../structures/tree/SproutingVolumetricCostEstimator.h"
#include"../structures/tree/AdimSproutingVolumetricCostEstimator.h"
#include"../structures/domain/AbstractDomain.h"
#include"../structures/domain/SimpleDomain2D.h"
#include"../structures/domain/SimpleDomain.h"
#include"../structures/domain/IntersectionVascularizedDomain.h"
#include"../structures/domain/DomainNVR.h"
#include"../structures/domain/PartiallyVascularizedDomain.h"
#include"../structures/domain/DummyDomain.h"
#include"../structures/domain/StagedDomain.h"
#include"../constrains/AbstractConstraintFunction.h"
#include"../constrains/ConstantConstraintFunction.h"
#include"../constrains/ConstantPiecewiseConstraintFunction.h"
#include"../structures/tree/AbstractObjectCCOTree.h"
#include"../structures/tree/SingleVesselCCOOTree.h"
#include"../constrains/AbstractConstraintFunction.h"

// This class interface
#include"StagedFRROTreeGeneratorLogger.h"

void logDomainFiles(FILE *fp, AbstractDomain *domain) 
{   
    domain->logDomainFiles(fp);
}

void logCostEstimator(FILE *fp, AbstractCostEstimator *costEstimator)
{
    costEstimator->logCostEstimator(fp);
}

void logGenData(FILE *fp, GeneratorData *data)
{   
    logCostEstimator(fp, data->costEstimator);
    fprintf(fp, "n_level_test = %d.\n", data->nLevelTest);
    fprintf(fp, "n_terminal_trial = %d.\n", data->nTerminalTrial);
    fprintf(fp, "d_lim_reduction_factor = %f.\n", data->dLimReductionFactor);
    fprintf(fp, "perfusion_area_factor = %f.\n", data->perfusionAreaFactor);
    fprintf(fp, "close_neighborhood factor = %f.\n", data->closeNeighborhoodFactor);
    fprintf(fp, "mid_point_d_lim_factor = %f.\n", data->midPointDlimFactor);
    fprintf(fp, "n_bifurcation_test = %d.\n", data->nBifurcationTest);
    fprintf(fp, "vessel_function = %d.\n", data->vesselFunction);
    fprintf(fp, "reset_d_lim = %d.\n", (int) data->resetsDLim);
}

void logConstraint(FILE *fp, AbstractConstraintFunction<double, int> *constraint) {
    const type_info& type_constant = typeid(ConstantConstraintFunction<double, int>);
    const type_info& constraint_type = typeid(*constraint);
    // Is ConstantConstraintFunction
    if (type_constant.hash_code() == constraint_type.hash_code()) {
        fprintf(fp, "ConstantConstraintFunction = %lf\n", constraint->getValue(0));
    }
    // Is ConstantPiecewiseConstraintFunction
    else {
        ConstantPiecewiseConstraintFunction<double, int> *pieceConstraint = static_cast<ConstantPiecewiseConstraintFunction<double, int> *>(constraint);
        fprintf(fp, "ConstantPiecewiseConstraintFunction\n");
        vector<double> values = pieceConstraint->getValues();
        vector<int> conditions = pieceConstraint->getConditions();
        size_t size = values.size();
        fprintf(fp, "Values = ");
        for (size_t i = 0; i < size; ++i) {
            fprintf(fp, "%lf ", values[i]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "Conditions = ");
        for (size_t i = 0; i < size; ++i) {
            fprintf(fp, "%d ", conditions[i]);
        }
        fprintf(fp, "\n");
    }
}

void logConstraint(FILE *fp, AbstractConstraintFunction<double, double> *constraint) {
    const type_info& type_constant = typeid(ConstantConstraintFunction<double, double>);
    const type_info& constraint_type = typeid(*constraint);
    // Is ConstantConstraintFunction
    if (type_constant.hash_code() == constraint_type.hash_code()) {
        fprintf(fp, "ConstantConstraintFunction = %lf\n", constraint->getValue(0));
    }
    // Is ConstantPiecewiseConstraintFunction
    else {
        ConstantPiecewiseConstraintFunction<double, double> *pieceConstraint = static_cast<ConstantPiecewiseConstraintFunction<double, double> *>(constraint);
        fprintf(fp, "ConstantPiecewiseConstraintFunction\n");
        vector<double> values = pieceConstraint->getValues();
        vector<double> conditions = pieceConstraint->getConditions();
        size_t size = values.size();
        fprintf(fp, "Values = ");
        for (size_t i = 0; i < size; ++i) {
            fprintf(fp, "%lf ", values[i]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "Conditions = ");
        for (size_t i = 0; i < size; ++i) {
            fprintf(fp, "%lf ", conditions[i]);
        }
        fprintf(fp, "\n");
    }
}

void logDomain(FILE *fp, AbstractDomain *domain, long long int n_term, AbstractConstraintFunction<double, int> *gam,
    AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu)
{
    fprintf(fp, "n_term = %lld.\n", n_term);
    fprintf(fp, "gamma\n");
    logConstraint(fp, gam);
    fprintf(fp, "eps_lim\n");
    logConstraint(fp, epsLim);
    fprintf(fp, "nu\n");
    logConstraint(fp, nu);
    fprintf(fp, "n_draw = %d.\n", domain->getDraw());
    fprintf(fp, "random seed = %d.\n", domain->getSeed());
    fprintf(fp, "characteristic_lenght = %f.\n", domain->getCharacteristicLength());
    fprintf(fp, "is_convex_domain = %d.\n", domain->isIsConvexDomain());
    fprintf(fp, "min_bif_angle = %f.\n", domain->getMinBifurcationAngle());
    fprintf(fp, "is_bif_plane_constrained = %d.\n", (int) domain->isIsBifPlaneContrained());
    fprintf(fp, "min_plane_angle = %f.\n", domain->getMinPlaneAngle());
    logGenData(fp, domain->getInstanceData());
}

StagedFRROTreeGeneratorLogger::StagedFRROTreeGeneratorLogger(FILE *fp, StagedFRROTreeGenerator *treeGen)
{
    this->file = fp;
    this->treeGenerator = treeGen;
};

StagedFRROTreeGeneratorLogger::~StagedFRROTreeGeneratorLogger()
{

};

void StagedFRROTreeGeneratorLogger::write()
{   
    FILE *fp = this->file;
    StagedFRROTreeGenerator *generator = this->treeGenerator;
    SingleVesselCCOOTree *tree = (SingleVesselCCOOTree *) generator->getTree();
    double q0 = tree->getQProx();
    StagedDomain* stagedDomain = this->treeGenerator->getDomain();
    vector<AbstractDomain *> *domains = stagedDomain->getDomains();
    vector<long long int> *nTerms = stagedDomain->getNTerminals();
    vector<AbstractConstraintFunction<double, int> *> *gams = generator->getGams();
    vector<AbstractConstraintFunction<double, int> *> *epsLims = generator->getEpsLims();
    vector<AbstractConstraintFunction<double, int> *> *nus = generator->getNus();
    int size = domains->size();
    string filenameCCO = tree->getFilenameCCO();
    if(filenameCCO.empty()) {
        point x0 = tree->getXProx();
        double r0 = tree->getRootRadius();
        fprintf(fp, "Root position = (%f, %f, %f).\n", x0.p[0], x0.p[1], x0.p[2]);
        fprintf(fp, "Root radius = %f.\n", r0);
        fprintf(fp, "Root influx = %f.\n", q0);
    }
    else {
        fprintf(fp, "Input CCO filename = %s\n", filenameCCO.c_str());
        fprintf(fp, "Root influx = %f.\n", q0);
    }
    fprintf(fp, "isInCm = %d\n", (int) tree->isInCm);
    fprintf(fp, "isFL = %d\n", (int) tree->isFL);
    fprintf(fp, "isGammaStage = %d\n", (int) tree->isGammaStage);
    if (tree->gamRadius) {
        fprintf(fp, "Gamma radius:\n");
        logConstraint(fp, tree->gamRadius);
    }
    if (tree->gamFlow) {
        fprintf(fp, "Gamma flow:\n");
        logConstraint(fp, tree->gamFlow);
    }
    for (int i = 0; i < size; ++i) {
        fprintf(fp, "\n");
        fprintf(fp, "Stage[%d]\n", i);
        logDomainFiles(fp, (*domains)[i]);
        logDomain(fp, (*domains)[i], (*nTerms)[i], (*gams)[i], (*epsLims)[i], (*nus)[i]);
    }
    
    fprintf(fp, "\n");
    fprintf(fp, "Initial dLim = %f.\n", generator->getDLimInitial());
    fprintf(fp, "Last dLim = %f.\n", generator->getDLimLast());
    time_t begin_time = generator->getBeginTime();
    time_t end_time = generator->getEndTime();
    struct tm initial_tm = *localtime(&begin_time);
    struct tm last_tm = *localtime(&end_time);
    char time_initial_c_string[21];
    char time_last_c_string[21];
    strftime(time_initial_c_string, 20, "%d_%m_%Y_%H_%M_%S", &initial_tm);
    strftime(time_last_c_string, 20, "%d_%m_%Y_%H_%M_%S", &last_tm);
    fprintf(fp, "\n");
    fprintf(fp, "Beginning of generation time = %s\n", time_initial_c_string);
    fprintf(fp, "End of generation time = %s\n", time_last_c_string);
}