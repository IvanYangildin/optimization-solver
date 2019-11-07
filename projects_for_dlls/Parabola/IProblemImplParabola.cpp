#include "IProblem.h"
#include "IBrocker.h"
#include "ISolver.h"
#include "ILog.h"
#include "SHARED_EXPORT.h"

#include <math.h>
namespace{
class IProblemImplParabola: public IProblem {
public:
    int getId() const;

    int goalFunction(IVector const* args, IVector const* params, double& res) const;
    int goalFunctionByArgs(IVector const*  args, double& res) const;
    int goalFunctionByParams(IVector const*  params, double& res) const;
    int getArgsDim(size_t& dim) const;
    int getParamsDim(size_t& dim) const;

    int setParams(IVector const* params);
    int setArgs(IVector const* args);

    int derivativeGoalFunction(size_t order, size_t idx, DerivedType dr, double& value, IVector const* args, IVector const* params) const;
    int derivativeGoalFunctionByArgs(size_t order, size_t idx, DerivedType dr, double& value, IVector const* args) const;
    int derivativeGoalFunctionByParams(size_t order, size_t idx, DerivedType dr, double& value, IVector const* params) const;

    ~IProblemImplParabola() {}

    IProblemImplParabola(){
        double* vals_args = new double[1];
        vals_args[0] = 0;
        IVector* argsV = IVector::createVector(1, vals_args);
        setArgs(argsV);

        double* vals_params = new double[1];
        vals_params[0] = 0;
        IVector* paramsV = IVector::createVector(1, vals_params);
        setParams(paramsV);

        delete vals_args;
        delete vals_params;
    }

private:

    IVector* m_params;
    IVector* m_args;
};
}

int IProblemImplParabola::goalFunction(IVector const* args, IVector const* params, double& res) const {
    if (args && (args->getDim() != 1))
        return ERR_DIMENSIONS_MISMATCH;
    double x;

    if (args->getCoord(0, x) != ERR_OK)
        return ERR_ANY_OTHER;
    double a;

    if (params->getCoord(0, a) != ERR_OK)
        return ERR_ANY_OTHER;

    res = x * x - a;
    return ERR_OK;
}

int IProblemImplParabola::goalFunctionByArgs(IVector const* args, double& res) const {
    if (args && (args->getDim() != 1))
        return ERR_DIMENSIONS_MISMATCH;
    double x;
    if (args->getCoord(0, x) != ERR_OK)
        return ERR_ANY_OTHER;

    double a;
    if (m_params->getCoord(0, a) != ERR_OK)
        return ERR_ANY_OTHER;

    res = x * x - a;
    return ERR_OK;
}

int IProblemImplParabola::goalFunctionByParams(IVector const* params, double& res) const {
    if (m_args && (m_args->getDim() != 1))
        return ERR_DIMENSIONS_MISMATCH;

    double x;
    if (m_args->getCoord(0, x) != ERR_OK)
        return ERR_ANY_OTHER;

    double a;
    if (params->getCoord(0, a) != ERR_OK)
        return ERR_ANY_OTHER;

    res = x * x - a;
    return ERR_OK;
}

int IProblemImplParabola::derivativeGoalFunction(size_t order, size_t idx, DerivedType dr, double& value, IVector const* args, IVector const* params) const {
    if (args && (args->getDim() != 1))
        return ERR_DIMENSIONS_MISMATCH;

    double x;
    if (args->getCoord(0, x) != ERR_OK)
        return ERR_ANY_OTHER;

    double a;
    if (m_params->getCoord(0, a) != ERR_OK)
        return ERR_ANY_OTHER;
    switch(dr)
    {
    case BY_PARAMS:
        switch(idx)
        {
        case 0:
            switch (order)
            {
            case 0:
                value = x * x - a;
                break;
            case 1:
                value = -1;
                break;
            default:
                value = 0;
                break;
            }
            return ERR_OK;
        default:
            return ERR_WRONG_ARG;
        }
        break;
    case BY_ARGS:
        switch(idx)
        {
        case 0:
            switch (order){
            case 0:
                value = x * x - a;
                break;
            case 1:
                value = 2 * x;
                break;
            case 2:
                value = 2;
                break;
            default:
                value = 0;
                break;
            }
            return ERR_OK;
        default:
            return ERR_WRONG_ARG;
        }
    default:
        return ERR_WRONG_ARG;
    }
}

int IProblemImplParabola::derivativeGoalFunctionByArgs(size_t order, size_t idx, DerivedType dr, double& value, IVector const* args) const{
    return derivativeGoalFunction(order, idx, dr, value, args, m_params);
}

int IProblemImplParabola::derivativeGoalFunctionByParams(size_t order, size_t idx, DerivedType dr, double& value, IVector const* params) const{
    return derivativeGoalFunction(order, idx, dr, value, m_args, params);
}

int IProblemImplParabola::getId() const{
    return INTERFACE_0;
}

int IProblemImplParabola::setParams(IVector const* params){
    m_params = params->clone();
    if (m_params == NULL)
        return ERR_MEMORY_ALLOCATION;
    return ERR_OK;
}
int IProblemImplParabola::setArgs(IVector const* args){
    m_args = args->clone();
    if (m_args == NULL)
        return ERR_MEMORY_ALLOCATION;
    return ERR_OK;
}

int IProblemImplParabola::getArgsDim(size_t& dim) const{
    if (!m_args)
        return ERR_ANY_OTHER;
    dim = m_args->getDim();
    return ERR_OK;
}
int IProblemImplParabola::getParamsDim(size_t& dim) const{
    if (!m_params)
        return ERR_ANY_OTHER;
    dim = m_params->getDim();
    return ERR_OK;
}

//---------------------------BROCKER---------------------------------------------------------
class BrockerParabola: public virtual IBrocker{
public:
    int getId() const;
    bool canCastTo(Type type) const;
    void* getInterfaceImpl(Type type) const;
    int release();

    BrockerParabola(IProblemImplParabola *problem);

    ~BrockerParabola();

    BrockerParabola();

private:
   IProblemImplParabola * m_problem;
};

int BrockerParabola::getId() const{
   return INTERFACE_0;
}

BrockerParabola::BrockerParabola(IProblemImplParabola *problem): m_problem(problem) {}

BrockerParabola::~BrockerParabola() {
    delete m_problem;
}

bool BrockerParabola::canCastTo(Type type) const{
   return type == IBrocker::PROBLEM ? true : false;
}

void* BrockerParabola::getInterfaceImpl(Type type) const{
   return type == IBrocker::PROBLEM ? m_problem : NULL;
}

int BrockerParabola::release(){
    delete this->m_problem;
    delete this;
    return ERR_OK;
}

extern "C" {
    SHARED_EXPORT void* getBrocker() {
        IProblemImplParabola *problem = new IProblemImplParabola();

        if (!problem) {
            ILog::report("getBrocker: not enough memory\n");
            return NULL;
        }

        BrockerParabola *brocker = new BrockerParabola(problem);

        if (!brocker) {
            ILog::report("getBrocker: not enough memory\n");
            delete problem;
            return NULL;
        }

        return brocker;
    }
}
