#include "ICompact.h"
#include "ILog.h"
#include <QVector>
#include <cfloat>
#include <limits.h>
#include <math.h>
#include <iostream>

#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SIGN(x) (x >= 0) ? ((x>0)?(1):(0)):(-1)
#define MY_EPS 1e-6

namespace{
    class ICompactImpl: public virtual ICompact
    {
    public:
        virtual int getId() const;

        /*factories*/
        static ICompact* createCompact(IVector const* const begin, IVector const* const end, IVector const* const step = 0);

        static ICompact* Intersection(ICompact const* const left, ICompact const* const right);
        static ICompact* Union(ICompact const* const left, ICompact const* const right);
        static ICompact* Difference(ICompact const* const left, ICompact const* const right);
        static ICompact* SymDifference(ICompact const* const left, ICompact const* const right);

        virtual int deleteIterator(IIterator * pIter);
        virtual int getByIterator(IIterator const* pIter, IVector*& pItem) const;

        virtual IIterator* end(IVector const* const step = 0);
        virtual IIterator* begin(IVector const* const step = 0);

        virtual int isContains(IVector const* const vec, bool& result) const;
        virtual int numContains(IVector const* const vec, int& result) const;
        virtual int isSubSet(ICompact const* const other) const;


        virtual int getNearestNeighbor(IVector const* vec, IVector *& nn) const;

        virtual ICompact* clone() const;

        /*dtor*/
        virtual ~ICompactImpl();

        void vecToIdx(IVector const* const vec, unsigned int& res) const;
        void idxToVec(unsigned int idx, IVector*& res) const;

        bool checkStep(IVector const* const step);

        class IIteratorImpl: public virtual ICompact::IIterator
        {
        public:
            //adds step to current value in iterator
            virtual int doStep();

            //change step
            virtual int setStep(IVector const* const step);

            const ICompact* getCompact() const {return m_com;}
            unsigned int getPos() const {return m_pos;}

            /*dtor*/
            virtual ~IIteratorImpl();
            IIteratorImpl(ICompact const* const compact, int pos, IVector const* const step);
        protected:




        private:
            ICompact* m_com;
            unsigned int m_pos;
            unsigned int m_com_num;
            IVector* m_step;

            /*non default copyable*/
            IIteratorImpl(const IIteratorImpl& other) = delete;
            void operator=(const IIteratorImpl& other) = delete;
        };

    protected:
        ICompactImpl() = default;
        ICompactImpl(IVector* begin, IVector* end, IVector* step, unsigned int card);
        ICompact* SimpleIntersection(ICompact const * const left, ICompact const * const right, int& error) const;
        //composite has not empty vec_compact or they are simple both
        int MediumIntersection(ICompact const * const simple, bool is_anti, ICompactImpl const * const composite, ICompactImpl*& new_comp) const;
        int get_compact_edges(ICompact const* const comp, IVector*& begin, IVector*& end) const;
        ICompact* simple_clone() const;
    private:

        /*non default copyable*/
        ICompactImpl(const ICompactImpl& other) = delete;
        void operator=(const ICompactImpl& other) = delete;

        QVector <ICompact*> vec_compact;
        QVector <bool> anti_compact;

        unsigned int m_dim;
        unsigned int m_cardinality; // how many points in the compact's grid

        IVector const* m_begin;
        IVector const* m_end;
        IVector const* m_steps;

        QVector <IIteratorImpl*> iterators;
    };
}

int ICompactImpl::get_compact_edges(ICompact const* const comp, IVector*& begin, IVector*& end) const
{
    ICompactImpl const* impl_comp = dynamic_cast<ICompactImpl const*>(comp);
    begin = impl_comp->m_begin->clone();
    end = impl_comp->m_end->clone();

    return ERR_OK;
}

ICompact::IIterator::IIterator(const ICompact *const compact, int pos, const IVector *const step){}

int ICompactImpl::getId() const
{
    return INTERFACE_0;
}

ICompact* ICompact::createCompact(const IVector *const begin, const IVector *const end, const IVector *const step)
{
    return ICompactImpl::createCompact(begin, end, step);
}

ICompact* ICompactImpl::createCompact(const IVector *const begin, const IVector *const end, const IVector *const step)
{
    if (!begin || !end)
    {
        ILog::report("ICompact.createCompact: constructor missed begin or end");
        return 0;
    }

    if (begin->getDim() != end->getDim())
    {
        ILog::report("ICompact.createCompact: dimension mismatch in createCompact");
        return 0;
    }

    int dim = begin->getDim();

    double* t_begin = new double[dim];
    double* t_end = new double[dim];
    if(!t_begin || !t_end)
    {
        delete[] t_begin;
        t_begin = 0;
        delete[] t_end;
        t_end = 0;

        ILog::report("ICompact.createCompact: not enough memory");
        return 0;
    }

    double b,e;
    for (unsigned int i = 0; i < dim; i++)
    {
        begin->getCoord(i, b);
        end->getCoord(i, e);
        t_begin[i] = (b < e)?b:e;
        t_end[i] = (b < e)?e:b;
    }

    IVector* try_begin = IVector::createVector(dim, t_begin);
    IVector* try_end = IVector::createVector(dim, t_end);

    delete[] t_begin;
    t_begin = 0;
    delete[] t_end;
    t_end = 0;

    if(!try_begin || !try_end)
    {
        delete try_begin;
        try_begin = 0;
        delete try_end;
        try_end = 0;

        ILog::report("ICompact.createCompact: not enough memory");
        return 0;
    }

    IVector* try_steps = 0;
    unsigned int card = 1;

    if (step) //non-default step
    {
        if (step->getDim() != begin->getDim() || step->getDim() != end->getDim() || begin->getDim() != end->getDim())
        {
            delete try_begin;
            try_begin = 0;
            delete try_end;
            try_end = 0;

            ILog::report("ICompact.createCompact: dimension mismatch in createCompact");
            return 0;
        }
        /*double b,e;
        for(unsigned int i = 0; i<end->getDim(); i++)
        {
            begin->getCoord(i,b);
            end->getCoord(i,e);
            /*if(b == e)
            {
                ILog::report("ICompact.createCompact: degenerate segment in compact construction");
                return 0;
            }
        }*/
        double d;

        for (unsigned int i = 0; i < dim; i++) // check if step data match the unsigned int
        {
            step->getCoord(i,d);
            if(d == 0)
            {
                delete try_begin;
                try_begin = 0;
                delete try_end;
                try_end = 0;

                ILog::report("ICompact.createCompact: zero step is unacceptable");
                return 0;
            }
            if((ABS(card) > (unsigned int)((double)UINT_MAX/ABS(d+1))))//Ivan's noting
            {
                delete try_begin;
                try_begin = 0;
                delete try_end;
                try_end = 0;

                ILog::report("ICompact.createCompact: there are too many vertices in the compact grid");
                return 0;
            }
            else
                card *= (unsigned int)ABS(d+1);//Ivan's noting
        }
        try_steps = step->clone();
        if(!try_steps)
        {
            ILog::report("ICompact.createCompact: not enough memory");
            return 0;
        }
        return new ICompactImpl(try_begin, try_end, try_steps, card);

    }
    else // default step
    {
        double gridStep = floor(pow((double)UINT_MAX,1./dim));
        double* tmp = new double[dim];
        if(!tmp)
        {
            delete try_begin;
            try_begin = 0;
            delete try_end;
            try_end = 0;

            ILog::report("ICompact.createCompact: not enough memory");
            return 0;
        }
        unsigned int card = 1;
        for(unsigned int i = 0; i < dim; i++)
        {
            tmp[i] = gridStep;
            card *= (unsigned int)ABS(gridStep);
        }
        try_steps = IVector::createVector(dim, tmp);

        delete[] tmp;
        tmp = 0;

        if(!try_steps)
        {
            delete try_begin;
            try_begin = 0;
            delete try_end;
            try_end = 0;

            ILog::report("ICompact.createCompact: not enough memory");
            return 0;
        }
        /*double b,e;
        for(unsigned int i = 0; i<end->getDim(); i++)
        {
            begin->getCoord(i,b);
            end->getCoord(i,e);
            if(b == e)
            {
                ILog::report("ICompact.createCompact: degenerate segment in compact construction");
                return 0;
            }
        }*/
        return new ICompactImpl(try_begin, try_end, try_steps, card);
    }
}

ICompactImpl::ICompactImpl(IVector* begin, IVector* end, IVector* step, unsigned int card):m_dim(begin->getDim()), m_begin(begin), m_end(end), m_steps(step), m_cardinality(card){}

ICompact* ICompactImpl::clone() const
{
    ICompactImpl* comp = dynamic_cast<ICompactImpl*>(ICompact::createCompact(m_begin, m_end, m_steps));
    if (!comp)
    {
        return 0;
    }
    for (int i = 0; i < vec_compact.size(); i++)
    {
        comp->anti_compact.append(anti_compact[i]);
        ICompact* sub_comp = vec_compact[i]->clone();
        if (!sub_comp)
        {
            delete comp;
            ILog::report("ICompact.clone: not enogh memory");
            return 0;
        }
        comp->vec_compact.append(sub_comp);
    }
    return comp;
}

ICompact* ICompactImpl::simple_clone() const
{
    return ICompact::createCompact(m_begin, m_end, m_steps);
}

int ICompactImpl::deleteIterator(IIterator *pIter)
{
    if(!pIter)
    {
        ILog::report("ICompact.deleteIterator: got nullptr, expected ICompact::IITerator pointer");
        return ERR_WRONG_ARG;
    }

    const ICompact* its_comp = dynamic_cast<IIteratorImpl*>(pIter)->getCompact();
    if(this != its_comp)
    {
        for (int i = 0; i < vec_compact.size(); i++)
        {
            if ((!anti_compact[i]) && (vec_compact[i] == its_comp))
            {
                int last_err = vec_compact[i]->deleteIterator(pIter);
                return last_err;
            }
        }
        ILog::report("ICompact.deleteIterator: tries to delete iterator which doesn't belong to this instance of ICompact");
        return ERR_WRONG_ARG;
    }
    int idx = iterators.indexOf(dynamic_cast<IIteratorImpl*>(pIter));
    if (idx == -1)
    {
        ILog::report("ICompact.deleteIterator: tries to delete iterator which doesn't belong to this instance of ICompact");
        return ERR_WRONG_ARG;
    }
    delete iterators[idx];
    iterators.remove(idx);

    return ERR_OK;
}


ICompact* ICompact::Intersection(ICompact const* const left, ICompact const* const right)
{
    return ICompactImpl::Intersection(left, right);
}

ICompact* ICompact::Union(ICompact const* const left, ICompact const* const right)
{
    return ICompactImpl::Union(left, right);
}

ICompact* ICompact::Difference(ICompact const* const left, ICompact const* const right)
{
    return ICompactImpl::Difference(left, right);
}

ICompact* ICompact::SymDifference(ICompact const* const left, ICompact const* const right)
{
    return ICompactImpl::SymDifference(left, right);
}


ICompact* ICompactImpl::SimpleIntersection(ICompact const * const left, ICompact const * const right, int& error) const
{
    IVector *l_begin = 0, *r_begin = 0;
    IVector *l_end = 0, *r_end = 0;
    error = get_compact_edges(left, l_begin, l_end);
    if (error != ERR_OK)
    {
        return 0;
    }
    error = get_compact_edges(right, r_begin, r_end);
    if (error != ERR_OK)
    {
        delete l_begin;
        delete l_end;
        return 0;
    }
    if (r_begin->getDim() != l_begin->getDim()){
        delete l_begin;
        delete l_end;
        delete r_begin;
        delete r_end;
        error = ERR_DIMENSIONS_MISMATCH;
        return 0;
    }

    unsigned dim = r_begin->getDim();
    double* arr = new double[dim];
    double l_coord, r_coord;
    for (unsigned i = 0; i < dim; i++){
        l_begin->getCoord(i, l_coord);
        r_begin->getCoord(i, r_coord);
        arr[i] = l_coord > r_coord? l_coord : r_coord;
    }
    delete l_begin;
    delete r_begin;
    l_begin = 0;
    r_begin = 0;
    IVector* new_begin = IVector::createVector(dim, arr);
    if (!new_begin){
        delete l_end;
        delete r_end;
        error = ERR_ANY_OTHER;
        return 0;
    }

    for (unsigned i = 0; i < dim; i++){
        l_end->getCoord(i, l_coord);
        r_end->getCoord(i, r_coord);
        arr[i] = l_coord < r_coord? l_coord : r_coord;
    }
    delete l_end;
    delete r_end;
    l_end = 0;
    r_end = 0;
    IVector* new_end = IVector::createVector(dim, arr);
    if (!new_end){
        delete new_begin;
        error = ERR_ANY_OTHER;
        return 0;
    }

    double e, b;
    for (unsigned i = 0; i < dim; i++){
        new_end->getCoord(i, e);
        new_begin->getCoord(i, b);
        if (e <= b){
            delete new_end;
            delete new_begin;
            return 0;
        }
    }
    IVector* new_step = dynamic_cast<ICompactImpl const * const>(left)->m_steps->clone();
    if (!new_step){
        delete new_begin;
        delete new_end;
        error = ERR_MEMORY_ALLOCATION;
        return 0;
    }
    ICompact* new_comp = ICompact::createCompact(new_begin, new_end, new_step);
    if (!new_comp)
        error = ERR_MEMORY_ALLOCATION;
    delete new_begin;
    delete new_step;
    delete new_end;

    return new_comp;
}

int ICompactImpl::MediumIntersection(ICompact const * const simple, bool is_anti, ICompactImpl const * const composite, ICompactImpl*& new_comp) const
{
    int error;
    if (!new_comp && is_anti)
        return 0;
    ICompact* temp = ICompactImpl::SimpleIntersection(simple, composite, error);
    if (error != ERR_OK)
        return error;
    if (temp){
        if (!new_comp)
            new_comp = dynamic_cast<ICompactImpl*>(temp);
        else{
            new_comp->vec_compact.append(temp);
            new_comp->anti_compact.append(is_anti);
        }
    }

    for (int i = 0; i < composite->vec_compact.size(); i++){
        if ((is_anti || !new_comp) && composite->anti_compact[i])
            continue;
        temp = ICompactImpl::SimpleIntersection(simple, composite->vec_compact[i], error);
        if (error != ERR_OK)
            return error;
        if (temp){
            if (!new_comp)
                new_comp = dynamic_cast<ICompactImpl*>(temp);
            else{
                new_comp->vec_compact.append(temp);
                new_comp->anti_compact.append(is_anti);
            }
        }
    }

    return error;
}

ICompact* ICompactImpl::Intersection(ICompact const* const left, ICompact const* const right)
{
    if (!left || !right)
        return 0;

    int is_need = left->isSubSet(right);
    if (is_need == -1){
        ILog::report("ICompact.Intersection: unable to compare compacts");
        return 0;
    }
    ICompactImpl const* impl_left = dynamic_cast<ICompactImpl const* const>(left);
    if (is_need == 1){
        return impl_left->clone();
    }
    ICompactImpl const* impl_right = dynamic_cast<ICompactImpl const* const>(right);
    ICompactImpl* new_comp = 0;

    int error;
    error = impl_left->MediumIntersection(impl_left, false, impl_right, new_comp);
    if (error != ERR_OK){
        ILog::report("ICompact.Intersection: unable to intersect compacts");
        delete new_comp;
        return 0;
    }
    for (int i = 0; i < impl_left->vec_compact.size(); i++){
        error = impl_left->MediumIntersection(impl_left->vec_compact[i], impl_left->anti_compact[i], impl_right, new_comp);
        if (error != ERR_OK){
            ILog::report("ICompact.Intersection: unable to intersect compacts");
            delete new_comp;
            return 0;
        }
    }

    return new_comp;
}

ICompact* ICompactImpl::Union(ICompact const* const left, ICompact const* const right)
{
    if (!right && !left)
        return 0;
    if (!right)
        return left->clone();
    if (!left)
        return right->clone();

    int is_need = right->isSubSet(left);
    if (is_need == -1){
        ILog::report("ICompact.Union: unable to compare compacts");
        return 0;
    }
    ICompactImpl const* impl_left = dynamic_cast<ICompactImpl const*>(left);
    if (is_need == 1){
        return impl_left->clone();
    }
    ICompactImpl const* impl_right = dynamic_cast<ICompactImpl const*>(right);

    ICompact* new_right = impl_right->simple_clone();
    if (!new_right){
        ILog::report("ICompact.Union: not enough memmory");
        return 0;
    }
    ICompactImpl* result = dynamic_cast<ICompactImpl*>(impl_left->clone());
    if (!result){
        ILog::report("ICompact.Union: not enough memmory");
        delete new_right;
        return 0;
    }

    result->vec_compact.append(new_right);
    result->anti_compact.append(false);

    for (int i = 0; i < impl_right->vec_compact.size(); i++){
        ICompact* temp = impl_right->vec_compact[i]->clone();
        if (!temp){
            ILog::report("ICompact.Union: not enough memmory");
            delete result;
            return 0;
        }
        result->vec_compact.append(temp);
        result->anti_compact.append(impl_right->anti_compact[i]);
    }

    return result;
}

ICompact* ICompactImpl::Difference(const ICompact *const left, const ICompact *const right)
{
    if (!left)
        return 0;
    if (!right)
        return left->clone();

    int code = left->isSubSet(right);
    if (code == -1){
        ILog::report("ICompact.Difference: unable to compare compacts");
        return 0;
    }
    if (code == 1)
        return 0;
    ICompact* temp = Intersection(left, right);
    if (!temp){
        ILog::report("ICompact.Difference: error while subtraction");
        return 0;
    }
    ICompactImpl* result = dynamic_cast<ICompactImpl*>(Union(left, temp));
    delete temp;
    if (!result){
        ILog::report("ICompact.Difference: error while subtraction");
        return 0;
    }
    ICompactImpl const* impl_left = dynamic_cast<ICompactImpl const*>(left);
    int left_size = impl_left->anti_compact.size();
    int new_size = result->anti_compact.size();
    for (int i = left_size; i < new_size; i++){
        result->anti_compact[i] = !result->anti_compact[i];
    }
    return result;
}

ICompact* ICompactImpl::SymDifference(const ICompact *const left, const ICompact *const right)
{
    ICompact* inters_comp = Intersection(left, right);
    if (!inters_comp){
        return 0;
    }
    return Difference(left, inters_comp);
}


int ICompactImpl::getByIterator(const IIterator *pIter, IVector *&pItem) const
{
    if(!pIter)
    {
        ILog::report("ICompact.getByIterator: got nullptr, expected ICompact::IITerator pointer");
        return ERR_WRONG_ARG;
    }
    ICompact const* its_comp = dynamic_cast<const IIteratorImpl*>(pIter)->getCompact();
    if(this != its_comp)
    {
        for (int i = 0; i < vec_compact.size(); i++)
        {
            if ((!anti_compact[i]) && (vec_compact[i] == its_comp))
            {
                int last_err = its_comp->getByIterator(pIter, pItem);
                return last_err;
            }
        }
        ILog::report("ICompact.getByIterator: got iterator which doesn't belong to this instance of ICompact");
        return ERR_WRONG_ARG;
    }

    if(pItem)
        delete pItem;

    unsigned int pos = dynamic_cast<const IIteratorImpl*>(pIter)->getPos();
    IVector* tmp = 0;
    idxToVec(pos,tmp);

    pItem = tmp->clone();
    delete tmp;

    return ERR_OK;
}

int ICompactImpl::getNearestNeighbor(const IVector *vec, IVector *&nn) const
{
    if(!vec)
    {
        ILog::report("Icompact.getNearestNeighbor: got nullptrm expected IVector pointer");
        return ERR_WRONG_ARG;
    }
    if(vec->getDim() != m_dim)
    {
        ILog::report("Icompact.getNearestNeighbor: dimension mismatch in getNearestNeighbor");
        return ERR_WRONG_ARG;
    }
    bool contains;
    isContains(vec, contains);
    if(!contains)
    {
        ILog::report("Icompact.getNearestNeighbor: can't find nearest for vector outside the compact");
        return ERR_WRONG_ARG;
    }
    double* tmp = new double[m_dim];
    double x,b,e,s;
    for(unsigned int i = 0; i < m_dim; i++)
    {
        vec->getCoord(i,x);
        m_begin->getCoord(i,b);
        m_end->getCoord(i,e);
        m_steps->getCoord(i,s);
        unsigned int numStepBeforeVec = (unsigned int)floor(((x-b)/(e-b))*s);
        if(ABS(x-numStepBeforeVec*(e-b)/s - b) > ABS(x- (numStepBeforeVec+1)*(e-b)/s - b))
        {
            tmp[i] = (numStepBeforeVec+1)*(e-b)/s + b;
        }
        else
        {
            tmp[i] = numStepBeforeVec*(e-b)/s + b;
        }
    }

    if(nn)
        delete nn;

    nn = IVector::createVector(m_dim, tmp);
    delete[] tmp;
    tmp = 0;
    return ERR_OK;
}

ICompactImpl::IIterator* ICompactImpl::end(IVector const* const step)
{
    IIterator* it = new IIteratorImpl(this,m_cardinality-1,step);

    if(!checkStep(step))
    {
        delete it;
        return 0;
    }

    iterators.push_back(dynamic_cast<IIteratorImpl*>(it));
    return it;
}

ICompactImpl::IIterator* ICompactImpl::begin(IVector const* const step)
{
    IIterator* it = new IIteratorImpl(this,0,step);


    if(!checkStep(step))
    {
        delete it;
        return 0;
    }

    iterators.push_back(dynamic_cast<IIteratorImpl*>(it));
    return it;
}

//-2 - nowhere
//-1 - base compact
int ICompactImpl::numContains(IVector const* const vec, int& result) const
{
    if(!vec)
    {
        ILog::report("IComapct: got nullptrm expected IVector pointer");
        return ERR_WRONG_ARG;
    }
    if(vec->getDim() != m_dim)
    {
        ILog::report("IComapct: dimension mismatch in numContains");
        return ERR_WRONG_ARG;
    }
    result = -2;
    double v,b,e;
    bool litmus;
    for(unsigned int i = 0; i < m_dim; i++)
    {
        m_begin->getCoord(i,b);
        m_end->getCoord(i,e);
        vec->getCoord(i,v);

        if (v - e > MY_EPS || b - v > MY_EPS)
        {
            int last_err;
            for (unsigned int j = 0; j < vec_compact.size(); j++)
            {
                if ((last_err = vec_compact[j]->isContains(vec, litmus)) != ERR_OK)
                {
                    ILog::report("IComapct.isContains: error in subset");
                    return last_err;
                }
                if (litmus && anti_compact[j])
                    result = -2;
                else if (litmus)
                    result = j;

            }
            return ERR_OK;
        }
    }
    bool test;
    int last_err;
    for (unsigned int j = 0; j < vec_compact.size(); j++)
    {
        if ((last_err = vec_compact[j]->isContains(vec, litmus)) != ERR_OK)
        {
            ILog::report("IComapct.isContains: error in subset");
            return last_err;
        }
        test = (litmus && anti_compact[j]) || (result && (!litmus || anti_compact[j]));
    }
    result = !test? -1: -2;
    return ERR_OK;

}

int ICompactImpl::isContains(IVector const* const vec, bool& result) const
{
    int temp;
    int error;
    if ((error = numContains(vec, temp)) != ERR_OK){
        ILog::report("IComapct.isContains: error while isContains");
        return error;
    }
    else{
        result = (temp != -2);
        return ERR_OK;
    }
}

//returning value:
//-1 - error
//0 - not subset
//1 - is subset
int ICompactImpl::isSubSet(ICompact const* const other) const//WHERE IS BOOL RESULT?????
{
    ICompactImpl* cpy = dynamic_cast<ICompactImpl*>(clone());
    if (!cpy){
        ILog::report("ICompact.isSubSet: not enough memory");
        return -1;
    }
    IIterator* it = cpy->begin();
    IVector* vec = 0;
    int code_end = ERR_OK, error;
    bool litmus;
    while (code_end == ERR_OK){
        error = getByIterator(it, vec);
        if (error != ERR_OK){
            cpy->deleteIterator(it);
            delete vec;
            delete cpy;
            ILog::report("ICompact.isSubSet: unable to get info by iterator");
            return -1;
        }
        other->isContains(vec, litmus);
        if (!litmus){
            cpy->deleteIterator(it);
            delete vec;
            delete cpy;
            return 0;
        }
        code_end = it->doStep();
    }
    cpy->deleteIterator(it);
    delete vec;
    delete cpy;
    return 0;
}

ICompactImpl::~ICompactImpl()
{
    delete m_begin;
    delete m_end;
    delete m_steps;
    for(int i = 0; i < iterators.size(); i++)
    {
        delete iterators[i];
    }
    for(int i = 0; i < vec_compact.size(); i++)
    {
        delete vec_compact[i];
    }
    iterators.clear();
    vec_compact.clear();
    anti_compact.clear();
}

int ICompactImpl::IIteratorImpl::doStep()
{
    int lastError = ERR_OK;
    ICompactImpl const* impl_m_com = dynamic_cast<const ICompactImpl*>(m_com);
    ICompactImpl const* impl_curr = impl_m_com;
    bool is_anti = false;
    if (m_com_num >= 0){
        impl_curr = dynamic_cast<const ICompactImpl*>(impl_m_com->vec_compact[m_com_num]);
        is_anti = impl_m_com->anti_compact[m_com_num];
    }

    if(m_step && !is_anti)
    {
        IVector* posVec = 0;
        impl_curr->idxToVec(m_pos,posVec);
        IVector* tmp = IVector::add(m_step,posVec);
        if(!tmp)
        {
            ILog::report("IIterator.doStep: error while make step");
            return ERR_ANY_OTHER;
        }
        int num_comp;
        lastError = impl_m_com->numContains(tmp, num_comp);
        if(lastError!= ERR_OK)
        {
            ILog::report("IIterator.doStep: error while make step ");
            return lastError;
        }
        if(num_comp == -2)
        {
            ILog::report("IIterator.doStep: compact iterator out of range");
            return ERR_OUT_OF_RANGE;
        }
        m_com_num = num_comp;
        is_anti = false;
        if (m_com_num >= 0){
            impl_curr = dynamic_cast<const ICompactImpl*>(impl_m_com->vec_compact[m_com_num]);
            is_anti = impl_m_com->anti_compact[m_com_num];
        }
        else
            impl_curr = impl_m_com;

        IVector* ttmp = 0;
        lastError = impl_curr->getNearestNeighbor(tmp, ttmp);
        if(lastError!= ERR_OK)
        {
            ILog::report("IIterator.doStep: error while make step");
            return lastError;
        }
        impl_curr->vecToIdx(ttmp,m_pos);
        delete ttmp;
        delete tmp;
        delete posVec;
    }
    else
    {
        if ((m_pos + 1 == impl_curr->m_cardinality) || !is_anti)
        {
            if ((impl_m_com->vec_compact.size() == m_com_num + 1)){
                ILog::report("IIterator.doStep: compact iterator out of range");
                return ERR_OUT_OF_RANGE;
            }
            else{
                m_pos = 0;
                m_com_num++;
                while (impl_m_com->anti_compact[m_com_num]){
                    if ((impl_m_com->vec_compact.size() == m_com_num + 1)){
                        ILog::report("IIterator.doStep: compact iterator out of range");
                        return ERR_OUT_OF_RANGE;
                    }
                    m_com_num++;
                }
            }
        }
        else
            m_pos++;
    }
    return lastError;
}

int ICompactImpl::IIteratorImpl::setStep(IVector const* const step)
{
    if(!step)
    {
        ILog::report("ICompact.setStep: got nullptr, expected IVector pointer");
        return ERR_WRONG_ARG;
    }
    if(step->getDim() != dynamic_cast<const ICompactImpl*>(m_com)->m_dim)
    {
        ILog::report("ICompact.setStep: dimension mismatch in setStep");
        return ERR_WRONG_ARG;
    }

    /*if(m_com->checkStep(step))
    {
        return ERR_INVALID_STEP;
    }*/

    IVector* tmp_step = step->clone();
    if(!tmp_step)
    {
        ILog::report("ICompact.IIterator.setStep: not enough memory");
        return ERR_MEMORY_ALLOCATION;
    }
    delete m_step;
    m_step = tmp_step;
}

ICompactImpl::IIteratorImpl::IIteratorImpl(ICompact const* const compact, int pos, IVector const* const step):ICompact::IIterator::IIterator(compact, pos, step)
{
    m_com = const_cast<ICompact*>(compact);
    m_pos = pos;
    m_com_num = -1;
    if(step)
        m_step = step->clone();
    else
        m_step = 0;
}

ICompactImpl::IIteratorImpl::~IIteratorImpl()
{
    delete m_step;
}

void ICompactImpl::vecToIdx(IVector const* const vec, unsigned int& res) const
{
    if(!vec)
        return;

    double result = 0;
    double x,b,e,s;
    double* ind = new double[m_dim];
    for(unsigned int i = 0; i < m_dim; i++)
    {
        m_begin->getCoord(i,b);
        m_end->getCoord(i,e);
        m_steps->getCoord(i,s);
        vec->getCoord(i,x);
        ind[i] = (unsigned int) (((x - b) / (e - b)) * (s + 1));

    }
    for (unsigned int i = m_dim - 1; i > 0; i--){
        m_steps->getCoord(i - 1, s);
        result += ind[i];
        result *= s + 1;
    }
    result += ind[0];
    res =(unsigned int)result;
    delete[] ind;
    ind = 0;
}

void ICompactImpl::idxToVec(unsigned int idx, IVector*& res) const
{
    double* vec = new double[m_dim];
    double b,e,s;
    for(unsigned int i = 0; i < m_dim; i++)
    {
        m_steps->getCoord(i,s);
        vec[i] = idx % (unsigned int)(s+1);
        idx /= (unsigned int)(s+1);
    }

    for(unsigned int i = 0; i < m_dim; i++)
    {
        m_begin->getCoord(i,b);
        m_end->getCoord(i,e);
        m_steps->getCoord(i,s);
        vec[i] *= (e - b)/(s+1);
    }

    for(unsigned int i = 0; i < m_dim; i++)
    {
        m_begin->getCoord(i,b);
        vec[i] += b;
    }
    res = IVector::createVector(m_dim, vec);
    delete[] vec;
    vec = 0;
}

bool ICompactImpl::checkStep(IVector const* const step)
{
    IIterator* it = new IIteratorImpl(this, 0, step);

    if(step)
    {
        IVector* it_pos = 0;
        int lastError = ERR_OK;
        if((lastError = getByIterator(it, it_pos)) != ERR_OK)
        {
            ILog::report("ICompact.begin: error while step verification");
            return false;
        }
        it->doStep();

        IVector* it_new_pos = 0;
        if((lastError = getByIterator(it, it_new_pos)) != ERR_OK)
        {
            ILog::report("ICompact.begin: error while step verification");
            return false;
        }

        double a, b, c = 0;
        for(unsigned int i = 0; i < m_dim; i++)
        {
            it_pos->getCoord(i, a);
            it_new_pos->getCoord(i, b);
            if (a == b)
                c++;
        }
        if(c == m_dim)
        {
            ILog::report("ICompact.begin: given step not suitable for iteration");
            return false;
        }

    }
    return true;
}

