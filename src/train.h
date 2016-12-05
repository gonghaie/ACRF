#ifndef _TRAIN_H
#define _TRAIN_H
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <assert.h>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include "lbfgs.h"
#include "calcfunc.h"
#include "feature.h"
#include "utils.h"

extern int singnum;
extern int honum;

extern string trainfn;
extern File file;
extern int thread_num;
extern vector<string> fname;
extern string label, ss, aa, uy;

extern vector < vector <vector <double> > > singletval;
extern vector< vector<double> > doubletval_ho;

struct Allparm
{
    vector<double> wei_singlet;
    vector<double> wei_doublet_ho;
    vector<double> wei_doublet_adj;
    int start;
    int end;
    vector<double> grad;
    double obj;
    vector<double> ratevec;

    Allparm(
            vector<double> & wei_singlet1,
            vector<double> & wei_doublet_ho1,
            vector<double> & wei_doublet_adj1,
            int & start1,
            int & end1,
            vector<double> & grad1,
            double & obj1,
            vector<double> & ratevec1
           ) : 
        wei_singlet(wei_singlet1), wei_doublet_ho(wei_doublet_ho1), wei_doublet_adj(wei_doublet_adj1),
        start(start1), end(end1), grad(grad1), obj(obj1), ratevec(ratevec1) {};
};

class Obj_func
{
protected:
    lbfgsfloatval_t *m_x;
public:
    Obj_func(
            );
    virtual ~Obj_func();
    int run(
            int N, 
            vector<double> & parm
           );
protected:
    static void * _calc_grad_obj( void * allparm);
    void * calc_grad_obj( void * allparm);
    
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        );
    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        );
    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        );
    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        );
};
#endif
