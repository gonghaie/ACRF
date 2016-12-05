#ifndef _CALCFUNC_H
#define _CALCFUNC_H
#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iterator>
#include <ctime>
#include <cfloat>
#include <iomanip>
#include "utils.h"
#include "feature.h"
using namespace std;
extern int singnum;
extern int honum;
extern string hostr;
class Calcfunc
{
public:
    vector<double> wei_singlet;
    vector<double> wei_doublet_ho;
    vector<double> wei_doublet_adj;
    vector < vector <vector <double> > > singletval;
    vector< vector<double> > doubletval_ho;
    int beg;
    vector< vector<int> > segment1;
    vector< vector<int> > segment2;
    string hostr;
    string label;
    string uy;
    string seqss;

    vector< vector<double> > fnewscore; // save the outscore_f of each segment for margins' caclulating
    vector<int> y_id;
    int T;
    int N;
    utility utils;
    double partition;
    
    vector< vector< vector< vector< vector< vector<double> > > > > > alpha_all;
    vector< vector< vector< vector< vector< vector<double> > > > > > beta_all;
    
    vector< vector<double> > table1;
    vector< vector< vector<double> > > table2;
    vector< vector< vector<double> > > table3;
    vector< vector< vector<double> > > table4;
    vector< vector< vector<double> > > table5;
    vector<double> sigma1;
    vector<double> sigma2;
    vector<double> sigma3;
    vector<double> sigma4;
    vector<double> sigma5;
    
public:
    virtual ~Calcfunc();
    Calcfunc(
            const vector<double> & wei_singlet1, 
            const vector<double> & wei_doublet_ho1, 
            const vector<double> & wei_doublet_adj1, 
            const vector < vector <vector <double> > > & singletval1, 
            const vector< vector<double> > & doubletval_ho1,
            const int & beg1,
            const vector< vector<int> > & segment1_1,
            const vector< vector<int> > & segment2_1,
            const string & hostr1,
            const string & label1,
            const string & uy1,
            const string & seqss
            );
    void calc_yid();
    double singlet2doublet(
        int newt,
        int si,
        int sj,
        int feat_index
        ) const;

    void calc_alpha_c(
        vector< vector<double> > & alpha_c,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        );

    void calc_beta_c(
        vector< vector<double> > & beta_c,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        );
    
    void calc_alpha_e(
        vector< vector< vector<double> > > & alpha_e,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        );

    void calc_beta_e(
        vector< vector< vector<double> > > & beta_e,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        );

    void calc_alpha_h(
        vector< vector< vector< vector< vector<double> > > > > & alpha_h,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        );
    void calc_beta_h(
        vector< vector< vector< vector< vector<double> > > > > & beta_h,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        );
    void margin_h(const int & cur, const int & start_f, const int & len );
    void margin_e(const int & cur, const int & start_f, const int & len );
    void margin_c(const int & cur, const int & start_f, const int & len );
    void calc_grad1(vector<double> & gradient) const;
    double calc_obj() const;
    //double psudolikelihood( const string & seqname ) const; // return precision of current protein
    void psudolikelihood( 
        const string & path_sing_file, 
        const string & seqname, 
        double & srate,
        const string & seqss, 
        const string & seqaa
        ) const;
};
#endif
