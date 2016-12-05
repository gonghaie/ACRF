#include <fstream>
#include <iostream>
#include <numeric>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <iterator>
#include "utils.h"
#include "calcfunc.h"
#include "feature.h"
using namespace std;
#define epsilon 1e-3
int singnum = 0;
int honum = 9;
string hostr;
vector< vector< vector<int> > > allsegment1, allsegment2;

static const double negDmax = -DBL_MAX/100.0;

void readparm(vector<double> & wei_singlet, vector<double> & wei_doublet_adj, vector<double> & wei_doublet_ho, const string & parmfn)
{
    FILE *fp = fopen(parmfn.c_str(), "r"); 
    vector<double> parm;
    while( !feof(fp) )
    {
        char tmp[10];
        fscanf(fp, "%s", tmp);
        string str(tmp);
        if ( !feof(fp))
        {
            parm.push_back( atof(str.c_str()) );
        }
    }
    wei_singlet.resize( 2*singnum );
    wei_doublet_adj.resize( 4*singnum );
    wei_doublet_ho.resize(honum);
    copy( parm.begin(), parm.begin()+2*singnum, wei_singlet.begin() );
    copy( parm.begin()+2*singnum, parm.begin()+6*singnum, wei_doublet_adj.begin() );
    copy( parm.begin()+6*singnum, parm.begin()+6*singnum+honum, wei_doublet_ho.begin() );
    fclose(fp);
}

int main(int argc, const char ** argv)
{
    string trainfn(argv[1]);
    string parmfn(argv[2]);
    string spathfn(argv[3]);
    string sratef(argv[4]);
    string singstr(argv[5]);
    string ho(argv[6]);
    hostr = ho;
    singnum = atoi( singstr.c_str() );
    
    File file;
    vector<string> fname, path_singlet;
    string label, ss, aa, seqname, seqss, seqaa, seqlabel;
    vector < vector <vector <double> > > singletval;
    vector< vector<double> > doubletval_ho;

    string uy = "be";
    Feat feat(trainfn);
    feat.train_info(file, label, ss, aa, fname);
    feat.read_feat(file.row, singletval, doubletval_ho);
    vector<double> wei_singlet, wei_doublet_ho, wei_doublet_adj;
    readparm(wei_singlet, wei_doublet_adj, wei_doublet_ho, parmfn);
    
    double srate;
    ofstream srateout(sratef.c_str());
    
    vector< vector<int> > segment1, segment2;
    vector< vector< vector< vector< vector<double> > > > > alpha_h, beta_h;
    vector< vector< vector<double> > > alpha_e, beta_e;
    vector< vector<double> > alpha_c, beta_c;
    int file_num = file.range.size();
    for(int f_id=0; f_id<file_num; f_id++)
    {
        int beg = file.range.at(f_id).beg;
        int end = file.range.at(f_id).end;
        seqlabel = label.substr(beg, end-beg+1);
        seqss = ss.substr(beg, end-beg+1);
        seqaa = aa.substr(beg, end-beg+1);
        seqname = fname.at(f_id);
        segment1 = allsegment1.at(f_id);
        segment2 = allsegment2.at(f_id);
        Calcfunc calcfunc(wei_singlet, wei_doublet_ho, wei_doublet_adj, singletval, doubletval_ho, file.range.at(f_id).beg, segment1, segment2, hostr, seqlabel, uy, seqss);

        int start_f = 0, len_f = 0;
        int start_b = 0, len_b = 0;
        int curss_f = 0;
        int curss_b = 0;

        vector<double> inscore_f(2, negDmax), outscore_f = inscore_f, inscore_b = inscore_f, outscore_b = inscore_f;
        calcfunc.alpha_all.clear();
        calcfunc.beta_all.clear();
        
        int segments_num = segment1.size();
        calcfunc.fnewscore.resize(segments_num);
        
        cout<<"\nf_id forward: "<<f_id<<endl;
        for(int i=0; i<segments_num; i++)
        {
            start_f = segment1.at(i).at(0);
            len_f = segment1.at(i).at(1);
            curss_f = segment1.at(i).at(2);
            if (start_f != 0)
                inscore_f = outscore_f; //上一次的输出作为当前输入

            vector< vector< vector< vector< vector<double> > > > > alpha_tmp(0);
            if (curss_f == 1)
            {
                calcfunc.calc_alpha_h(alpha_h, start_f, len_f, inscore_f, outscore_f);
                calcfunc.alpha_all.push_back( alpha_h );
            }
            if (curss_f == 2)
            {
                calcfunc.calc_alpha_e(alpha_e, start_f, len_f, inscore_f, outscore_f);
                alpha_tmp.resize(1);
                alpha_tmp.at(0).push_back( alpha_e );
                calcfunc.alpha_all.push_back( alpha_tmp );
            }
            if (curss_f == 0)
            {
                calcfunc.calc_alpha_c(alpha_c, start_f, len_f, inscore_f, outscore_f);
                alpha_tmp.resize(1);
                alpha_tmp.at(0).resize(1);
                alpha_tmp.at(0).at(0).push_back( alpha_c );
                calcfunc.alpha_all.push_back( alpha_tmp );
            }
            calcfunc.fnewscore.at(i) = outscore_f;
        }

        cout<<"\nf_id backward: "<<f_id<<endl;
        for(int i=0; i<segments_num; i++)
        {
            start_b = segment2.at(i).at(0);
            len_b = segment2.at(i).at(1);
            curss_b = segment2.at(i).at(2);
            if (start_b != calcfunc.T-1)
                inscore_b = outscore_b; //上一次的输出作为当前输入
            
            vector< vector< vector< vector< vector<double> > > > > beta_tmp(0);

            if (curss_b == 1)
            {
                calcfunc.calc_beta_h(beta_h, start_b, len_b, inscore_b, outscore_b);
                calcfunc.beta_all.push_back( beta_h ); 
            }
            if (curss_b == 2)
            {
                calcfunc.calc_beta_e(beta_e, start_b, len_b, inscore_b, outscore_b);
                beta_tmp.resize(1);
                beta_tmp.at(0).push_back( beta_e );
                calcfunc.beta_all.push_back( beta_tmp ); 
            }
            if (curss_b == 0)
            {
                calcfunc.calc_beta_c(beta_c, start_b, len_b, inscore_b, outscore_b);
                beta_tmp.resize(1);
                beta_tmp.at(0).resize(1);
                beta_tmp.at(0).at(0).push_back( beta_c );
                calcfunc.beta_all.push_back( beta_tmp ); 
            }
        }

        for(int i=0; i<segments_num; i++)
            if ( i < (segments_num>>1) )
                swap( calcfunc.beta_all.at(i), calcfunc.beta_all.at(segments_num-i-1) );
        
        double outscore_f_sum = negDmax;
        double outscore_b_sum = negDmax;
        utility utils;
        int out_num = outscore_f.size();
        for(int i=0; i<out_num; i++)
        {
            outscore_f_sum = utils.logadd(outscore_f_sum, outscore_f.at(i));
            outscore_b_sum = utils.logadd(outscore_b_sum, outscore_b.at(i));
        }
        //cout<<"final partition: "<<setprecision(3)<<outscore_f_sum<<" "<<outscore_b_sum<<endl;
        //assert( fabs( outscore_f_sum - outscore_b_sum ) < epsilon );
        calcfunc.partition = outscore_f_sum;
        
        for(int cur=0; cur<segments_num; cur++)
        {
            start_f = segment1.at(cur).at(0);
            len_f = segment1.at(cur).at(1);
            curss_f = segment1.at(cur).at(2);
            
            if (curss_f == 1)
                calcfunc.margin_h( cur, start_f, len_f );
            if (curss_f == 2)
                calcfunc.margin_e( cur, start_f, len_f );
            if (curss_f == 0)
                calcfunc.margin_c( cur, start_f, len_f );
        }
    
        calcfunc.psudolikelihood(spathfn, seqname, srate, seqss, seqaa);
        srateout<<seqname<<" "<<srate<<endl;
    }
    srateout.close();
    return 0;
}


