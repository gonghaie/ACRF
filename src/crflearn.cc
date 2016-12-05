#include <fstream>
#include <iostream>
#include <numeric>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <time.h>
#include <iterator>
#include "utils.h"
#include "feature.h"
#include "calcfunc.h"
#include "train.h"
using namespace std;

int singnum = 0;
int honum = 9;

string trainfn;
File file;
int thread_num;
vector<string> fname;
string label, ss, aa, uy, hostr;
vector < vector <vector <double> > > singletval;
vector< vector<double> > doubletval_ho;
vector< vector< vector<int> > > allsegment1, allsegment2;

int main(int argc, const char ** argv)
{
    clock_t mainstart = clock();
    string trainfn1(argv[1]);
    string parmfn(argv[2]);
    string thread_str(argv[3]);
    string singstr(argv[4]);
    string ho(argv[5]);
    hostr = ho;
    singnum = atoi( singstr.c_str() );

    thread_num = atoi(thread_str.c_str());

    trainfn = trainfn1;
    string uy_tmp("be");
    uy = uy_tmp;
    Feat feat(trainfn);
    feat.train_info(file, label, ss, aa, fname);
    feat.read_feat(file.row, singletval, doubletval_ho);
    int N = singnum*6 + honum;
    Obj_func obj;
    vector<double> parm;
    int ret = obj.run(N, parm);
    
    if (ret)
        cout<<"\n"<<ret<<" failed!\n";
    else
        cout<<"\n"<<ret<<" passed!\n";
    ofstream outf(parmfn.c_str());
    for(int i=0; i<N; i++)
        outf<<parm.at(i)<<" ";
    outf.close();
    clock_t mainfinish = clock();
    /*
    cout<<"whole training time (s): "<<double(mainfinish-mainstart)/1000000.0<<endl;
    cout<<"whole training time (m): "<<double(mainfinish-mainstart)/60000000.0<<endl;
    */
    cout<<"whole training time (s): "<<double(mainfinish-mainstart)/CLOCKS_PER_SEC<<endl;
    cout<<"whole training time (m): "<<double(mainfinish-mainstart)/CLOCKS_PER_SEC/60<<endl;
    cout<<"whole training time (h): "<<double(mainfinish-mainstart)/CLOCKS_PER_SEC/3600<<endl;
    cout<<"whole training time (d): "<<double(mainfinish-mainstart)/CLOCKS_PER_SEC/3600/24<<endl;
    return 0;
}
