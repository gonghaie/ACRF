#ifndef _FEATURE_H
#define _FEATURE_H
#include <map>
#include <set>
#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cstdio>
#include <cmath>
using namespace std;

extern int singnum;
extern int honum;
extern string hostr;
extern vector< vector< vector<int> > > allsegment1, allsegment2;

class CurInfo
{
public:
    string ssstr;
    string fn;
};

class Couple
{
public:
    int beg;
    int end;
};

class File
{
public:
    int row;
    vector<Couple> range;
};

class Feat
{
public:
    string trainfn;
public:
    Feat(string trainfile1);
    ~Feat();

    void train_info(
            File & file, 
            string & label,
            string & ss,
            string & aa,
            vector<string> & fname
            );
    
    void read_feat(
        int row, 
        vector < vector <vector <double> > > & singletval, 
        vector <vector<double> > & doubletval_ho
        );
};
#endif
