#ifndef _UTILS_H
#define _UTILS_H
#include <vector>
#include <string>
using namespace std;
class utility
{
public:
    ~utility();
    double logadd(const double & a, const double & b);
    //void getid(vector<string> y, vector<string> uy, vector<int> & y_id);
    void getid(const string & y, const string & uy, vector<int> & y_id);
};



#endif
