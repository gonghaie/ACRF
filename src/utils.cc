#include <iostream>
#include <cmath>
#include "utils.h"
using namespace std;
utility::~utility()
{
}
double utility::logadd(const double & a, const double & b)
{
    double res;
    if(b<a)
        res = ( a + log( exp(b-a) + 1 ) );
    else
        res = ( b + log( exp(a-b) + 1 ) );
    return res;
}
/*
 * calculate the index sequence of the real state sequence; initialize the margins
 * input: y
 * output: y_id
 */
void utility::getid(const string & y, const string & uy, vector<int> & y_id)
{
    int T = y.size();
    y_id.resize(T, 0);
    for(int t=0; t<T; t++)
    {
        if( y.at(t) == uy.at(0) )
            y_id.at(t) = 0;
        else y_id.at(t) = 1;
    }
}
