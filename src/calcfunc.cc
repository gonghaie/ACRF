#include "calcfunc.h"
#define epsilon 1e-3
static const double nearzero = 1e-10;
static const double negDmax = -DBL_MAX/100.0;
using namespace std;
Calcfunc::Calcfunc(
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
        const string & seqss1
        )
{
    T = label1.size();
    N = uy1.size();
    wei_singlet = wei_singlet1;
    wei_doublet_ho = wei_doublet_ho1;
    wei_doublet_adj = wei_doublet_adj1;
    singletval = singletval1;
    doubletval_ho = doubletval_ho1;
    beg = beg1;
    segment1 = segment1_1;
    segment2 = segment2_1;
    hostr = hostr1;
    label = label1;
    uy = uy1;
    seqss = seqss1;

    vector< vector<double> > empty1;
    vector< vector<double> > tmp1( T, vector<double>(N, 0.5) );
    table1 = tmp1;
    swap(tmp1, empty1);
    
    vector< vector< vector<double> > > empty2;
    vector< vector< vector<double> > > tmp2( T, vector< vector<double> >(N, vector<double> (N, 0.25) ) );
    table2 = tmp2;
    swap(tmp2, empty2);

    table3 = table2;
    table4 = table2;
    table5 = table2;

    sigma1.resize(T, 0);
    sigma2.resize(T, 0);
    sigma3.resize(T, 0);
    sigma4.resize(T, 0);
    sigma5.resize(T, 0);

}
Calcfunc::~Calcfunc()
{
}

void Calcfunc::calc_yid()
{
    utils.getid(label, uy, y_id);
}

double Calcfunc::singlet2doublet(
        int newt,
        int si,
        int sj,
        int feat_index
        ) const
{
    if (hostr.at(0)=='0') return 0;
    int id = feat_index % singnum;
    double tmpval = singletval.at(newt).at(0).at(id);
    if ( feat_index >= ( si*2+sj ) * singnum && feat_index < ( si*2+sj+1 ) * singnum )
        return tmpval;
    else return 0;
}

/*
 * objective: calculate the forward matrix of a given segment whose SS is C
 * input: inscore, a segment
 * output: outscore, alpha_c
 */
void Calcfunc::calc_alpha_c(
        vector< vector<double> > & alpha_c,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        )
{
    vector< vector<double> > empty_c;
    vector< vector<double> > tmp_c( len, vector<double> (N, negDmax) );
    alpha_c = tmp_c;
    swap(tmp_c, empty_c);
    

    int t = 0;
    int newt = t+beg+start;
    
    vector<double> newscore2(N, 0);//inscore + adjacent doublet
    if(start!=0)
    {
        for(int i=0; i<N; i++)
        {
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum4 = 0;
                
                for(int k1=0; k1<4*singnum; k1++) 
                {
                    sum4 += wei_doublet_adj.at(k1) * singlet2doublet(newt, j, i, k1);
                }

                double add = inscore.at(j) + sum4;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            } 
            newscore2.at(i) += dd;
        }
    }

    for(int i=0; i<N; i++)
    {
        double sum1 = 0; 
        for(int k1=0; k1<2*singnum; k1++)
        {
            sum1 += wei_singlet.at(k1) * ( 
                    singletval.at(newt).at(i).at(k1) 
                    );
        }
        alpha_c.at(t).at(i) = sum1 + newscore2.at(i);
    }
    
    for(int t=1; t<len; t++)
    {
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            double sum1 = 0;
            for(int k1=0; k1<2*singnum; k1++)
            {
                sum1 += wei_singlet.at(k1) * singletval.at(newt).at(i).at(k1);
            }
            alpha_c.at(t).at(i) = sum1;
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum3 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(newt, j, i, k1);
                }
                double add = alpha_c.at(t-1).at(j)+sum3;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            }
            alpha_c.at(t).at(i) += dd;
        }
    }

    outscore.resize(N, negDmax);
    for(int i=0; i<N; i++)
    {
        outscore.at(i) = alpha_c.at(len-1).at(i);
    }
}

/*
 * objective: calculate the backward matrix of a given segment whose SS is C
 * input: inscore, a segment
 * output: outscore, beta_c
 */
void Calcfunc::calc_beta_c(
        vector< vector<double> > & beta_c,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        )
{
    
    
    vector< vector<double> > empty_c;
    vector< vector<double> > tmp_c( len, vector<double> (N, negDmax) );
    beta_c = tmp_c;
    swap(tmp_c, empty_c);
    
    int t = len-1;
    int newt = beg+start-(len-1-t);
    vector<double> newscore2(N, 0); //inscore + adjacent doublet
    if(start!=T-1)
    {
        for(int i=0; i<N; i++)
        {
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum4 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum4 += wei_doublet_adj.at(k1) * singlet2doublet(newt+1, i, j, k1);
                }
                double add = inscore.at(j) + sum4;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            } 
            newscore2.at(i) += dd;
        }
    }
    
    for(int i=0; i<N; i++)
    {
        double sum1 = 0; 
        for(int k1=0; k1<2*singnum; k1++)
        {
            sum1 += wei_singlet.at(k1) * ( 
                    singletval.at(newt).at(i).at(k1) 
                    );
        }
        beta_c.at(t).at(i) = sum1 + newscore2.at(i);
    }

    for(int t=len-2; t>=0; t--)
    {
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            double sum1 = 0;
            for(int k1=0; k1<2*singnum; k1++)
            {
                sum1 += wei_singlet.at(k1) * singletval.at(newt).at(i).at(k1);
            }
            beta_c.at(t).at(i) = sum1;
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {

                double sum3 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(newt+1, i, j, k1);
                }

                double add = beta_c.at(t+1).at(j)+sum3;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            }
            beta_c.at(t).at(i) += dd;
        }
        
    }
    
    outscore.resize(N, negDmax);
    for(int i=0; i<N; i++)
    {
        outscore.at(i) = beta_c.at(0).at(i);
    }
}

/*
 * objective: calculate the forward matrix of a given segment whose SS is E
 * input: inscore, a segment
 * output: outscore, alpha_e
 */
void Calcfunc::calc_alpha_e(
        vector< vector< vector<double> > >  & alpha_e,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        )
{
    vector< vector< vector<double> > > empty_e;
    vector< vector< vector<double> > > tmp_e( len, vector< vector<double> > ( N, vector<double> ( N, negDmax) ) );
    alpha_e = tmp_e;
    swap(tmp_e, empty_e);
    
    int t = 0;
    int newt = t+beg+start;
    
    double newscore1 = 0;
    vector<double> newscore2(N, 0);
    if(start!=0)
    {
        for(int i=0; i<N; i++)
        {
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum4 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                    sum4 += wei_doublet_adj.at(k1) * singlet2doublet(newt, j, i, k1);
                double add = inscore.at(j) + sum4;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            } 
            newscore2.at(i) += dd;
        }
    }

    
    if (len>=1)
    {
        int t = 0;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double sum1 = 0; 
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(newt).at(j).at(k1) 
                            );
                }
                alpha_e.at(t).at(i).at(j) = sum1 + newscore2.at(j);
            }
        }
    }
    
    if (len>=2)
    {
        int t = 1;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double sum1 = 0; 
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(newt).at(j).at(k1) 
                            + singletval.at(newt-1).at(i).at(k1) 
                            );
                }

                double sum3 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(newt, i, j, k1);
                }
                alpha_e.at(t).at(i).at(j) = sum1 + sum3 + newscore2.at(i);
            }
        }
    }
    
    if (len>=3)
    {
        int t = 2;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double dd = negDmax;
                for(int b=0; b<N; b++)
                {
                    double sum1=0, sum2=0;
                    for(int k1=0; k1<2*singnum; k1++)
                    {
                        sum1 += wei_singlet.at(k1) * ( 
                                singletval.at(newt).at(j).at(k1) 
                                + singletval.at(newt-1).at(i).at(k1) 
                                + singletval.at(newt-2).at(b).at(k1) 
                                );
                    }
                    for(int p=0; p<3; p++)
                        sum2 += wei_doublet_ho.at(p) * ( (b==j) ? doubletval_ho.at(newt).at(p) : 0 );

                    double sum3 = 0;
                    for(int k1=0; k1<4*singnum; k1++)
                    {
                        sum3 += wei_doublet_adj.at(k1) * ( 
                                singlet2doublet(newt, i, j, k1)
                                + singlet2doublet(newt-1, b, i, k1)
                                );
                    }
                    
                    double add = sum1+sum2+sum3 + newscore2.at(b);
                    if ( fabs(add) < nearzero ) add = negDmax;
                    dd = utils.logadd( dd, add );
                }
                alpha_e.at(t).at(i).at(j) = dd;
            }
        }
    }

    for(int t=3; t<len; t++)
    {
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double sum1 = 0; 
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(newt).at(j).at(k1) 
                            + singletval.at(newt-1).at(i).at(k1) 
                            );
                }
                alpha_e.at(t).at(i).at(j) = sum1;
                double dd = negDmax;
                for(int a=0; a<N; a++)
                {
                    for(int b=0; b<N; b++)
                    {
                        double sum2 = 0;
                        for(int p=0; p<3; p++)
                            sum2 += wei_doublet_ho.at(p) * ( 
                                    ( (b==j) ? doubletval_ho.at(newt).at(p) : 0 )
                                    + ( (a==i) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                    );

                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * ( 
                                    singlet2doublet(newt, i, j, k1)
                                    + singlet2doublet(newt-1, b, i, k1)
                                    );
                        }
                        double add = alpha_e.at(t-2).at(a).at(b)+sum2+sum3;
                        if ( fabs(add) < nearzero ) add = negDmax;
                        dd = utils.logadd( dd, add );
                    }
                }
                alpha_e.at(t).at(i).at(j) += dd;
            }
        }
    }

    outscore.resize(N, negDmax);
    for(int i=0; i<N; i++)
    {
        outscore.at(i) = utils.logadd(alpha_e.at( len-1 ).at(0).at(i), alpha_e.at( len-1 ).at(1).at(i) );
    }
}

/*
 * objective: calculate the backward matrix of a given segment whose SS is E
 * input: inscore, a segment
 * output: outscore, beta_e
 */
void Calcfunc::calc_beta_e(
        vector< vector< vector<double> > > & beta_e,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        )
{
    
    
    vector< vector< vector<double> > > empty_e;
    vector< vector< vector<double> > > tmp_e( len, vector< vector<double> > ( N, vector<double> ( N, negDmax) ) );
    beta_e = tmp_e;
    swap(tmp_e, empty_e);
    
    int t = len-1;
    int newt = beg+start-(len-1-t);
    
    double newscore1 = 0;
    vector<double> newscore2(N, 0);
    if(start!=T-1)
    {
        for(int i=0; i<N; i++)
        {
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum4 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                    sum4 += wei_doublet_adj.at(k1) * singlet2doublet(newt+1, i, j, k1);
                double add = inscore.at(j) + sum4;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            } 
            newscore2.at(i) += dd;
        }
    }

    if (len>=1)
    {
        int t = len-1;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double sum1 = 0; 
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(newt).at(j).at(k1) 
                            );
                }
                beta_e.at(t).at(i).at(j) = sum1 + newscore2.at(j);
            }
        }
    }
    
    if (len>=2)
    {
        int t = len-2;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double sum1 = 0; 
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(newt).at(j).at(k1) 
                            + singletval.at(newt+1).at(i).at(k1) 
                            );
                }

                double sum3 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(newt+1, j, i, k1);
                }

                beta_e.at(t).at(i).at(j) = sum1 + sum3 + newscore2.at(i);
            }
        }
    }
    
    if (len>=3)
    {
        int t = len-3;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double dd = negDmax;
                for(int b=0; b<N; b++)
                {
                    double sum1=0, sum2=0;
                    for(int k1=0; k1<2*singnum; k1++)
                    {
                        sum1 += wei_singlet.at(k1) * ( 
                                singletval.at(newt).at(j).at(k1) 
                                + singletval.at(newt+1).at(i).at(k1) 
                                + singletval.at(newt+2).at(b).at(k1) 
                                );
                    }
                    for(int p=0; p<3; p++)
                        sum2 += wei_doublet_ho.at(p) * ( (j==b) ? doubletval_ho.at(newt+2).at(p) : 0 );

                    double sum3 = 0;
                    for(int k1=0; k1<4*singnum; k1++)
                    {
                        sum3 += wei_doublet_adj.at(k1) * ( 
                                singlet2doublet(newt+1, j, i, k1)
                                + singlet2doublet(newt+2, i, b, k1)
                                );
                    }

                    double add = sum1+sum2+sum3 + newscore2.at(b);
                    if ( fabs(add) < nearzero ) add = negDmax;
                    dd = utils.logadd( dd, add );
                }
                beta_e.at(t).at(i).at(j) = dd;
            }
        }
    }
    for(int t=len-4; t>=0; t--)
    {
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                double sum1 = 0; 
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(newt).at(j).at(k1) 
                            + singletval.at(newt+1).at(i).at(k1) 
                            );
                }
                beta_e.at(t).at(i).at(j) = sum1;
                double dd = negDmax;
                for(int a=0; a<N; a++)
                {
                    for(int b=0; b<N; b++)
                    {
                        double sum2 = 0;
                        for(int p=0; p<3; p++)
                            sum2 += wei_doublet_ho.at(p) * ( 
                                ( (j==b) ? doubletval_ho.at(newt+2).at(p) : 0 )
                                + ( (i==a) ? doubletval_ho.at(newt+3).at(p) : 0 )
                                );
                        
                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * ( 
                                    singlet2doublet(newt+1, j, i, k1)
                                    + singlet2doublet(newt+2, i, b, k1)
                                    );
                        }
                        double add = beta_e.at(t+2).at(a).at(b)+sum2+sum3;
                        if ( fabs(add) < nearzero ) add = negDmax;
                        dd = utils.logadd( dd, add );
                    }
                }
                beta_e.at(t).at(i).at(j) += dd;
            }
        }
    }

    outscore.resize(N, negDmax);
    for(int i=0; i<N; i++)
    {
        outscore.at(i) = utils.logadd(beta_e.at(0).at(0).at(i), beta_e.at(0).at(1).at(i) );
    }
}

/*
 * objective: calculate the forward matrix of a given segment whose SS is H
 * input: inscore, a segment
 * output: outscore, alpha_h
 */
void Calcfunc::calc_alpha_h(
        vector< vector< vector< vector< vector<double> > > > > & alpha_h,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        )
{
    
    
    vector< vector< vector< vector< vector<double> > > > > empty_h;
    vector< vector< vector< vector< vector<double> > > > > tmp_h( len, vector< vector< vector< vector<double> > > > ( N, vector< vector< vector<double> > > ( N, vector< vector<double> > (N, vector<double> (N, negDmax) ) )  ) );
    alpha_h = tmp_h;
    swap(tmp_h, empty_h);
    
    int t = 0;
    int newt = t+beg+start;
    vector<double> newscore2(N, 0);
    if(start!=0)
    {
        for(int i=0; i<N; i++)
        {
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum4 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                    sum4 += wei_doublet_adj.at(k1) * singlet2doublet(newt, j, i, k1);
                double add = inscore.at(j) + sum4;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            } 
            newscore2.at(i) += dd;
        }
    }
    
    if (len>=1)
    {
        int t = 0;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    );
                        }
                        alpha_h.at(t).at(i).at(j).at(m).at(n) = sum1 + newscore2.at(n);
                    }
                }
            }
        }
    }
    
    if (len>=2)
    {
        int t = 1;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt-1).at(m).at(k1) 
                                    );
                        }
                        
                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * singlet2doublet(newt, m, n, k1);
                        }

                        alpha_h.at(t).at(i).at(j).at(m).at(n) = sum1 + sum3 + newscore2.at(m);
                    }
                }
            }
        }
    }
    
    if (len>=3)
    {
        int t = 2;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt-1).at(m).at(k1) 
                                    + singletval.at(newt-2).at(j).at(k1) 
                                    );
                        }

                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * ( 
                                    singlet2doublet(newt, m, n, k1)
                                    + singlet2doublet(newt-1, j, m, k1)
                                    );
                        }

                        alpha_h.at(t).at(i).at(j).at(m).at(n) = sum1 + sum3 + newscore2.at(j);
                    }
                }
            }
        }
    }

    if (len>=4)
    {
        int t = 3;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0, sum2_hh=0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt-1).at(m).at(k1) 
                                    + singletval.at(newt-2).at(j).at(k1) 
                                    + singletval.at(newt-3).at(i).at(k1) 
                                    );
                        }
                        for(int p=3; p<6; p++)
                            sum2_hh += wei_doublet_ho.at(p) * ( (i==n) ? doubletval_ho.at(newt).at(p) : 0 );
                        
                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * ( 
                                    singlet2doublet(newt, m, n, k1)
                                    + singlet2doublet(newt-1, j, m, k1)
                                    + singlet2doublet(newt-2, i, j, k1)
                                    );
                        }

                        alpha_h.at(t).at(i).at(j).at(m).at(n) = sum1 + sum3 + sum2_hh + newscore2.at(i);
                    }
                }
            }
        }
    }

    if (len>=5)
    {
        int t = 4;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double dd = negDmax;
                        for(int c=0; c<N; c++)
                        {
                            double sum1=0, sum2=0, sum2_hh=0;
                            for(int k1=0; k1<2*singnum; k1++)
                            {
                                sum1 += wei_singlet.at(k1) * ( 
                                        singletval.at(newt).at(n).at(k1) 
                                        + singletval.at(newt-1).at(m).at(k1) 
                                        + singletval.at(newt-2).at(j).at(k1) 
                                        + singletval.at(newt-3).at(i).at(k1) 
                                        + singletval.at(newt-4).at(c).at(k1) 
                                        );
                            }
                            for(int p=6; p<9; p++)
                                sum2 += wei_doublet_ho.at(p) * ( (c==n) ? doubletval_ho.at(newt).at(p) : 0 );
                            
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) *( 
                                    ( (i==n) ? doubletval_ho.at(newt).at(p) : 0 )
                                    + ( (c==m) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                    );
                            
                            double sum3 = 0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * (
                                        singlet2doublet(newt, m, n, k1)
                                        + singlet2doublet(newt-1, j, m, k1)
                                        + singlet2doublet(newt-2, i, j, k1)
                                        + singlet2doublet(newt-3, c, i, k1)
                                        );
                            }

                            double add = sum1+sum2+sum2_hh+sum3 + newscore2.at(c);
                            if ( fabs(add) < nearzero ) add = negDmax;
                            dd = utils.logadd( dd, add );
                        }
                        alpha_h.at(t).at(i).at(j).at(m).at(n) = dd;
                    }
                }
            }
        }
    }

    if (len>=6)
    {
        int t = 5;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double dd = negDmax;
                        for(int b=0; b<N; b++)
                        {
                            for(int c=0; c<N; c++)
                            {
                                double sum1=0, sum2=0, sum2_hh=0;
                                for(int k1=0; k1<2*singnum; k1++)
                                {
                                    sum1 += wei_singlet.at(k1) * ( 
                                            singletval.at(newt).at(n).at(k1) 
                                            + singletval.at(newt-1).at(m).at(k1) 
                                            + singletval.at(newt-2).at(j).at(k1) 
                                            + singletval.at(newt-3).at(i).at(k1) 
                                            + singletval.at(newt-4).at(c).at(k1) 
                                            + singletval.at(newt-5).at(b).at(k1) 
                                            );
                                }
                                for(int p=6; p<9; p++)
                                    sum2 += wei_doublet_ho.at(p) * ( 
                                            ( ( c==n ) ? doubletval_ho.at(newt).at(p) : 0 )
                                        + ( ( b==m ) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                        );
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) *( 
                                        ( ( i==n ) ? doubletval_ho.at(newt).at(p) : 0)
                                        + ( ( c==m ) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                        + ( ( b==j ) ? doubletval_ho.at(newt-2).at(p) : 0 )
                                        );

                                double sum3 = 0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * (
                                            singlet2doublet(newt, m, n, k1)
                                            + singlet2doublet(newt-1, j, m, k1)
                                            + singlet2doublet(newt-2, i, j, k1)
                                            + singlet2doublet(newt-3, c, i, k1)
                                            + singlet2doublet(newt-4, b, c, k1)
                                            );
                                }

                            double add = sum1+sum2+sum2_hh+sum3 + newscore2.at(b);
                            if ( fabs(add) < nearzero ) add = negDmax;
                            dd = utils.logadd( dd, add );
                            }
                        }
                        alpha_h.at(t).at(i).at(j).at(m).at(n) = dd;
                    }
                }       
            }
        }
    }

    if (len>=7)
    {
        int t = 6;
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double dd = negDmax;
                        for(int a=0; a<N; a++)
                        {
                            for(int b=0; b<N; b++)
                            {
                                for(int c=0; c<N; c++)
                                {
                                    double sum1=0, sum2=0, sum2_hh=0;
                                    for(int k1=0; k1<2*singnum; k1++)
                                    {
                                        sum1 += wei_singlet.at(k1) * ( 
                                                singletval.at(newt).at(n).at(k1) 
                                                + singletval.at(newt-1).at(m).at(k1) 
                                                + singletval.at(newt-2).at(j).at(k1) 
                                                + singletval.at(newt-3).at(i).at(k1) 
                                                + singletval.at(newt-4).at(c).at(k1) 
                                                + singletval.at(newt-5).at(b).at(k1) 
                                                + singletval.at(newt-6).at(a).at(k1) 
                                                );
                                    }
                                    for(int p=6; p<9; p++)
                                        sum2 += wei_doublet_ho.at(p) * ( 
                                            ( (c==n) ? doubletval_ho.at(newt).at(p) : 0 )
                                            + ( ( b==m ) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                            + ( ( a==j ) ? doubletval_ho.at(newt-2).at(p) : 0 )
                                            );
                                    for(int p=3; p<6; p++)
                                        sum2_hh += wei_doublet_ho.at(p) *( 
                                            ( (i==n) ? doubletval_ho.at(newt).at(p) : 0 )
                                            + ( ( c==m ) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                            + ( ( b==j ) ? doubletval_ho.at(newt-2).at(p) : 0 )
                                            + ( ( a==i ) ? doubletval_ho.at(newt-3).at(p) : 0 )
                                            );

                                    double sum3 = 0;
                                    for(int k1=0; k1<4*singnum; k1++)
                                    {
                                        sum3 += wei_doublet_adj.at(k1) * (
                                                singlet2doublet(newt, m, n, k1)
                                                + singlet2doublet(newt-1, j, m, k1)
                                                + singlet2doublet(newt-2, i, j, k1)
                                                + singlet2doublet(newt-3, c, i, k1)
                                                + singlet2doublet(newt-4, b, c, k1)
                                                + singlet2doublet(newt-5, a, b, k1)
                                                );
                                    }

                                    double add = sum1+sum2+sum2_hh+sum3 + newscore2.at(a);
                                    if ( fabs(add) < nearzero ) add = negDmax;
                                    dd = utils.logadd( dd, add );
                                }
                            }
                            alpha_h.at(t).at(i).at(j).at(m).at(n) = dd;
                        }
                    }
                }       
            }
        }
    }

    for(int t=7; t<len; t++)
    {
        int newt = t+beg+start;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt-1).at(m).at(k1) 
                                    + singletval.at(newt-2).at(j).at(k1) 
                                    + singletval.at(newt-3).at(i).at(k1) 
                                    );
                        }
                        alpha_h.at(t).at(i).at(j).at(m).at(n) = sum1;
                        double dd = negDmax;
                        for(int a=0; a<N; a++)
                        {
                            for(int b=0; b<N; b++)
                            {
                                for(int c=0; c<N; c++)
                                {
                                    for(int d=0; d<N; d++)
                                    {
                                        double sum2 = 0, sum2_hh = 0;
                                        for(int p=6; p<9; p++)
                                            sum2 += wei_doublet_ho.at(p) * ( 
                                                + ( ( d==n ) ? doubletval_ho.at(newt).at(p) : 0 )
                                                + ( ( c==m ) ? doubletval_ho.at(newt-1).at(p) : 0 ) 
                                                + ( ( b==j ) ? doubletval_ho.at(newt-2).at(p) : 0 ) 
                                                + ( ( a==i ) ? doubletval_ho.at(newt-3).at(p) : 0 ) 
                                                );
                                        for(int p=3; p<6; p++)
                                            sum2_hh += wei_doublet_ho.at(p) *( 
                                                + ( ( i==n ) ? doubletval_ho.at(newt).at(p): 0 )
                                                + ( ( d==m ) ? doubletval_ho.at(newt-1).at(p) : 0 )
                                                + ( ( c==j ) ? doubletval_ho.at(newt-2).at(p) : 0 )
                                                + ( ( b==i ) ? doubletval_ho.at(newt-3).at(p) : 0 )
                                                );
                                        
                                        double sum3 = 0;
                                        for(int k1=0; k1<4*singnum; k1++)
                                        {
                                            sum3 += wei_doublet_adj.at(k1) * ( 
                                                    singlet2doublet(newt, m, n, k1)
                                                    + singlet2doublet(newt-1, j, m, k1)
                                                    + singlet2doublet(newt-2, i, j, k1)
                                                    + singlet2doublet(newt-3, d, i, k1)
                                                    );
                                        }

                                    double add = alpha_h.at(t-4).at(a).at(b).at(c).at(d) + sum2 + sum2_hh+sum3;
                                    if ( fabs(add) < nearzero ) add = negDmax;
                                    dd = utils.logadd( dd, add );
                                    }
                                }
                            }
                        }
                        alpha_h.at(t).at(i).at(j).at(m).at(n) += dd;
                    }
                }
            }
        }
    }
    outscore.resize(N, negDmax);
    for(int n=0; n<N; n++)
    {
        double dd = negDmax;
        for(int m=0; m<N; m++)
            for(int j=0; j<N; j++)
                for(int i=0; i<N; i++)
                { 
                    double add = alpha_h.at( len-1 ).at(i).at(j).at(m).at(n);
                    if ( fabs(add) < nearzero ) add = negDmax;
                    dd = utils.logadd( dd, add);
                }
        outscore.at(n) = dd;
    }
}

/*
 * objective: calculate the backward matrix of a given segment whose SS is H
 * input: inscore, a segment
 * output: outscore, beta_h
 */
void Calcfunc::calc_beta_h(
        vector< vector< vector< vector< vector<double> > > > > & beta_h,
        int start,
        int len,
        const vector<double> & inscore,
        vector<double> & outscore
        )
{
    vector< vector< vector< vector< vector<double> > > > > empty_h;
    vector< vector< vector< vector< vector<double> > > > > tmp_h( len, vector< vector< vector< vector<double> > > > ( N, vector< vector< vector<double> > > ( N, vector< vector<double> > (N, vector<double> (N, negDmax) ) )  ) );
    beta_h = tmp_h;
    swap(tmp_h, empty_h);

    int t = len-1;
    int newt = beg+start-(len-1-t);
    
    vector<double> newscore2(N, 0);
    if(start!=T-1)
    {
        for(int i=0; i<N; i++)
        {
            double dd = negDmax;
            for(int j=0; j<N; j++)
            {
                double sum4 = 0;
                for(int k1=0; k1<4*singnum; k1++)
                    sum4 += wei_doublet_adj.at(k1) * singlet2doublet(newt+1, i, j, k1);
                double add = inscore.at(j) + sum4;
                if ( fabs(add) < nearzero ) add = negDmax;
                dd = utils.logadd( dd, add );
            } 
            newscore2.at(i) += dd;
        }
    }
    
    if (len>=1)
    {
        int t = len-1;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    );
                        }
                        beta_h.at(t).at(i).at(j).at(m).at(n) = sum1 + newscore2.at(n);
                    }
                }
            }
        }
    }
    
    if (len>=2)
    {
        int t = len-2;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt+1).at(m).at(k1) 
                                    );
                        }

                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * singlet2doublet(newt+1, n, m, k1);
                        }

                        beta_h.at(t).at(i).at(j).at(m).at(n) = sum1 + sum3 + newscore2.at(m);
                    }
                }
            }
        }
    }
    
    if (len>=3)
    {
        int t = len-3;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt+1).at(m).at(k1) 
                                    + singletval.at(newt+2).at(j).at(k1) 
                                    );
                        }
                        
                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * ( 
                                    singlet2doublet(newt+1, n, m, k1)
                                    + singlet2doublet(newt+2, m, j, k1)
                                    );
                        }

                        beta_h.at(t).at(i).at(j).at(m).at(n) = sum1 + sum3 + newscore2.at(j);
                    }
                }
            }
        }
    }

    if (len>=4)
    {
        int t = len-4;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0, sum2_hh = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt+1).at(m).at(k1) 
                                    + singletval.at(newt+2).at(j).at(k1) 
                                    + singletval.at(newt+3).at(i).at(k1) 
                                    );
                        }
                        for(int p=3; p<6; p++)
                            sum2_hh += wei_doublet_ho.at(p) * ( (n==i) ? doubletval_ho.at(newt+3).at(p): 0 );
                        
                        double sum3 = 0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * ( 
                                    singlet2doublet(newt+1, n, m, k1)
                                    + singlet2doublet(newt+2, m, j, k1)
                                    + singlet2doublet(newt+3, j, i, k1)
                                    );
                        }

                        beta_h.at(t).at(i).at(j).at(m).at(n) = sum1 + sum3 + sum2_hh + newscore2.at(i);
                    }
                }
            }
        }
    }
    if (len>=5)
    {
        int t = len-5;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double dd = negDmax;
                        for(int c=0; c<N; c++)
                        {
                            double sum1=0, sum2=0, sum2_hh=0;
                            for(int k1=0; k1<2*singnum; k1++)
                            {
                                sum1 += wei_singlet.at(k1) * ( 
                                        singletval.at(newt).at(n).at(k1) 
                                        + singletval.at(newt+1).at(m).at(k1) 
                                        + singletval.at(newt+2).at(j).at(k1) 
                                        + singletval.at(newt+3).at(i).at(k1) 
                                        + singletval.at(newt+4).at(c).at(k1) 
                                        );
                            }
                            for(int p=6; p<9; p++)
                                sum2 += wei_doublet_ho.at(p) * ( ( n==c ) ? doubletval_ho.at(newt+4).at(p) : 0 );
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) *( 
                                        ( ( n==i ) ? doubletval_ho.at(newt+3).at(p) : 0 )
                                        + ( ( m==c ) ? doubletval_ho.at(newt+4).at(p) : 0 )
                                        );

                            double sum3 = 0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * ( 
                                        singlet2doublet(newt+1, n, m, k1)
                                        + singlet2doublet(newt+2, m, j, k1)
                                        + singlet2doublet(newt+3, j, i, k1)
                                        + singlet2doublet(newt+4, i, c, k1)
                                        );
                            }

                            double add = sum1+sum2+sum2_hh+sum3 + newscore2.at(c);
                            if ( fabs(add) < nearzero ) add = negDmax;
                            dd = utils.logadd( dd, add );
                        }
                        beta_h.at(t).at(i).at(j).at(m).at(n) = dd;
                    }
                }
            }
        }
    }

    if (len>=6)
    {
        int t = len-6;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double dd = negDmax;
                        for(int b=0; b<N; b++)
                        {
                            for(int c=0; c<N; c++)
                            {
                                double sum1=0, sum2=0, sum2_hh=0;
                                for(int k1=0; k1<2*singnum; k1++)
                                {
                                    sum1 += wei_singlet.at(k1) * ( 
                                            singletval.at(newt).at(n).at(k1) 
                                            + singletval.at(newt+1).at(m).at(k1) 
                                            + singletval.at(newt+2).at(j).at(k1) 
                                            + singletval.at(newt+3).at(i).at(k1) 
                                            + singletval.at(newt+4).at(c).at(k1) 
                                            + singletval.at(newt+5).at(b).at(k1) 
                                            );
                                }
                                for(int p=6; p<9; p++)
                                    sum2 += wei_doublet_ho.at(p) * ( 
                                        ( ( n==c ) ? doubletval_ho.at(newt+4).at(p) : 0 )
                                        + ( ( m==b ) ? doubletval_ho.at(newt+5).at(p) : 0 )
                                        );
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) *( 
                                       + ( ( n==i ) ? doubletval_ho.at(newt+3).at(p) : 0 )
                                       + ( ( m==c ) ? doubletval_ho.at(newt+4).at(p) : 0 )
                                       + ( ( j==b ) ? doubletval_ho.at(newt+5).at(p) : 0 )
                                        );

                                double sum3 = 0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * ( 
                                            singlet2doublet(newt+1, n, m, k1)
                                            + singlet2doublet(newt+2, m, j, k1)
                                            + singlet2doublet(newt+3, j, i, k1)
                                            + singlet2doublet(newt+4, i, c, k1)
                                            + singlet2doublet(newt+5, c, b, k1)
                                            );
                                }

                                double add = sum1+sum2+sum2_hh+sum3 + newscore2.at(b);
                                if ( fabs(add) < nearzero ) add = negDmax;
                                dd = utils.logadd( dd, add );
                            }
                        }
                        beta_h.at(t).at(i).at(j).at(m).at(n) = dd;
                    }
                }       
            }
        }
    }

    if (len>=7)
    {
        int t = len-7;
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double dd = negDmax;
                        for(int a=0; a<N; a++)
                        {
                            for(int b=0; b<N; b++)
                            {
                                for(int c=0; c<N; c++)
                                {
                                    double sum1=0, sum2=0, sum2_hh=0;
                                    for(int k1=0; k1<2*singnum; k1++)
                                    {
                                        sum1 += wei_singlet.at(k1) * ( 
                                                singletval.at(newt).at(n).at(k1) 
                                                + singletval.at(newt+1).at(m).at(k1) 
                                                + singletval.at(newt+2).at(j).at(k1) 
                                                + singletval.at(newt+3).at(i).at(k1) 
                                                + singletval.at(newt+4).at(c).at(k1) 
                                                + singletval.at(newt+5).at(b).at(k1) 
                                                + singletval.at(newt+6).at(a).at(k1) 
                                                );
                                    }
                                    for(int p=6; p<9; p++)
                                        sum2 += wei_doublet_ho.at(p) * ( 
                                            + ( ( n==c ) ? doubletval_ho.at(newt+4).at(p) : 0 )
                                            + ( ( m==b ) ? doubletval_ho.at(newt+5).at(p) : 0 )
                                            + ( ( j==a ) ? doubletval_ho.at(newt+6).at(p) : 0 )
                                            );
                                    for(int p=3; p<6; p++)
                                        sum2_hh += wei_doublet_ho.at(p) *( 
                                            + ( ( n==i ) ? doubletval_ho.at(newt+3).at(p) : 0 )
                                            + ( ( m==c ) ? doubletval_ho.at(newt+4).at(p) : 0 )
                                            + ( ( j==b ) ? doubletval_ho.at(newt+5).at(p) : 0 )
                                            + ( ( i==a ) ? doubletval_ho.at(newt+6).at(p) : 0 )
                                            );

                                    double sum3 = 0;
                                    for(int k1=0; k1<4*singnum; k1++)
                                    {
                                        sum3 += wei_doublet_adj.at(k1) * ( 
                                                singlet2doublet(newt+1, n, m, k1)
                                                + singlet2doublet(newt+2, m, j, k1)
                                                + singlet2doublet(newt+3, j, i, k1)
                                                + singlet2doublet(newt+4, i, c, k1)
                                                + singlet2doublet(newt+5, c, b, k1)
                                                + singlet2doublet(newt+6, b, a, k1)
                                                );
                                    }

                                    double add = sum1+sum2+sum2_hh+sum3 + newscore2.at(a);
                                    if ( fabs(add) < nearzero ) add = negDmax;
                                    dd = utils.logadd( dd, add );
                                }
                            }
                            beta_h.at(t).at(i).at(j).at(m).at(n) = dd;
                        }
                    }
                }       
            }
        }
    }
    for(int t=len-8; t>=0; t--)
    {
        int newt = beg+start-(len-1-t);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                for(int m=0; m<N; m++)
                {
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0; 
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * ( 
                                    singletval.at(newt).at(n).at(k1) 
                                    + singletval.at(newt+1).at(m).at(k1) 
                                    + singletval.at(newt+2).at(j).at(k1) 
                                    + singletval.at(newt+3).at(i).at(k1) 
                                    );
                        }
                        beta_h.at(t).at(i).at(j).at(m).at(n) = sum1;
                        double dd = negDmax;
                        for(int a=0; a<N; a++)
                        {
                            for(int b=0; b<N; b++)
                            {
                                for(int c=0; c<N; c++)
                                {
                                    for(int d=0; d<N; d++)
                                    {
                                        double sum2 = 0, sum2_hh = 0;
                                        for(int p=6; p<9; p++)
                                            sum2 += wei_doublet_ho.at(p) * ( 
                                                + ( ( n==d ) ?  doubletval_ho.at(newt+4).at(p) : 0 )
                                                + ( ( m==c ) ?  doubletval_ho.at(newt+5).at(p) : 0 )
                                                + ( ( j==b ) ?  doubletval_ho.at(newt+6).at(p) : 0 )
                                                + ( ( i==a ) ?  doubletval_ho.at(newt+7).at(p) : 0 )
                                                );
                                        for(int p=3; p<6; p++)
                                            sum2_hh += wei_doublet_ho.at(p) *( 
                                            + ( ( n==i ) ? doubletval_ho.at(newt+3).at(p) : 0 )
                                            + ( ( m==d ) ? doubletval_ho.at(newt+4).at(p) : 0 )
                                            + ( ( j==c ) ? doubletval_ho.at(newt+5).at(p) : 0 )
                                            + ( ( i==b ) ? doubletval_ho.at(newt+6).at(p) : 0 )
                                            );

                                        double sum3 = 0;
                                        for(int k1=0; k1<4*singnum; k1++)
                                        {
                                            sum3 += wei_doublet_adj.at(k1) * ( 
                                                    singlet2doublet(newt+1, n, m, k1)
                                                    + singlet2doublet(newt+2, m, j, k1)
                                                    + singlet2doublet(newt+3, j, i, k1)
                                                    + singlet2doublet(newt+4, i, d, k1)
                                                    );
                                        }

                                    double add = beta_h.at(t+4).at(a).at(b).at(c).at(d)+sum2+sum2_hh+sum3;
                                    if ( fabs(add) < nearzero ) add = negDmax;
                                    dd = utils.logadd( dd, add );
                                    }
                                }
                            }
                        }
                        beta_h.at(t).at(i).at(j).at(m).at(n) += dd;
                    }
                }
            }
        }
    }
    outscore.resize(N, negDmax);
    for(int n=0; n<N; n++)
    {
        double dd = negDmax;
        for(int m=0; m<N; m++)
            for(int j=0; j<N; j++)
                for(int i=0; i<N; i++)
                {
                    double add = beta_h.at(0).at(i).at(j).at(m).at(n);
                    if ( fabs(add) < nearzero ) add = negDmax;
                    dd = utils.logadd( dd, add);
                }
        outscore.at(n) = dd;
    }
}

/*
 * calculate the margins(table1 and table5) of a segment whose SS is 'C'
 */
void Calcfunc::margin_c(const int & cur, const int & start_f, const int & len )
{
    //cout<<"\ntable1 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    
    int t;
    for(int t=0; t<len; t++)
    {
        //cout<<"t:"<<t<<" ";
        for(int n=0; n<N; n++)
        {
            double dd = negDmax;
            double sum1=0;
            for(int k1=0; k1<2*singnum; k1++)
            {
                sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
            }
            dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(0).at(t).at(n) + beta_all.at(cur).at(0).at(0).at(0).at(t).at(n) - sum1 );
            table1.at(t+start_f).at(n) = exp( dd - partition );
            //cout<<table1.at(t+start_f).at(n)<<" ";
            //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
        }
        //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
        //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
    }
    
    //cout<<"\ntable5 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    t=0;
    if (cur>0)
    {
        //cout<<"t:"<<t<<" ";
        for(int n=0; n<N; n++)
        {
            for(int m=0; m<N; m++)
            {
                double dd = negDmax;
                double sum3=0;
                for(int k1=0; k1<4*singnum; k1++)
                { 
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, m, n, k1);
                }
                dd = utils.logadd( dd, fnewscore.at(cur-1).at(m) + beta_all.at(cur).at(0).at(0).at(0).at(t).at(n) + sum3 );
                table5.at(t+start_f).at(m).at(n) = exp( dd - partition );
                //cout<<table5.at(t+start_f).at(m).at(n)<<" ";
                //sigma5.at(t+start_f) += table5.at(t+start_f).at(m).at(n);
            }
        }
        //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
        //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
    }

    for(int t=1; t<len; t++)
    {
        //cout<<"t:"<<t<<" ";
        for(int n=0; n<N; n++)
        {
            for(int m=0; m<N; m++)
            {
                double dd = negDmax;
                double sum3=0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, m, n, k1);
                }
                dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(0).at(t-1).at(m) + beta_all.at(cur).at(0).at(0).at(0).at(t).at(n) + sum3 );
                table5.at(t+start_f).at(m).at(n) = exp( dd - partition );
                //cout<<table5.at(t+start_f).at(m).at(n)<<" ";
                //sigma5.at(t+start_f) += table5.at(t+start_f).at(m).at(n);
            }
        }
        //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
        //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
        continue;
    }
    
}

/*
 * calculate the margins(table1, table5, table3) of a segment whose SS is 'E'
 */
void Calcfunc::margin_e(const int & cur, const int & start_f, const int & len )
{
    //cout<<"\ntable1 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;

    for(int t=0; t<len; t++)
    {
        //cout<<"t:"<<t<<" ";

        if(t==0)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int m=0; m<N; m++)
                {
                    dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(t).at(m).at(n) + beta_all.at(cur).at(0).at(0).at(t).at(m).at(n) -sum1 );
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==len-1)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int m=0; m<N; m++)
                {
                    dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(t).at(m).at(n) + beta_all.at(cur).at(0).at(0).at(t).at(m).at(n) -sum1 );
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if ( t>0 && t<len-1)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1 = 0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int c=0; c<N; c++)
                {
                    for(int m=0; m<N; m++)
                    {
                        double sum2 = 0;
                        for(int p=0; p<3; p++)
                            sum2 += wei_doublet_ho.at(p) * ( 
                                ( ( m==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                ); 
                        dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(t).at(m).at(n) + beta_all.at(cur).at(0).at(0).at(t).at(c).at(n) + sum2 - sum1);
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
    }

    //cout<<"\ntable3 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    for(int t=0; t<len; t++)
    {
        if (t<2)
        {
            continue;
        }
        //cout<<"t:"<<t<<" ";
        if (t>=2)
        {
            for(int c=0; c<N; c++)
            {
                for(int m=0; m<N; m++)
                {
                    double dd = negDmax;
                    for(int n=0; n<N; n++)
                    {
                        double sum1 = 0;
                        for(int k1=0; k1<2*singnum; k1++)
                        {
                            sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg-1).at(n).at(k1);
                        }
                        double sum2 = 0;
                        for(int p=0; p<3; p++)
                            sum2 += wei_doublet_ho.at(p) * ( ( m==c ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 );
                        dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(t-1).at(m).at(n) + beta_all.at(cur).at(0).at(0).at(t-1).at(c).at(n) + sum2 - sum1);
                    }
                    table3.at(t+start_f).at(m).at(c) =  exp( dd - partition );
                    //cout<<table3.at(t+start_f).at(m).at(c) << " ";
                    //sigma3.at(t+start_f) += table3.at(t+start_f).at(m).at(c);
                }
            }
            //cout<<setprecision(10)<<//sigma3.at(t+start_f)<<endl;
            //assert( fabs(//sigma3.at(t+start_f)-1) < epsilon );
            continue;
        }
    }

    //cout<<"\ntable5 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    int t=0;
    if (cur>0)
    {
        //cout<<"t:"<<t<<" ";
        for(int n=0; n<N; n++)
        {
            for(int m=0; m<N; m++)
            {
                double dd = negDmax;
                double sum3=0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, m, n, k1);
                }
                for(int j=0; j<N; j++)
                {
                    dd = utils.logadd( dd, fnewscore.at(cur-1).at(m) + beta_all.at(cur).at(0).at(0).at(t).at(j).at(n) + sum3 );
                }
                table5.at(t+start_f).at(m).at(n) = exp( dd - partition );
                //cout<<table5.at(t+start_f).at(m).at(n)<<" ";
                //sigma5.at(t+start_f) += table5.at(t+start_f).at(m).at(n);
            }
        }
        //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
        //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
    }

    for(int t=1; t<len; t++)
    {
        //cout<<"t:"<<t<<" ";
        for(int m=0; m<N; m++)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1 = 0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * ( 
                            singletval.at(t+start_f+beg).at(n).at(k1)
                            + singletval.at(t+start_f+beg-1).at(m).at(k1)
                            );
                }
                double sum3=0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, m, n, k1);
                }
                dd = utils.logadd( dd, alpha_all.at(cur).at(0).at(0).at(t).at(m).at(n) + beta_all.at(cur).at(0).at(0).at(t-1).at(n).at(m) - sum3 - sum1);
                table5.at(t+start_f).at(m).at(n) =  exp( dd - partition );
                //cout<<table5.at(t+start_f).at(m).at(n) << " ";
                //sigma5.at(t+start_f) += table5.at(t+start_f).at(m).at(n);
            }
        }
        //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
        //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
        continue;
    }
}

/*
 * calculate the margins(table1, table5, table4, table2) of a segment whose SS is 'H'
 */
void Calcfunc::margin_h(const int & cur, const int & start_f, const int & len )
{
    //cout<<"\ntable1 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    for(int t=0; t<len && len<6; t++)
    {
        //cout<<"t:"<<t<<" ";

        if(t==0)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int i=0; i<N; i++)
                {
                    for(int j=0; j<N; j++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) -sum1 );
                        }
                    }
                }

                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==1 && len==4)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int j=0; j<N; j++)
                    {
                        double sum2_hh = 0;
                        for(int p=3; p<6; p++)
                            sum2_hh += wei_doublet_ho.at(p) * ( ( c==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 );
                        double sum3=0;
                        for(int k1=0; k1<4*singnum; k1++)
                        {
                            sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, c, n, k1);
                        }
                        for(int m=0; m<N; m++)
                        {
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t-1).at(0).at(0).at(0).at(c) + beta_all.at(cur).at(t).at(0).at(j).at(m).at(n) + sum2_hh + sum3 );
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
        if(t==1 && len==5)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int i=0; i<N; i++)
                    {
                        double sum2 = 0;
                        for(int p=6; p<9; p++)
                            sum2 += wei_doublet_ho.at(p) * ( ( c==i ) ? doubletval_ho.at(t+start_f+beg+3).at(p) : 0 );
                        for(int j=0; j<N; j++)
                        {
                            double sum2_hh = 0;
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) * (
                                    ( ( c==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                    );
                            double sum3=0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, c, n, k1);
                            }
                            for(int m=0; m<N; m++)
                            {
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t-1).at(0).at(0).at(0).at(c) + beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
        if(t==2 && len==4)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int b=0; b<N; b++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            double sum2_hh = 0;
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) * ( ( c==m ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 );
                            double sum3=0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, b, n, k1);
                            }
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t-1).at(0).at(0).at(c).at(b) + beta_all.at(cur).at(t).at(0).at(0).at(m).at(n ) + sum2_hh + sum3);
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
        if(t==2 && len==5)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int j=0; j<N; j++)
                    {
                        double sum2 = 0;
                        for(int p=6; p<9; p++)
                            sum2 += wei_doublet_ho.at(p) * ( ( c==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 );
                        for(int b=0; b<N; b++)
                        {
                            for(int m=0; m<N; m++)
                            {
                                double sum2_hh = 0;
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) * (
                                        ( ( c==m ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                        + ( ( b==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                        );
                                double sum3=0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, b, n, k1);
                                }
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t-1).at(0).at(0).at(c).at(b) + beta_all.at(cur).at(t).at(0).at(j).at(m).at(n) + sum2 + sum2_hh + sum3);
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==3 && len==5)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int a=0; a<N; a++)
                {
                    for(int m=0; m<N; m++)
                    {
                        double sum2 = 0;
                        for(int p=6; p<9; p++)
                            sum2 += wei_doublet_ho.at(p) * ( ( a==m ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 );
                        for(int b=0; b<N; b++)
                        {
                            for(int c=0; c<N; c++)
                            {
                                double sum2_hh = 0;
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) * (
                                        ( ( a==n ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                        + ( ( b==m ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                        );
                                double sum3=0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, c, n, k1);
                                }
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t-1).at(0).at(a).at(b).at(c) + beta_all.at(cur).at(t).at(0).at(0).at(m).at(n) + sum2 + sum2_hh + sum3 );
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==len-1 || len==4)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int i=0; i<N; i++)
                {
                    for(int j=0; j<N; j++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) - sum1);
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
    }

    for(int t=0; t<len && len>=6; t++)
    {
        //cout<<"t:"<<t<<" ";
        if(t==0)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int i=0; i<N; i++)
                {
                    for(int j=0; j<N; j++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) -sum1 );
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==1)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int i=0; i<N; i++)
                    {
                        double sum2 = 0;
                        if (len>=5)
                        {
                            for(int p=6; p<9; p++)
                                sum2 += wei_doublet_ho.at(p) * ( ( c==i ) ? doubletval_ho.at(t+start_f+beg+3).at(p) : 0 );
                        }
                        for(int j=0; j<N; j++)
                        {
                            double sum2_hh = 0;
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) * (
                                    ( ( c==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                    );
                            double sum3=0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, c, n, k1);
                            }
                            for(int m=0; m<N; m++)
                            {
                                dd = utils.logadd( dd, beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) + alpha_all.at(cur).at(t-1).at(0).at(0).at(0).at(c) + sum2 + sum2_hh + sum3 );
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
        if(t==2)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int b=0; b<N; b++)
                    {
                        for(int j=0; j<N; j++)
                        {
                            for(int i=0; i<N; i++)
                            {
                                double sum2 = 0;
                                if (t+2<len)
                                {
                                    for(int p=6; p<9; p++)
                                        sum2 += wei_doublet_ho.at(p) * (
                                            ( ( b==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                            );
                                }
                                if (t+3<len)
                                {
                                    for(int p=6; p<9; p++)
                                        sum2 += wei_doublet_ho.at(p) * (
                                            ( ( c==i ) ? doubletval_ho.at(t+start_f+beg+3).at(p) : 0 )
                                            );
                                }
                                for(int m=0; m<N; m++)
                                {
                                    double sum2_hh = 0;
                                    for(int p=3; p<6; p++)
                                        sum2_hh += wei_doublet_ho.at(p) * (
                                                ( ( b==m ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                                );
                                    if (t+2<len)
                                    {
                                        for(int p=3; p<6; p++)
                                            sum2_hh += wei_doublet_ho.at(p) * (
                                                ( ( c==j ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                                );
                                    }
                                    double sum3=0;
                                    for(int k1=0; k1<4*singnum; k1++)
                                    {
                                        sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, c, n, k1);
                                    }
                                    dd = utils.logadd( dd, beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) + alpha_all.at(cur).at(t-1).at(0).at(0).at(b).at(c) + sum2 + sum2_hh + sum3 );
                                }
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }


        if(t==len-1)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int i=0; i<N; i++)
                {
                    for(int j=0; j<N; j++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t).at(i).at(j).at(m).at(n) -sum1 );
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==len-2)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int i=0; i<N; i++)
                    {
                        double sum2 = 0;
                        if (len>=5)
                        {
                            for(int p=6; p<9; p++)
                                sum2 += wei_doublet_ho.at(p) * ( ( i==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 );
                        }
                        for(int j=0; j<N; j++)
                        {
                            for(int m=0; m<N; m++)
                            {
                                double sum2_hh = 0;
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) * (
                                        ( ( j==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                        );
                                double sum3=0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg+1, n, c, k1);
                                }
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t+1).at(0).at(0).at(0).at(c) + sum2 + sum2_hh + sum3 );
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if(t==len-3)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                for(int c=0; c<N; c++)
                {
                    for(int b=0; b<N; b++)
                    {
                        for(int j=0; j<N; j++)
                        {
                            for(int i=0; i<N; i++)
                            {
                                double sum2 = 0;
                                for(int p=6; p<9; p++)
                                    sum2 += wei_doublet_ho.at(p) * (
                                        + ( ( i==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                        + ( ( j==b ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                        );
                                for(int m=0; m<N; m++)
                                {
                                    double sum2_hh = 0;
                                    for(int p=3; p<6; p++)
                                        sum2_hh += wei_doublet_ho.at(p) * (
                                            + ( ( j==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                            + ( ( m==b ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                            );
                                    double sum3=0;
                                    for(int k1=0; k1<4*singnum; k1++)
                                    {
                                        sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg+1, n, c, k1);
                                    }
                                    dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t+1).at(0).at(0).at(b).at(c) + sum2 + sum2_hh + sum3 );
                                }
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }

        if ( t>2 && t<len-3)
        {
            for(int n=0; n<N; n++)
            {
                double dd = negDmax;
                double sum1 = 0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg).at(n).at(k1);
                }
                for(int b=0; b<N; b++)
                {
                    for(int c=0; c<N; c++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            for(int i=0; i<N; i++)
                            {
                                for(int j=0; j<N; j++)
                                {
                                    for(int a=0; a<N; a++)
                                    {
                                        double sum2 = 0;
                                        for(int p=6; p<9; p++)
                                            sum2 += wei_doublet_ho.at(p) * ( 
                                                + ( ( i==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                                + ( ( j==b ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                                + ( ( m==a ) ? doubletval_ho.at(t+start_f+beg+3).at(p) : 0 )
                                                ); 
                                        double sum2_hh = 0;
                                        for(int p=3; p<6; p++)
                                                sum2_hh += wei_doublet_ho.at(p) * (
                                                ( ( m==b ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                                + ( ( j==c ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                                );
                                        dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t).at(a).at(b).at(c).at(n) + sum2 - sum1 + sum2_hh );
                                    }
                                }
                            }
                        }
                    }
                }
                table1.at(t+start_f).at(n) = exp( dd - partition );
                //cout<<table1.at(t+start_f).at(n)<<" ";
                //sigma1.at(t+start_f) += table1.at(t+start_f).at(n);
            }
            //cout<<setprecision(10)<<//sigma1.at(t+start_f)<<endl;
            //assert( fabs(//sigma1.at(t+start_f)-1) < epsilon );
            continue;
        }
    }
    
    //cout<<"\ntable2 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    for(int t=0; t<len; t++)
    {
        if (t<4)
        {
            continue;
        }
        //cout<<t<<" ";

        if (t==4)
        {
            for(int d=0; d<N; d++)
            {
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int m=0; m<N; m++)
                    {
                        for(int n=0; n<N; n++)
                        {
                            for(int j=0; j<N; j++)
                            {
                                double sum2 = 0;
                                for(int p=6; p<9; p++)
                                    sum2 += wei_doublet_ho.at(p) * ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 );
                                double sum2_hh = 0;
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) * ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p): 0 );
                                double sum3=0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg-3, d, n, k1);
                                }
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t-4).at(0).at(0).at(0).at(d) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                            }
                        }
                    }
                    table2.at(t+start_f).at(d).at(i) =  exp( dd - partition );
                    //cout<<table2.at(t+start_f).at(d).at(i) << " ";
                    //sigma2.at(t+start_f) += table2.at(t+start_f).at(d).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma2.at(t+start_f)<<endl;
            //assert( fabs(//sigma2.at(t+start_f)-1) < epsilon );
            continue;
        }
        if (t==5)
        {
            for(int d=0; d<N; d++)
            {
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int c=0; c<N; c++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            for(int n=0; n<N; n++)
                            {
                                for(int j=0; j<N; j++)
                                {
                                    double sum2 = 0;
                                    for(int p=6; p<9; p++)
                                        sum2 += wei_doublet_ho.at(p) * ( 
                                                ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                                + ( (c==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                );
                                    double sum2_hh = 0;
                                    for(int p=3; p<6; p++)
                                        sum2_hh += wei_doublet_ho.at(p) * (
                                            ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                            + ( ( c==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                            );
                                    double sum3=0;
                                    for(int k1=0; k1<4*singnum; k1++)
                                    {
                                        sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg-3, d, n, k1);
                                    }
                                    dd = utils.logadd( dd, alpha_all.at(cur).at(t-4).at(0).at(0).at(c).at(d) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                                }
                            }
                        }
                    }
                    table2.at(t+start_f).at(d).at(i) =  exp( dd - partition );
                    //cout<<table2.at(t+start_f).at(d).at(i) << " ";
                    //sigma2.at(t+start_f) += table2.at(t+start_f).at(d).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma2.at(t+start_f)<<endl;
            //assert( fabs(//sigma2.at(t+start_f)-1) < epsilon );
            continue;
        }

        if (t==6)
        {
            for(int d=0; d<N; d++)
            {
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int c=0; c<N; c++)
                    {
                        for(int b=0; b<N; b++)
                        {
                            for(int m=0; m<N; m++)
                            {
                                for(int n=0; n<N; n++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        double sum2 = 0;
                                        for(int p=6; p<9; p++)
                                            sum2 += wei_doublet_ho.at(p) * ( 
                                                    + ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                                    + ( ( c==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                    + ( ( b==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                                    );
                                        double sum2_hh = 0;
                                        for(int p=3; p<6; p++)
                                            sum2_hh += wei_doublet_ho.at(p) * (
                                                    + ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                    + ( ( c==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                                    + ( ( b==n ) ? doubletval_ho.at(t+start_f+beg-3).at(p) : 0 )
                                            );
                                        double sum3=0;
                                        for(int k1=0; k1<4*singnum; k1++)
                                        {
                                            sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg-3, d, n, k1);
                                        }
                                        dd = utils.logadd( dd, alpha_all.at(cur).at(t-4).at(0).at(b).at(c).at(d) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                                    }
                                }
                            }
                        }
                    }
                    table2.at(t+start_f).at(d).at(i) =  exp( dd - partition );
                    //cout<<table2.at(t+start_f).at(d).at(i) << " ";
                    //sigma2.at(t+start_f) += table2.at(t+start_f).at(d).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma2.at(t+start_f)<<endl;
            //assert( fabs(//sigma2.at(t+start_f)-1) < epsilon );
            continue;
        }

        if (t>6)
        {
            for(int d=0; d<N; d++)
            {
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int c=0; c<N; c++)
                    {
                        for(int a=0; a<N; a++)
                        {
                            for(int b=0; b<N; b++)
                            {
                                for(int m=0; m<N; m++)
                                {
                                    for(int n=0; n<N; n++)
                                    {
                                        for(int j=0; j<N; j++)
                                        {
                                            double sum2 = 0;
                                            for(int p=6; p<9; p++)
                                                sum2 += wei_doublet_ho.at(p) * ( 
                                                    + ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                                    + ( ( c==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                    + ( ( b==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                                    + ( ( a==n ) ? doubletval_ho.at(t+start_f+beg-3).at(p) : 0 )
                                                    );
                                            double sum2_hh = 0;
                                            for(int p=3; p<6; p++)
                                                sum2_hh += wei_doublet_ho.at(p) * (
                                                    + ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                    + ( ( c==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                                    + ( ( b==n ) ? doubletval_ho.at(t+start_f+beg-3).at(p) : 0 )
                                                    );
                                            double sum3=0;
                                            for(int k1=0; k1<4*singnum; k1++)
                                            {
                                                sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg-3, d, n, k1);
                                            }
                                            dd = utils.logadd( dd, alpha_all.at(cur).at(t-4).at(a).at(b).at(c).at(d) 
                                                    + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    table2.at(t+start_f).at(d).at(i) =  exp( dd - partition );
                    //cout<<table2.at(t+start_f).at(d).at(i) << " ";
                    //sigma2.at(t+start_f) += table2.at(t+start_f).at(d).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma2.at(t+start_f)<<endl;
            //assert( fabs(//sigma2.at(t+start_f)-1) < epsilon );
            continue;
        }
    }

    //cout<<"\ntable4 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl; //可以与table2和table1合并
    for(int t=0; t<len; t++)
    {
        if (t<3)
        {
            continue;
        }
        //cout<<t<<" ";
        
        if (t==3)
        {
            for(int n=0; n<N; n++)
            {
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg-3).at(n).at(k1);
                }
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int m=0; m<N; m++)
                    {
                        for(int j=0; j<N; j++)
                        {
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t-3).at(0).at(0).at(0).at(n) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) - sum1 );
                        }
                    }
                    table4.at(t+start_f).at(n).at(i) =  exp( dd - partition );
                    //cout<<table4.at(t+start_f).at(n).at(i) << " ";
                    //sigma4.at(t+start_f) += table4.at(t+start_f).at(n).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma4.at(t+start_f)<<endl;
            //assert( fabs(//sigma4.at(t+start_f)-1) < epsilon );
            continue;
        }

        if (t==4)
        {
            for(int n=0; n<N; n++)
            {
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int m=0; m<N; m++)
                    {
                        for(int d=0; d<N; d++)
                        {
                            for(int j=0; j<N; j++)
                            {
                                double sum2 = 0;
                                for(int p=6; p<9; p++)
                                    sum2 += wei_doublet_ho.at(p) * ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 );
                                double sum2_hh = 0;
                                for(int p=3; p<6; p++)
                                    sum2_hh += wei_doublet_ho.at(p) * ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 );
                                double sum3=0;
                                for(int k1=0; k1<4*singnum; k1++)
                                {
                                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg-3, d, n, k1);
                                }
                                dd = utils.logadd( dd, alpha_all.at(cur).at(t-4).at(0).at(0).at(0).at(d) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                            }
                        }
                    }
                    table4.at(t+start_f).at(n).at(i) =  exp( dd - partition );                     
                    //cout<<table4.at(t+start_f).at(n).at(i) << " ";
                    //sigma4.at(t+start_f) += table4.at(t+start_f).at(n).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma4.at(t+start_f)<<endl;
            //assert( fabs(//sigma4.at(t+start_f)-1) < epsilon );
            continue;
        }
        if (t==5)
        {
            for(int n=0; n<N; n++)
            {
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int c=0; c<N; c++)
                    {
                        for(int m=0; m<N; m++)
                        {
                            for(int d=0; d<N; d++)
                            {
                                for(int j=0; j<N; j++)
                                {
                                    double sum2 = 0;
                                    for(int p=6; p<9; p++)
                                        sum2 += wei_doublet_ho.at(p) * ( 
                                            ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                            + ( ( c==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                            );
                                    double sum2_hh = 0;
                                    for(int p=3; p<6; p++)
                                        sum2_hh += wei_doublet_ho.at(p) * (
                                            ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                            + ( ( c==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                            );
                                    double sum3=0;
                                    for(int k1=0; k1<4*singnum; k1++)
                                    {
                                        sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg-3, d, n, k1);
                                    }
                                    dd = utils.logadd( dd, alpha_all.at(cur).at(t-4).at(0).at(0).at(c).at(d) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh + sum3 );
                                }
                            }
                        }
                    }
                    table4.at(t+start_f).at(n).at(i) =  exp( dd - partition );
                    //cout<<table4.at(t+start_f).at(n).at(i) << " ";
                    //sigma4.at(t+start_f) += table4.at(t+start_f).at(n).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma4.at(t+start_f)<<endl;
            //assert( fabs(//sigma4.at(t+start_f)-1) < epsilon );
            continue;
        }

        if (t>=6)
        {
            for(int n=0; n<N; n++)
            {
                double sum1=0;
                for(int k1=0; k1<2*singnum; k1++)
                {
                    sum1 += wei_singlet.at(k1) * singletval.at(t+start_f+beg-3).at(n).at(k1);
                }
                for(int i=0; i<N; i++)
                {
                    double dd = negDmax;
                    for(int c=0; c<N; c++)
                    {
                        for(int b=0; b<N; b++)
                        {
                            for(int m=0; m<N; m++)
                            {
                                for(int d=0; d<N; d++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        double sum2 = 0;
                                        for(int p=6; p<9; p++)
                                            sum2 += wei_doublet_ho.at(p) * ( 
                                                + ( ( d==i ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                                + ( ( c==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                + ( ( b==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                                );
                                        double sum2_hh = 0;
                                        for(int p=3; p<6; p++)
                                            sum2_hh += wei_doublet_ho.at(p) * (
                                                + ( ( d==j ) ? doubletval_ho.at(t+start_f+beg-1).at(p) : 0 )
                                                + ( ( c==m ) ? doubletval_ho.at(t+start_f+beg-2).at(p) : 0 )
                                                );
                                        dd = utils.logadd( dd, alpha_all.at(cur).at(t-3).at(b).at(c).at(d).at(n) + beta_all.at(cur).at(t-3).at(i).at(j).at(m).at(n) + sum2 + sum2_hh - sum1 );
                                    }
                                }
                            }
                        }
                    }
                    table4.at(t+start_f).at(n).at(i) =  exp( dd - partition );
                    //cout<<table4.at(t+start_f).at(n).at(i) << " ";
                    //sigma4.at(t+start_f) += table4.at(t+start_f).at(n).at(i);
                }
            }
            //cout<<setprecision(10)<<//sigma4.at(t+start_f)<<endl;
            //assert( fabs(//sigma4.at(t+start_f)-1) < epsilon );
            continue;
        }
    }

    //cout<<"\ntable5 start len cur: "<<start_f<<" "<<len<<" "<<cur<<endl;
    int t=0;
    if (cur>0)
    {
        //cout<<"t:"<<t<<" ";
        for(int n=0; n<N; n++)
        {
            for(int m=0; m<N; m++)
            {
                double dd = negDmax;
                double sum3=0;
                for(int k1=0; k1<4*singnum; k1++)
                {
                    sum3 += wei_doublet_adj.at(k1) * singlet2doublet(t+start_f+beg, m, n, k1);
                }
                for(int j=0; j<N; j++)
                {
                    for(int i=0; i<N; i++)
                    {
                        for(int a=0; a<N; a++)
                            dd = utils.logadd( dd, fnewscore.at(cur-1).at(m) + beta_all.at(cur).at(t).at(i).at(j).at(a).at(n) + sum3 );
                    }
                }
                table5.at(t+start_f).at(m).at(n) = exp( dd - partition );
                //cout<<table5.at(t+start_f).at(m).at(n)<<" ";
                //sigma5.at(t+start_f) += table5.at(t+start_f).at(m).at(n);
            }
        }
        //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
        //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
    }
    for(int t=1; t<len; t++)
    {
        //cout<<"t:"<<t<<" ";
        if (t==1)
        {
            for(int i=0; i<N; i++)
            {
                for(int j=0; j<N; j++)
                {
                    double dd = negDmax;
                    for(int m=0; m<N; m++)
                    {
                        for(int n=0; n<N; n++)
                        {
                            double sum1=0;
                            for(int k1=0; k1<2*singnum; k1++)
                            {
                                sum1 += wei_singlet.at(k1) * (
                                        singletval.at(t+start_f+beg-1).at(i).at(k1)
                                        + singletval.at(t+start_f+beg).at(j).at(k1)
                                        + singletval.at(t+start_f+beg+1).at(m).at(k1)
                                        + singletval.at(t+start_f+beg+2).at(n).at(k1)
                                        );
                            }
                            double sum3=0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * (
                                        singlet2doublet(t+start_f+beg, i, j, k1)
                                        + singlet2doublet(t+start_f+beg+1, j, m, k1)
                                        + singlet2doublet(t+start_f+beg+2, m, n, k1)
                                        );
                            }
                            double sum2_hh = 0;
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) * (
                                        ( ( i==n ) ? doubletval_ho.at(t+start_f+beg+2).at(p) : 0 )
                                        );
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t+2).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t-1).at(n).at(m).at(j).at(i) - sum3 - sum2_hh - sum1 );

                        }
                    }
                    table5.at(t+start_f).at(i).at(j) =  exp( dd - partition );
                    //cout<<table5.at(t+start_f).at(i).at(j) << " ";
                    //sigma5.at(t+start_f) += table5.at(t+start_f).at(i).at(j);
                }
            }
            //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
            //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
            continue;
        }
        if (t==2)
        {
            for(int m=0; m<N; m++)
            {
                for(int j=0; j<N; j++)
                {
                    double dd = negDmax;
                    for(int i=0; i<N; i++)
                    {
                        for(int n=0; n<N; n++)
                        {
                            double sum1=0;
                            for(int k1=0; k1<2*singnum; k1++)
                            {
                                sum1 += wei_singlet.at(k1) * (
                                        singletval.at(t+start_f+beg-2).at(i).at(k1)
                                        + singletval.at(t+start_f+beg-1).at(j).at(k1)
                                        + singletval.at(t+start_f+beg).at(m).at(k1)
                                        + singletval.at(t+start_f+beg+1).at(n).at(k1)
                                        );
                            }
                            double sum3=0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * (
                                        singlet2doublet(t+start_f+beg-1, i, j, k1)
                                        + singlet2doublet(t+start_f+beg, j, m, k1)
                                        + singlet2doublet(t+start_f+beg+1, m, n, k1)
                                        );
                            }
                            double sum2_hh = 0;
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) * (
                                        ( ( i==n ) ? doubletval_ho.at(t+start_f+beg+1).at(p) : 0 )
                                        );
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t+1).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t-2).at(n).at(m).at(j).at(i) - sum3 - sum2_hh - sum1 );

                        }
                    }
                    table5.at(t+start_f).at(j).at(m) =  exp( dd - partition );
                    //cout<<table5.at(t+start_f).at(j).at(m) << " ";
                    //sigma5.at(t+start_f) += table5.at(t+start_f).at(j).at(m);
                }
            }
            //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
            //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
            continue;
        }
        if (t>=3)
        {
            for(int n=0; n<N; n++)
            {
                for(int m=0; m<N; m++)
                {
                    double dd = negDmax;
                    for(int i=0; i<N; i++)
                    {
                        for(int j=0; j<N; j++)
                        {
                            double sum1=0;
                            for(int k1=0; k1<2*singnum; k1++)
                            {
                                sum1 += wei_singlet.at(k1) * (
                                        singletval.at(t+start_f+beg).at(n).at(k1)
                                        + singletval.at(t+start_f+beg-1).at(m).at(k1)
                                        + singletval.at(t+start_f+beg-2).at(j).at(k1)
                                        + singletval.at(t+start_f+beg-3).at(i).at(k1)
                                        );
                            }
                            double sum3=0;
                            for(int k1=0; k1<4*singnum; k1++)
                            {
                                sum3 += wei_doublet_adj.at(k1) * (
                                        singlet2doublet(t+start_f+beg, m, n, k1)
                                        + singlet2doublet(t+start_f+beg-1, j, m, k1)
                                        + singlet2doublet(t+start_f+beg-2, i, j, k1)
                                        );
                            }
                            double sum2_hh = 0;
                            for(int p=3; p<6; p++)
                                sum2_hh += wei_doublet_ho.at(p) * (
                                        ( ( i==n ) ? doubletval_ho.at(t+start_f+beg).at(p) : 0 )
                                        );
                            dd = utils.logadd( dd, alpha_all.at(cur).at(t).at(i).at(j).at(m).at(n) + beta_all.at(cur).at(t-3).at(n).at(m).at(j).at(i) - sum3 - sum2_hh - sum1 );
                        }
                    }
                    table5.at(t+start_f).at(m).at(n) =  exp( dd - partition );
                    //cout<<table5.at(t+start_f).at(m).at(n) << " ";
                    //sigma5.at(t+start_f) += table5.at(t+start_f).at(m).at(n);
                }
            }
            //cout<<setprecision(10)<<//sigma5.at(t+start_f)<<endl;
            //assert( fabs(//sigma5.at(t+start_f)-1) < epsilon );
            continue;
        }
    }
}

//a customized functor
template<class T>
struct negative
{
    void operator()( T & x) const
    {
        x = -x;
    }
};
/*
 * calculate the gradient for given margins, singlet values, doublet values
 * debugging to check whether the real gradient is equal to the calculated one
 */
void Calcfunc::calc_grad1( vector<double> & gradient ) const 
{
    int K1 = 2*singnum;
    int K2 = 4*singnum;;
    int K3 = honum;
    gradient.resize(K1 + K2 + K3, 0);
    for(int k=0; k<K1; k++)
    {
        for(int t=0; t<T; t++)
        {
            for(int i=0; i<N; i++)
            {
                gradient.at(k) += table1.at(t).at(i) * singletval.at(t+beg).at(i).at(k);
            }
        }
    }
    
    for(int k=0; k<K2; k++)
    {  
        for(int t=1; t<T; t++)
        {
            for(int i=0; i<N; i++)
            {
                for(int j=0; j<N; j++)
                {
                    gradient.at(K1+k) += table5.at(t).at(i).at(j) * singlet2doublet(t+beg, i, j, k);
                }
            }
        }
    }
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            for(int t=0; t<T; t++)
            {
                for(int k=0; k<3; k++)
                    gradient.at(K1+K2+k) += table3.at(t).at(i).at(j) * ( ( i==j ) ? doubletval_ho.at(t+beg).at(k) : 0 );
                for(int k=3; k<6; k++)
                    gradient.at(K1+K2+k) += table4.at(t).at(i).at(j) * ( ( i==j ) ? doubletval_ho.at(t+beg).at(k) : 0 );
                for(int k=6; k<9; k++)
                    gradient.at(K1+K2+k) += table2.at(t).at(i).at(j) * ( ( i==j ) ? doubletval_ho.at(t+beg).at(k) : 0 );
            }
        }
    }

    for(int k=0; k<K1; k++)
    {
        double sum = 0;
        for(int t=0; t<T; t++)
        {
            sum += singletval.at(t+beg).at( y_id.at(t) ).at(k);
        }
        gradient.at(k) = sum - gradient.at(k);
    }

    for(int k=0; k<K2; k++)
    {
        double sum_adj = 0;
        for(int t=1; t<T; t++)
        {
            sum_adj += singlet2doublet( t+beg, y_id.at(t-1), y_id.at(t), k);
        }
        gradient.at(K1+k) = sum_adj - gradient.at(K1+k);
    }
    

    for(int k=0; k<3; k++)
    {
        double sum_ho = 0;
        for(int t=2; t<T; t++)
        {
            sum_ho += ( ( y_id.at(t-2) == y_id.at(t) ) ? doubletval_ho.at(t+beg).at(k) : 0 );
        }
        gradient.at(K1+K2+k) = sum_ho - gradient.at(K1+K2+k);
    }
    for(int k=3; k<6; k++)
    {
        double sum_ho = 0;
        for(int t=3; t<T; t++)
        {
            sum_ho += ( ( y_id.at(t-3) == y_id.at(t) ) ? doubletval_ho.at(t+beg).at(k): 0 );
        }
        gradient.at(K1+K2+k) = sum_ho - gradient.at(K1+K2+k);
    }
    for(int k=6; k<9; k++)
    {
        double sum_ho = 0;
        for(int t=4; t<T; t++)
        {
            sum_ho += ( ( y_id.at(t-4) == y_id.at(t) ) ? doubletval_ho.at(t+beg).at(k) : 0 );
        }
        gradient.at(K1+K2+k) = sum_ho - gradient.at(K1+K2+k);
    }
    for_each(gradient.begin(), gradient.end(), negative<double>());
}

/*
 * calculate the value of the objective function with given weights, singlet values, doublet values
 */
double Calcfunc::calc_obj() const
{
    int K1 = 2*singnum;
    int K2 = 4*singnum;;
    int K3 = honum;
    double sum = 0;
    for(int t=0; t<T; t++)
    {
        for(int k=0; k<K1; k++)
            sum += wei_singlet.at(k) * singletval.at(t+beg).at( y_id.at(t) ).at(k);
        if (t>=1)
            for(int k=0; k<K2; k++)
                sum += wei_doublet_adj.at(k) * singlet2doublet( t+beg, y_id.at(t-1), y_id.at(t), k);
        if (t>=2)
            for(int k=0; k<3; k++)
                sum += wei_doublet_ho.at(k) * ( ( y_id.at(t-2) == y_id.at(t) ) ? doubletval_ho.at(t+beg).at(k) : 0 );
        if (t>=3)
            for(int k=3; k<6; k++)
                sum += wei_doublet_ho.at(k) * ( ( y_id.at(t-3) == y_id.at(t) ) ? doubletval_ho.at(t+beg).at(k) : 0 );
        if (t>=4)
            for(int k=6; k<9; k++)
                sum += wei_doublet_ho.at(k) * ( ( y_id.at(t-4) == y_id.at(t) ) ? doubletval_ho.at(t+beg).at(k) : 0 );
    }
    sum -= partition;
    sum = 0 - sum;

    //cout<<"tmp_fx: "<<sum<<endl;

    return sum;
}
/*
double Calcfunc::psudolikelihood( const string & seqname ) const
{
    string path;
    path.resize(T);
    for(int t=0; t<T; t++)
    {
        double max = negDmax;
        for(int i=0; i<N; i++)
        {
            if(table1.at(t).at(i) > max)
            {
                max = table1.at(t).at(i);
                path.at(t)= uy.at(i);
            }
        }
    }
    int right = 0;
    for(int t=0; t<T; t++)
    {
        if ( path.at(t) == label.at(t) )
            right++;
    }
    double srate = (double) right / (double) T;
    return srate;
}
*/

void Calcfunc::psudolikelihood( 
        const string & path_sing_file, 
        const string & seqname, 
        double & srate,
        const string & seqss, 
        const string & seqaa
        ) const
{
    ofstream path_sing(path_sing_file.c_str(), ios::app);
    string path;
    path.resize(T);

    for(int t=0; t<T; t++)
    {
        double max = negDmax;
        for(int i=0; i<N; i++)
        {
            if(table1.at(t).at(i) > max)
            {
                max = table1.at(t).at(i);
                path.at(t) = uy.at(i);
            }
        }
        path_sing<<seqname<<" "<<path.at(t)<<" "<<seqss.at(t)<<" "<<seqaa.at(t)<<" "<<setprecision(2)<<table1.at(t).at(0) <<endl;
    }
    path_sing<<endl;
    path_sing.close();
    int right = 0;
    for(int t=0; t<T; t++)
        if ( path.at(t) == label.at(t) )
            right++;
    srate = (double) right / (double) (T);
}
