#include "train.h"
using namespace std;

static const double negDmax = -DBL_MAX/100.0;
#define epsilon 1e-3

Obj_func::Obj_func(
        )
{
    m_x = NULL;
}

Obj_func::~Obj_func()
{
    if (m_x != NULL) 
    {
        lbfgs_free(m_x);
        m_x = NULL;
    }
}

int Obj_func::run(
        int N, 
        vector<double> & parm
        )
{
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_x = lbfgs_malloc(N);
    if (m_x == NULL) 
    {
        cout<<"ERROR: Failed to allocate a memory block for variables.\n";
        return -1;
    }
    lbfgs_parameter_t param;
    param.delta = 1e-2;
    param.past = 1;
    lbfgs_parameter_init(&param);
    srand((unsigned)time(NULL));
    for(int i=0; i<N; i++)
        m_x[i] = rand()%21+1;
    cout<<"initial weights: \n";
    for(int i=0; i<N; i++)
        cout<<m_x[i]<<" ";
    cout<<endl;
    
    int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, & param);
    parm.resize(N);
    for(int i=0; i<N; i++)
        parm.at(i) = m_x[i];
    return ret;
}

/*
 * In the function lbfgs(...),  _evaluate and _progress should passed with their function pointers. But _evaluate and _progress are
 * member functions, so they should be converted to function pointers using reinterpret_cast and use member functions evaluate or progress 
 * as the auxiliary. So are _calc_grad_obj and calc_grad_obj. 
 */

lbfgsfloatval_t Obj_func::_evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
{
    return reinterpret_cast<Obj_func*>(instance)->evaluate(x, g, n, step);
}

int Obj_func::_progress(
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
        )
{
    return reinterpret_cast<Obj_func*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}


void * Obj_func::_calc_grad_obj( void * allparm)
{
    pid_t pid = getpid();       //获取当前进程id
    pthread_t tid = pthread_self(); //获取当前线程id
    printf("pid: %u tid: %u (0x%x)\n",
            (unsigned int)pid,
            (unsigned int)tid );
    void *instance;
    reinterpret_cast<Obj_func*>(instance)->calc_grad_obj(allparm); 
}

void * Obj_func::calc_grad_obj( void * allparm)
{
    vector<double> wei_singlet = ( (Allparm *)allparm)->wei_singlet;
    vector<double> wei_doublet_ho = ( ( Allparm *)allparm)->wei_doublet_ho;
    vector<double> wei_doublet_adj = ( ( Allparm *) allparm)->wei_doublet_adj;
    int start = ( ( Allparm *) allparm)->start;
    int end = ( ( Allparm *) allparm)->end;
    
    vector< vector< vector< vector< vector<double> > > > > alpha_h, beta_h;
    vector< vector< vector<double> > > alpha_e, beta_e;
    vector< vector<double> > alpha_c, beta_c;
    string seqname, seqlabel, seqss;
    vector< vector<int> > segment1, segment2;

    for(int f_id=start; f_id<=end; f_id++)
    {
        Feat feat(trainfn);
        int fbeg = file.range.at(f_id).beg;
        int fend = file.range.at(f_id).end;
        seqlabel = label.substr(fbeg, fend-fbeg+1);
        seqss = ss.substr(fbeg, fend-fbeg+1);
        seqname = fname.at(f_id);
        segment1 = allsegment1.at(f_id);
        segment2 = allsegment2.at(f_id);
        Calcfunc calcfunc( wei_singlet, wei_doublet_ho, wei_doublet_adj, singletval, doubletval_ho, fbeg, segment1, segment2, hostr, seqlabel, uy, seqss );

        int start_f = 0, len_f = 0;
        int start_b = 0, len_b = 0;
        int curss_f = 0;
        int curss_b = 0;

        vector<double> inscore_f(2, negDmax), outscore_f = inscore_f, inscore_b = inscore_f, outscore_b = inscore_f;
        calcfunc.alpha_all.clear();
        calcfunc.beta_all.clear();
        
        int segments_num = calcfunc.segment1.size();
        calcfunc.fnewscore.resize(segments_num);

        cout<<"\nf_id forward: "<<f_id<<endl;
        for(int i=0; i<segments_num; i++)
        {
            start_f = calcfunc.segment1.at(i).at(0);
            len_f = calcfunc.segment1.at(i).at(1);
            curss_f = calcfunc.segment1.at(i).at(2);
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
            //calcfunc.realpartition(start_f, len_f, inscore_f); //debugging
        }
        cout<<"\nf_id backward: "<<f_id<<endl;
        for(int i=0; i<segments_num; i++)
        {
            start_b = calcfunc.segment2.at(i).at(0);
            len_b = calcfunc.segment2.at(i).at(1);
            curss_b = calcfunc.segment2.at(i).at(2);
            if (start_b != calcfunc.T-1)
                inscore_b = outscore_b;
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
        //cout<<"f_id final_partition: "<<f_id<<" "<<setprecision(3)<<outscore_f_sum<<" "<<outscore_b_sum<<endl;
        //assert( fabs( outscore_f_sum - outscore_b_sum ) < epsilon );
        calcfunc.partition = outscore_f_sum;
        for(int cur=0; cur<segments_num; cur++)
        {
            start_f = calcfunc.segment1.at(cur).at(0);
            len_f = calcfunc.segment1.at(cur).at(1);
            curss_f = calcfunc.segment1.at(cur).at(2);
            if (curss_f == 0)
                calcfunc.margin_c( cur, start_f, len_f );
            if (curss_f == 1)
                calcfunc.margin_h( cur, start_f, len_f );
            if (curss_f == 2)
                calcfunc.margin_e( cur, start_f, len_f );
        }
        calcfunc.calc_yid();
        /*
        for(int i=0; i<seqlabel.size(); i++)
            cout<<calcfunc.y_id.at(i);
        cout<<endl;
        */
        vector<double> grad(0);
        calcfunc.calc_grad1(grad);
        double obj = calcfunc.calc_obj();
        
        ( ( Allparm *) allparm )->obj += obj;
        for(int i=0; i<6*singnum+honum; i++)
            ( ( Allparm *) allparm )->grad.at(i) += grad.at(i);
    }
}

lbfgsfloatval_t Obj_func::evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
{

    vector<double> wei_singlet, wei_doublet_ho, wei_doublet_adj;
    wei_singlet.resize(2*singnum);
    wei_doublet_adj.resize(4*singnum);
    wei_doublet_ho.resize(honum);
    copy(x, x+2*singnum, wei_singlet.begin());
    copy(x+2*singnum, x+6*singnum, wei_doublet_adj.begin());
    copy(x+6*singnum, x+6*singnum+honum, wei_doublet_ho.begin());

    int K = file.range.size();

    int num = thread_num; //线程的个数num
    int fnum = (K-1)/num + 1; //每个线程至多处理的文件个数
    vector< pair<int, int> > fileid( num );
    int tid = 0;
    for( int i=0; i<K; i +=fnum )
    {
        pair<int, int> newpair = make_pair(i, min(i+fnum-1, K-1) );
        fileid.at(tid++) = newpair;
    }

    if(num > K) 
    {
        num = K;
        fileid.resize(num);
        tid = 0;
        for( int i=0; i<num; i++)
        {
            pair<int, int> newpair = make_pair(i, i);
            fileid.at(tid++) = newpair;
        }
    }
    
    /*
    cout<<"fileid: \n";
    for(int i=0; i<fileid.size(); i++)
        cout<<fileid.at(i).first<<" "<<fileid.at(i).second<<endl;
    */

    vector<pthread_t> pth(num); //线程的个数num
    vector<Allparm *> allparm(num);

    for(int tid=0; tid<num; tid++)
    {
        vector<double> grad(n, 0);
        double obj = 0;
        vector<double> ratevec(0);
        allparm.at(tid) = new Allparm(wei_singlet, wei_doublet_ho, wei_doublet_adj, fileid.at(tid).first, fileid.at(tid).second, grad, obj, ratevec);
        
        int err = pthread_create( &(pth.at(tid)), NULL, _calc_grad_obj, static_cast<void*> (allparm.at(tid)) );

        if(err != 0)
        {
            printf("create thread error: %s\n", strerror(err));
            return 0;
        }
    }

    vector<double> gradval(n, 0);
    double objval = 0; 
    
    for(int tid=0; tid<num; tid++)
    {
        //if ( pth.at(tid) == '0' )
        pthread_join( pth.at(tid), NULL);
    }
    for(int tid=0; tid<num; tid++)
    {
        for(int k=0; k<n; k++)
            gradval.at(k) += ( ( Allparm *) (allparm.at(tid)) )->grad.at(k);
        objval += ( ( Allparm *) (allparm.at(tid)) )->obj;
    }
    for(int i=0; i<n; i++)
        g[i] = gradval.at(i) / (double)K;
    lbfgsfloatval_t fx = (lbfgsfloatval_t) objval / (double)K;
    return fx;
}

int Obj_func::progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
{
    cout<<"\nfx: "<<fx<<endl;
    printf("Iteration %d: ", k);
    for(int i=0; i<n; i++)
        cout<<x[i]<<" ";   
    printf("\n gradient:\n");
    for(int i=0; i<n; i++)
        cout<<g[i]<<" ";
    printf("\n");
    return 0;
}
