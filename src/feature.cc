#include "feature.h"
using namespace std;
static const double epsilon = 1e-5;
static const int dist_h = 5; //A1A5
static const int dist_hh = 4; //A1A4
static const int dist_e = 3; //B1B3

Feat::Feat(string trainfile1)
{
    trainfn = trainfile1;
}

Feat::~Feat()
{
}
/*
 * row: number of lines of the whole data set.
 * singletval: T*N*(2*singnum), singlet features.
 * doubletval_adj: T*N*N*(4*singnum), singlet features.
 * doubletval_ho: honum*T*N*N, singlet features. B1B3( bin + sim + mi ), A1A4( bin + sim + mi ), A1A5( bin + sim + mi )
 * uy: state dictionary
 * objective: read feature values from train data.
 * input: train data smoothed according to SS.
 * output: singletval, doubletval_adj, doubletval_ho
 */
// head ( 7: d1a0aa_ 1 M C 0.875 e 63 )
// pos (1) + edge(1) + aa (20) + ss (3) + pssm (20) + contactnum (1) + disorder (1) + bindsite (1) + endpoint (11) + phychem (7) + ss3 (3) + ss8 (8) + ccpc (20) + cc80 (20) + seq_consv(1) + stru_consv(1) + mi (len) + cos(len) + contact (len)
// 7+1+1+20+3+20+1+1+1+11+7+3+8+20+20+1+1+3*len everyline
// 1+1+20+3+20+1+1+1+11+7+3+8+20+20+2 = 117+2 = 119 singlet
// 1, 2, 22, 25, 45, 46, 47, 48, 59, 66, 69, 77, 97, 117, 118, 119
void Feat::read_feat(
        int row, 
        vector < vector <vector <double> > > & singletval, 
        vector <vector<double> > & doubletval_ho
        )
{
    vector< vector< vector<double> > > empty_sing;
    vector< vector< vector<double> > > tmp1( row, vector< vector<double> >( 2, vector<double> (singnum*2, 0) ) );
    singletval = tmp1;
    swap(empty_sing, tmp1);
    vector< vector<double> > empty_ho;
    vector< vector<double> > tmp3( row, vector<double> ( honum, 0 ) );
    doubletval_ho = tmp3;
    swap(empty_ho, tmp3);

    ifstream inf(trainfn.c_str());
    string line;
    long t=0;
    string lasttag("-1"), curtag("");
    int ct = 0;
    string curss("");
    while(getline(inf, line))
    {
        stringstream sstr(line);
        string curstr(""), fname;
        sstr>>fname;
        sstr>>curtag>>curstr>>curstr;
        if (lasttag=="-1")
        {
            ct = 0;
            curss = curstr;
        }
        else
        {
            ct++;
            curss += curstr;
        }
        const int datacol = 6;
        for(int i=0; i<datacol-4; i++)
            sstr>>curstr;
        int len = 0;
        sstr >> len;
        
        
        
        // hostr: ho + adj + 
        // len (2) + aa (20) + ss (3) + pssm (20) + contactnum (1) + disorder (1) + bindsite (1) + endpoint (11) + phychem (7) + ss3 (3) + ss8 (8) + ccpc (20) + cc80 (20) + mi (len) + cos(len) + contact (len)
        double curval = 0;
        vector<double> singlet( singnum, 0);
        for(int i=0; i<singnum; i++)
        {
            sstr>>curval;
            singlet.at(i) = curval;
        }
        if(hostr.at(1)=='0') 
            for(int ind=0; ind<2; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(2)=='0') 
            for(int ind=2; ind<2+20; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(3)=='0') 
            for(int ind=2+20; ind<2+20+3; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(4)=='0') 
            for(int ind=2+20+3; ind<2+20+3+20; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(5)=='0') 
            for(int ind=2+20+3+20; ind<2+20+3+20+1; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(6)=='0') 
            for(int ind=2+20+3+20+1; ind<2+20+3+20+1+1; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(7)=='0') 
            for(int ind=2+20+3+20+1+1; ind<2+20+3+20+1+1+1; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(8)=='0') 
            for(int ind=2+20+3+20+1+1+1; ind<2+20+3+20+1+1+1+11; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(9)=='0') 
            for(int ind=2+20+3+20+1+1+1+11; ind<2+20+3+20+1+1+1+11+7; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(10)=='0') 
            for(int ind=2+20+3+20+1+1+1+11+7; ind<2+20+3+20+1+1+1+11+7+3; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(11)=='0') 
            for(int ind=2+20+3+20+1+1+1+11+7+3; ind<2+20+3+20+1+1+1+11+7+3+8; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(12)=='0') 
            for(int ind=2+20+3+20+1+1+1+11+7+3+8; ind<2+20+3+20+1+1+1+11+7+3+8+20; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(13)=='0') 
            for(int ind=2+20+3+20+1+1+1+11+7+3+8+20; ind<2+20+3+20+1+1+1+11+7+3+8+20+20; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(14)=='0') 
            for(int ind=2+20+3+20+1+1+1+11+7+3+8+20+20; ind<2+20+3+20+1+1+1+11+7+3+8+20+20+1; ind++)
                singlet.at(ind) = 0;
        if(hostr.at(15)=='0') 
            for(int ind=2+20+3+20+1+1+1+11+7+3+8+20+20+1; ind<2+20+3+20+1+1+1+11+7+3+8+20+20+1+1; ind++)
                singlet.at(ind) = 0;

        copy( singlet.begin(), singlet.end(), singletval.at(t).at(0).begin() );
        copy( singlet.begin(), singlet.end(), singletval.at(t).at(1).begin() + singnum );

        if(hostr.substr(16, 3) != "000")
        {
            vector<double> mi_cos_cont(0);
            while( sstr>>curstr )
            {
                if (curstr=="nan") curstr = "0";
                mi_cos_cont.push_back( atof(curstr.c_str()) );
            }
            int newlen = mi_cos_cont.size() / 3;
            assert( newlen == len );

            vector<double> mi( mi_cos_cont.begin(), mi_cos_cont.begin()+len );
            vector<double> cos( mi_cos_cont.begin()+len, mi_cos_cont.begin()+2*len );
            vector<double> cont( mi_cos_cont.begin()+2*len, mi_cos_cont.begin()+3*len );
            if ( hostr.at(16) == '0' ) mi.resize( len, 0 );
            if ( hostr.at(17) == '0' ) cos.resize( len, 0 );
            if ( hostr.at(18) == '0' ) cont.resize( len, 0 );
            if(ct>=2)
            {
                if( curss.substr(ct-2, 3) == "EEE" )
                {
                    doubletval_ho.at(t).at(0) = mi.at( ct-2 );
                    doubletval_ho.at(t).at(1) = cos.at( ct-2 );
                    doubletval_ho.at(t).at(2) = cont.at( ct-2 );
                }
            }
            if(ct>=3)
            {
                if( curss.substr(ct-3, 4) == "HHHH" )
                {
                    doubletval_ho.at(t).at(3) = mi.at( ct-3 );
                    doubletval_ho.at(t).at(4) = cos.at( ct-3 );
                    doubletval_ho.at(t).at(5) = cont.at( ct-3 );
                }
            }
            if(ct>=4)
            {
                if( curss.substr(ct-4, 5) == "HHHHH" )
                {
                    doubletval_ho.at(t).at(6) = mi.at( ct-4 );
                    doubletval_ho.at(t).at(7) = cos.at( ct-4 );
                    doubletval_ho.at(t).at(8) = cont.at( ct-4 );
                }
            }
        }
        lasttag = curtag;
        t++;
    }
    inf.close();


    /*
       cout<<"singletval: \n";
       for(int t=0; t<row; t++)
       {
       cout<<" t: "<<t+1<<endl;
       for(int i=0; i<2; i++)
       {
       for(int k=0; k<2*singnum; k++)
       {
       cout<<singletval.at(t).at(i).at(k)<<" ";
       }
       cout<<endl;
       }
       cout<<endl;
       }

       cout<<"\ndoubletval_ho: "<<row<<endl;
       for(int t=0; t<row; t++)
       {
    //cout<<"curt:"<<t+1<<endl;
    for(int p=0; p<honum; p++)
    //cout<<"curp: "<<p<<" "<<doubletval_ho.at(t).at(p)<<endl;
    cout<<doubletval_ho.at(t).at(p)<<" ";
    cout<<endl;
    }
    */
}

/*
 * fold && testtag: if testtag=='==', this function will extract all protein infos whose fold numbers are 'fold'; if testtag=='!=', this function will extract all protein infos whose fold numbers are not 'fold';  
 * file: beg and end indices of each protein whose fold number is 'fold' in the original data set
 * label: the whole label sequence of all proteins whose fold number is 'fold' after concatenation
 * ss: the whole ss sequence of all proteins whose fold number is 'fold' after concatenation
 * aa: the whole aa sequence of all proteins whose fold number is 'fold' after concatenation
 * fname: all protein names wholse fold numbers are 'fold'
 * objective: read info of each protein from original train data
 * input: original train data, fold, testtag
 * output: file, label, ss, aa, fname
 */
void Feat::train_info(
        File & file, 
        string & label, 
        string & ss,
        string & aa,
        vector<string> & fname
        )
{
    string str("");
    ifstream infile(trainfn.c_str());
    label.clear();
    ss.clear();
    aa.clear();
    Couple couple;
    int ct = 0; //record the current line number
    couple.beg = 0;
    string lastseq = "";
    while(getline(infile, str))
    {
        //d1a0ca_ 1 F E 0.327103 e 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.327 0.030 0.702 0 0.13 0 0 0 0.01 0.01 0 0.29 0 0 0 0 0.23 0 0 0 0.07 0.25 0.01 1.70552 6 1 0.0640834 0.714826 0 0 0 0 0 0
        stringstream sstr(str);
        string seqname(""), endf("");
        char curaa, curss, be;
        double rsa;
        sstr>>seqname>>endf>>curaa>>curss>>rsa>>be;
        label += be;
        if (endf != "1")
        {
            couple.end = ct;
            file.range.push_back(couple);
            couple.beg = ct+1;
        }
        aa += curaa;
        ss += curss;
        if (seqname != lastseq)
        {
            fname.push_back(seqname);
            lastseq = seqname;
        }
        ct++;
    }
    file.row = ct;
    infile.close();

    int file_num = file.range.size();
    vector<int> beg(file_num);
    vector<int> end(file_num);
    if(!file_num)
        return;

    cout<<"allfname: \n";

    int fsize = fname.size();
    for(int i=0; i<fsize; i++)
        cout<<fname.at(i)<<" ";
    cout<<endl;
    cout<<"all couples: "<<file_num<<endl;
    for(int i=0; i<file_num; i++)
        cout<<file.range.at(i).beg<<" "<<file.range.at(i).end<<endl;

    for(int i=0; i<file_num; i++)
    {
        int beg = file.range.at(i).beg;
        int end = file.range.at(i).end;
        string seqlabel("");   
        string seqss("");
        seqlabel = label.substr(beg, end-beg+1);
        seqss = ss.substr(beg, end-beg+1);
        vector< vector<int> > segment1(0), segment2(0);
        int start = 0, len = 0;
        int curss = ( ( seqss.at(0)=='C' ) ? 0 : ( (seqss.at(0)=='H') ? 1 : 2 ) ) ;
        int curlen = seqss.size();
        for(int pos=1; pos<curlen; pos++)
        {

            len++;
            if( seqss.at(pos) == seqss.at(pos-1) ) continue;
            int arr1[] = {start, len, curss};
            int arr2[] = {start+len-1, len, curss};
            vector<int> triple1(arr1, arr1+3);
            vector<int> triple2(arr2, arr2+3);
            segment1.push_back(triple1);
            segment2.push_back(triple2);
            start = pos;
            len = 0;
            curss = ( ( seqss.at(pos)=='C' ) ? 0 : ( (seqss.at(pos)=='H') ? 1 : 2 ) ) ;
        }
        len++;
        int arr1[] = {start, len, curss};
        int arr2[] = {start+len-1, len, curss};
        vector<int> triple1(arr1, arr1+3);
        vector<int> triple2(arr2, arr2+3);
        segment1.push_back(triple1);
        segment2.push_back(triple2);
        int seg_num = segment2.size();
        for(int j=0; j<seg_num; j++)
            if (j< (seg_num >> 1) )
                swap( segment2.at(seg_num-j-1), segment2.at(j) );
        allsegment1.push_back(segment1);
        allsegment2.push_back(segment2);
    }
}

