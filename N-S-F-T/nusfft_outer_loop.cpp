
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <complex.h>
#include <map>
#include <vector>
#include "fftw3.h"
#include "nusfft_outer_loop.h"
//#include "mex.h" //mx函数，mex函数用到的头文件
#include "mex.hpp"
#include "mexAdapter.hpp"
#pragma comment(lib, "D:/fftw-3.3.5-dll64/libfftw3-3.lib");
#pragma comment(lib, "D:/fftw-3.3.5-dll64/libfftw3f-3.lib");
#pragma comment(lib, "D:/fftw-3.3.5-dll64/libfftw3l-3.lib");

//#include <ctime>
//#include <math.h>
//#include "stdio.h"
//#include <string>
//#include "utils.h"
//#include "timer.h"
//#include "computefourier.h"
//#include "filters.h"
//#include "fft.h"
//#include "fftw.h"

struct Filter {
    complex_t* time;
    int sizet;
    complex_t* freq;
};

std::map<int, fftw_plan> fftw_plans;
using namespace matlab::data;
using matlab::mex::ArgumentList;
constexpr auto PI = 3.14159265358979323846;
//using namespace matlab::engine;
//timespec start_time;
//bool ALGORITHM1 = false;
int num_candidates;
double* pst;
complex_t* out_I, *out;

//void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
class MexFunction : public matlab::mex::Function {
   // ArrayFactory factory;
   //std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
  // std::ostringstream stream;
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
       // checkArguments(outputs, inputs);       
        //******************************************************传入参数************************************* 
        struct Filter filter;
        struct Filter filter_Est;
        const int n = (int)inputs[1][0];
        filter.sizet = (int)inputs[3][0];
        filter_Est.sizet = (int)inputs[6][0];
        const int B_est = (int)inputs[8][0];
        const int B_thresh = (int)inputs[9][0];
        const int B_loc = (int)inputs[10][0];
        const int loops_thresh = (int)inputs[11][0];
        const int loops_loc = (int)inputs[12][0];
        const int loops = (int)inputs[13][0];
       //stream << "outputs[0]" << std::endl;
       // displayOnMATLAB(stream);       
        complex_t* origx = (complex_t*)calloc( n,sizeof(complex_t));
       // complex_t* f = (complex_t*)calloc(sizeof(complex_t) * n);
        filter.time = (complex_t*)calloc(n, sizeof(complex_t));
        filter.freq = (complex_t*)calloc(n, sizeof(complex_t));
        filter_Est.time = (complex_t*)calloc(n, sizeof(complex_t));
        filter_Est.freq = (complex_t*)calloc(n, sizeof(complex_t));
       // memset(x, 0, n * sizeof(complex_t));
           
        complex_t** pe_major = new complex_t * [n];//定义二维数组指针
        complex_t** inv_pe_major = new complex_t * [n];
        for (int i = 0; i < n; i++)
        {
            pe_major[i] = new complex_t[n];
            inv_pe_major[i] = new complex_t[n];
        }
        complex_t** sa_major = new complex_t * [B_loc];//定义二维数组指针
        for (int i = 0; i < B_loc; i++)
        {
            sa_major[i] = new complex_t[n];
        }
        //传入数组值
        for (int i = 0; i < n; i++) {
            origx[i] = (complex_t)inputs[0][i];
            filter.time[i] = (complex_t)inputs[2][i];
            filter.freq[i] = (complex_t)inputs[4][i];
            //filter_Est.time[i] = (complex_t)inputs[5][i];free(filter_Est.time);
            filter_Est.freq[i] = (complex_t)inputs[7][i];
            //f[i] = (complex_t)inputs[17][i];
            for (int j = 0;  j< n; j++) {
                pe_major[i][j] = (complex_t)inputs[14][i][j];
                inv_pe_major[i][j] = (complex_t)inputs[15][i][j];
            }
        }  
        for (int i = 0; i < B_loc; i++)
        {
            for (int j = 0; j < n; j++) {
                sa_major[i][j] = (complex_t)inputs[16][i][j];
            }           
        }       
        /*
        TypedArray<complex_t*>  x = std::move(inputs[0]); 
        TypedArray<complex_t*>  filtertime = std::move(inputs[2]);
        TypedArray<complex_t*>  filterfreq = std::move(inputs[4]);
        TypedArray<complex_t*>  filter_Esttime = std::move(inputs[5]);
        TypedArray<complex_t*>  filter_Estfreq = std::move(inputs[7]);
        TypedArray<complex_t*>  pe_major = std::move(inputs[14]);      
        TypedArray<complex_t*>  inv_pe_major = std::move(inputs[15]);
        TypedArray<complex_t*>  sa_major = std::move(inputs[16]);
        TypedArray<complex_t*>  f = std::move(inputs[17]);
        */
        //****************************************************执行操作*********************************************
/*#ifdef DEBUG
        complex_t* tmp = (complex_t*)malloc(n* sizeof(*tmp));
        memcpy(tmp, filter.time, w_loc * sizeof(*tmp));
        free(filter.time); filter.time = tmp;
        //fftw_dft(filterf, n, filtert);

        tmp = (complex_t*)malloc(n* sizeof(*tmp));
        memcpy(tmp, filter_Est.time, w_est * sizeof(*tmp));
        free(filter_Est.time); filter_Est.time = tmp;
        //fftw_dft(filterf_est, n, filtert_est);

        //plot("filter time series", map_abs(Vec(filtert, n)));
        plot("filter fourier series", map_abs(Vec(filter.freq, 4 * n / B_loc)));

        //plot("filterest time series", map_abs(Vec(filtert_est, n)));
        plot("filterest fourier series", map_abs(Vec(filter_Est.freq, 4 * n / B_loc)));

        complex_t* out = (complex_t*)malloc(n * sizeof(*out));
        fftw_dft(out, n, x);
        for (int i = 0; i < n; i++)
            out[i] /= n;

        //plot("Original time series", map_real(Vec(x, n)));
        plot("Original fourier series", map_abs(Vec(out, n)));
        free(out);
#endif*/
       // std::map<int, complex_t> ans;
        outer_loop(origx, n, filter, filter_Est, B_est, B_thresh, B_loc, loops_thresh, loops_loc, loops, pe_major, inv_pe_major, sa_major);
        free(origx);
        free(filter.time);
        free(filter.freq);
        free(filter_Est.time);
        free(filter_Est.freq);
        for (int i = 0; i < n; i++)
        {
            delete [] pe_major[i];
            delete [] inv_pe_major[i];
        }
        
        for (int i = 0; i < B_loc; i++)
        {
            delete [] sa_major[i] ;
        }

        delete [] pe_major;
        delete [] inv_pe_major;
        delete [] sa_major;
        //****************************************************传出矩阵*********************************************      
       
        for (int io = 0; io < num_candidates; io++) {
            outputs[0][io] = pst[io];
            outputs[1][io] = out_I[io];
        }
        for (int iou = 0; iou < n; iou++) {          
            outputs[2][iou] = out[iou];
        }
        
    /*
        int n, B_est, B_thresh, B_loc, loops_thresh, loops_loc, loops;
        complex_t* x, * f;
        complex_t* pe_major, * inv_pe_major, *sa_major;
        //******************************************************传入参数************************************* 
        n = mxGetScalar(prhs[1]);
        struct Filter filter, filter_Est;
        filter.sizet = mxGetScalar(prhs[3]);
        filter_Est.sizet = mxGetScalar(prhs[6]);
        B_est = mxGetScalar(prhs[8]);
        B_thresh = mxGetScalar(prhs[9]);
        B_loc = mxGetScalar(prhs[10]);
        loops_thresh = mxGetScalar(prhs[11]);
        loops_loc = mxGetScalar(prhs[12]);
        loops = mxGetScalar(prhs[13]);
        //mxComplexSingle*
        x = (complex_t*)mxGetComplexSingles(prhs[0]);
        filter.time = (complex_t*)mxGetComplexSingles(prhs[2]);
        filter.freq = (complex_t*)mxGetComplexSingles(prhs[4]);
        filter_Est.time = (complex_t*)mxGetComplexSingles(prhs[5]);
        filter_Est.freq = (complex_t*)mxGetComplexSingles(prhs[7]);
        pe_major = (complex_t*)mxGetComplexSingles(prhs[14]);
        inv_pe_major = (complex_t*)mxGetComplexSingles(prhs[15]);
        sa_major = (complex_t*)mxGetComplexSingles(prhs[16]);
        f = (complex_t*)mxGetComplexSingles(prhs[17]); 
       
         //****************************************************执行操作*********************************************
        #ifdef DEBUG
        complex_t* tmp = (complex_t*)calloc(n, sizeof(*tmp));
        memcpy(tmp, filter.time, w_loc * sizeof(*tmp));
        free(filter.time); filter.time = tmp;
        //fftw_dft(filterf, n, filtert);

        tmp = (complex_t*)calloc(n, sizeof(*tmp));
        memcpy(tmp, filter_Est.time, w_est * sizeof(*tmp));
        free(filter_Est.time); filter_Est.time = tmp;
        //fftw_dft(filterf_est, n, filtert_est);

        //plot("filter time series", map_abs(Vec(filtert, n)));
        plot("filter fourier series", map_abs(Vec(filter.freq, 4 * n / B_loc)));

        //plot("filterest time series", map_abs(Vec(filtert_est, n)));
        plot("filterest fourier series", map_abs(Vec(filter_Est.freq, 4 * n / B_loc)));

        complex_t* out = (complex_t*)calloc(n * sizeof(*out));
        fftw_dft(out, n, x);
        for (int i = 0; i < n; i++)
            out[i] /= n;

        //plot("Original time series", map_real(Vec(x, n)));
        plot("Original fourier series", map_abs(Vec(out, n)));
        free(out);
        #endif

        //std::map<int, complex_t> ans;
        //ans = outer_loop(x, n, filter, filter_Est, B_est, B_thresh, B_loc, loops_thresh, loops_loc, loops, pe_major, inv_pe_major, sa_major, f);
        outer_loop(x, n, filter, filter_Est, B_est, B_thresh, B_loc, loops_thresh, loops_loc, loops, pe_major, inv_pe_major, sa_major, f);
        //****************************************************传出矩阵*********************************************
        
        plhs[0] = mxCreateDoubleMatrix(1, num_candidates, mxREAL);//创建传出的矩阵
        double* p0 = mxGetPr(plhs[0]);//获取矩阵指针
        p0 = pst;
        plhs[1] = mxCreateDoubleMatrix(1, num_candidates, mxCOMPLEX);//创建传出的矩阵      
        complex<double>* p1 = (complex<double>*)mxGetComplexSingles(plhs[1]);//获取矩阵指针
        p1 = (complex<double>*) out_I;
        plhs[2] = mxCreateDoubleMatrix(1, n, mxCOMPLEX);//创建传出的矩阵
        complex<double>* p2 = (complex<double>*)mxGetComplexSingles(plhs[2]);//获取矩阵指针
        p2 = (complex<double>*) out;   
        */    
    }


    inline const int timesmod(const int& x, const int& a, const int& n) {
    return int((((long long int)x) * a) % n);
}

    int fftw_dft(complex_t* out, int n, complex_t* x, int backwards) {
    fftw_plan p;
    if (OPTIMIZE_FFTW) {  //measure the best plan the first time
        if (fftw_plans.find(n) == fftw_plans.end()) { // no plan made yet
            fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(*in) * n);
            fftw_complex* out2 = (fftw_complex*)fftw_malloc(sizeof(*out2) * n);
            p = fftw_plan_dft_1d(n, in, out2,
                backwards ? FFTW_BACKWARD : FFTW_FORWARD,
                FFTW_MEASURE);
            fftw_plans.insert(std::make_pair(n, p));
            fftw_free(in);
            fftw_free(out2);
        }
    }
    fftw_complex* xf = (fftw_complex*)x;
    fftw_complex* outf = (fftw_complex*)x;
    p = fftw_plan_dft_1d(n, xf, outf,
        backwards ? FFTW_BACKWARD : FFTW_FORWARD,
        FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return 0;
}
    
    double cabs2(complex_t x) {
        return (creal(x) * creal(x) + cimag(x) * cimag(x));
    }
   
    int gcd(int a, int b) {
        if (a % b == 0) return b;
        return gcd(b, a % b);
    }

    int mod_inverse(int a, int n) {
        int i = n, v = 0, d = 1;
        while (a > 0) {
            int t = i / a, x = a;
            a = i % x;
            i = x;
            x = d;
            d = v - t * x;
            v = x;
        }
        v %= n;
        if (v < 0) v = (v + n) % n;
        return v;
    }

    real_t nth_element_immutable(real_t* input, int n, int num) {
        real_t* x = (real_t*)malloc(n * sizeof(*x));
      
        if(x){
             memcpy(x, input, n * sizeof(*x));                
             std::nth_element(x, x + num, x + n);
             real_t ans; ans = x[num];
             free(x);
             return ans;
        }       
    }

    void find_largest_indices(int* output, int num, real_t* samples, int n) {
        assert(n >= num + 1);
        //use num+1 so we can use > cutoff and probably get exactly num.
        //if we get fewer, the second pass uses == cutoff.
        real_t cutoff = nth_element_immutable(samples, n, n - num - 1);

        int count = 0;
        for (int i = 0; i < n; i++)
            if (samples[i] > cutoff)
                output[count++] = i;
        if (count < num) {
            for (int i = 0; i < n; i++) {
                if (samples[i] == cutoff) {
                    output[count++] = i;
                    if (count >= num)
                        break;
                }
            }
            std::sort(output, output + count);
        }
        assert(count == num);
    }

    int inner_loop_locate(complex_t* origx, int n, const Filter& filter, int num, int B,
        int a, int ai, int b, complex_t* x_samp, int* J, double& PF_T, double& BC_T,
        complex_t** pe_major, complex_t** inv_pe_major, complex_t** sa_major) {

        if (n % B)
            fprintf(stderr, "Warning: n is not divisible by B, which algorithm expects.\n");
 
        complex_t* x1 = (complex_t*)calloc(n,sizeof(x1) );
        complex_t* x2 = (complex_t*)calloc(n,sizeof(x2) );
        complex_t* x3 = (complex_t*)calloc( n,sizeof(x3) );
        complex_t* x4 = (complex_t*)calloc(n,sizeof(complex_t)  );
        complex_t* x_sampt = (complex_t*)calloc(n , sizeof(*x_sampt));
        //memset(x_sampt, 0, B * sizeof(x_sampt[0]));

        //Permute, dot, collate all in one loop.
        //for (int i = 0; i < n; i++) {
          //  x1[i] = 0;
            //x3[i] = 0;
        //}
        //随机重排
        for (int i1 = 0; i1 < n; i1++) {
            for (int j1 = 0; j1 < n; j1++) {
                x1[i1] += inv_pe_major[i1][j1] * origx[j1];
            }
        }

        int index = b;
        for (int ii = 0; ii < n; ii++) {
            x2[ii] = x1[index];
            index = (index + ai) % n;
        }
        free(x1);

        for (int i3 = 0; i3 < n; i3++) {
            for (int j3 = 0; j3 < n; j3++) {
                x3[i3] += pe_major[i3][j3] * x2[j3];
            }
        }
        free(x2);

        //滤波
        for (int i4 = 0; i4 < filter.sizet; i4++) {
            x4[i4] = x3[i4] * filter.time[i4];
        }
        free(x3);

        //频域下采样
        for (int i5 = 0; i5 < B; i5++) {
            for (int j5 = 0; j5 < n; j5++) {
                x_sampt[i5] += (complex_t) (n / B) * sa_major[i5][j5] * x4[j5];
            }
        }
        free(x4);
        //fft
        fftw_dft(x_samp, B, x_sampt,0);
        free(x_sampt);
        
        //寻找大桶位置集
        real_t* samples = (real_t*)calloc(B , sizeof(*samples));
        for (int i = 0; i < B; i++)
            if (x_samp) {
                 samples[i] = cabs2(x_samp[i]); //降采样后的功率谱
            }
            

        find_largest_indices(J, num, samples, B);

        /* #ifdef DEBUG
            debug_inner_loop(origx, n, filter, num, B, a, ai, b, J, samples);
            #endif
        */
        free(samples);
        return 0;
    }


    int inner_loop_filter_regular(int* J, int n, int num, int B, int a, int ai, int b, int loop_threshold,
        int* score, int* hits, int& hits_found, double& G_T) {
 
        // Given the set of large samples, find the locations in [n] that map there
        // and output them

        for (int i = 0; i < num; i++) {
            int low, high;
            low = (int(ceil((J[i] - 0.5) * n / B)) + n) % n;
            high = (int(ceil((J[i] + 0.5) * n / B)) + n) % n;
            int loc = timesmod(low, a, n);
            for (int j = low; j != high; j = (j + 1) % n) {
                score[loc]++;
                if (score[loc] == loop_threshold)
                    hits[hits_found++] = loc;
                loc = (loc + a) % n;
            }
        }

        return 0;
    }


    std::map<int, complex_t>
        estimate_values(const int* hits, const int& hits_found, complex_t** x_samp, const int& loops, int n, const int* permute, const int B, const int B2, const Filter& filter, const Filter& filter_Est, int location_loops) {
        std::map<int, complex_t> ans;
        real_t* values[2];
        for (int a = 0; a < 2; a++)
            values[a] = (real_t*)malloc(loops * sizeof(*values[a]));

        for (int i = 0; i < hits_found; i++) {
            int position = 0;

            for (int j = 0; j < loops; j++) {
                int cur_B = (j < location_loops) ? B : B2;
                const Filter& cur_filter = (j < location_loops) ? filter : filter_Est;
                int permuted_index = timesmod(permute[j], hits[i], n);
                int hashed_to = permuted_index / (n / cur_B);
                int dist = permuted_index % (n / cur_B);
                if (dist > (n / cur_B) / 2) {
                    hashed_to = (hashed_to + 1) % cur_B;
                    dist -= n / cur_B;
                }
                dist = (n - dist) % n;
                complex_t filter_value = cur_filter.freq[dist];// * cexp(2*PI * I * timesmod(permuteb[j], hits[i], n) / n);
                values[0][position] = creal(x_samp[j][hashed_to] / filter_value);
                values[1][position] = cimag(x_samp[j][hashed_to] / filter_value);
                position++;
                //printf("MOO %d %lf %lf: %lf %d %lf %lf+%lfj\n", hits[i], permuted_index * 1./n, hashed_to * 1./B, hashed_to * (n * 1. /B), dist, cabs(filter_value), values[0][position-1], values[1][position-1]);
            }
            int location = (loops - 1) / 2;

            for (int a = 0; a < 2; a++)
                std::nth_element(values[a], values[a] + location, values[a] + position);
            real_t realv = values[0][location];
            real_t imagv = values[1][location];
            ans[hits[i]]= (realv ,  imagv);
        }       
        //free(filter.freq);
        //free(filter_Est.freq);
        for (int a = 0; a < 2; a++)
            free(values[a]);
      
        return ans;
    }


     //std::map<int, complex_t>
    void outer_loop(complex_t* origx, const int n, const Filter& filter, const Filter& filter_Est, const int B2,
        const int num, const int B, const int loop_threshold, const int location_loops, const int loops,
        complex_t** pe_major, complex_t** inv_pe_major, complex_t** sa_major)
    {
        int* permute = (int*)calloc(loops , sizeof(*permute));
        //int* permuteb = (int*)malloc(loops * sizeof(*permuteb));
        /*
        complex_t* x_samp[loops];
        for (int i = 0; i < loops; i++) {
            if (i < location_loops)
                x_samp[i] = (complex_t*)calloc(B, sizeof(*x_samp[i]));
            else
                x_samp[i] = (complex_t*)calloc(B2, sizeof(*x_samp[i]));
        }*/
        complex_t** x_samp = new complex_t* [loops];
        for (int i = 0; i < loops; i++)
        {
            if (i < location_loops)
                x_samp[i] = new complex_t [B];
            else
                x_samp[i] = new complex_t [B2];
        }

        int hits_found = 0;

        //Variables used for timing
        double SCORE_T = 0;
        double PF_T = 0;
        double G_T = 0;
        double BC_T = 0;
        double PF_ALL = 0;
        double G_ALL = 0;
        double BC_ALL = 0;

        // calloc is faster if few pages are hit, while malloc/memset is
        // faster if most pages are hit.
        double pages_hit = double (num) * double(n / B) * double(location_loops);
        int* score;
        if (pages_hit > n / 1024) {
            score = (int*)malloc(n * sizeof(*score));
            memset(score, 0, n * sizeof(*score));
        }
        else {
            score = (int*)calloc(n, sizeof(*score));
        }
      
        // printf("Created score array: %lf\n", get_time() - DDD);

        int* hits = (int*)malloc(n * sizeof(*hits));

        double PF_LOC = 0;
        double G_LOC = 0;

        //BEGIN INNER LOOPS

        for (int i = 0; i < loops; i++) {
            int a = 0;
            int b = 0;//random() % n;
            while (gcd(a, n) != 1) {
                a = int(rand() % n);
            }
            int ai = mod_inverse(a, n);

            permute[i] = ai;
            //permuteb[i] = b;

            int perform_location = (i < location_loops);
            //assert(ALGORITHM1 || !perform_location);
            Filter cur_filter = perform_location ? filter : filter_Est;
            int cur_B = perform_location ? B : B2;

            int* J = (int*)malloc(num * sizeof(*J));

            inner_loop_locate(origx, n, cur_filter, num, cur_B, a, ai, b, x_samp[i], J, PF_T, BC_T, pe_major, inv_pe_major, sa_major);

            if (perform_location) {
                inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold, score, hits, hits_found, G_T);
            }
            free(J);

            PF_ALL += PF_T;
            BC_ALL += BC_T;
            if (perform_location) {
                PF_LOC += PF_T;
                G_LOC += G_T;
                G_ALL += G_T;
            }
        }
        //free(filter.time);
        //free(origx);

        //END INNER LOOPS

        //BEGIN ESTIMATION
        std::map<int, complex_t> ans = estimate_values(hits, hits_found, x_samp, loops, n, permute, B, B2, filter, filter_Est, location_loops);

        int num_candidates = (int)ans.size();
        //double* pst = (double*)calloc(num_candidates , sizeof(*pst));
        //complex_t* out_I = (complex_t*)calloc(num_candidates,  sizeof(*out_I));
        //complex_t* out = (complex_t*)calloc(n , sizeof(*out));
        for (int p2 = 0; p2 < n; p2++)//初始化
        {
            out[p2] = (0,0);
        }
        int counter = 0;
        for (std::map<int, complex_t>::iterator it = ans.begin(); it != ans.end(); it++) {
            if(pst){ pst[counter] = (double)it->first; }//t.first就是在迭代器中获取map键值
          
            complex_t value = (complex_t)it->second;// second指向数值// value=out_I
            if (out_I) { out_I[counter] = value; }// second指向数值      
            if (out) { out[it->first] = value; }
        }

        //END ESTIMATION
#ifdef DEBUG
        real_t* x = (real_t*)malloc(n * sizeof(*x));
        memset(x, 0, n * sizeof(*x));
        for (std::map<int, complex_t>::iterator it = ans.begin();
            it != ans.end(); it++) {
            x[it->first] = cabs(it->second);
        }

        real_t* xc = (real_t*)malloc(n * sizeof(*xc));
        memset(xc, 0, n * sizeof(*xc));
        for (int i = 0; i < n; i++)
            xc[i] = score[i] * 1. / loops;

        complex_t* xf = (complex_t*)malloc(n * sizeof(*xf));
        fftw_dft(xf, n, origx);
        for (int i = 0; i < n; i++)
            xf[i] /= n;
        plot("outer loop", "counts\ntrue signal\nreconstruction",
            Vec(xc, n), map_abs(Vec(xf, n)), Vec(x, n));

        free(x);
        free(xc);
        free(xf);
#endif //DEBUG

        //for (int i = 0; i < loops; i++)
        for (int i = 0; i < loops; i++)
            delete[] x_samp[i];
        delete [] x_samp;
        free(permute);
        free(score);
        free(hits);
        return;
        // return ans;
    }


/*
    void displayOnMATLAB(std::ostringstream& stream) {
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) })
        );
            stream.str("");
    }


    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        // Get array factory
        ArrayFactory factory;

        // Check first input argument
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE ||
            inputs[0].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("First input must be scalar double") }));
        }

        // Check second input argument
        if (inputs[1].getType() != ArrayType::DOUBLE ||
            inputs[1].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be double array") }));
        }
        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }*/
};