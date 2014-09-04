/*****************************************************************************************
    hotspot_kernal.h

    (c) 2014 - Nikhil R Podduturi, J. Seth Strattan
    J. Michael Cherry Lab, Department of Genetics, Stanford University School of Medicine

    Licensed under the GNU General Public License 2.0 license.
******************************************************************************************/
#include "alignment.h"
#include "hotspot.h"

#include <thrust/sort.h>
#include <thrust/random.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

using namespace nvbio;
using namespace hotspot;

struct sort_aln
{
    __host__ __device__
    bool operator()(const Alignment& a, const Alignment& b)
    {
        return a.pos < b.pos;
    }
};

struct gen_rand
{
    __device__
    double operator () (int idx)
    {
        thrust::default_random_engine rand_eng(1);
        thrust::uniform_int_distribution<int> uni_dist;
        rand_eng.discard(idx);
        int rand = uni_dist(rand_eng);
        return (double) rand/RAND_MAX;
    }
};

__global__ void compute_hs(Alignment *aln, int w_size, int n, bool use_fuzzy,
    double thresh, double *random, int count, double prob, double mean, double sd)
{
    int tid = threadIdx.x+blockIdx.x*blockDim.x;
    if( tid < n )
    {
        int contained = 1;
        float pos_avg = aln[tid].pos;
        if(tid < (n-1))
        {
            int temp = tid + 1;
            while ((aln[tid].pos + w_size >= aln[temp].pos ) && (temp <= n) )
            {
                contained = contained + 1;
                pos_avg = pos_avg + aln[temp].pos;
                temp = temp + 1;
            }
        }
        if (tid != 0)
        {
            int temp2 = tid - 1;
            while ((aln[tid].pos - w_size <= aln[temp2].pos) && (temp2 > 0))
            {
                contained = contained + 1;
                pos_avg = pos_avg + aln[temp2].pos;
                temp2 = temp2 - 1;
            }
        }
        pos_avg = pos_avg/contained;
        if (use_fuzzy)
        {
            if (fabs(contained - thresh) <= 0.5)
            {
                thresh = thresh + (random[tid] - 0.5);
            }
        }
        double cont_frac = contained/(double)count;
        double diff = fabs(cont_frac - prob);
        if(contained > thresh)
        {
            double cur_SD = (contained - 1 - mean) / sd;
        }
    }
}

template <typename V, typename T>
void compute_hotspots(V& d_chr, T& hotspots, int int_low, int int_high, int int_inc,
    float genome_size, int totaltagcount, int num_sd, bool use_fuzzy, int fuzzy_seed)
{
    thrust::sort(thrust::device, d_chr.begin(), d_chr.end(), sort_aln());
    thrust::host_vector<Alignment> h_chr = d_chr;

    double disc     = 0.0;
    int wincount    = 0;
    int n           = d_chr.size();

    //random number generation code
    thrust::device_vector<double> random_vec(n);
    if(use_fuzzy)
    {
        thrust::transform(
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(n),
            random_vec.begin(),
            gen_rand()
        );
    }

    for (int winsize = int_low; winsize <= int_high; winsize += int_inc, wincount++ )
    {
        double prob = winsize / genome_size;
        double mean = prob * totaltagcount;
        double sd = std::sqrt(prob*(1-prob)*totaltagcount);
        double detectThresh = 1 + mean + num_sd * sd;

        int gridDim = nvbio::round( (float) n/256 ) + 1;
        int w_size = winsize/2;

        Alignment* d_chr_ptr = thrust::raw_pointer_cast(&d_chr[0]);
        double* random_vec_ptr = thrust::raw_pointer_cast(&random_vec[0]);
        compute_hs<<<gridDim, 256, 1>>> (d_chr_ptr, w_size, n, use_fuzzy,
            detectThresh, random_vec_ptr, totaltagcount, prob, mean, sd);
    }
    log_info(stderr, "  computed for tags - %lu\n", n);
    random_vec.clear();
    random_vec.shrink_to_fit();
}
