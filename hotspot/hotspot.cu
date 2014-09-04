/*****************************************************************************
  hotspot.cu

  (c) 2014 - Nikhil R Podduturi, J. Seth Strattan
  J. Michael Cherry Lab, Department of Genetics, Stanford University School of Medicine

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <ostream>
#include <fstream>
#include <nvbio/basic/console.h>
#include <nvbio/basic/timer.h>
#include <nvbio/basic/shared_pointer.h>
#include <cuda_runtime_api.h>
#include <thrust/partition.h>

#include "hotspot_kernal.h"

void crcInit();

using namespace nvbio;
using namespace hotspot;

struct comp_chr
{
    int chr;
    comp_chr(int _chr) {chr=_chr;}

    __host__ __device__
    bool operator()(const Alignment x)
    {
        return x.ref_id != chr;
    }
};

int main(int argc, char* argv[])
{

    Timer timer;
    timer.start();

    cudaSetDeviceFlags( cudaDeviceMapHost | cudaDeviceLmemResizeToMax );

    crcInit();

    int cuda_device  = -1;

    const char* library_file    =   NULL;
    const char* density_file    =   NULL;
    const char* output_file     =   NULL;

    // hotspot defaults for now
    int low_int         = 250;
    int high_int        = 1000;
    int int_increment   = 100;
    int fuzzy_seed      = 1;
    double genome_size  = 2.55E9;
    double min_SD       = 3.0;
    bool use_fuzzy      = false;

    int arg = 1;
    while (arg < argc)
    {
        if (strcmp( argv[arg], "-device" ) == 0)
        {
            cuda_device = atoi(argv[++arg]);
            ++arg;
        }
        else if (strcmp( argv[arg], "-i" ) == 0)
        {
            library_file = argv[++arg];
            ++arg;
        }
        else if (strcmp( argv[arg], "-k" ) == 0)
        {
            density_file = argv[++arg];
            ++arg;
        }
        else if (strcmp( argv[arg], "-o" ) == 0)
        {
            output_file = argv[++arg];
            ++arg;
        }
        else if (strcmp( argv[arg], "-range") == 0)
        {
            low_int         = atoi(argv[arg + 1]);
            high_int        = atoi(argv[arg + 2]);
            int_increment   = atoi(argv[arg + 3]);
            arg = arg + 3;
        }
        else if (strcmp( argv[arg], "-bckgnmsize") == 0)
        {
            genome_size = atof(argv[++arg]);
            ++arg;
        }
        else if (strcmp( argv[arg], "-minsd") == 0)
        {
            min_SD = atof(argv[++arg]);
            ++arg;
        }
        else if (strcmp( argv[arg], "-fuzzy") == 0)
        {
            use_fuzzy = true;
            ++arg;
        }
        else if (strcmp( argv[arg], "-fuzzy-seed") == 0)
        {
            fuzzy_seed = atoi(argv[++arg]);
            ++arg;
        }
        else
            break;
    }

    // inspect and select cuda devices
    int device_count;
    cudaGetDeviceCount(&device_count);
    log_verbose(stderr, "  cuda devices : %d\n", device_count);

    if (device_count)
    {
        if (cuda_device == -1)
        {
            int best_device = 0;
            cudaDeviceProp best_device_prop;
            cudaGetDeviceProperties( &best_device_prop, best_device );
            for (int device = 0; device < device_count; ++device)
            {
                cudaDeviceProp device_prop;
                cudaGetDeviceProperties( &device_prop, device );
                if (device_prop.major >= best_device_prop.major &&
                        device_prop.minor >= best_device_prop.minor)
                {
                    best_device_prop = device_prop;
                    best_device      = device;
                }
            }
            cuda_device = best_device;
        }
        log_verbose(stderr, "  chosen device %d\n", cuda_device);
        {
            cudaDeviceProp device_prop;
            cudaGetDeviceProperties( &device_prop, cuda_device );
            log_verbose(stderr, "    device name        : %s\n", device_prop.name);
            log_verbose(stderr, "    compute capability : %d.%d\n", device_prop.major, device_prop.minor);
        }
        cudaSetDevice( cuda_device );
    }

    SharedPointer<AlignmentStream> library_stream = SharedPointer<AlignmentStream>( open_alignment_file( library_file ) );
    if (library_stream == NULL || library_stream->is_ok() == false)
    {
        log_error(stderr, "failed opening \"%s\"\n", library_file);
        exit(1);
    }

    thrust::host_vector<Alignment> h_alignments;
    const uint32 total_tag_count = library_stream->read( &h_alignments );

    log_info(stderr, "Total tags found = %d\n", total_tag_count);
    log_info(stderr, "Identifying hotspots started\n");
    thrust::device_vector<Alignment> d_alignments( h_alignments );

    thrust::host_vector<Hotspot> hotspots;
    for(int i=0; i<24; ++i)
    {
        thrust::device_vector<Alignment>::iterator iter = thrust::stable_partition(
            thrust::device,
            d_alignments.begin(),
            d_alignments.end(),
            comp_chr(i)
        );
        thrust::device_vector<Alignment> d_chr(iter, d_alignments.end());
        d_alignments.erase(iter, d_alignments.end());

        log_info(stderr, "Calculating hotspots for chr%d\n", i+1);
        compute_hotspots(d_chr, hotspots, low_int, high_int, int_increment,
            genome_size, total_tag_count, min_SD, use_fuzzy, fuzzy_seed);
    }
    timer.stop();
    log_info(stderr, "Time taken - %um:%us\n",
        uint32(timer.seconds()/60),
        uint32(timer.seconds())%60);
    return 0;
}
