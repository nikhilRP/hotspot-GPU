/*****************************************************************************
  alignment.h

  (c) 2014 - Nikhil R Podduturi, J. Seth Strattan
  J. Michael Cherry Lab, Department of Genetics, Stanford University School of Medicine

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#pragma once

#include <nvbio/basic/types.h>
#include <nvbio/basic/numbers.h>
#include <thrust/host_vector.h>

namespace nvbio {
namespace hotspot {

struct Alignment
{
    Alignment() :
        read_len( 0 ),
        pos( 0 ),
        ref_id( uint8(-1) ) {}

    uint32 read_len;
    uint32 pos;
    uint8 ref_id;
};

struct AlignmentStream
{
    /// virtual destructor
    ///
    virtual ~AlignmentStream() {}

    /// return if the stream is ok
    ///
    virtual bool is_ok() { return false; }

    /// read the data from bam file
    ///
    virtual uint32 read(thrust::host_vector<Alignment> *batch) { return 0; }
};

/// open an alignment file
///
AlignmentStream* open_alignment_file(const char* file_name);

} // hotspot namespace
} // nvbio namespace
