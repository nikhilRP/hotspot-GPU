/*
 * nvbio
 * Copyright (C) 2012-2014, NVIDIA Corporation
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * version 2 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#pragma once

#include <nvbio/basic/types.h>
#include <nvbio/basic/numbers.h>
#include <thrust/host_vector.h>

namespace nvbio {
namespace hotspot {

struct Alignment
{
    Alignment() :
        read_id( uint32(-1) ),
        read_len( 0 ),
        pos( 0 ),
        ref_id( uint32(-1) ) {}

    uint32 read_id;
    uint32 read_len;
    uint32 pos;
    uint32 ref_id;
};

struct AlignmentStream
{
    /// virtual destructor
    ///
    virtual ~AlignmentStream() {}

    /// return if the stream is ok
    ///
    virtual bool is_ok() { return false; }

    /// get the next batch
    ///
    virtual uint32 read(thrust::host_vector<Alignment> *batch) { return 0; }
};

/// open an alignment file
///
AlignmentStream* open_alignment_file(const char* file_name);

} // hotspot namespace
} // nvbio namespace
