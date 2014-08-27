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

#include "alignment.h"

#include <contrib/bamtools/BamReader.h>
#include <nvbio/basic/console.h>
#include <crc/crc.h>

namespace nvbio {
namespace hotspot {

struct BAMAlignmentStream : public AlignmentStream
{
    BAMAlignmentStream(const char* file_name)
    {
        log_verbose(stderr, "opening BAM file \"%s\"... started\n", file_name);
        m_bam_reader.Open( file_name );
        m_offset = 0;
        log_verbose(stderr, "opening BAM file \"%s\"... done\n", file_name);
    }

    // return if the stream is ok
    //
    bool is_ok() { return true; } // TODO: add a mechanism to bamtools to know whether the file opened correctly

    // get the next batch
    //
    uint32 read(thrust::host_vector<Alignment> *batch)
    {
        log_info(stderr, "Loading bam file started\n");
        uint32 n_read = 0;
        BamTools::BamAlignment bam_aln;
        while(m_bam_reader.GetNextAlignment(bam_aln))
        {
            Alignment aln;

            aln.read_id    =   uint32( crcCalc( bam_aln.Name.c_str(), uint32(bam_aln.Name.length()) ) );
            aln.read_len   =   bam_aln.Length;
            aln.pos        =   bam_aln.Position;
            aln.ref_id     =   bam_aln.RefID;
            batch->push_back(aln);
            n_read++;
        }
        log_info(stderr, "  Done\n");
        return n_read;
    }

    BamTools::BamReader m_bam_reader;
    uint32              m_offset;
};

AlignmentStream* open_bam_file(const char* file_name)
{
    return new BAMAlignmentStream( file_name );
}

} // hotspot namespace
} // nvbio namespace


