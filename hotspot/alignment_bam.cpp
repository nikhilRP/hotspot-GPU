/*****************************************************************************
  alignment_bam.cpp

  (c) 2014 - Nikhil R Podduturi, J. Seth Strattan
  J. Michael Cherry Lab, Department of Genetics, Stanford University School of Medicine

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
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

    bool is_ok() { return true; }
    uint32 read(thrust::host_vector<Alignment> *batch)
    {
        log_info(stderr, "Loading bam file started\n");
        uint32 n_read = 0;
        BamTools::BamAlignment bam_aln;
        while(m_bam_reader.GetNextAlignment(bam_aln))
        {
            Alignment aln;

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


