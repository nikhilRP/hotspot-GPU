/*****************************************************************************
  alignment.cpp

  (c) 2014 - Nikhil R Podduturi, J. Seth Strattan
  J. Michael Cherry Lab, Department of Genetics, Stanford University School of Medicine

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include "alignment.h"

#include <zlib/zlib.h>
#include <string.h>

namespace nvbio {
namespace hotspot {

AlignmentStream* open_bam_file(const char* file_name);

AlignmentStream* open_alignment_file(const char* file_name)
{
    if (strcmp( file_name + strlen(file_name) - 4u, ".bam" ) == 0)
        return open_bam_file( file_name );

    return NULL;
}

} // hotspot namespace
} // nvbio namespace
