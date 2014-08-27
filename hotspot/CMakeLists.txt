nvbio_module(hotspot)

addsources(
    alignment.h
    alignment.cpp
    alignment_bam.cpp
    hotspot.cu
)

cuda_add_executable(hotspot ${hotspot_srcs})
target_link_libraries(hotspot bamtools nvbio zlibstatic crcstatic ${SYSTEM_LINK_LIBRARIES})
