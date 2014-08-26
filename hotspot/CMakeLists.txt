nvbio_module(hotspot)

addsources(
hotspot.cu
)

cuda_add_executable(hotspot ${idr_srcs})
target_link_libraries(hotspot nvbio zlibstatic crcstatic ${SYSTEM_LINK_LIBRARIES})
