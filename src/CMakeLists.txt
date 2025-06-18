add_library(DWR 
testBottleneck.cpp
bottleneckGraph.cpp
timer.cpp 
)

target_include_directories(DWR PUBLIC ${CMAKE_SOURCE_DIR}/include)