
add_executable(bench_dynamic_DWR bench_dynamic_DWR.cpp)
target_link_libraries(bench_dynamic_DWR DWR gtest gtest_main TBB::tbb)

add_executable(bench_static_query bench_static_query.cpp)
target_link_libraries(bench_static_query gtest gtest_main)

add_executable(bench_dynamic_HIndex bench_dynamic_HIndex.cpp)
target_link_libraries(bench_dynamic_HIndex gtest gtest_main TBB::tbb)

set_target_properties(${TARGET} PROPERTIES
COMPILE_FLAGS "-pthread"
LINK_FLAGS "-pthread")


enable_testing()