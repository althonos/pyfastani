set(
    FASTANI_HEADERS
    src/cgi/include/cgid_types.hpp
    src/cgi/include/computeCoreIdentity.hpp
    src/common/clipp.h
    src/common/kseq.h
    src/common/murmur3.h
    src/common/prettyprint.hpp
    src/map/include/base_types.hpp
    src/map/include/commonFunc.hpp
    src/map/include/computeMap.hpp
    src/map/include/map_parameters.hpp
    src/map/include/map_stats.hpp
    src/map/include/MIIteratorL2.hpp
    src/map/include/parseCmdArgs.hpp
    src/map/include/slidingMap.hpp
    src/map/include/winSketch.hpp
)
foreach(_header ${FASTANI_HEADERS})
    if(_header STREQUAL "src/common/clipp.h")
        file(COPY ${_header} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/src/common)
    else()
        add_custom_command(
            OUTPUT  
                ${CMAKE_CURRENT_BINARY_DIR}/include/${_header}
            COMMAND 
                ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/patch_header.py -i ${CMAKE_CURRENT_SOURCE_DIR}/${_header} -o ${CMAKE_CURRENT_BINARY_DIR}/include/${_header}
            MAIN_DEPENDENCY
                ${_header}
            VERBATIM
        )
    endif()
    set(FASTANI_GENERATED_HEADERS ${FASTANI_GENERATED_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/include/${_header})
endforeach()

file(COPY omp.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/src)
add_library(fastani src/cgi/core_genome_identity.cpp omp.cpp ${FASTANI_GENERATED_HEADERS})
target_compile_definitions(fastani PUBLIC USE_BOOST=1 BOOST_NO_EXCEPTIONS=1)
target_include_directories(fastani PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include/src)

if(Boost_FOUND)
    target_include_directories(fastani PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(fastani PUBLIC ${Boost_LIBRARIES})
    target_link_directories(fastani PUBLIC ${Boost_LIBRARY_DIRS})
else()
    target_link_libraries(fastani PUBLIC boost_math)
endif()
