# coding: utf-8
# cython: language_level=3, linetrace=True, language=cpp

cimport libcpp11.chrono
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, realloc, free
from libcpp.functional cimport function
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

cimport fastani.map.map_stats
from fastani.cgi.compute_core_identity cimport (
    computeCGI,
    correctRefGenomeIds,
    computeGenomeLengths,
    outputCGI
)
from fastani.cgi.cgid_types cimport CGI_Results
from fastani.map cimport base_types
from fastani.map.compute_map cimport Map
from fastani.map.map_parameters cimport Parameters
from fastani.map.map_stats cimport recommendedWindowSize
from fastani.map.win_sketch cimport Sketch
from fastani.map.base_types cimport (
    MappingResult,
    MappingResultsVector_t,
)

from _utils cimport new_map_with_result_vector

import os



cdef extern from "<chrono>" namespace "std" nogil:
    pass


# cdef class MinimizerInfo:
#     cdef base_types.MinimizerInfo* _mi
#
#     def __cinit__(self):
#         self._mi = new base_types.MinimizerInfo()
#
#     def __dealloc__(self):
#         free(self._mi)


def main(query_file, ref_file, output_file):

    cdef Parameters parameters
    parameters.kmerSize = 16;
    parameters.minReadLength = 3000;
    parameters.alphabetSize = 4;
    parameters.minFraction = 0.2;
    parameters.threads = 1;
    parameters.p_value = 1e-03;
    parameters.percentageIdentity = 80;
    parameters.visualize = False;
    parameters.matrixOutput = False;
    parameters.referenceSize = 5000000;
    parameters.reportAll = True;

    # Add each query as an individual genome
    parameters.querySequences.push_back(os.fsencode(query_file))
    parameters.refSequences.push_back(os.fsencode(ref_file))

    # Add output file
    parameters.outFileName = os.fsencode(output_file)

    # Compute optimal window size
    parameters.windowSize = fastani.map.map_stats.recommendedWindowSize(
        parameters.p_value,
        parameters.kmerSize,
        parameters.alphabetSize,
        parameters.percentageIdentity,
        parameters.minReadLength,
        parameters.referenceSize
    )

    # Build the sketch for the reference
    cdef Sketch* reference_sketch = new Sketch(parameters)

    # Final output vector
    cdef vector[CGI_Results] finalResults
    cdef vector[CGI_Results] finalResults_local


    # query/ref output vector
    cdef MappingResultsVector_t mapResults
    cdef uint64_t totalQueryFragments = 0

    # mapper
    cdef Map* mapper = new_map_with_result_vector(parameters, reference_sketch[0], totalQueryFragments, 0, mapResults)

    #
    computeCGI(parameters, mapResults, mapper[0], reference_sketch[0], totalQueryFragments, 0, parameters.outFileName, finalResults_local);
    # correctRefGenomeIds(finalResults_local) <-- only useful for OpenMP multithreaded

    finalResults.insert(finalResults.end(), finalResults_local.begin(), finalResults_local.end());

    # report
    cdef unordered_map[string, uint64_t] genomeLengths
    computeGenomeLengths(parameters, genomeLengths)
    outputCGI(parameters, genomeLengths, finalResults, parameters.outFileName);
