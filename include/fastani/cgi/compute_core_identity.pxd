from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from fastani.cgi.cgid_types cimport MappingResult_CGI, CGI_Results
from fastani.map.base_types cimport MappingResultsVector_t
from fastani.map.compute_map cimport Map
from fastani.map.map_parameters cimport Parameters
from fastani.map.win_sketch cimport Sketch


cdef extern from "cgi/include/computeCoreIdentity.hpp" namespace "cgi" nogil:


    void reviseRefIdToGenomeId(vector[MappingResult_CGI] &shortResults, Sketch &refSketch)
    void computeGenomeLengths(Parameters &parameters, unordered_map[string, uint64_t] &genomeLengths)

    void outputVisualizationFile(
        Parameters &parameters,
        vector[MappingResult_CGI] &mappings_2way,
        Map &mapper,
        Sketch &refSketch,
        uint64_t queryFileNo,
        string &fileName
    )

    void computeCGI(
        Parameters &parameters,
        MappingResultsVector_t &results,
        Map &mapper,
        Sketch &refSketch,
        uint64_t totalQueryFragments,
        uint64_t queryFileNo,
        string &fileName,
        vector[CGI_Results] &CGI_ResultsVector
    )

    void outputCGI(
        Parameters &parameters,
        unordered_map[string, uint64_t] &genomeLengths,
        vector[CGI_Results] &CGI_ResultsVector,
        string &fileName
    )

    void outputPhylip(
        Parameters &parameters,
        unordered_map[string, uint64_t] &genomeLengths,
        vector[CGI_Results] &CGI_ResultsVector,
        string &fileName
    )

    void splitReferenceGenomes(
        Parameters &parameters,
        vector[Parameters] &parameters_split
    )

    void correctRefGenomeIds(vector[CGI_Results] &CGI_ResultsVector)
