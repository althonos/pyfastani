add_subdirectory(_sequtils)

cython_extension(_fasta   LINKS sequtils         )
cython_extension(_fastani LINKS fastani sequtils SOURCES _utils.cpp)