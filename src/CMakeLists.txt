
# compiler options
add_compile_options(-std=c++14)
add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-Werror)
add_compile_options(-Wfatal-errors)
add_compile_options(-pedantic)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

add_definitions (${SEQAN_DEFINITIONS})

# include file location
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${SEQAN_INCLUDE_DIRS})
include_directories("${scATACseq_SOURCE_DIR}/src/Matrix")
include_directories("${scATACseq_SOURCE_DIR}/src/Clustering")
include_directories("${scATACseq_SOURCE_DIR}/src/Random")
include_directories("${scATACseq_SOURCE_DIR}/src/Parallel")
include_directories("${scATACseq_SOURCE_DIR}/src/Statistics")
include_directories("${scATACseq_SOURCE_DIR}/src/GUI")
include_directories("${scATACseq_SOURCE_DIR}/src/Applications")
include_directories("${scATACseq_SOURCE_DIR}/src/Matrix")
include_directories("${scATACseq_SOURCE_DIR}/src/GenomicTools")
include_directories("${scATACseq_SOURCE_DIR}/src/Utility")
# compile modules into static libraries
## set output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/lib")
## build instructions
add_library(Clustering "Clustering/DataLayer.cpp"
                       "Clustering/Data2DLayer.cpp"
                       "Clustering/Data3DLayer.cpp"
                       "Clustering/ReadLayer.cpp"
                       "Clustering/SequenceLayer.cpp"
                       "Clustering/ConsensusSequenceLayer.cpp"
                       "Clustering/ModelComputer.cpp"
                       "Clustering/ReadModelComputer.cpp"
                       "Clustering/SequenceModelComputer.cpp"
                       "Clustering/ConsensusSequenceModelComputer.cpp"
                       "Clustering/EMBase.cpp"
                       "Clustering/EMRead.cpp"
                       "Clustering/EMSequence.cpp"
                       "Clustering/EMConsensusSequence.cpp"
                       "Clustering/EMJoint.cpp")

add_library(Random       "Random/Random.cpp"
                         "Random/RandomNumberGenerator.cpp")
add_library(Parallel     "Parallel/ThreadPool.cpp")
add_library(Statistics   "Statistics/Statistics.cpp"
                         "Statistics/KmersStatistics.cpp") 
add_library(GUI          "GUI/ConsoleProgressBar.cpp"
                         "GUI/Diplayable.cpp"
                         "GUI/Updatable.cpp")
add_library(GenomicTools "GenomicTools/MatrixCreator.cpp"
                         "GenomicTools/ReadMatrixCreator.cpp"
                         "GenomicTools/CorrelationMatrixCreator.cpp"
                         "GenomicTools/ClassReadDataCreator.cpp"
                         "GenomicTools/ClassSequenceDataCreator.cpp"
                         "GenomicTools/SequenceMatrixCreator.cpp"
                         "GenomicTools/GenomeRegion.cpp")
add_library(Utility      "Utility/matrix_utility.cpp"
                         "Utility/dna_utility.cpp")

## resolve dependencies
target_link_libraries(Utility      ${SEQAN_LIBRARIES})
target_link_libraries(Clustering   Utility Random Statistics GUI Parallel ${SEQAN_LIBRARIES})
target_link_libraries(Parallel     Threads::Threads)
target_link_libraries(GenomicTools Utility ${SEQAN_LIBRARIES})

# executables
## a toy for seqan
set(EXE_MAIN_TEST "main_test")
add_executable(${EXE_MAIN_TEST} "main_test.cpp")
target_link_libraries(${EXE_MAIN_TEST} GenomicTools Clustering)
set_target_properties(${EXE_MAIN_TEST} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## a toy for correlation matrix
set(EXE_MAIN_CORMAT "main_cormat")
add_executable(${EXE_MAIN_CORMAT} "main_cormat.cpp")
target_link_libraries(${EXE_MAIN_CORMAT} ${SEQAN_LIBRARIES} Utility GenomicTools Random)
set_target_properties(${EXE_MAIN_CORMAT} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an application to read binary matrices and display them in txt
set(EXE_MAIN_MATBIN2TXT "MatrixBinToTxt")
add_executable(${EXE_MAIN_MATBIN2TXT} "Applications/MatrixBinToTxtApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_MAIN_MATBIN2TXT} Utility Boost::program_options)
set_target_properties(${EXE_MAIN_MATBIN2TXT} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an application to know the indices of 0-signal rows in integer 2D matrices
set(EXE_MAIN_WHICHNULLROWS "WhichNullRows")
add_executable(${EXE_MAIN_WHICHNULLROWS} "Applications/WhichNullRowsApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_MAIN_WHICHNULLROWS} Boost::program_options)
set_target_properties(${EXE_MAIN_WHICHNULLROWS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an application to create a matrix from BED and a BAM file
set(EXE_MAIN_BAMMATRIX "CorrelationMatrixCreator")
add_executable(${EXE_MAIN_BAMMATRIX} "Applications/CorrelationMatrixCreatorApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_MAIN_BAMMATRIX} GenomicTools Utility Boost::program_options)
set_target_properties(${EXE_MAIN_BAMMATRIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an application to create a sequence matrix from BED and a fasta file
set(EXE_MAIN_SEQMATRIX "SequenceMatrixCreator")
add_executable(${EXE_MAIN_SEQMATRIX} "Applications/SequenceMatrixCreatorApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_MAIN_SEQMATRIX} GenomicTools Utility Boost::program_options)
set_target_properties(${EXE_MAIN_SEQMATRIX} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an EMRead standalone
set(EXE_EMREAD "EMRead")
add_executable(${EXE_EMREAD}   "Applications/EMReadApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_EMREAD} Clustering Boost::program_options)
set_target_properties(${EXE_EMREAD} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an EMSequence standalone
set(EXE_EMSEQ "EMSequence")
add_executable(${EXE_EMSEQ}   "Applications/EMSequenceApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_EMSEQ} Clustering Utility Boost::program_options)
set_target_properties(${EXE_EMSEQ} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an EMConsensusSequence standalone
set(EXE_EMCONSSEQ "EMConsensusSequence")
add_executable(${EXE_EMCONSSEQ}   "Applications/EMConsensusSequenceApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_EMCONSSEQ} Clustering Utility Boost::program_options)
set_target_properties(${EXE_EMCONSSEQ} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an EMJoint standalone
set(EXE_EMJOINT "EMJoint")
add_executable(${EXE_EMJOINT}   "Applications/EMJointApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_EMJOINT} Clustering Utility Boost::program_options)
set_target_properties(${EXE_EMJOINT} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an executable to compute data models from the data and the post prob of an EM classification
set(EXE_PROB2REF "ProbToModel")
add_executable(${EXE_PROB2REF}   "Applications/ProbToModelApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_PROB2REF} Clustering Utility Boost::program_options)
set_target_properties(${EXE_PROB2REF} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
## an executable to extend read models from an EM classification
set(EXE_READMODELEXTENDER "ReadModelExtender")
add_executable(${EXE_READMODELEXTENDER}   "Applications/ReadModelExtenderApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_READMODELEXTENDER} Clustering GenomicTools Utility Boost::program_options)
set_target_properties(${EXE_READMODELEXTENDER} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")

## an executable to extend read models from an EM classification
set(EXE_SEQUENCEMODELEXTENDER "SequenceModelExtender")
add_executable(${EXE_SEQUENCEMODELEXTENDER}   "Applications/SequenceModelExtenderApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_SEQUENCEMODELEXTENDER} Clustering GenomicTools Utility Boost::program_options)
set_target_properties(${EXE_SEQUENCEMODELEXTENDER} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")

## an executable to extract read data from an EM classification
set(EXE_CLASSREADDATACREATOR "ClassReadDataCreator")
add_executable(${EXE_CLASSREADDATACREATOR}   "Applications/ClassReadDataCreatorApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_CLASSREADDATACREATOR} Clustering GenomicTools Utility Boost::program_options)
set_target_properties(${EXE_CLASSREADDATACREATOR} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")

## an executable to extract sequence data from an EM classification
set(EXE_CLASSSEQDATACREATOR "ClassSequenceDataCreator")
add_executable(${EXE_CLASSSEQDATACREATOR}   "Applications/ClassSequenceDataCreatorApplication.cpp" "Applications/ApplicationInterface.cpp")
target_link_libraries(${EXE_CLASSSEQDATACREATOR} Clustering GenomicTools Utility Boost::program_options)
set_target_properties(${EXE_CLASSSEQDATACREATOR} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")


## a test suite
set(EXE_TESTS "unittests")
add_executable(${EXE_TESTS} "unittests.cpp"
                            "Unittests/unittests_matrix.cpp"
                            "Unittests/unittests_genomictools.cpp")
target_link_libraries(${EXE_TESTS} ${UNITTEST_LIB} ${SEQAN_LIBRARIES} GenomicTools)
set_target_properties(${EXE_TESTS} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${scATACseq_SOURCE_DIR}/bin")
