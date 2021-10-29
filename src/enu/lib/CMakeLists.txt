include_directories(${shared_path})
   
set (headers
     Alias.h
     AdjacencyMatrixEnu.h
     Enumerator.h
     EnumRun.h
     EnumSpecification.h
     Family.h
     FamilySpecifications.h
     HydrogenDistribution.h
     InputOutput.h
     InputFiles.h
     MaxFillAlg.h
     printInfo.h
     Pseudoatom.h
     readLambdas.h
     StereoGenerator.h
     Substructure.h
     matchingEnu.h
)

set (sources
     Alias.cpp
     Enumerator.cpp
     EnumRun.cpp
     Family.cpp
     FamilySpecifications.cpp
     HydrogenDistribution.cpp
     InputOutput.cpp
     MaxFillAlg.cpp
     printInfo.cpp
     Pseudoatom.cpp
     readLambdas.cpp
     StereoGenerator.cpp
     Substructure.cpp
     matchingEnu.cpp
)

add_library(enuLib STATIC ${sources} ${headers})
