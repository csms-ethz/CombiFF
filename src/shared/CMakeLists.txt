set(sources
     Atom.cpp
     Permutation.cpp
     smilesToMatrix.cpp
     XmlParser.cpp
     Range.cpp
     exceptions.cpp
     Matrix.cpp
     ContainerOperators.cpp)

set(headers
     XmlParser.h
     Atom.h
     ContainerOperators.h
     ElementTypes.h
     exceptions.h
     LambdaVector.h
     Range.h
     Matrix.h
     Permutation.h
     SmilesGenerator.h
     smilesToMatrix.h
     StringVector.h
     Time.h
     version.h
     )

add_library(sharedLib STATIC ${sources} ${headers})
