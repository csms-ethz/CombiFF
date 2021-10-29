include_directories(${shared_path})

set (headers
     AdjacencyMatrixCnv.h
     AdjacencyMatrixHandler.h
     FamilyIdentifierHandler.h
     FormulaHandler.h
     Handler.h
     InputOutput.h
     NameHandler.h
     SmilesHandler.h
)

set (sources
     AdjacencyMatrixHandler.cpp
     FamilyIdentifierHandler.cpp
     FormulaHandler.cpp
     Handler.cpp
     InputOutput.cpp
     NameHandler.cpp
     SmilesHandler.cpp
)

add_library(cnvLib STATIC ${sources} ${headers})
