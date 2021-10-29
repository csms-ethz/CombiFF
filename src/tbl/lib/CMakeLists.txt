include_directories(${shared_path})
   
set(headers
     MatrixTbl.h
     InputOutput.h
     Link.h
     MoleculeDecomposer.h
     TblFragment.h
     TblRun.h
     TopologicalProperty.h
     AtomTypes.h
     InputFiles.h
     TblAtom.h
     moleculeMacros.h
     Linkage.h
     mtb.h
     )

set(sources
     AtomTypes.cpp
     InputOutput.cpp
     TblFragment.cpp
     Link.cpp
     MoleculeDecomposer.cpp
     TblRun.cpp
     TopologicalProperty.cpp
     moleculeMacros.cpp
     mtb.cpp
     TblAtom.cpp)

add_library(tblLib STATIC ${sources} ${headers})
