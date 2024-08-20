set(BASICALGORITHMS_HEADERS 
  arrayRef.h 
  containerHelper.h 
  basicAlgorithms.h 
  averagingBuffer.h 
  fileIO.h 
  filterIterator.h
  kmeans.h 
  reversibleHeap.h 
  stringHelper.h 
  vectorStack.h 
  range.h 
  disjointSet.h
  graphSearchAlgorithms.h
)

set(BASICALGORITHMS_SOURCES 
  containerHelper.cpp 
  averagingBuffer.cpp 
  fileIO.cpp 
  stringHelper.cpp 
  disjointSet.cpp
)

add_libpgo_lib(basicAlgorithms "${BASICALGORITHMS_SOURCES}" "${BASICALGORITHMS_HEADERS}")