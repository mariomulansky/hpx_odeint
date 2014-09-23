CXXFLAGS = -O3 -fopenmp -ffast-math -std=c++0x -I${BOOST_ROOT}
LDFLAGS = -lgomp -lboost_chrono -lboost_timer -lboost_system

perf_omp: perf_omp.cpp
	${CXX} ${CXXFLAGS} perf_omp.cpp -o perf_omp ${LDFLAGS} 
