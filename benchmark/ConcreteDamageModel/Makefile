CXX=clang++
CXXFLAGS=-std=c++14 -g 
SOURCE_FILES= ConcreteDamageModel.cpp

compress_1D_cyclic: compress_1D_cyclic.cpp $(SOURCE_FILES)
	$(CXX) -o test_concrete.out compress_1D_cyclic.cpp $(SOURCE_FILES) $(CXXFLAGS)
compress_1D_monotonic: compress_1D_monotonic.cpp $(SOURCE_FILES)
	$(CXX) -o test_concrete.out compress_1D_monotonic.cpp $(SOURCE_FILES) $(CXXFLAGS)
test_eigenvalue3:
	$(CXX) -o eigenvalue3.out test_eigenvalue3.cpp

clean:
	rm -f eigenvalue3.out
	rm -f test_concrete.out