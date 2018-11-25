CXX=clang++
CXXFLAGS=-std=c++14 -g 
SOURCE_FILES= ConcreteDamageModel.cpp

monotonic_1D_compress: monotonic_1D_compress.cpp $(SOURCE_FILES)
	$(CXX) -o test_concrete.out monotonic_1D_compress.cpp $(SOURCE_FILES) $(CXXFLAGS)
test_eigenvalue3:
	$(CXX) -o eigenvalue3.out test_eigenvalue3.cpp

clean:
	rm -f eigenvalue3.out
	rm -f test_concrete.out