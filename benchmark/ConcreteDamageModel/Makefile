CXX=g++
CXXFLAGS=-std=c++14 -O2
SOURCE_FILES=test_concrete.cpp ConcreteDamageModel.cpp

test_concrete_material: $(SOURCE_FILES)
	$(CXX) -o test_concrete.out $(SOURCE_FILES) $(CXXFLAGS)
test_eigenvalue3:
	$(CXX) -o eigenvalue3 test_eigenvalue3.cpp

clean:
	rm -f eigenvalue3