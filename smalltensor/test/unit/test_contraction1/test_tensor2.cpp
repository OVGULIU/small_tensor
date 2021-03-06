#include "../include_small_tensor.h"

#include <iostream>
using namespace std; using namespace smalltensor;
int main(int argc, char const *argv[])
{
	tensor3<float, 3,3,3> obj1;
	tensor2<float, 3,3> obj2;
	tensor2<float, 3,3> obj3;
	tensor1<float, 3> obj4;
	tensor2<float, 3,3> obj5;

	tensor2<float, 3,3> obj6;
	tensor2<float, 3,3> obj7;

	tensor2<float, 2,2> obj8;
	tensor2<float, 2,2> obj9;
	tensor2<float, 2,2> obj10;

	obj2(2,2) = 2.;
	eindex<'i'> i;
	eindex<'j'> j;
	eindex<'k'> k;

	// Test 1
	obj1(2,2,2)=1;
	obj4(2) = 3;
	obj2(i,j)=obj1(i,j,k)*obj4(k);
	ASSERT_MSG(obj2(2,2)==3,"tensor2(_i,_j) contraction1  operator error");
	
	// Test 2
	obj3(1,1)=9; obj5(1,1) = 3;
	obj2(i,j)=obj3(i,k)*obj5(k,j);
	ASSERT_MSG(obj2(1,1)==27,"tensor2(_i,_j) contraction1  operator error");

	// Test 3
	obj1(0,0,0)=32;
	obj4(0) = 2;
	obj2(i,j)=obj4(k)*obj1(k,i,j);
	ASSERT_MSG(obj2(0,0)==64,"tensor2(_i,_j) contraction1  operator error");

	// Test 4
	obj3(0,0) = 1 ; obj3(0,1) = 2 ; obj3(0,2) = 3 ; 
	obj3(1,0) = 4 ; obj3(1,1) = 5 ; obj3(1,2) = 6 ; 
	obj3(2,0) = 7 ; obj3(2,1) = 8 ; obj3(2,2) = 9 ; 
	obj6(0,0) = 2 ; obj6(0,1) = 3 ; obj6(0,2) = 4 ; 
	obj6(1,0) = 5 ; obj6(1,1) = 6 ; obj6(1,2) = 7 ; 
	obj6(2,0) = 8 ; obj6(2,1) = 9 ; obj6(2,2) = 10 ; 
	obj7(i,j) = obj3(i,k) * obj6(k,j) ; 
	ASSERT_MSG(obj7(0,0)==36,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(0,1)==42,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(0,2)==48,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(1,0)==81,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(1,1)==96,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(1,2)==111,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(2,0)==126,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(2,1)==150,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(2,2)==174,"tensor2(_i,_j) contraction1  operator error");

	// Test 5
	obj8(0,0) = 1 ; obj8(0,1) = 2 ; 
	obj8(1,0) = 3 ; obj8(1,1) = 4 ; 
	obj9(0,0) = 11 ; obj9(0,1) = 12 ; 
	obj9(1,0) = 13 ; obj9(1,1) = 14 ; 
	obj10(i,j) = obj8(i,k) * obj9(k,j) ; 
	ASSERT_MSG(obj10(0,0)==37,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj10(0,1)==40,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj10(1,0)==85,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj10(1,1)==92,"tensor2(_i,_j) contraction1  operator error");

	// Test 5
	tensor2<float, 2,3> m1;
	tensor2<float, 3,2> m2;
	tensor2<float, 2,2> m3;
	m1(0,0) = 1 ; m1(0,1) = 2 ; m1(0,2) = 3 ; 
	m1(1,0) = 4 ; m1(1,1) = 5 ; m1(1,2) = 6 ; 
	m2(0,0) = 2 ; m2(0,1) = 3 ; 
	m2(1,0) = 4 ; m2(1,1) = 5 ; 
	m2(2,0) = 6 ; m2(2,1) = 7 ; 
	m3(i,j) = m1(i,k) * m2(k,j) ; 
	// std::cerr << "m3 =" << m3 << std::endl ;
	ASSERT_MSG(m3(0,0)==28,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(m3(0,1)==34,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(m3(1,0)==64,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(m3(1,1)==79,"tensor2(_i,_j) contraction1  operator error");

	// Test 6
	tensor2<float, 3, 3> t1 ; 
	tensor2<float, 3, 3> t2 ; 
	tensor2<float, 3, 3> t3 ; 

	t1(0,0) = 1 ; t1(0,1) = 4 ; t1(0,2) = 7 ; 
	t1(1,0) = 2 ; t1(1,1) = 5 ; t1(1,2) = 8 ; 
	t1(2,0) = 3 ; t1(2,1) = 6 ; t1(2,2) = 9 ; 
	t2(0,0) = 2 ; t2(0,1) = 3 ; t2(0,2) = 4 ; 
	t2(1,0) = 5 ; t2(1,1) = 6 ; t2(1,2) = 7 ; 
	t2(2,0) = 8 ; t2(2,1) = 9 ; t2(2,2) = 10 ; 
	t3(i,j) = t1(k,i) * t2(k,j) ; 
	ASSERT_MSG(t3(0,0)==36,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(0,1)==42,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(0,2)==48,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(1,0)==81,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(1,1)==96,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(1,2)==111,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(2,0)==126,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(2,1)==150,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t3(2,2)==174,"tensor2(_i,_j) contraction1  operator error");

	// Test 7
	tensor2<float, 3,3> t4 ; 
	t4(0,0) = 1 ; t4(0,1) = 4 ; t4(0,2) = 7 ; 
	t4(1,0) = 2 ; t4(1,1) = 5 ; t4(1,2) = 8 ; 
	t4(2,0) = 3 ; t4(2,1) = 6 ; t4(2,2) = 10 ; 
	auto t5 = t4.Inv() ; 
	ASSERT_MSG(t5(0,2)==1,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t5(1,2)==-2,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t5(2,0)==1,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t5(2,1)==-2,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t5(2,2)==1,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(t4.compute_Determinant()==-3,"tensor2(_i,_j) contraction1  operator error");

	// Test 8
	obj3(0,0) = 1 ; obj3(0,1) = 2 ; obj3(0,2) = 3 ; 
	obj3(1,0) = 4 ; obj3(1,1) = 5 ; obj3(1,2) = 6 ; 
	obj3(2,0) = 7 ; obj3(2,1) = 8 ; obj3(2,2) = 9 ; 
	obj6(0,0) = 2 ; obj6(0,1) = 5 ; obj6(0,2) = 8 ; 
	obj6(1,0) = 3 ; obj6(1,1) = 6 ; obj6(1,2) = 9 ; 
	obj6(2,0) = 4 ; obj6(2,1) = 7 ; obj6(2,2) = 10 ; 
	obj7(i,j) = obj3(i,k) * obj6(j,k) ; 
	ASSERT_MSG(obj7(0,0)==36,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(0,1)==42,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(0,2)==48,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(1,0)==81,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(1,1)==96,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(1,2)==111,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(2,0)==126,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(2,1)==150,"tensor2(_i,_j) contraction1  operator error");
	ASSERT_MSG(obj7(2,2)==174,"tensor2(_i,_j) contraction1  operator error");

	cout<<"Done execution. Exiting..." <<endl;

	return 0;
}