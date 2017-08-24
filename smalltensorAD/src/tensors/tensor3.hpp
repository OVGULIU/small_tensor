#pragma once
#include "../expressions/expressions.h"

template <typename __dat_t, std::size_t __d1, std::size_t __d2, std::size_t __d3>
class tensor3
{
public:
	__dat_t* __restrict__ _data;

	tensor3(): _data{new __dat_t[__d1*__d2*__d3]}{}
	tensor3(tensor3 const& rhs_): _data{new __dat_t[__d1*__d2*__d3]}{
		DEBUG_MSG("tensor3 copy constructor is called");
		std::memcpy(_data, rhs_._data, sizeof(__dat_t)*__d1*__d2*__d3);
	}
	tensor3& operator=(tensor3 const& rhs_){
		DEBUG_MSG("tensor3 copy assignment operator is called");
		if(this != &rhs_){
	        std::memcpy(_data, rhs_._data, sizeof(__dat_t)*__d1*__d2*__d3);		
		}
		return *this;
	}
    tensor3(tensor3&& rhs_) noexcept
    : _data{new __dat_t[__d1*__d2*__d3]}
    {
    	std::swap(_data, rhs_._data);
    	DEBUG_MSG("tensor3 move constructor is called");
    }
    tensor3& operator=(tensor3&& rhs_) noexcept{
    	DEBUG_MSG("tensor3 move assignment operator is called");
    	if(this != &rhs_){
    		std::swap(_data, rhs_._data);
    	}
    	return *this;
    }
	~tensor3(){
		if (_data!=nullptr){
		    delete[] _data;
		    _data=nullptr;
		}
	}

	template <typename val_type>
	tensor3(ad_graph<val_type>& graph_, val_type value_)
	: _data{new __dat_t[__d1*__d2*__d3]}
	{
		DEBUG_MSG("tensor3 constructor with Graph is called");
		for (std::size_t n1 = 0; n1 < __d1; ++n1){
			for (std::size_t n2 = 0; n2 < __d2; ++n2){
				for (std::size_t n3 = 0; n3 < __d3; ++n3){
					(*this)(n1,n2,n3) = __dat_t(graph_, value_) ; 
				}
			}
		}
	}

	template <typename val_type>
	tensor3(ad_graph<val_type>& graph_, int value_)
	: _data{new __dat_t[__d1*__d2*__d3]}
	{
		DEBUG_MSG("tensor3 constructor with Graph is called");
		for (std::size_t n1 = 0; n1 < __d1; ++n1){
			for (std::size_t n2 = 0; n2 < __d2; ++n2){
				for (std::size_t n3 = 0; n3 < __d3; ++n3){
					(*this)(n1,n2,n3) = __dat_t(graph_, value_) ; 
				}
			}
		}
	}

	inline __dat_t& operator()(std::size_t n1_, std::size_t n2_, std::size_t n3_){
		ASSERT_MSG(n1_< __d1 && n2_ < __d2 && n3_ < __d3, "tensor3() index out of bounds in lvalue. ");
		return _data[ n1_ * __d2 *__d3 + n2_ * __d3 + n3_];
	}
	inline __dat_t operator()(std::size_t n1_, std::size_t n2_, std::size_t n3_)const{
		ASSERT_MSG(n1_< __d1 && n2_ < __d2 && n3_ < __d3, "tensor3() index out of bounds in rvalue. ");
		return _data[ n1_ * __d2 *__d3 + n2_ * __d3 + n3_];
	}
	template <char i, char j, char k>
	inline expr3<__dat_t, __d1, __d2, __d3, i, j, k>& operator()(Ident<i> i_, Ident<j> j_, Ident<k> k_){
        return static_cast<expr3<__dat_t, __d1, __d2, __d3, i, j, k>&>(*this);
	}
	template <char i, char j, char k>
	inline expr3<__dat_t, __d1, __d2, __d3, i, j, k> const& operator()(Ident<i> i_, Ident<j> j_, Ident<k> k_)const{
        return static_cast<expr3<__dat_t, __d1, __d2, __d3, i, j, k>const&>(*this);
	}
	template <char i, char j>
	inline expr1<__dat_t, __d1, i> operator()(Ident<i> i_, Ident<j> j_, Ident<j/*same*/> k_){
		ASSERT_MSG(__d2 == __d3, "Dimension size should be equal for dummy indices. ");
		typedef expr1<__dat_t, __d1, i> ret_type;
		ret_type ret_i(*((*this)(0,0,0).get_graph()), 0.) ;
		for (std::size_t n1 = 0; n1 < __d1; ++n1){
			for (std::size_t n2 = 0; n2 < __d2; ++n2)
			{
				ret_i(n1) = ret_i(n1) + (*this)(n1,n2,n2);
			}
		}
        return ret_i;
	}

	template<typename scalar_type>
	inline tensor3& operator*=(scalar_type const& scalar_){
		for (std::size_t n1 = 0; n1 < __d1; ++n1){
			for (std::size_t n2 = 0; n2 < __d2; ++n2){
				for (std::size_t n3 = 0; n3 < __d3; ++n3){
					(*this)(n1,n2,n3) = (*this)(n1,n2,n3) * scalar_ ;
				}
			}
		}
		return (*this);
	}
};

