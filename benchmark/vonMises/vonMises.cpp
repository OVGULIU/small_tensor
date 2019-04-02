#include "vonMises.h"


namespace openus{

const Mat33 vonMises::_zero_strain ;
const Mat33 vonMises::_zero_stress ;
const Mat33 vonMises::kronecker_delta("identity") ;

vonMises::vonMises(entt::registry<>* registry)
: FEMaterialBase(registry)
{

}

void 
vonMises::Initialize() {
	_stiff_tensor.setZero() ; 
	auto mu = _E / ( 2 * ( 1 + _nu ) );
	auto lambda = ( _nu * _E ) / ( ( 1 + _nu ) * ( 1 - 2 * _nu ) );
	_stiff_tensor( 0, 0, 0, 0 ) = lambda + 2 * mu;
	_stiff_tensor( 0, 0, 1, 1 ) = lambda;
	_stiff_tensor( 0, 0, 2, 2 ) = lambda;
	_stiff_tensor( 0, 1, 0, 1 ) = mu;
	_stiff_tensor( 0, 1, 1, 0 ) = mu;
	_stiff_tensor( 0, 2, 0, 2 ) = mu;
	_stiff_tensor( 0, 2, 2, 0 ) = mu;
	_stiff_tensor( 1, 0, 0, 1 ) = mu;
	_stiff_tensor( 1, 0, 1, 0 ) = mu;
	_stiff_tensor( 1, 1, 0, 0 ) = lambda;
	_stiff_tensor( 1, 1, 1, 1 ) = lambda + 2 * mu;
	_stiff_tensor( 1, 1, 2, 2 ) = lambda;
	_stiff_tensor( 1, 2, 1, 2 ) = mu;
	_stiff_tensor( 1, 2, 2, 1 ) = mu;
	_stiff_tensor( 2, 0, 0, 2 ) = mu;
	_stiff_tensor( 2, 0, 2, 0 ) = mu;
	_stiff_tensor( 2, 1, 1, 2 ) = mu;
	_stiff_tensor( 2, 1, 2, 1 ) = mu;
	_stiff_tensor( 2, 2, 0, 0 ) = lambda;
	_stiff_tensor( 2, 2, 1, 1 ) = lambda;
	_stiff_tensor( 2, 2, 2, 2 ) = lambda + 2 * mu;
}

void 
vonMises::Update() {

}

Tensor4 const&
vonMises::GetStiffnessTensor() const {
	return _stiff_tensor ;
}




} // namespace openus

