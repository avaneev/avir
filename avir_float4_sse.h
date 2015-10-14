//$ nocpp

#ifndef AVIR_FLOAT4_SSE_INCLUDED
#define AVIR_FLOAT4_SSE_INCLUDED

#include <xmmintrin.h>

namespace avir {

/**
 * @brief SIMD packed 4-float type.
 *
 * This class implements a packed 4-float type that can be used to perform
 * parallel computation using SIMD instructions on SSE-enabled processors.
 * This class can be used as the "fptype" argument of the avir::CImageResizer
 * class.
 */

class float4
{
public:
	float4()
	{
	}

	float4( const float4& s )
		: value( s.value )
	{
	}

	float4( const __m128 s )
		: value( s )
	{
	}

	float4( const float s )
		: value( _mm_set1_ps( s ))
	{
	}

	float4& operator = ( const float4& s )
	{
		value = s.value;
		return( *this );
	}

	float4& operator = ( const __m128 s )
	{
		value = s;
		return( *this );
	}

	float4& operator = ( const float s )
	{
		value = _mm_set1_ps( s );
		return( *this );
	}

	float4& operator += ( const float4& s )
	{
		value = _mm_add_ps( value, s.value );
		return( *this );
	}

	float4& operator -= ( const float4& s )
	{
		value = _mm_sub_ps( value, s.value );
		return( *this );
	}

	float4& operator *= ( const float4& s )
	{
		value = _mm_mul_ps( value, s.value );
		return( *this );
	}

	float4& operator /= ( const float4& s )
	{
		value = _mm_div_ps( value, s.value );
		return( *this );
	}

	float4 operator + ( const float4& s ) const
	{
		return( _mm_add_ps( value, s.value ));
	}

	float4 operator - ( const float4& s ) const
	{
		return( _mm_sub_ps( value, s.value ));
	}

	float4 operator * ( const float4& s ) const
	{
		return( _mm_mul_ps( value, s.value ));
	}

	float4 operator / ( const float4& s ) const
	{
		return( _mm_div_ps( value, s.value ));
	}

	operator float () const
	{
		float v;
		_mm_store_ss( &v, value );
		return( v );
	}

	__m128 value; ///< Packed value of 4 floats.
		///<
};

/**
 * SIMD rounding function, based on the trunc() function. Biased result.
 *
 * @param v Value to round.
 * @return Rounded SIMD value. Some bias may be introduced.
 */

inline float4 round( const float4& v )
{
	const int SignMask = _mm_movemask_ps( v.value );

	static const float Mults[ 16 ][ 4 ] = {
		{  1.0,  1.0,  1.0,  1.0 },
		{ -1.0,  1.0,  1.0,  1.0 },
		{  1.0, -1.0,  1.0,  1.0 },
		{ -1.0, -1.0,  1.0,  1.0 },
		{  1.0,  1.0, -1.0,  1.0 },
		{ -1.0,  1.0, -1.0,  1.0 },
		{  1.0, -1.0, -1.0,  1.0 },
		{ -1.0, -1.0, -1.0,  1.0 },
		{  1.0,  1.0,  1.0, -1.0 },
		{ -1.0,  1.0,  1.0, -1.0 },
		{  1.0, -1.0,  1.0, -1.0 },
		{ -1.0, -1.0,  1.0, -1.0 },
		{  1.0,  1.0, -1.0, -1.0 },
		{ -1.0,  1.0, -1.0, -1.0 },
		{  1.0, -1.0, -1.0, -1.0 },
		{ -1.0, -1.0, -1.0, -1.0 }
	};

	const __m128 msign = _mm_loadu_ps( Mults[ SignMask ]);
	const __m128 bias = _mm_add_ps( _mm_set1_ps( 0.5 ),
		_mm_mul_ps( v.value, msign ));

	const __m64 lo = _mm_cvttps_pi32( bias );
	const __m64 hi = _mm_cvttps_pi32( _mm_movehl_ps( bias, bias ));

	const __m128 res = _mm_cvtpi32x2_ps( lo, hi );

	return( _mm_mul_ps( res, msign ));
}

/**
 * SIMD function "clamps" (clips) the specified packed values so that they are
 * not lesser than "minv", and not greater than "maxv".
 *
 * @param Value Value to clamp.
 * @param minv Minimal allowed value.
 * @param maxv Maximal allowed value.
 * @return The clamped value.
 */

inline float4 clamp( const float4& Value, const float4 minv,
	const float4 maxv )
{
	return( _mm_min_ps( _mm_max_ps( Value.value, minv.value ), maxv.value ));
}

typedef fpclass_def< avir :: float4, float > fpclass_float4; ///<
	///< Class that can be used as the "fpclass" template parameter of the
	///< CImageResizer class to perform calculation using non-SIMD algorithms,
	///< but using SIMD float4 type.
	///<

} // namespace avir

#endif // AVIR_FLOAT4_SSE_INCLUDED
