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

/*	float4( const float* const p )
		: value( _mm_loadu_ps( p ))
	{
	}

	float4( const float* const p, int lim )
	{
		if( lim > 2 )
		{
			if( lim > 3 )
			{
				value = _mm_loadu_ps( p );
			}
			else
			{
				value = _mm_set_ps( 0.0f, p[ 2 ], p[ 1 ], p[ 0 ]);
			}
		}
		else
		{
			if( lim == 2 )
			{
				value = _mm_set_ps( 0.0f, 0.0f, p[ 1 ], p[ 0 ]);
			}
			else
			{
				value = _mm_load_ss( p );
			}
		}
	}
*/
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

	operator float () const
	{
		return( _mm_cvtss_f32( value ));
	}

/*	void storeu( float* const p ) const
	{
		_mm_storeu_ps( p, value );
	}

	void storeu( float* const p, int lim ) const
	{
		if( lim > 2 )
		{
			if( lim > 3 )
			{
				_mm_storeu_ps( p, value );
			}
			else
			{
				_mm_storel_pi( (__m64*) p, value );
				_mm_store_ss( p + 2, _mm_movehl_ps( value, value ));
			}
		}
		else
		{
			if( lim == 2 )
			{
				_mm_storel_pi( (__m64*) p, value );
			}
			else
			{
				_mm_store_ss( p, value );
			}
		}
	}
*/
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

	__m128 value; ///< Packed value of 4 floats.
		///<
};

/**
 * SIMD rounding function, exact result.
 *
 * @param v Value to round.
 * @return Rounded SIMD value.
 */

inline float4 round( const float4& v )
{
	unsigned int prevrm = _MM_GET_ROUNDING_MODE();
	_MM_SET_ROUNDING_MODE( _MM_ROUND_NEAREST );

	const __m128 res = _mm_cvtpi32x2_ps( _mm_cvtps_pi32( v.value ),
		_mm_cvtps_pi32( _mm_movehl_ps( v.value, v.value )));

	_MM_SET_ROUNDING_MODE( prevrm );

	return( res );
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
	///< avir::CImageResizer class to perform calculation using default
	///< interleaved algorithm, using SIMD float4 type.
	///<

} // namespace avir

#endif // AVIR_FLOAT4_SSE_INCLUDED
