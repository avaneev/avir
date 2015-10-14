//$ nobt
//$ nocpp

/**
 * @file avir.h
 *
 * @brief The "main" and only inclusion file with all required classes and
 * functions.
 *
 * This is the "main" inclusion file for the "AVIR" image resizer. This
 * inclusion file contains implementation of the AVIR image resizing algorithm
 * in its entirety. Also includes several classes and functions that can be
 * useful elsewhere.
 *
 * AVIR Copyright (c) 2015 Aleksey Vaneev
 *
 * @mainpage
 *
 * @section intro_sec Introduction
 *
 * Me, Aleksey Vaneev, is happy to offer you an open source image resizing
 * library which has reached a production level of quality, and is ready to be
 * incorporated into any project. This library features routines for both
 * down- and upsizing of 8- and 16-bit, 1 to 4-channel images. Image resizing
 * routines were implemented in multi-platform C++ code, and have a high level
 * of optimality. Beside resizing, this library offers a sub-pixel shift
 * operation.
 *
 * The resizing algorithm at first produces 2X upsized image (relative to the
 * source image size, or relative to the destination image size if downsizing
 * is performed) and then performs interpolation using a bank of sinc
 * function-based fractional delay filters. On the last stage a correction
 * filter is applied which fixes smoothing introduced on previous steps.
 *
 * An important element utilized by this algorithm is the so called Peaked
 * Cosine window function, which is applied over sinc function in all filters.
 * Please consult the documentation for more details.
 *
 * AVIR is devoted to women. Your digital photos can look good at any size!
 *
 * @section license License
 *
 * AVIR License Agreement
 *
 * AVIR Copyright (c) 2015 Aleksey Vaneev
 *
 * 1. AVIR image resizing software library and its parts and associated
 * documentation files (collectively, "AVIR") is licensed, not sold. AVIR and
 * all rights are owned by AVIR copyright holder and are protected by
 * copyright law and international copyright treaties. You (individual person
 * or a company) acknowledge that no title to the intellectual property in
 * AVIR is transferred to you. You further acknowledge that title and full
 * ownership rights to AVIR will remain the exclusive property of AVIR
 * copyright holder and you will not acquire any rights to AVIR except as
 * expressly set forth in this License Agreement.
 *
 * 2. This License Agreement is only applicable to non-commercial and
 * low-income commercial use of AVIR. This means that you have to sign a
 * special license agreement with AVIR copyright holder if you or your
 * company sells and/or licenses (directly or via resellers) software products
 * and/or non-software products with firmware where AVIR is included or may be
 * included, and if or when the combined sales of such products reach at your
 * end a sum of more than 40000 (forty thousand) US dollars per year (or
 * equivalent value in your native currency). Sales also include subscriptions
 * where a subscribed end user may get a functioning, and/or feature-extended,
 * and/or advertisement-free product where AVIR is included or may be
 * included. Software products also include user-operated server-based image
 * and video editing applications where AVIR is included or may be included.
 * Donation-based, but otherwise free from fees, software products where AVIR
 * is included or may be included are treated herein as non-commercial.
 *
 * 3. The execution of AVIR image resizing function in non-commercial,
 * commercial and any other environments (operations) of any size and for any
 * purpose is allowed without a fee.
 *
 * 4. If you have a considerable development budget, and/or make a
 * considerable profit from donations (grants), and/or make advertising profit
 * associated with software where AVIR is included, or if AVIR image resizing
 * function is an important function (e.g. on a photo hosting, a social web
 * site) which is performed using proprietary software product where AVIR is
 * included, then consider voluntarily signing a special license agreement
 * with AVIR copyright holder.
 *
 * 5. If you sublicense to a non-end user third party, or transfer rights to
 * your software where AVIR or its parts are included, to a third party, this
 * third party should agree to this License Agreement or sign a special
 * license agreement with AVIR copyright holder, or remove AVIR dependence
 * from the received software.
 *
 * 6. Permission is hereby granted, free of charge, to any person or company
 * (subject to the restrictions outlined above) obtaining a copy of AVIR, to
 * deal in AVIR without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense (except as
 * noted above), and/or sell copies of AVIR, and to permit persons to whom
 * AVIR is furnished to do so, subject to the following conditions:
 * The above copyright notice and this License Agreement shall be included in
 * all copies or portions of AVIR.
 *
 * 7. AVIR IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH AVIR OR THE USE OR OTHER DEALINGS IN
 * AVIR.
 *
 * 8. AVIR copyright holder's e-mail: aleksey.vaneev@gmail.com
 * AVIR copyright holder's registered address:
 * Aleksey Vaneev
 * 6-54 Lenin Street,
 * Syktyvkar, Komi Republic 167000
 * Russia
 *
 * Please credit the author of this library in your documentation in the
 * following way: "AVIR image resizing algorithm designed by Aleksey Vaneev"
 *
 * @version 1.3
 */

#ifndef AVIR_CIMAGERESIZER_INCLUDED
#define AVIR_CIMAGERESIZER_INCLUDED

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

namespace avir {

#if !defined( M_PI )
	/**
	 * The macro equals to "pi" constant, fits 53-bit floating point mantissa.
	 */

	#define M_PI 3.1415926535897932
#endif // M_PI

#if !defined( M_PId2 )
	/**
	 * The macro equals to "pi divided by 2" constant, fits 53-bit floating
	 * point mantissa.
	 */

	#define M_PId2 1.5707963267948966
#endif // M_PId2

/**
 * Rounding function, based on the floor() function. Biased result.
 *
 * @param d Value to round.
 * @return Rounded value. Some bias may be introduced.
 */

template< class T >
inline T round( const T d )
{
	return( d < 0.0 ? -floor( (T) 0.5 - d ) : floor( d + (T) 0.5 ));
}

/**
 * Template function "clamps" (clips) the specified value so that it is not
 * lesser than "minv", and not greater than "maxv".
 *
 * @param Value Value to clamp.
 * @param minv Minimal allowed value.
 * @param maxv Maximal allowed value.
 * @return The clamped value.
 */

template< class T >
inline T clamp( const T& Value, const T minv, const T maxv )
{
	if( Value < minv )
	{
		return( minv );
	}
	else
	if( Value > maxv )
	{
		return( maxv );
	}
	else
	{
		return( Value );
	}
}

/**
 * Function converts (via typecast) specified array of type T1 values of
 * length l into array of type T2 values. If T1 is the same as T2, copy
 * operation is performed. When copying data at overlapping address spaces,
 * "op" should be lower than "ip".
 *
 * @param ip Input buffer.
 * @param[out] op Output buffer.
 * @param l The number of elements to copy.
 * @param ip Input buffer pointer increment.
 * @param op Output buffer pointer increment.
 */

template< class T1, class T2 >
inline void copyArray( const T1* ip, T2* op, int l,
	const int ipinc = 1, const int opinc = 1 )
{
	while( l > 0 )
	{
		*op = (T2) *ip;
		op += opinc;
		ip += ipinc;
		l--;
	}
}

/**
 * Function adds values located in array "ip" to array "op".
 *
 * @param ip Input buffer.
 * @param[out] op Output buffer.
 * @param l The number of elements to add.
 * @param ip Input buffer pointer increment.
 * @param op Output buffer pointer increment.
 */

template< class T1, class T2 >
inline void addArray( const T1* ip, T2* op, int l,
	const int ipinc = 1, const int opinc = 1 )
{
	while( l > 0 )
	{
		*op += *ip;
		op += opinc;
		ip += ipinc;
		l--;
	}
}

/**
 * Function calculates frequency response of the specified FIR filter at the
 * specified circular frequency. Phase can be calculated as atan2( im, re ).
 *
 * @param flt FIR filter's coefficients.
 * @param fltlen Number of coefficients (taps) in the filter.
 * @param th Circular frequency [0; pi].
 * @param[out] re0 Resulting real part of the complex frequency response.
 * @param[out] im0 Resulting imaginary part of the complex frequency response.
 * @param fltlat Filter's latency in samples (taps).
 */

template< class T >
inline void calcFIRFilterResponse( const T* flt, int fltlen,
	const double th, double& re0, double& im0, const int fltlat = 0 )
{
	const double sincr = 2.0 * cos( th );
	double cvalue1;
	double svalue1;

	if( fltlat == 0 )
	{
		cvalue1 = 1.0;
		svalue1 = 0.0;
	}
	else
	{
		cvalue1 = cos( -fltlat * th );
		svalue1 = sin( -fltlat * th );
	}

	double cvalue2 = cos( -( fltlat + 1 ) * th );
	double svalue2 = sin( -( fltlat + 1 ) * th );

	double re = 0.0;
	double im = 0.0;

	while( fltlen > 0 )
	{
		re += cvalue1 * flt[ 0 ];
		im += svalue1 * flt[ 0 ];
		flt++;
		fltlen--;

		double tmp = cvalue1;
		cvalue1 = sincr * cvalue1 - cvalue2;
		cvalue2 = tmp;

		tmp = svalue1;
		svalue1 = sincr * svalue1 - svalue2;
		svalue2 = tmp;
	}

	re0 = re;
	im0 = im;
}

/**
 * Function normalizes FIR filter so that its frequency response at DC is
 * equal to DCGain.
 *
 * @param[in,out] p Filter coefficients.
 * @param l Filter length.
 * @param DCGain Filter's gain at DC.
 * @param pstep "p" array step.
 */

template< class T >
inline void normalizeFIRFilter( T* const p, const int l, const double DCGain,
	const int pstep = 1 )
{
	double s = 0.0;
	T* pp = p;
	int i = l;

	while( i > 0 )
	{
		s += *pp;
		pp += pstep;
		i--;
	}

	s = DCGain / s;
	pp = p;
	i = l;

	while( i > 0 )
	{
		*pp = (T) ( *pp * s );
		pp += pstep;
		i--;
	}
}

/**
 * @brief Memory buffer class for element array storage, with capacity
 * tracking.
 *
 * Allows easier handling of memory blocks allocation and automatic
 * deallocation for arrays (buffers) consisting of elements of specified
 * class. Tracks buffer's capacity in "int" variable; unsuitable for
 * allocation of very large memory blocks (with more than 2 billion elements).
 *
 * This class manages memory space only - it does not perform element class
 * construction (initialization) operations. Buffer's required memory address
 * alignment specification is supported.
 *
 * Uses standard library to allocate and deallocate memory.
 *
 * @tparam T Buffer element's type.
 */

template< class T >
class CBuffer
{
public:
	CBuffer()
		: Data( NULL )
		, DataAligned( NULL )
		, Capacity( 0 )
		, Alignment( 0 )
	{
	}

	/**
	 * Constructor creates the buffer with the specified capacity.
	 *
	 * @param aCapacity Buffer's capacity.
	 * @param aAlignment Buffer's required memory address alignment. 0 - use
	 * stdlib's default alignment.
	 */

	CBuffer( const int aCapacity, const int aAlignment = 0 )
	{
		allocinit( aCapacity, aAlignment );
	}

	CBuffer( const CBuffer& Source )
	{
		allocinit( Source.Capacity, Source.Alignment );
		memcpy( DataAligned, Source.DataAligned, Capacity * sizeof( T ));
	}

	~CBuffer()
	{
		freeData();
	}

	CBuffer& operator = ( const CBuffer& Source )
	{
		alloc( Source.Capacity, Source.Alignment );
		memcpy( DataAligned, Source.DataAligned, Capacity * sizeof( T ));
		return( *this );
	}

	/**
	 * Function allocates memory so that the specified number of elements
	 * can be stored in *this buffer object.
	 *
	 * @param aCapacity Storage for this number of elements to allocate.
	 * @param aAlignment Buffer's required memory address alignment,
	 * power-of-2 values only. 0 - use stdlib's default alignment.
	 */

	void alloc( const int aCapacity, const int aAlignment = 0 )
	{
		freeData();
		allocinit( aCapacity, aAlignment );
	}

	/**
	 * Function deallocates any previously allocated buffer.
	 */

	void free()
	{
		freeData();
		Data = NULL;
		Capacity = 0;
		Alignment = 0;
	}

	/**
	 * @return The capacity of the element buffer.
	 */

	int getCapacity() const
	{
		return( Capacity );
	}

	/**
	 * Function reallocates *this buffer to a larger size so that it will be
	 * able to hold the specified number of elements. Downsizing is not
	 * performed. Alignment is not changed.
	 *
	 * @param NewCapacity New (increased) capacity.
	 */

	void increaseCapacity( const int NewCapacity )
	{
		if( NewCapacity < Capacity )
		{
			return;
		}

		const int PrevCapacity = Capacity;
		T* const PrevData = Data;
		T* const PrevDataAligned = DataAligned;

		allocinit( NewCapacity, Alignment );
		memcpy( DataAligned, PrevDataAligned, PrevCapacity * sizeof( T ));

		:: free( PrevData );
	}

	/**
	 * Function "truncates" (reduces) capacity of the buffer without
	 * reallocating it. Alignment is not changed.
	 *
	 * @param NewCapacity New required capacity.
	 */

	void truncateCapacity( const int NewCapacity )
	{
		if( NewCapacity >= Capacity )
		{
			return;
		}

		Capacity = NewCapacity;
	}

	/**
	 * Function increases capacity so that the specified number of
	 * elements can be stored. This function increases the previous capacity
	 * value by third the current capacity value until space for the required
	 * number of elements is available. Alignment is not changed.
	 *
	 * @param ReqCapacity Required capacity.
	 */

	void updateCapacity( const int ReqCapacity )
	{
		if( ReqCapacity <= Capacity )
		{
			return;
		}

		int NewCapacity = Capacity;

		while( NewCapacity < ReqCapacity )
		{
			NewCapacity += NewCapacity / 3 + 1;
		}

		increaseCapacity( NewCapacity );
	}

	operator T* () const
	{
		return( DataAligned );
	}

private:
	T* Data; ///< Element buffer pointer.
		///<
	T* DataAligned; ///< Memory address-aligned element buffer pointer.
		///<
	int Capacity; ///< Element buffer capacity.
		///<
	int Alignment; ///< Memory address alignment in use. 0 - use stdlib's
		///< default alignment.
		///<

	/**
	 * Internal element buffer allocation function used during object
	 * construction.
	 *
	 * @param aCapacity Storage for this number of elements to allocate.
	 * @param aAlignment Buffer's required memory address alignment. 0 - use
	 * stdlib's default alignment.
	 */

	void allocinit( const int aCapacity, const int aAlignment )
	{
		if( aAlignment == 0 )
		{
			Data = (T*) :: malloc( aCapacity * sizeof( T ));
			DataAligned = Data;
			Alignment = 0;
		}
		else
		{
			Data = (T*) :: malloc( aCapacity * sizeof( T ) + aAlignment );
			DataAligned = alignptr( Data, aAlignment );
			Alignment = aAlignment;
		}

		Capacity = aCapacity;
	}

	/**
	 * Function frees a previously allocated Data buffer.
	 */

	void freeData()
	{
		:: free( Data );
	}

	/**
	 * Function modifies the specified pointer so that it becomes memory
	 * address-aligned.
	 *
	 * @param ptr Pointer to align.
	 * @param align Alignment in bytes to apply.
	 * @return Pointer aligned to align bytes. Works with power-of-2
	 * alignments only. If no alignment is necessary, "align" bytes will be
	 * added to the pointer value.
	 */

	template< class Tp >
	inline Tp alignptr( const Tp ptr, const uintptr_t align )
	{
		return( (Tp) ( (uintptr_t) ptr + align -
			( (uintptr_t) ptr & ( align - 1 ))) );
	}
};

/**
 * Function optimizes the length of the symmetric FIR filter by removing left-
 * and right-most elements that are below specific threshold.
 *
 * @param[in,out] Flt Buffer that contains filter being optimized.
 * @param[in,out] FltLatency Variable that holds the current latency of the
 * filter. May be adjusted on function return.
 * @param Threshold Threshold level.
 */

template< class T >
inline void optimizeFIRFilter( CBuffer< T >& Flt, int& FltLatency,
	T const Threshold = (T) 0.00001 )
{
	int i;

	for( i = 0; i <= FltLatency; i++ )
	{
		if( fabs( Flt[ i ]) >= Threshold || i == FltLatency )
		{
			if( i > 0 )
			{
				const int NewCapacity = Flt.getCapacity() - i * 2;
				copyArray( &Flt[ i ], &Flt[ 0 ], NewCapacity );
				Flt.truncateCapacity( NewCapacity );
				FltLatency -= i;
			}

			break;
		}
	}
}

/**
 * @brief Array of structured objects.
 *
 * Implements allocation of a linear array of objects of class T (which are
 * initialized), addressable via operator[]. Each object is created via the
 * "operator new". New object insertions are quick since implementation uses
 * prior space allocation (capacity), thus not requiring frequent memory block
 * reallocations.
 *
 * @tparam T Array element's type.
 */

template< class T >
class CStructArray
{
public:
	CStructArray()
		: ItemCount( 0 )
	{
	}

	CStructArray( const CStructArray& Source )
		: ItemCount( 0 )
		, Items( Source.getItemCount() )
	{
		while( ItemCount < Source.getItemCount() )
		{
			Items[ ItemCount ] = new T( Source[ ItemCount ]);
			ItemCount++;
		}
	}

	~CStructArray()
	{
		clear();
	}

	CStructArray& operator = ( const CStructArray& Source )
	{
		clear();

		const int NewCount = Source.ItemCount;
		Items.updateCapacity( NewCount );

		while( ItemCount < NewCount )
		{
			Items[ ItemCount ] = new T( Source[ ItemCount ]);
			ItemCount++;
		}

		return( *this );
	}

	T& operator []( const int Index )
	{
		return( *Items[ Index ]);
	}

	const T& operator []( const int Index ) const
	{
		return( *Items[ Index ]);
	}

	/**
	 * Function creates a new object of type T with the default constructor
	 * and adds this object to the array.
	 *
	 * @return Reference to a newly added object.
	 */

	T& add()
	{
		if( ItemCount == Items.getCapacity() )
		{
			Items.increaseCapacity( ItemCount * 3 / 2 + 1 );
		}

		Items[ ItemCount ] = new T();
		ItemCount++;

		return( (*this)[ ItemCount - 1 ]);
	}

	/**
	 * Function erases all items of *this array.
	 */

	void clear()
	{
		while( ItemCount > 0 )
		{
			ItemCount--;
			delete Items[ ItemCount ];
		}
	}

	/**
	 * @return The number of allocated items.
	 */

	int getItemCount() const
	{
		return( ItemCount );
	}

private:
	int ItemCount; ///< The number of items available in the array.
		///<
	CBuffer< T* > Items; ///< Element buffer.
		///<
};

/**
 * @brief Sine signal generator class.
 *
 * Class implements sine signal generator without biasing, with
 * constructor-based initalization only.
 */

class CSineGen
{
public:
	/**
	 * Constructor initializes *this sine signal generator.
	 *
	 * @param si Sine function increment, in radians.
	 * @param ph Starting phase, in radians. Add 0.5 * M_PI for cosine
	 * function.
	 */

	CSineGen( const double si, const double ph )
		: svalue1( sin( ph ))
		, svalue2( sin( ph - si ))
		, sincr( 2.0 * cos( si ))
	{
	}

	/**
	 * @return The next value of the sine function, without biasing.
	 */

	double generate()
	{
		const double res = svalue1;

		svalue1 = sincr * res - svalue2;
		svalue2 = res;

		return( res );
	}

private:
	double svalue1; ///< Current sine value.
		///<
	double svalue2; ///< Previous sine value.
		///<
	double sincr; ///< Sine value increment.
		///<
};

/**
 * @brief Peaked Cosine window function generator class.
 *
 * Class implements Peaked Cosine window function generator. Generates the
 * right-handed half of the window function. The Alpha parameter of this
 * window function offers the control of the balance between the early and
 * later taps of the filter. E.g. at Alpha=1 both early and later taps are
 * attenuated, but at Alpha=4 mostly later taps are attenuated. This offers a
 * great control over ringing artifacts produced by a low-pass filter in image
 * processing, without compromising achieved image sharpness.
 */

class CDSPWindowGenPeakedCosine
{
public:
	/**
	 * Constructor initializes *this window function generator.
	 *
	 * @param aAlpha Alpha parameter, affects the peak shape (peak
	 * augmentation) of the window function. Should be >= 1.0.
	 * @param aLen2 Half filter kernel's length (non-truncated).
	 */

	CDSPWindowGenPeakedCosine( const double aAlpha, const double aLen2 )
		: Alpha( aAlpha )
		, Len2( aLen2 )
		, wn( 0 )
		, w1( M_PId2 / Len2, M_PI * 0.5 )
	{
	}

	/**
	 * @return The next Peaked Cosine window function coefficient.
	 */

	double generate()
	{
		const double h = pow( wn / Len2, Alpha );
		wn++;

		return( w1.generate() * ( 1.0 - h ));
	}

private:
	double Alpha; ///< Alpha parameter, affects the peak shape of window.
		///<
	double Len2; ///< Half length of the window function.
		///<
	int wn; ///< Window function integer position. 0 - center of the
		///< window function.
		///<
	CSineGen w1; ///< Sine-wave generator.
		///<
};

/**
 * @brief FIR filter-based equalizer generator.
 *
 * Class implements an object used to generate FIR filter kernels with the
 * specified frequency response (aka paragraphic equalizer). The calculated
 * filter kernel is windowed by the Peaked Cosine window function.
 *
 * In image processing, due to short length of filters being used (6-8 taps)
 * the resulting frequency response of the filter is approximate and may be
 * mathematically imperfect, but still adequate to the visual requirements.
 *
 * On a side note, this equalizer generator can be successfully used for audio
 * signal equalization as well: for example, it is used in almost the same
 * form in Voxengo Marvel GEQ equalizer plug-in.
 */

class CDSPFIREQ
{
public:
	/**
	 * Function initializes *this object with the required parameters. The
	 * gain of frequencies beyond the MinFreq..MaxFreq range are controlled by
	 * the first and the last band's gain.
	 *
	 * @param SampleRate Processing sample rate (use 2 for image processing).
	 * @param aKernelLength Required kernel length in samples (taps). The
	 * actual kernel length is truncated to an integer value.
	 * @param aBandCount Number of band crossover points required to control,
	 * including bands at MinFreq and MaxFreq.
	 * @param MinFreq Minimal frequency that should be controlled.
	 * @param MaxFreq Maximal frequency that should be controlled.
	 * @param IsLogBands "True" if the bands should be spaced logarithmically.
	 * @param WFAlpha Peaked Cosine window function's Alpha parameter.
	 */

	void init( const double SampleRate, const double aKernelLength,
		const int aBandCount, const double MinFreq, const double MaxFreq,
		const bool IsLogBands, const double WFAlpha )
	{
		KernelLength = aKernelLength;
		BandCount = aBandCount;

		CenterFreqs.alloc( BandCount );

		z = (int) ceil( KernelLength * 0.5 );
		zi = z + ( z & 1 );
		z2 = z * 2;

		CBuffer< double > oscbuf( z2 );
		initOscBuf( oscbuf );

		CBuffer< double > winbuf( z );
		initWinBuf( winbuf, WFAlpha );

		UseFirstVirtBand = ( MinFreq > 0.0 );
		const int k = zi * ( BandCount + ( UseFirstVirtBand ? 1 : 0 ));
		Kernels1.alloc( k );
		Kernels2.alloc( k );

		double m;
		double mo;

		if( IsLogBands )
		{
			m = exp( log( MaxFreq / MinFreq ) / ( BandCount - 1 ));
			mo = 0.0;
		}
		else
		{
			m = 1.0;
			mo = ( MaxFreq - MinFreq ) / ( BandCount - 1 );
		}

		double f = MinFreq;
		double x1 = 0.0;
		double x2;
		int si;

		if( UseFirstVirtBand )
		{
			si = 0;
		}
		else
		{
			si = 1;
			CenterFreqs[ 0 ] = 0.0;
			f = f * m + mo;
		}

		double* kernbuf1 = &Kernels1[ 0 ];
		double* kernbuf2 = &Kernels2[ 0 ];
		int i;

		for( i = si; i < BandCount; i++ )
		{
			x2 = f * 2.0 / SampleRate;
			CenterFreqs[ i ] = x2;

			fillKernelBand( x1, x2, kernbuf1, kernbuf2, oscbuf, winbuf );

			kernbuf1 += zi;
			kernbuf2 += zi;
			x1 = x2;
			f = f * m + mo;
		}

		if( x1 < 1.0 )
		{
			UseLastVirtBand = true;
			fillKernelBand( x1, 1.0, kernbuf1, kernbuf2, oscbuf, winbuf );
		}
		else
		{
			UseLastVirtBand = false;
		}
	}

	/**
	 * @return Kernel's length, in samples (taps).
	 */

	int getKernelLength() const
	{
		return( z2 - 1 );
	}

	/**
	 * @return Kernel's output latency, in samples (taps).
	 */

	int getKernelLatency() const
	{
		return( z - 1 );
	}

	/**
	 * Function creates FIR kernel with the specified gain levels at band
	 * crossover points.
	 *
	 * @param BandGains Array of linear gain levels, count=BandCount specified
	 * in the init() function.
	 * @param[out] Kernel Output kernel buffer, length = getKernelLength().
	 */

	void buildKernel( const double* const BandGains, double* const Kernel )
	{
		const double* kernbuf1 = &Kernels1[ 0 ];
		const double* kernbuf2 = &Kernels2[ 0 ];
		double x1 = 0.0;
		double y1 = BandGains[ 0 ];
		double x2;
		double y2;

		int i;
		int si;

		if( UseFirstVirtBand )
		{
			si = 1;
			x2 = CenterFreqs[ 0 ];
			y2 = y1;
		}
		else
		{
			si = 2;
			x2 = CenterFreqs[ 1 ];
			y2 = BandGains[ 1 ];
		}

		copyKernelBand( Kernel, kernbuf1, kernbuf2, y1 - y2,
			x1 * y2 - x2 * y1 );

		kernbuf1 += zi;
		kernbuf2 += zi;
		x1 = x2;
		y1 = y2;

		for( i = si; i < BandCount; i++ )
		{
			x2 = CenterFreqs[ i ];
			y2 = BandGains[ i ];

			addKernelBand( Kernel, kernbuf1, kernbuf2, y1 - y2,
				x1 * y2 - x2 * y1 );

			kernbuf1 += zi;
			kernbuf2 += zi;
			x1 = x2;
			y1 = y2;
		}

		if( UseLastVirtBand )
		{
			addKernelBand( Kernel, kernbuf1, kernbuf2, y1 - y2,
				x1 * y2 - y1 );
		}

		for( i = 0; i < z - 1; i++ )
		{
			Kernel[ z + i ] = Kernel[ z - 2 - i ];
		}
	}

private:
	double KernelLength; ///< Length of kernel.
		///<
	int z; ///< Equals (int) floor( KernelLength * 0.5 ).
		///<
	int zi; ///< Equals "z" if z is even, or z + 1 if z is odd. Used as a
		///< Kernels1 and Kernels2 size multiplier and kernel buffer increment
		///< to make sure each kernel buffer is 16-byte aligned.
		///<
	int z2; ///< Equals z * 2.
		///<
	int BandCount; ///< Number of controllable bands.
		///<
	CBuffer< double > CenterFreqs; ///< Center frequencies for all bands,
		///< normalized to 0.0-1.0 range.
		///<
	CBuffer< double > Kernels1; ///< Half-length kernel buffers for each
		///< spectral band (linear part).
		///<
	CBuffer< double > Kernels2; ///< Half-length kernel buffers for each
		///< spectral band (ramp part).
		///<
	bool UseFirstVirtBand; ///< "True" if the first virtual band
		///< (between 0.0 and MinFreq) should be used. The first virtual band
		///< won't be used if MinFreq equals 0.0.
		///<
	bool UseLastVirtBand; ///< "True" if the last virtual band (between
		///< MaxFreq and SampleRate * 0.5) should be used. The last virtual
		///< band won't be used if MaxFreq * 2.0 equals SampleRate.
		///<

	/**
	 * Function initializes the "oscbuf" used in the fillKernelBand function.
	 *
	 * @param oscbuf Oscillator buffer, length = z * 2.
	 */

	void initOscBuf( double* oscbuf ) const
	{
		int i = z;

		while( i > 0 )
		{
			oscbuf[ 0 ] = 0.0;
			oscbuf[ 1 ] = 1.0;
			oscbuf += 2;
			i--;
		}
	}

	/**
	 * Function initializes window function buffer. This function generates
	 * Peaked Cosine window function.
	 *
	 * @param winbuf Windowing buffer.
	 * @param Alpha Peaked Cosine alpha parameter.
	 */

	void initWinBuf( double* winbuf, const double Alpha ) const
	{
		CDSPWindowGenPeakedCosine wf( Alpha, KernelLength * 0.5 );
		int i;

		for( i = 1; i <= z; i++ )
		{
			winbuf[ z - i ] = wf.generate();
		}
	}

	/**
	 * Function fills first half of symmetric FIR kernel for the band. This
	 * function should be called successively for adjacent bands. Previous
	 * band's x2 should be equal to current band's x1. A band kernel consists
	 * of 2 elements: linear kernel and ramp kernel.
	 *
	 * @param x1 Band's left corner frequency (0..1).
	 * @param x2 Band's right corner frequency (0..1).
	 * @param kernbuf1 Kernel band buffer 1 (linear part), length = z.
	 * @param kernbuf2 Kernel band buffer 2 (ramp part), length = z.
	 * @param oscbuf Oscillation buffer. Before the first call of the
	 * fillKernelBand() should be initialized with the call of the
	 * initOscBuf() function.
	 * @param winbuf Buffer that contains windowing function.
	 */

	void fillKernelBand( const double x1, const double x2, double* kernbuf1,
		double* kernbuf2, double* oscbuf, const double* const winbuf )
	{
		const double s2_incr = M_PI * x2;
		const double s2_coeff = 2.0 * cos( s2_incr );

		double s2_value1 = sin( s2_incr * ( -z + 1 ));
		double c2_value1 = sin( s2_incr * ( -z + 1 ) + M_PI * 0.5 );
		oscbuf[ 0 ] = sin( s2_incr * -z );
		oscbuf[ 1 ] = sin( s2_incr * -z + M_PI * 0.5 );

		int ks;

		for( ks = 1; ks < z; ks++ )
		{
			const int ks2 = ks * 2;
			const double s1_value1 = oscbuf[ ks2 ];
			const double c1_value1 = oscbuf[ ks2 + 1 ];
			oscbuf[ ks2 ] = s2_value1;
			oscbuf[ ks2 + 1 ] = c2_value1;

			const double x = M_PI * ( ks - z );
			const double v0 = winbuf[ ks - 1 ] / (( x1 - x2 ) * x );

			kernbuf1[ ks - 1 ] = ( x2 * s2_value1 - x1 * s1_value1 +
				( c2_value1 - c1_value1 ) / x ) * v0;

			kernbuf2[ ks - 1 ] = ( s2_value1 - s1_value1 ) * v0;

			s2_value1 = s2_coeff * s2_value1 - oscbuf[ ks2 - 2 ];
			c2_value1 = s2_coeff * c2_value1 - oscbuf[ ks2 - 1 ];
		}

		kernbuf1[ z - 1 ] = ( x2 * x2 - x1 * x1 ) / ( x1 - x2 ) * 0.5;
		kernbuf2[ z - 1 ] = -1.0;
	}

	/**
	 * Function copies kernel band's elements to the output buffer.
	 *
	 * @param outbuf Output buffer.
	 * @param kernbuf1 Kernel buffer 1 (linear part).
	 * @param kernbuf2 Kernel buffer 2 (ramp part).
	 * @param c Multiplier for linear kernel element.
	 * @param d Multiplier for ramp kernel element.
	 */

	void copyKernelBand( double* outbuf, const double* const kernbuf1,
		const double* const kernbuf2, const double c, const double d ) const
	{
		int ks;

		for( ks = 0; ks < z; ks++ )
		{
			outbuf[ ks ] = c * kernbuf1[ ks ] + d * kernbuf2[ ks ];
		}
	}

	/**
	 * Function adds kernel band's elements to the output buffer.
	 *
	 * @param outbuf Output buffer.
	 * @param kernbuf1 Kernel buffer 1 (linear part).
	 * @param kernbuf2 Kernel buffer 2 (ramp part).
	 * @param c Multiplier for linear kernel element.
	 * @param d Multiplier for ramp kernel element.
	 */

	void addKernelBand( double* outbuf, const double* const kernbuf1,
		const double* const kernbuf2, const double c, const double d ) const
	{
		int ks;

		for( ks = 0; ks < z; ks++ )
		{
			outbuf[ ks ] += c * kernbuf1[ ks ] + d * kernbuf2[ ks ];
		}
	}
};

/**
 * @brief Low-pass filter windowed by Peaked Cosine window function.
 *
 * This class implements calculation of linear-phase FIR low-pass filter
 * kernel windowed by the Peaked Cosine window function, for image processing
 * applications.
 */

class CDSPPeakedCosineLPF
{
public:
	int fl2; ///< Half filter's length, excluding the peak value. This value
		///< can be also used as filter's latency in samples (taps).
		///<
	int KernelLen; ///< Filter kernel's length in samples (taps).
		///<

	/**
	 * Constructor initalizes *this object.
	 *
	 * @param aLen2 Half-length (non-truncated) of low-pass filter's kernel,
	 * in samples (taps).
	 * @param aFreq2 Low-pass filter's corner frequency [0; pi].
	 * @param aAlpha Peaked Cosine window function Alpha parameter.
	 */

	CDSPPeakedCosineLPF( const double aLen2, const double aFreq2,
		const double aAlpha )
		: fl2( (int) ceil( aLen2 ) - 1 )
		, KernelLen( fl2 + fl2 + 1 )
		, Len2( aLen2 )
		, Freq2( aFreq2 )
		, Alpha( aAlpha )
	{
	}

	/**
	 * Function generates a linear-phase low-pass filter kernel windowed by
	 * Peaked Cosine window function.
	 *
	 * @param[out] op Output buffer, length = KernelLen (fl2 * 2 + 1).
	 * @param DCGain Required gain at DC. The resulting filter will be
	 * normalized to achieve this DC gain.
	 */

	template< class T >
	void generateLPF( T* op, const double DCGain )
	{
		CDSPWindowGenPeakedCosine wf( Alpha, Len2 );
		CSineGen f2( Freq2, 0.0 );

		op += fl2;
		T* op2 = op;
		f2.generate();
		int t = 1;

		*op = (T) ( Freq2 * wf.generate() / M_PI );
		double s = *op;

		while( t <= fl2 )
		{
			const double v = f2.generate() * wf.generate() / t / M_PI;
			op++;
			op2--;
			*op = (T) v;
			*op2 = (T) v;
			s += *op + *op2;
			t++;
		}

		t = KernelLen;
		s = DCGain / s;

		while( t > 0 )
		{
			*op2 = (T) ( *op2 * s );
			op2++;
			t--;
		}
	}

private:
	double Len2; ///< Half-length (non-truncated) of low-pass filter, in
		///< samples (taps).
		///<
	double Freq2; ///< Low-pass filter's corner frequency.
		///<
	double Alpha; ///< Peaked Cosine window function Alpha parameter.
		///<
};

/**
 * @brief Sinc function-based fractional delay filter bank.
 *
 * Class implements storage and initialization of a bank of sinc-based
 * fractional delay filters, expressed as 1st order polynomial interpolation
 * coefficients. The filters are produced from a single "long" windowed
 * low-pass filter kernel. Also supports zero-order bank ("nearest neighbor"
 * interpolation).
 *
 * @tparam fptype Specifies storage type of the filter coefficients bank. The
 * filters are initially calculated using the "double" precision.
 */

template< class fptype >
class CDSPFracFilterBankLin
{
public:
	/**
	 * Function initializes (builds) the filter bank based on the supplied
	 * parameters.
	 *
	 * @param ReqFracCount Required number of fractional delays in the filter
	 * bank. The minimal value is 2.
	 * @param ReqOrder Required order of the interpolation polynomial
	 * (0 or 1).
	 * @param BaseLen Low-pass filter's base length, in samples (taps).
	 * Affects the actual length of the filter and its overall steepness.
	 * @param Cutoff Low-pass filter's normalized cutoff frequency [0; 1].
	 * @param WFAlpha Peaked Cosine window function's Alpha parameter.
	 * @param Alignment Memory alignment of the filter bank, power-of-2 value.
	 * 0 - use default stdlib alignment.
	 * @param FltLenAlign Filter's length alignment, power-of-2 value.
	 */

	void init( const int ReqFracCount, const int ReqOrder,
		const double BaseLen, const double Cutoff, const double WFAlpha,
		const int Alignment = 0, const int FltLenAlign = 1 )
	{
		CDSPPeakedCosineLPF p( 0.5 * BaseLen * ReqFracCount,
			M_PI * Cutoff / ReqFracCount, WFAlpha );

		Order = ReqOrder;
		const int ElementSize = ReqOrder + 1;
		FilterLen = ( p.fl2 / ReqFracCount + 1 ) * 2;
		FilterLen = ( FilterLen + FltLenAlign - 1 ) & ~( FltLenAlign - 1 );
		FilterSize = FilterLen * ElementSize;

		const int BufLen = FilterLen * ReqFracCount + InterpPoints - 1;
		const int BufOffs = InterpPoints / 2 - 1;
		const int BufCenter = FilterLen * ReqFracCount / 2 + BufOffs;

		CBuffer< double > Buf( BufLen );
		memset( Buf, 0, ( BufCenter - p.fl2 ) * sizeof( double ));
		int i = BufLen - BufCenter - p.fl2 - 1;
		memset( &Buf[ BufLen - i ], 0, i * sizeof( double ));

		p.generateLPF( &Buf[ BufCenter - p.fl2 ], ReqFracCount );

		Table.alloc( ReqFracCount * FilterSize, Alignment );
		fptype* op = Table;
		int j;

		if( ElementSize == 2 )
		{
			for( i = ReqFracCount; i > 0; i-- )
			{
				double* p = Buf + BufOffs + i;

				for( j = 0; j < FilterLen; j++ )
				{
					op[ 0 ] = (fptype) p[ 0 ];
					op[ FilterLen ] = (fptype) ( p[ -1 ] - p[ 0 ]);
					op++;
					p += ReqFracCount;
				}

				op += FilterLen;
			}
		}
		else
		{
			for( i = ReqFracCount; i > 0; i-- )
			{
				double* p = Buf + BufOffs + i;

				for( j = 0; j < FilterLen; j++ )
				{
					*op = (fptype) *p;
					op++;
					p += ReqFracCount;
				}
			}
		}
	}

	/**
	 * Operator returns the pointer to the specified interpolation table
	 * filter.
	 *
	 * @param i Filter (fractional delay) index, in the range 0 to
	 * ReqFracCount - 1, inclusive.
	 * @return Reference to the filter elements. Higher order polynomial
	 * coefficients are stored after after previous order coefficients.
	 */

	const fptype& operator []( const int i ) const
	{
		return( Table[ i * FilterSize ]);
	}

	/**
	 * @return The length of each fractional delay filter, in samples (taps).
	 * Always an even value.
	 */

	int getFilterLen() const
	{
		return( FilterLen );
	}

	/**
	 * @return The order of the interpolation polynomial.
	 */

	int getOrder() const
	{
		return( Order );
	}

private:
	static const int InterpPoints = 2; ///< The maximal number of points the
		///< interpolation is based on.
		///<
	int Order; ///< The order of the interpolation polynomial.
		///<
	int FilterLen; ///< Specifies the number of samples (taps) each fractional
		///< delay filter has. This is an even value, depends on the
		///< ReqNormFreq and ReqAtten parameters.
		///<
	int FilterSize; ///< The size of a single filter element, equals
		///< FilterLen * ElementSize.
		///<
	CBuffer< fptype > Table; ///< Interpolation table, size equals to
		///< ReqFracCount * FilterLen * ElementSize.
		///<
};

/**
 * @brief Thread pool for multi-threaded image resizing operation.
 *
 * This base class is used to organize a multi-threaded image resizing
 * operation. The thread pool should consist of threads that initially wait
 * for a signal. Upon receiving a signal (via the startAllWorkloads()
 * function) each previously added thread should execute its workload's
 * process() function once, and return to the wait signal state again. The
 * thread pool should be also able to efficiently wait for all workloads to
 * finish via the waitAllWorkloadsToFinish() function.
 *
 * The image resizing algorithm makes calls to functions of this class.
 */

class CImageResizerThreadPool
{
public:
	CImageResizerThreadPool()
	{
	}

	virtual ~CImageResizerThreadPool()
	{
	}

	/**
	 * @brief Thread pool's workload object.
	 *
	 * This class should be used as a base class for objects that perform the
	 * actual work spread over several threads.
	 */

	class CWorkload
	{
	public:
		virtual ~CWorkload()
		{
		}

		/**
		 * Function that gets called from the thread when thread pool's
		 * startAllWorkloads() function is called.
		 */

		virtual void process() = 0;
	};

	/**
	 * @return The suggested number of workloads (and their associated
	 * threads) to add. The minimal value this function can return is 1. The
	 * usual value may depend on the number of physical and virtual cores
	 * present in the system, and on other considerations.
	 */

	virtual int getSuggestedWorkloadCount() const
	{
		return( 1 );
	}

	/**
	 * Function adds a new workload (and possibly thread) to the thread pool.
	 * The caller decides how many parallel workloads (and threads) it
	 * requires, but this number will not exceed the value returned by the
	 * getSuggestedWorkloadCount() function. It is implementation-specific how
	 * many workloads to associate with a single thread. But for efficiency
	 * reasons each workload should be associated with its own thread.
	 *
	 * Note that the same set of workload objects will be processed each time
	 * the startAllWorkloads() function is called. This means that workload
	 * objects are added only once. The caller changes the state of the
	 * workload objects and then calls the startAllWorkloads() function to
	 * process them.
	 *
	 * @param Workload Workload object whose process() function will be called
	 * from within the thread when the startAllWorkloads() function is called.
	 */

	virtual void addWorkload( CWorkload* const Workload )
	{
	}

	/**
	 * Function starts all workloads associated with threads previously added
	 * via the addWorkload() function. It is assumed that this function
	 * performs the necessary "memory barrier" (or "cache sync") kind of
	 * operation so that all threads catch up the prior changes made to the
	 * workload objects.
	 */

	virtual void startAllWorkloads()
	{
	}

	/**
	 * Function waits for all workloads to finish.
	 */

	virtual void waitAllWorkloadsToFinish()
	{
	}

	/**
	 * Function removes all workloads previously added via the addWorkload()
	 * function. This function gets called only after the
	 * waitAllWorkloadsToFinish() function call.
	 */

	virtual void removeAllWorkloads()
	{
	}
};

/**
 * @brief Resizing algorithm parameters structure.
 *
 * This structure holds all selectable parameters used by the resizing
 * algorithm on various stages, for both downsizing and upsizing. There are no
 * other parameters exist that can optimize the performance of the resizing
 * algorithm. Filter length parameters can take fractional values.
 *
 * Beside quality, these parameters (except Alpha parameters) directly affect
 * the computative cost of the resizing algorithm. It is possible to trade
 * the visual quality for computative cost.
 *
 * Anti-alias filtering during downsizing can be defined as a considerable
 * reduction of contrast of smallest features of an image. Unfortunately, such
 * de-contrasting partially affects features of all sizes thus producing a
 * non-linearity of frequency response. All pre-defined parameter sets are
 * described by 3 values separated by slashes. The first value is the
 * de-contrasting factor of small features (which are being removed) while
 * the second value is the de-contrasting factor of large features (which
 * should remain intact), with value of 1 equating to "no contrast change".
 * The third value is the optimization score (see below), with value of 0
 * equating to the "perfect" linearity of frequency response.
 *
 * The pre-defined parameter sets offered by this library were auto-optimized
 * for the given LPFltBaseLen, IntFltLen and CorrFltAlpha values. The
 * optimization goal was to minimize the score: the sum of squares of the
 * difference between original and processed images (which was not actually
 * resized, k=1). The original image was a 0.5 megapixel uniformly-distributed
 * white-noise image with pixel intensities in the 0-1 range. Such goal
 * converges very well and produces filtering system with the flattest
 * frequency response possible for the given constraints. With this goal,
 * increasing the LPFltBaseLen value reduces the general amount of aliasing
 * artifacts.
 */

struct CImageResizerParams
{
	double CorrFltAlpha; ///< Alpha parameter of the Peaked Cosine window
		///< function used on the correction filter. The "usable" values are
		///< in the narrow range 1.0 to 1.5.
		///<
	double CorrFltLen; ///< Correction filter's length in samples (taps). The
		///< "usable" range is narrow, 5.5 to 8, as to minimize the
		///< "overcorrection" which is mathematically precise, but visually
		///< unacceptable.
		///<
	double IntFltAlpha; ///< Alpha parameter of the Peaked Cosine window
		///< function used on the interpolation low-pass filter. The "usable"
		///< values are in the range 1.5 to 2.5.
		///<
	double IntFltCutoff; ///< Interpolation low-pass filter's cutoff frequency
		///< (normalized, [0; 1]). The "usable" range is 0.6 to 0.8.
		///<
	double IntFltLen; ///< Interpolation low-pass filter's length in samples
		///< (taps). The length value should be at least 18 or otherwise a
		///< "dark grid" artifact will be introduced if a further sharpening
		///< is applied. IntFltLen together with other IntFlt parameters
		///< should be tuned in a way that produces the flattest frequency
		///< response in 0-0.5 normalized frequency range (this range is due
		///< to 2X upsampling).
		///<
	double LPFltAlpha; ///< Alpha parameter of the Peaked Cosine window
		///< function used on the low-pass filter. The "usable" values are
		///< in the range 1.5 to 6.5.
		///<
	double LPFltBaseLen; ///< Base length of the low-pass (aka anti-aliasing
		///< or reconstruction) filter, in samples (taps), further adjusted by
		///< the actual cutoff frequency, upsampling and downsampling factors.
		///< The "usable" range is between 6 and 9.
		///<
	double LPFltCutoffMult; ///< Low-pass filter's cutoff frequency
		///< multiplier. This value can be both below and above 1.0 as
		///< low-pass filters are inserted on downsampling and upsampling
		///< steps and always have corner frequency equal to or below 0.5pi.
		///< This multiplier shifts low-pass filter's corner frequency towards
		///< lower (if below 1.0) or higher (if above 1.0) frequencies. This
		///< multiplier can be way below 1.0 since any additional
		///< high-frequency damping will be partially corrected by the
		///< correction filter. The "usable" range is 0.3 to 1.0.
		///<

	CImageResizerParams()
		: HBFltAlpha( 1.95296 )
		, HBFltCutoff( 0.51776 )
		, HBFltLen( 20.00000 )
	{
	}

	double HBFltAlpha; ///< Half-band filter's Alpha. Assigned internally.
		///<
	double HBFltCutoff; ///< Half-band filter's cutoff point [0; 1]. Assigned
		///< internally.
		///<
	double HBFltLen; ///< Length of the half-band low-pass filter. Assigned
		///< internally. Internally used to perform 2X downsampling when
		///< downsizing is considerable (more than 8 times). These filter
		///< parameters should be treated as "technical" and do not require
		///< adjustment as they were tuned to suit all combinations of other
		///< parameters. This half-band filter provides a wide transition
		///< band (for minimal ringing artifacts) and a high stop-band
		///< attenuation (for minimal aliasing).
		///<
};

/**
 * @brief The default set of resizing algorithm parameters
 * (10.01/1.029/0.019169).
 *
 * This is the default set of resizing parameters that was designed to deliver
 * a sharp image while still providing a low amount of ringing artifacts, and
 * having a reasonable computational cost.
 */

struct CImageResizerParamsDef : public CImageResizerParams
{
	CImageResizerParamsDef()
	{
		CorrFltAlpha = 1.0;//10.01/1.88/1.029(522.43)/0.019169:258648,446808
		CorrFltLen = 6.30770;
		IntFltAlpha = 2.27825;
		IntFltCutoff = 0.75493;
		IntFltLen = 18.0;
		LPFltAlpha = 3.40127;
		LPFltBaseLen = 7.78;
		LPFltCutoffMult = 0.78797;
	}
};

/**
 * @brief Set of resizing algorithm parameters for low-ringing performance
 * (7.86/1.065/0.000106).
 *
 * This set of resizing algorithm parameters offers a very low-ringing
 * performance at the expense of higher aliasing artifacts and a slightly
 * reduced contrast.
 */

struct CImageResizerParamsLR : public CImageResizerParams
{
	CImageResizerParamsLR()
	{
		CorrFltAlpha = 1.0;//7.86/1.96/1.065(73865.02)/0.000106:258636,437381
		CorrFltLen = 5.87671;
		IntFltAlpha = 2.25322;
		IntFltCutoff = 0.74090;
		IntFltLen = 18.0;
		LPFltAlpha = 1.79306;
		LPFltBaseLen = 7.00;
		LPFltCutoffMult = 0.68881;
	}
};

/**
 * @brief Set of resizing algorithm parameters for ultra low-aliasing
 * resizing (13.62/1.001/0.000557).
 *
 * This set of resizing algorithm parameters offers a very considerable
 * anti-aliasing performance with a good frequency response linearity (and
 * contrast). This set of parameters is computationally expensive and may
 * produce visible ringing artifacts on sharp features.
 */

struct CImageResizerParamsUltra : public CImageResizerParams
{
	CImageResizerParamsUltra()
	{
		CorrFltAlpha = 1.0;//13.62/1.79/1.001(24448.76)/0.000557:258654,457910
		CorrFltLen = 7.47096;
		IntFltAlpha = 1.94183;
		IntFltCutoff = 0.75473;
		IntFltLen = 18.0;
		LPFltAlpha = 5.51118;
		LPFltBaseLen = 8.34;
		LPFltCutoffMult = 0.78020;
	}
};

/**
 * @brief Image resizing variables class.
 * 
 * This is an utility "catch all" class that defines various variables used
 * during image resizing. Several variables that are explicitly initialized in
 * this class' constructor are also used as "input" variables to the image
 * resizing function.
 */

class CImageResizerVars
{
public:
	int ElCount; ///< The number of "fptype" elements used to store 1 pixel.
		///<
	int ElCountIO; ///< The number of source and destination image's elements
		///< used to store 1 pixel.
		///<
	int fppack; ///< The number of atomic types stored in a single "fptype"
		///< element.
		///<
	int fpalign; ///< Suggested alignment size in bytes. This is not a
		///< required alignment, because image resizing algorithm cannot be
		///< made to have a strictly aligned data access at all steps (e.g.
		///< interpolation cannot perform aligned accesses).
		///<
	int elalign; ///< Length alignment of arrays of elements. This applies to
		///< filters and intermediate buffers: this constant forces filters
		///< and scanlines to have a length which is a multiple of this value,
		///< for more efficient SIMD implementation. Value different to 1
		///< also means image pixels are de-interleaved during processing.
	int IntFltOrder; ///< Interpolation filter's order.
		///<
	int IntFltLen; ///< Interpolation filter's length in samples (taps).
		///<
	int BufLen; ///< Intermediate buffer's length in "fptype" elements.
	int BufOffset; ///< Offset into the intermediate buffer, used to provide
		///< prefix elements required during processing so that no "out of
		///< range" access happens. This offset is a multiple of ElCount if
		///< pixels are stored in interleaved form.
		///<
	int BufIncr; ///< Intermediate buffer increment in "fptype" elements. Used
		///< by de-interleaved processing algorithms: in this case each
		///< image's channel is stored independently, BufIncr elements apart.
		///<
	int RndSeed; ///< Random seed parameter. This parameter is incremented by
		///< 1 after each random generator initialization. The use of this
		///< variable depends on the ditherer implementation.
		///<
	CImageResizerThreadPool* ThreadPool; ///< Thread pool to be used by the
		///< image resizing function. Set to NULL to use single-threaded
		///< processing.
		///<

	CImageResizerVars()
		: RndSeed( 0 )
		, ThreadPool( NULL )
	{
	}
};

/**
 * @brief Image resizer's filtering step class.
 *
 * Class defines data to perform a single filtering step over a whole
 * horizontal or vertical scanline. Resizing consists of 1 or more steps that
 * may be performed before the actual resizing takes place. Filtering may also
 * follow a resizing step. Each step must ensure that scanline data contains
 * enough pixels to perform the next step (which may be resizing) without
 * exceeding scanline's bounds.
 *
 * A derived class must implement several "const" and "static" functions that
 * are used to perform the actual filtering in interleaved or de-interleaved
 * mode.
 *
 * @tparam fptype Floating point type to use for storing pixel data. SIMD
 * types can be used: in this case each element may hold a whole pixel.
 * @tparam fptypeatom The atomic type the "fptype" consists of.
 */

template< class fptype, class fptypeatom >
class CImageResizerFilterStep
{
public:
	bool IsUpsample; ///< "True" if this step is an upsampling step, "false"
		///< if downsampling step. Should be set to "false" if ResampleFactor
		///< equals 0.
		///<
	int ResampleFactor; ///< Resample factor (>=1). If 0, this is a resizing
		///< step. This value should be >1 if IsUpsample equals "true".
		///<
	CBuffer< fptype > Flt; ///< Filter kernel to use at this step.
		///<
	double DCGain; ///< DC gain which was applied to the filter. Not defined
		///< if ResampleFactor = 0.
		///<
	int FltLatency; ///< Filter's latency (shift) in pixels.
		///<
	const CImageResizerVars* Vars; ///< Image resizing-related variables.
		///<
	int InLen; ///< Input scanline's length in pixels.
		///<
	int InBuf; ///< Input buffer index, 0 or 1.
		///<
	int InPrefix; ///< Required input prefix pixels. These prefix pixels will
		///< be filled with source scanline's first pixel value. If IsUpsample
		///< is "true", this is the additional number of times the first pixel
		///< will be filtered before processing scanline, this number is also
		///< reflected in the OutPrefix.
		///<
	int InSuffix; ///< Required input suffix pixels. These suffix pixels will
		///< be filled with source scanline's last pixel value. If IsUpsample
		///< is "true", this is the additional number of times the last pixel
		///< will be filtered before processing scanline, this number is also
		///< reflected in the OutSuffix.
		///<
	int OutLen; ///< Length of the resulting scanline.
		///<
	int OutBuf; ///< Output buffer index. 0 or 1; 2 for the last step.
		///<
	int OutPrefix; ///< Required output prefix pixels. These prefix pixels
		///< will not be pre-filled with any values. Value is valid only if
		///< IsUpsample equals "true".
		///<
	int OutSuffix; ///< Required input suffix pixels. These suffix pixels will
		///< not be pre-filled with any values. Value is valid only if
		///< IsUpsample equals "true".
		///<
	CBuffer< fptype > PrefixDC; ///< DC component fluctuations added at the
		///< start of the resulting scanline, used when IsUpsample equals
		///< "true".
		///<
	CBuffer< fptype > SuffixDC; ///< DC component fluctuations added at the
		///< end of the resulting scanline, used when IsUpsample equals
		///< "true".
		///<

	/**
	 * @brief Resizing position structure.
	 *
	 * Structure holds resizing position and pointer to fractional delay
	 * filter.
	 */

	struct CResizePos
	{
		fptypeatom x; ///< "X" interpolation coefficient.
			///<
		const fptype* ftp; ///< Fractional delay filter pointer.
			///<
		int SrcOffs; ///< Source scanline offset.
			///<
	};

	CBuffer< CResizePos > RPosBuf; ///< Resizing positions buffer. Used when
		///< ResampleFactor equals 0 (resizing step).
		///<
};

/**
 * @brief Non-interleaved filtering steps implementation class.
 *
 * This class implements scanline filtering functions in non-interleaved mode.
 * This means that each pixel is processed independently.
 *
 * @tparam fptype Floating point type to use for storing pixel data. SIMD
 * types can be used: in this case each element may hold a whole pixel.
 * @tparam fptypeatom The atomic type the "fptype" consists of.
 */

template< class fptype, class fptypeatom >
class CImageResizerFilterStepNI :
	public CImageResizerFilterStep< fptype, fptypeatom >
{
public:
	using CImageResizerFilterStep< fptype, fptypeatom > :: IsUpsample;
	using CImageResizerFilterStep< fptype, fptypeatom > :: ResampleFactor;
	using CImageResizerFilterStep< fptype, fptypeatom > :: Flt;
	using CImageResizerFilterStep< fptype, fptypeatom > :: FltLatency;
	using CImageResizerFilterStep< fptype, fptypeatom > :: Vars;
	using CImageResizerFilterStep< fptype, fptypeatom > :: InLen;
	using CImageResizerFilterStep< fptype, fptypeatom > :: InPrefix;
	using CImageResizerFilterStep< fptype, fptypeatom > :: InSuffix;
	using CImageResizerFilterStep< fptype, fptypeatom > :: OutLen;
	using CImageResizerFilterStep< fptype, fptypeatom > :: OutPrefix;
	using CImageResizerFilterStep< fptype, fptypeatom > :: OutSuffix;
	using CImageResizerFilterStep< fptype, fptypeatom > :: PrefixDC;
	using CImageResizerFilterStep< fptype, fptypeatom > :: SuffixDC;
	using CImageResizerFilterStep< fptype, fptypeatom > :: RPosBuf;

	/**
	 * Function performs "packing" of a scanline and type conversion
	 * Scanline, depending on the "fptype" can be potentially stored as a
	 * packed SIMD values having a certain atomic type.
	 *
	 * @param ip Input scanline.
	 * @param op0 Output scanline.
	 * @param l0 The number of pixels to "unpack".
	 * @param Vars Image resizing-related variables.
	 */

	template< class Tin >
	static void packScanline( const Tin* ip, fptype* const op0, const int l0,
		const CImageResizerVars& Vars0 )
	{
		const int ElCount = Vars0.ElCount;
		const int ElCountIO = Vars0.ElCountIO;
		fptype* op = op0;
		int l = l0;

		if( ElCountIO == 1 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op;
				v[ 0 ] = (fptype) ip[ 0 ];
				op += ElCount;
				ip++;
				l--;
			}
		}
		else
		if( ElCountIO == 4 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op;
				v[ 0 ] = (fptype) ip[ 0 ];
				v[ 1 ] = (fptype) ip[ 1 ];
				v[ 2 ] = (fptype) ip[ 2 ];
				v[ 3 ] = (fptype) ip[ 3 ];
				op += ElCount;
				ip += 4;
				l--;
			}
		}
		else
		if( ElCountIO == 3 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op;
				v[ 0 ] = (fptype) ip[ 0 ];
				v[ 1 ] = (fptype) ip[ 1 ];
				v[ 2 ] = (fptype) ip[ 2 ];
				op += ElCount;
				ip += 3;
				l--;
			}
		}
		else
		if( ElCountIO == 2 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op;
				v[ 0 ] = (fptype) ip[ 0 ];
				v[ 1 ] = (fptype) ip[ 1 ];
				op += ElCount;
				ip += 2;
				l--;
			}
		}

		const int ZeroCount = ElCount * Vars0.fppack - ElCountIO;
		op = op0;
		l = l0;

		if( ZeroCount == 1 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op + ElCountIO;
				v[ 0 ] = (fptype) 0.0;
				op += ElCount;
				l--;
			}
		}
		else
		if( ZeroCount == 2 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op + ElCountIO;
				v[ 0 ] = (fptype) 0.0;
				v[ 1 ] = (fptype) 0.0;
				op += ElCount;
				l--;
			}
		}
		else
		if( ZeroCount == 3 )
		{
			while( l > 0 )
			{
				fptypeatom* v = (fptypeatom*) op + ElCountIO;
				v[ 0 ] = (fptype) 0.0;
				v[ 1 ] = (fptype) 0.0;
				v[ 2 ] = (fptype) 0.0;
				op += ElCount;
				l--;
			}
		}
	}

	/**
	 * Function converts vertical scanline to horizontal scanline. This
	 * function is called by the image resizer when image is resized
	 * vertically. This means that the vertical scanline is stored in the
	 * same format produced by the packScanline() function.
	 *
	 * @param ip Input vertical scanline.
	 * @param op Output buffer (temporary buffer used during resizing).
	 * @param SrcLen The number of pixels in the input scanline, also used to
	 * calculate input buffer increment.
	 * @param SrcIncr Input buffer increment to the next vertical pixel.
	 * @param Vars0 Image resizing-related variables.
	 */

	static void convertVtoH( const fptype* ip, fptype* op, const int SrcLen,
		const int SrcIncr, const CImageResizerVars& Vars0 )
	{
		const int ElCount = Vars0.ElCount;
		int j;

		if( ElCount == 1 )
		{
			for( j = 0; j < SrcLen; j++ )
			{
				op[ 0 ] = ip[ 0 ];
				ip += SrcIncr;
				op++;
			}
		}
		else
		if( ElCount == 4 )
		{
			for( j = 0; j < SrcLen; j++ )
			{
				op[ 0 ] = ip[ 0 ];
				op[ 1 ] = ip[ 1 ];
				op[ 2 ] = ip[ 2 ];
				op[ 3 ] = ip[ 3 ];
				ip += SrcIncr;
				op += 4;
			}
		}
		else
		if( ElCount == 3 )
		{
			for( j = 0; j < SrcLen; j++ )
			{
				op[ 0 ] = ip[ 0 ];
				op[ 1 ] = ip[ 1 ];
				op[ 2 ] = ip[ 2 ];
				ip += SrcIncr;
				op += 3;
			}
		}
		else
		if( ElCount == 2 )
		{
			for( j = 0; j < SrcLen; j++ )
			{
				op[ 0 ] = ip[ 0 ];
				op[ 1 ] = ip[ 1 ];
				ip += SrcIncr;
				op += 2;
			}
		}
	}

	/**
	 * Function performs "unpacking" of a scanline and type conversion
	 * (truncation is used when floating point is converted to integer).
	 * Scanline, depending on the "fptype" can be potentially stored as a
	 * packed SIMD values having a certain atomic type. The unpacking function
	 * assumes that scanline is stored in the style produced by the
	 * packScanline() function.
	 *
	 * @param ip Input scanline.
	 * @param op Output scanline.
	 * @param l The number of pixels to "unpack".
	 * @param Vars0 Image resizing-related variables.
	 */

	template< class Tout >
	static void unpackScanline( const fptype* ip, Tout* op, int l,
		const CImageResizerVars& Vars0 )
	{
		const int ElCount = Vars0.ElCount;
		const int ElCountIO = Vars0.ElCountIO;

		if( ElCountIO == 1 )
		{
			while( l > 0 )
			{
				const fptypeatom* v = (const fptypeatom*) ip;
				op[ 0 ] = (Tout) v[ 0 ];
				ip += ElCount;
				op++;
				l--;
			}
		}
		else
		if( ElCountIO == 4 )
		{
			while( l > 0 )
			{
				const fptypeatom* v = (const fptypeatom*) ip;
				op[ 0 ] = (Tout) v[ 0 ];
				op[ 1 ] = (Tout) v[ 1 ];
				op[ 2 ] = (Tout) v[ 2 ];
				op[ 3 ] = (Tout) v[ 3 ];
				ip += ElCount;
				op += 4;
				l--;
			}
		}
		else
		if( ElCountIO == 3 )
		{
			while( l > 0 )
			{
				const fptypeatom* v = (const fptypeatom*) ip;
				op[ 0 ] = (Tout) v[ 0 ];
				op[ 1 ] = (Tout) v[ 1 ];
				op[ 2 ] = (Tout) v[ 2 ];
				ip += ElCount;
				op += 3;
				l--;
			}
		}
		else
		if( ElCountIO == 2 )
		{
			while( l > 0 )
			{
				const fptypeatom* v = (const fptypeatom*) ip;
				op[ 0 ] = (Tout) v[ 0 ];
				op[ 1 ] = (Tout) v[ 1 ];
				ip += ElCount;
				op += 2;
				l--;
			}
		}
	}

	/**
	 * Function that replicates a set of adjacent elements several times in a
	 * row. This operation is usually used to replicate pixels at the start or
	 * end of image's scanline.
	 *
	 * @param ip Source array.
	 * @param ipl Source array length (usually 1..4, but can be any number).
	 * @param[out] op Destination buffer.
	 * @param l Number of times the source array should be replicated (the
	 * destination buffer should be able to hold ipl * l number of elements).
	 * @param opinc Destination buffer position increment after replicating
	 * the source array. This value should be equal to at least ipl.
	 */

	template< class T1, class T2 >
	static void replicateArray( const T1* const ip, const int ipl, T2* op,
		int l, const int opinc )
	{
		if( ipl == 1 )
		{
			while( l > 0 )
			{
				op[ 0 ] = ip[ 0 ];
				op += opinc;
				l--;
			}
		}
		else
		if( ipl == 4 )
		{
			while( l > 0 )
			{
				op[ 0 ] = ip[ 0 ];
				op[ 1 ] = ip[ 1 ];
				op[ 2 ] = ip[ 2 ];
				op[ 3 ] = ip[ 3 ];
				op += opinc;
				l--;
			}
		}
		else
		if( ipl == 3 )
		{
			while( l > 0 )
			{
				op[ 0 ] = ip[ 0 ];
				op[ 1 ] = ip[ 1 ];
				op[ 2 ] = ip[ 2 ];
				op += opinc;
				l--;
			}
		}
		else
		if( ipl == 2 )
		{
			while( l > 0 )
			{
				op[ 0 ] = ip[ 0 ];
				op[ 1 ] = ip[ 1 ];
				op += opinc;
				l--;
			}
		}
		else
		{
			while( l > 0 )
			{
				int i;

				for( i = 0; i < ipl; i++ )
				{
					op[ i ] = ip[ i ];
				}

				op += opinc;
				l--;
			}
		}
	}

	/**
	 * Function prepares input scanline buffer for *this filtering step.
	 * Left- and right-most pixels are replicated to make sure no buffer
	 * overrun happens. Such approach also allows to bypass any pointer
	 * range checks.
	 *
	 * @param Src Source buffer.
	 */

	void prepareInBuf( fptype* Src ) const
	{
		if( IsUpsample || InPrefix + InSuffix == 0 )
		{
			return;
		}

		const int ElCount = Vars -> ElCount;
		replicateArray( Src, ElCount, Src - ElCount, InPrefix, -ElCount );

		Src += ( InLen - 1 ) * ElCount;
		replicateArray( Src, ElCount, Src + ElCount, InSuffix, ElCount );
	}

	/**
	 * Function peforms scanline upsampling with filtering.
	 *
	 * @param Src Source scanline buffer (length = this -> InLen). Source
	 * scanline increment will be equal to ElCount.
	 * @param Dst Destination scanline buffer.
	 */

	void doUpsample( const fptype* const Src, fptype* const Dst ) const
	{
		const int ElCount = Vars -> ElCount;
		fptype* op0 = &Dst[ -OutPrefix * ElCount ];
		memset( op0, 0, ( OutPrefix + OutLen + OutSuffix ) * ElCount *
			sizeof( fptype ));

		const fptype* const f = Flt;
		const int flen = Flt.getCapacity();
		const fptype* ip = Src;
		fptype* op;
		const int opstep = ElCount * ResampleFactor;
		int l;
		int i;

		if( ElCount == 1 )
		{
			l = InPrefix;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op++;
				}

				op0 += opstep;
				l--;
			}

			l = InLen - 1;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op++;
				}

				ip += ElCount;
				op0 += opstep;
				l--;
			}

			l = InSuffix;

			while( l >= 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op++;
				}

				op0 += opstep;
				l--;
			}
		}
		else
		if( ElCount == 4 )
		{
			l = InPrefix;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op[ 2 ] += f[ i ] * ip[ 2 ];
					op[ 3 ] += f[ i ] * ip[ 3 ];
					op += 4;
				}

				op0 += opstep;
				l--;
			}

			l = InLen - 1;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op[ 2 ] += f[ i ] * ip[ 2 ];
					op[ 3 ] += f[ i ] * ip[ 3 ];
					op += 4;
				}

				ip += ElCount;
				op0 += opstep;
				l--;
			}

			l = InSuffix;

			while( l >= 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op[ 2 ] += f[ i ] * ip[ 2 ];
					op[ 3 ] += f[ i ] * ip[ 3 ];
					op += 4;
				}

				op0 += opstep;
				l--;
			}
		}
		else
		if( ElCount == 3 )
		{
			l = InPrefix;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op[ 2 ] += f[ i ] * ip[ 2 ];
					op += 3;
				}

				op0 += opstep;
				l--;
			}

			l = InLen - 1;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op[ 2 ] += f[ i ] * ip[ 2 ];
					op += 3;
				}

				ip += ElCount;
				op0 += opstep;
				l--;
			}

			l = InSuffix;

			while( l >= 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op[ 2 ] += f[ i ] * ip[ 2 ];
					op += 3;
				}

				op0 += opstep;
				l--;
			}
		}
		else
		if( ElCount == 2 )
		{
			l = InPrefix;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op += 2;
				}

				op0 += opstep;
				l--;
			}

			l = InLen - 1;

			while( l > 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op += 2;
				}

				ip += ElCount;
				op0 += opstep;
				l--;
			}

			l = InSuffix;

			while( l >= 0 )
			{
				op = op0;

				for( i = 0; i < flen; i++ )
				{
					op[ 0 ] += f[ i ] * ip[ 0 ];
					op[ 1 ] += f[ i ] * ip[ 1 ];
					op += 2;
				}

				op0 += opstep;
				l--;
			}
		}

		op = op0;
		const fptype* dc = SuffixDC;
		l = SuffixDC.getCapacity();

		if( ElCount == 1 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				dc++;
				op++;
				l--;
			}
		}
		else
		if( ElCount == 4 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				op[ 1 ] += ip[ 1 ] * dc[ 0 ];
				op[ 2 ] += ip[ 2 ] * dc[ 0 ];
				op[ 3 ] += ip[ 3 ] * dc[ 0 ];
				dc++;
				op += 4;
				l--;
			}
		}
		else
		if( ElCount == 3 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				op[ 1 ] += ip[ 1 ] * dc[ 0 ];
				op[ 2 ] += ip[ 2 ] * dc[ 0 ];
				dc++;
				op += 3;
				l--;
			}
		}
		else
		if( ElCount == 2 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				op[ 1 ] += ip[ 1 ] * dc[ 0 ];
				dc++;
				op += 2;
				l--;
			}
		}

		ip = Src;
		op = Dst - InPrefix * opstep;
		dc = PrefixDC;
		l = PrefixDC.getCapacity();

		if( ElCount == 1 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				dc++;
				op++;
				l--;
			}
		}
		else
		if( ElCount == 4 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				op[ 1 ] += ip[ 1 ] * dc[ 0 ];
				op[ 2 ] += ip[ 2 ] * dc[ 0 ];
				op[ 3 ] += ip[ 3 ] * dc[ 0 ];
				dc++;
				op += 4;
				l--;
			}
		}
		else
		if( ElCount == 3 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				op[ 1 ] += ip[ 1 ] * dc[ 0 ];
				op[ 2 ] += ip[ 2 ] * dc[ 0 ];
				dc++;
				op += 3;
				l--;
			}
		}
		else
		if( ElCount == 2 )
		{
			while( l > 0 )
			{
				op[ 0 ] += ip[ 0 ] * dc[ 0 ];
				op[ 1 ] += ip[ 1 ] * dc[ 0 ];
				dc++;
				op += 2;
				l--;
			}
		}
	}

	/**
	 * Function extends *this upsampling step so that it produces more
	 * upsampled pixels that cover the prefix and suffix needs of the next
	 * step. After the call to this function the InPrefix and InSuffix
	 * variables of the next step will be set to zero.
	 */

	void extendUpsample( CImageResizerFilterStepNI& NextStep )
	{
		InPrefix = ( NextStep.InPrefix + ResampleFactor - 1 ) /
			ResampleFactor;

		OutPrefix += InPrefix * ResampleFactor;
		NextStep.InPrefix = 0;

		InSuffix = ( NextStep.InSuffix + ResampleFactor - 1 ) /
			ResampleFactor;

		OutSuffix += InSuffix * ResampleFactor;
		NextStep.InSuffix = 0;
	}

	/**
	 * Function peforms scanline filtering with optional downsampling.
	 * Function makes use of the symmetry of the filter kernel.
	 *
	 * @param Src Source scanline buffer (length = this -> InLen). Source
	 * scanline increment will be equal to ElCount.
	 * @param Dst Destination scanline buffer.
	 * @param DstIncr Destination scanline buffer increment.
	 */

	void doFilter( const fptype* const Src, fptype* Dst,
		const int DstIncr ) const
	{
		const int ElCount = Vars -> ElCount;
		const fptype* const f = &Flt[ FltLatency ];
		const int flen = FltLatency + 1;
		const fptype* ip = Src;
		const fptype* ip1;
		const fptype* ip2;
		const int ipstep = ElCount * ResampleFactor;
		int l = OutLen;
		int i;

		if( ElCount == 1 )
		{
			while( l > 0 )
			{
				fptype s = f[ 0 ] * ip[ 0 ];
				ip1 = ip;
				ip2 = ip;

				for( i = 1; i < flen; i++ )
				{
					ip1++;
					ip2--;
					s += f[ i ] * ( ip1[ 0 ] + ip2[ 0 ]);
				}

				Dst[ 0 ] = s;
				Dst += DstIncr;
				ip += ipstep;
				l--;
			}
		}
		else
		if( ElCount == 4 )
		{
			while( l > 0 )
			{
				fptype s1 = f[ 0 ] * ip[ 0 ];
				fptype s2 = f[ 0 ] * ip[ 1 ];
				fptype s3 = f[ 0 ] * ip[ 2 ];
				fptype s4 = f[ 0 ] * ip[ 3 ];
				ip1 = ip;
				ip2 = ip;

				for( i = 1; i < flen; i++ )
				{
					ip1 += 4;
					ip2 -= 4;
					s1 += f[ i ] * ( ip1[ 0 ] + ip2[ 0 ]);
					s2 += f[ i ] * ( ip1[ 1 ] + ip2[ 1 ]);
					s3 += f[ i ] * ( ip1[ 2 ] + ip2[ 2 ]);
					s4 += f[ i ] * ( ip1[ 3 ] + ip2[ 3 ]);
				}

				Dst[ 0 ] = s1;
				Dst[ 1 ] = s2;
				Dst[ 2 ] = s3;
				Dst[ 3 ] = s4;
				Dst += DstIncr;
				ip += ipstep;
				l--;
			}
		}
		else
		if( ElCount == 3 )
		{
			while( l > 0 )
			{
				fptype s1 = f[ 0 ] * ip[ 0 ];
				fptype s2 = f[ 0 ] * ip[ 1 ];
				fptype s3 = f[ 0 ] * ip[ 2 ];
				ip1 = ip;
				ip2 = ip;

				for( i = 1; i < flen; i++ )
				{
					ip1 += 3;
					ip2 -= 3;
					s1 += f[ i ] * ( ip1[ 0 ] + ip2[ 0 ]);
					s2 += f[ i ] * ( ip1[ 1 ] + ip2[ 1 ]);
					s3 += f[ i ] * ( ip1[ 2 ] + ip2[ 2 ]);
				}

				Dst[ 0 ] = s1;
				Dst[ 1 ] = s2;
				Dst[ 2 ] = s3;
				Dst += DstIncr;
				ip += ipstep;
				l--;
			}
		}
		else
		if( ElCount == 2 )
		{
			while( l > 0 )
			{
				fptype s1 = f[ 0 ] * ip[ 0 ];
				fptype s2 = f[ 0 ] * ip[ 1 ];
				ip1 = ip;
				ip2 = ip;

				for( i = 1; i < flen; i++ )
				{
					ip1 += 2;
					ip2 -= 2;
					s1 += f[ i ] * ( ip1[ 0 ] + ip2[ 0 ]);
					s2 += f[ i ] * ( ip1[ 1 ] + ip2[ 1 ]);
				}

				Dst[ 0 ] = s1;
				Dst[ 1 ] = s2;
				Dst += DstIncr;
				ip += ipstep;
				l--;
			}
		}
	}

	/**
	 * Function performs resizing of a single scanline. This function does
	 * not "know" about the length of the source scanline buffer. This buffer
	 * should be padded with enough pixels so that ( SrcPos - FilterLenD2 ) is
	 * always >= 0 and ( SrcPos + ( DstLineLen - 1 ) * k + FilterLenD2 + 1 )
	 * does not exceed source scanline's buffer length. SrcLine's increment is
	 * assumed to be equal to ElCount.
	 *
	 * @param SrcLine Source scanline buffer.
	 * @param DstLine Destination (resized) scanline buffer.
	 * @param DstLineInc Destination scanline position increment, should be
	 * divisible by ElCount.
	 */

	void doResize( const fptype* SrcLine, fptype* DstLine,
		const int DstLineInc ) const
	{
		const int IntFltLen = Vars -> IntFltLen;
		const int ElCount = Vars -> ElCount;
		const typename CImageResizerFilterStep< fptype, fptypeatom > ::
			CResizePos* rpos = &RPosBuf[ 0 ];

		int DstLineLen = OutLen;

#define AVIR_RESIZE_PART1 \
			while( DstLineLen > 0 ) \
			{ \
				const fptype x = (fptype) rpos -> x; \
				const fptype* ftp = rpos -> ftp; \
				const fptype* Src = SrcLine + rpos -> SrcOffs; \
				int l = IntFltLen;

#define AVIR_RESIZE_PART1nx \
			while( DstLineLen > 0 ) \
			{ \
				const fptype* ftp = rpos -> ftp; \
				const fptype* Src = SrcLine + rpos -> SrcOffs; \
				int l = IntFltLen;

#define AVIR_RESIZE_PART2 \
				DstLineLen--; \
				DstLine += DstLineInc; \
				rpos++; \
			}

		if( Vars -> IntFltOrder == 1 )
		{
			if( ElCount == 1 )
			{
				AVIR_RESIZE_PART1

				fptype sum = 0.0;

				while( l > 0 )
				{
					sum += ( ftp[ 0 ] + ftp[ IntFltLen ] * x ) * Src[ 0 ];
					ftp++;
					Src++;
					l--;
				}

				DstLine[ 0 ] = sum;

				AVIR_RESIZE_PART2
			}
			else
			if( ElCount == 4 )
			{
				AVIR_RESIZE_PART1

				fptype sum[ 4 ];
				sum[ 0 ] = 0.0;
				sum[ 1 ] = 0.0;
				sum[ 2 ] = 0.0;
				sum[ 3 ] = 0.0;

				while( l > 0 )
				{
					const fptype xx = ftp[ 0 ] + ftp[ IntFltLen ] * x;
					sum[ 0 ] += xx * Src[ 0 ];
					sum[ 1 ] += xx * Src[ 1 ];
					sum[ 2 ] += xx * Src[ 2 ];
					sum[ 3 ] += xx * Src[ 3 ];
					ftp++;
					Src += 4;
					l--;
				}

				DstLine[ 0 ] = sum[ 0 ];
				DstLine[ 1 ] = sum[ 1 ];
				DstLine[ 2 ] = sum[ 2 ];
				DstLine[ 3 ] = sum[ 3 ];

				AVIR_RESIZE_PART2
			}
			else
			if( ElCount == 3 )
			{
				AVIR_RESIZE_PART1

				fptype sum[ 3 ];
				sum[ 0 ] = 0.0;
				sum[ 1 ] = 0.0;
				sum[ 2 ] = 0.0;

				while( l > 0 )
				{
					const fptype xx = ftp[ 0 ] + ftp[ IntFltLen ] * x;
					sum[ 0 ] += xx * Src[ 0 ];
					sum[ 1 ] += xx * Src[ 1 ];
					sum[ 2 ] += xx * Src[ 2 ];
					ftp++;
					Src += 3;
					l--;
				}

				DstLine[ 0 ] = sum[ 0 ];
				DstLine[ 1 ] = sum[ 1 ];
				DstLine[ 2 ] = sum[ 2 ];

				AVIR_RESIZE_PART2
			}
			else
			if( ElCount == 2 )
			{
				AVIR_RESIZE_PART1

				fptype sum[ 2 ];
				sum[ 0 ] = 0.0;
				sum[ 1 ] = 0.0;

				while( l > 0 )
				{
					const fptype xx = ftp[ 0 ] + ftp[ IntFltLen ] * x;
					sum[ 0 ] += xx * Src[ 0 ];
					sum[ 1 ] += xx * Src[ 1 ];
					ftp++;
					Src += 2;
					l--;
				}

				DstLine[ 0 ] = sum[ 0 ];
				DstLine[ 1 ] = sum[ 1 ];

				AVIR_RESIZE_PART2
			}
		}
		else
		{
			if( ElCount == 1 )
			{
				AVIR_RESIZE_PART1nx

				fptype sum = 0.0;

				while( l > 0 )
				{
					sum += ftp[ 0 ] * Src[ 0 ];
					ftp++;
					Src++;
					l--;
				}

				DstLine[ 0 ] = sum;

				AVIR_RESIZE_PART2
			}
			else
			if( ElCount == 4 )
			{
				AVIR_RESIZE_PART1nx

				fptype sum[ 4 ];
				sum[ 0 ] = 0.0;
				sum[ 1 ] = 0.0;
				sum[ 2 ] = 0.0;
				sum[ 3 ] = 0.0;

				while( l > 0 )
				{
					const fptype xx = ftp[ 0 ];
					sum[ 0 ] += xx * Src[ 0 ];
					sum[ 1 ] += xx * Src[ 1 ];
					sum[ 2 ] += xx * Src[ 2 ];
					sum[ 3 ] += xx * Src[ 3 ];
					ftp++;
					Src += 4;
					l--;
				}

				DstLine[ 0 ] = sum[ 0 ];
				DstLine[ 1 ] = sum[ 1 ];
				DstLine[ 2 ] = sum[ 2 ];
				DstLine[ 3 ] = sum[ 3 ];

				AVIR_RESIZE_PART2
			}
			else
			if( ElCount == 3 )
			{
				AVIR_RESIZE_PART1nx

				fptype sum[ 3 ];
				sum[ 0 ] = 0.0;
				sum[ 1 ] = 0.0;
				sum[ 2 ] = 0.0;

				while( l > 0 )
				{
					const fptype xx = ftp[ 0 ];
					sum[ 0 ] += xx * Src[ 0 ];
					sum[ 1 ] += xx * Src[ 1 ];
					sum[ 2 ] += xx * Src[ 2 ];
					ftp++;
					Src += 3;
					l--;
				}

				DstLine[ 0 ] = sum[ 0 ];
				DstLine[ 1 ] = sum[ 1 ];
				DstLine[ 2 ] = sum[ 2 ];

				AVIR_RESIZE_PART2
			}
			else
			if( ElCount == 2 )
			{
				AVIR_RESIZE_PART1nx

				fptype sum[ 2 ];
				sum[ 0 ] = 0.0;
				sum[ 1 ] = 0.0;

				while( l > 0 )
				{
					const fptype xx = ftp[ 0 ];
					sum[ 0 ] += xx * Src[ 0 ];
					sum[ 1 ] += xx * Src[ 1 ];
					ftp++;
					Src += 2;
					l--;
				}

				DstLine[ 0 ] = sum[ 0 ];
				DstLine[ 1 ] = sum[ 1 ];

				AVIR_RESIZE_PART2
			}
		}
	}
#undef AVIR_RESIZE_PART1
#undef AVIR_RESIZE_PART2
};

/**
 * @brief Image resizer's default dithering class.
 *
 * This class defines an object that performs rounding, clipping and dithering
 * operations over horizontal scanline pixels before scanline is stored in the
 * output buffer.
 *
 * The ditherer should expect the same storage order of the pixels in a
 * scanline as used in the "filter step" class. So, a separate ditherer class
 * should be defined for each scanline pixel storage style. The default
 * ditherer implements a simple rounding without dithering: it can be used for
 * an efficient dithering method which can be multi-threaded.
 *
 * @tparam fptype Floating point type to use for storing pixel data. SIMD
 * types can be used.
 */

template< class fptype >
class CImageResizerDithererDef
{
public:
	/**
	 * Function initializes the ditherer object.
	 *
	 * @param aLen Scanline length in pixels to process.
	 * @param aVars Image resizing-related variables.
	 * @param aTrMul Bit-depth truncation multiplier. 1 - no additional
	 * truncation.
	 * @param aPkOut Peak output value allowed.
	 */

	void init( const int aLen, const CImageResizerVars& aVars,
		const fptype aTrMul, const fptype aPkOut )
	{
		Len = aLen;
		Vars = &aVars;
		LenE = aLen * Vars -> ElCount;
		TrMul = aTrMul;
		PkOut = aPkOut;
	}

	/**
	 * @return "True" if dithering is recursive relative to scanlines meaning
	 * multi-threaded execution is not supported by this dithering method.
	 */

	static bool isRecursive()
	{
		return( false );
	}

	/**
	 * Function performs rounding and clipping operations.
	 */

	void dither( fptype* const ResScanline ) const
	{
		int j;

		for( j = 0; j < LenE; j++ )
		{
			const fptype z0 = round( ResScanline[ j ] / TrMul ) * TrMul;
			ResScanline[ j ] = clamp( z0, (fptype) 0.0, PkOut );
		}
	}

protected:
	int Len; ///< Scanline's length in pixels.
		///<
	const CImageResizerVars* Vars; ///< Image resizing-related variables.
		///<
	int LenE; ///< = LenE * ElCount.
		///<
	fptype TrMul; ///< Bit-depth truncation multiplier.
		///<
	fptype PkOut; ///< Peak output value allowed.
		///<
};

/**
 * @brief Image resizer's quasi-random dithering class.
 *
 * This ditherer implements a classic quasi-random dithering which looks OK
 * and whose results are compressed by PNG well.
 *
 * @tparam fptype Floating point type to use for storing pixel data. SIMD
 * types can be used.
 */

template< class fptype >
class CImageResizerDithererQRnd : public CImageResizerDithererDef< fptype >
{
public:
	/**
	 * Function initializes the ditherer object.
	 *
	 * @param aLen Scanline length in pixels to process.
	 * @param aVars Image resizing-related variables.
	 * @param aTrMul Bit-depth truncation multiplier. 1 - no additional
	 * truncation.
	 * @param aPkOut Peak output value allowed.
	 */

	void init( const int aLen, const CImageResizerVars& aVars,
		const fptype aTrMul, const fptype aPkOut )
	{
		CImageResizerDithererDef< fptype > :: init( aLen, aVars, aTrMul,
			aPkOut );

		ResScanlineDith0.alloc( LenE + Vars -> ElCount, sizeof( fptype ));
		ResScanlineDith = ResScanlineDith0 + Vars -> ElCount;
		int i;

		for( i = 0; i < LenE + Vars -> ElCount; i++ )
		{
			ResScanlineDith0[ i ] = 0.0;
		}
	}

	static bool isRecursive()
	{
		return( true );
	}

	void dither( fptype* const ResScanline )
	{
		const int ElCount = Vars -> ElCount;
		int j;

		for( j = 0; j < LenE; j++ )
		{
			ResScanline[ j ] += ResScanlineDith[ j ];
			ResScanlineDith[ j ] = 0.0;
		}

		for( j = 0; j < LenE - ElCount; j++ )
		{
			// Perform rounding, noise estimation and saturation.

			const fptype z0 = round( ResScanline[ j ] / TrMul ) * TrMul;
			const fptype Noise = ResScanline[ j ] - z0;
			ResScanline[ j ] = clamp( z0, (fptype) 0.0, PkOut );

			ResScanline[ j + ElCount ] += Noise * (fptype) 0.4375;
			ResScanlineDith[ j + ElCount ] += Noise * (fptype) 0.0625;
			ResScanlineDith[ j ] += Noise * (fptype) 0.3125;
			ResScanlineDith[ j - ElCount ] += Noise * (fptype) 0.1875;
		}

		while( j < LenE )
		{
			const fptype z0 = round( ResScanline[ j ] / TrMul ) * TrMul;
			const fptype Noise = ResScanline[ j ] - z0;
			ResScanline[ j ] = clamp( z0, (fptype) 0.0, PkOut );

			ResScanlineDith[ j ] += Noise * (fptype) 0.3125;
			ResScanlineDith[ j - ElCount ] += Noise * (fptype) 0.1875;
			j++;
		}
	}

protected:
	using CImageResizerDithererDef< fptype > :: Len;
	using CImageResizerDithererDef< fptype > :: Vars;
	using CImageResizerDithererDef< fptype > :: LenE;
	using CImageResizerDithererDef< fptype > :: TrMul;
	using CImageResizerDithererDef< fptype > :: PkOut;

	CBuffer< fptype > ResScanlineDith0; ///< Error propagation buffer for
		///< dithering, first pixel unused.
		///<
	fptype* ResScanlineDith; ///< Error propagation buffer pointer which skips
		///< the first ElCount elements.
		///<
};

/**
 * @brief Floating-point processing definition and abstraction class.
 *
 * This class defines several constants and typedefs that point to classes
 * that should be used by the image resizing algorithm. Such "definition
 * class" can be used to define alternative scanline processing algorithms
 * (e.g. SIMD) and image scanline packing styles used during processing. This
 * class also offers an abstraction layer for dithering, rounding and
 * saturation operation.
 * 
 * The fpclass_def class can be used to define processing using both SIMD and
 * non-SIMD types, but using algorithms that are non-SIMD optimized
 * themselves.
 *
 * @tparam afptype Floating point type to use for storing intermediate data
 * and variables. For variables that are not used in intensive calculations
 * the "double" type is always used. On the latest Intel processors (like
 * i7-4770K) there is almost no performance difference between "double" and
 * "float". Image quality differences between "double" and "float" are not
 * apparent on 8-bit images. At the same time the "float" uses half amount of
 * working memory the "double" type uses. SIMD types can be used. The
 * functions round() and clamp() in the "avir" or other visible namespace
 * should be available for the specified type.
 * @tparam afptypeatom The atomic type the "afptype" consists of.
 */

template< class afptype, class afptypeatom = afptype >
class fpclass_def
{
public:
	typedef afptype fptype; ///< Floating-point type to use during processing.
		///<
	typedef afptypeatom fptypeatom; ///< Atomic type.
		///<
	static const int fppack = sizeof( fptype ) / sizeof( fptypeatom ); ///<
		///< The number of atomic types stored in a single "fptype" element.
		///<
	static const int fpalign = sizeof( fptype ); ///< Suggested alignment size
		///< in bytes. This is not a required alignment, because image
		///< resizing algorithm cannot be made to have a strictly aligned data
		///< access at all steps (e.g. interpolation cannot perform aligned
		///< accesses).
		///<
	static const int elalign = 1; ///< Length alignment of arrays of elements.
		///< This applies to filters and intermediate buffers: this constant
		///< forces filters and scanlines to have a length which is a multiple
		///< of this value, for more efficient SIMD implementation. Value
		///< different to 1 also means image pixels are de-interleaved during
		///< processing.
		///<
	typedef CImageResizerFilterStepNI< fptype, fptypeatom > CFilterStep; ///<
		///< Filter step class to use during processing.
		///<
	typedef CImageResizerDithererQRnd< fptype > CDitherer; ///<
		///< Ditherer class to use during processing.
		///<
};

/**
 * @brief Image resizer class.
 *
 * The object of this class can be used to resize 1-4 channel images to any
 * required size. Resizing is performed by utilizing interpolated sinc
 * fractional delay filters plus (if necessary) a cascade of built-in
 * sinc function-based 2X upsampling or 2X downsampling stages, followed by a
 * correction filtering.
 *
 * Object of this class can be allocated on stack.
 *
 * @tparam fpclass Floating-point processing definition class to use. See
 * avir::fpclass_def for more details.
 */

template< class fpclass = fpclass_def< float > >
class CImageResizer
{
public:
	/**
	 * Constructor initializes the resizer.
	 *
	 * @param aResBitDepth Required bit depth of resulting image (4-16). If
	 * integer value output is used (e.g. uint8_t), the bit depth also affects
	 * rounding: for example, if aResBitDepth=6 and "Tout" is uint8_t, the
	 * result will be rounded to 6 most significant bits (2 least significant
	 * bits truncated, with dithering applied). The source image may have any
	 * real bit-depth: if this image was correctly dithered, during downsizing
	 * its bit-depth will increase proportionally to the downsizing factor.
	 * @param aSrcBitDepth Source image's real bit-depth. Set to 0 to use
	 * aResBitDepth.
	 * @param aParams Resizing algorithm's parameters to use. Leave out for
	 * default values. Can be useful when performing automatic optimization of
	 * parameters.
	 */

	CImageResizer( const int aResBitDepth = 8, const int aSrcBitDepth = 0,
		const CImageResizerParams& aParams = CImageResizerParamsDef() )
		: Params( aParams )
		, ResBitDepth( aResBitDepth )
	{
		const int SrcBitDepth = ( aSrcBitDepth == 0 ? ResBitDepth :
			aSrcBitDepth );

		const int IntBitDepth = ( ResBitDepth > SrcBitDepth ? ResBitDepth :
			SrcBitDepth );

		const double SNR = -6.02 * ( IntBitDepth + 3 );
		int UseOrder;

		if( IntBitDepth > 8 )
		{
			UseOrder = 1; // -146 dB max
			FilterFracs = (int) ceil( 0.23134052 * exp( -0.058062929 * SNR ));
		}
		else
		{
			UseOrder = 0; // -72 dB max
			FilterFracs = (int) ceil( 0.33287686 * exp( -0.11334583 * SNR ));
		}

		if( FilterFracs < 2 )
		{
			FilterFracs = 2;
		}

		FilterBank.init( FilterFracs, UseOrder, Params.IntFltLen,
			Params.IntFltCutoff, Params.IntFltAlpha, fpclass :: fpalign,
			fpclass :: elalign );
	}

	/**
	 * Function resizes image.
	 *
	 * @param SrcBuf Source image buffer.
	 * @param SrcWidth Source image width.
	 * @param SrcHeight Source image height.
	 * @param SrcScanlineSize Physical size of source scanline in elements. If
	 * this value is below 1, SrcWidth * ElCountIO will be used as the
	 * physical source scanline size.
	 * @param NewBuf Buffer to accept the resized image. Can be equal to
	 * SrcBuf if the size of the resized image is smaller or equal to source
	 * image in size.
	 * @param NewWidth New image width.
	 * @param NewHeight New image height.
	 * @param ElCountIO The number of elements (channels) used to store each
	 * source and destination pixel (1-4).
	 * @param k Resizing step (one output pixel corresponds to "k" input
	 * pixels). A downsizing factor if > 1.0; upsizing factor if <= 1.0.
	 * Multiply by -1 if you would like to bypass "ox" and "oy" adjustment
	 * which is done by default to produce a centered image. If step value
	 * equals 0, the step size will be chosen automatically and independently
	 * for horizontal and vertical resizing.
	 * @param ox Start X pixel offset within source image (can be negative).
	 * Positive offset moves image to the left.
	 * @param oy Start Y pixel offset within source image (can be negative).
	 * Positive offset moves image to the top.
	 * @param aVars Pointer to variables to be passed to the image resizing
	 * function. Can be NULL. Only variables that are initialized in default
	 * constructor are accepted by this function. Any variables in this object
	 * can be modified by this function. The access to this object is not
	 * thread-safe, each concurrent instance of this function should use a
	 * separate aVars object.
	 * @tparam Tin Input buffer element's type. Can be uint8_t (0-255 value
	 * range), uint16_t (0-65535 value range), float (0.0-1.0 value range),
	 * double (0.0-1.0 value range). Larger integer types are treated as
	 * uint16_t. Signed integer types are unsupported.
	 * @tparam Tout Output buffer element's type. Can be uint8_t (0-255 value
	 * range), uint16_t (0-65535 value range), float (0.0-1.0 value range),
	 * double (0.0-1.0 value range). Larger integer types are treated as
	 * uint16_t. Signed integer types are unsupported.
	 */

	template< class Tin, class Tout >
	void resizeImage( const Tin* const SrcBuf, const int SrcWidth,
		const int SrcHeight, int SrcScanlineSize, Tout* const NewBuf,
		const int NewWidth, const int NewHeight, const int ElCountIO,
		const double k, double ox = 0.0, double oy = 0.0,
		CImageResizerVars* const aVars = NULL ) const
	{
		if( SrcWidth == 0 || SrcHeight == 0 )
		{
			memset( NewBuf, 0, NewWidth * NewHeight * sizeof( Tout ));
			return;
		}
		else
		if( NewWidth == 0 || NewHeight == 0 )
		{
			return;
		}

		CImageResizerVars DefVars;
		CImageResizerVars& Vars = ( aVars == NULL ? DefVars : *aVars );

		CImageResizerThreadPool DefThreadPool;
		CImageResizerThreadPool& ThreadPool = ( Vars.ThreadPool == NULL ?
			DefThreadPool : *Vars.ThreadPool );

		double kx;
		double ky;

		if( k == 0.0 )
		{
			if( NewWidth > SrcWidth )
			{
				kx = (double) SrcWidth /
					( NewWidth + ( (double) NewWidth / SrcWidth - 1.0 ));
			}
			else
			{
				kx = (double) SrcWidth / NewWidth;
				ox += ( kx - 1.0 ) * 0.5;
			}

			if( NewHeight > SrcHeight )
			{
				ky = (double) SrcHeight /
					( NewHeight + ( (double) NewHeight / SrcHeight - 1.0 ));
			}
			else
			{
				ky = (double) SrcHeight / NewHeight;
				oy += ( ky - 1.0 ) * 0.5;
			}
		}
		else
		if( k > 0.0 )
		{
			kx = k;
			ky = k;
			const double ko = ( k >= 1.0 ? ( k - 1.0 ) : k ) * 0.5;
			ox += ko;
			oy += ko;
		}
		else
		{
			kx = -k;
			ky = -k;
		}

		// Evaluate pre-multipliers used on the output stage.

		const bool IsInFloat = ( (Tin) 0.4 != 0 );
		const bool IsOutFloat = ( (Tout) 0.4 != 0 );
		fptype OutMul; // Output multiplier.

		if( IsOutFloat )
		{
			OutMul = 1.0;
		}
		else
		{
			OutMul = ( sizeof( Tout ) == 1 ? 255.0 : 65535.0 );
		}

		if( !IsInFloat )
		{
			OutMul /= ( sizeof( Tin ) == 1 ? 255.0 : 65535.0 );
		}

		const int ElCount = ( ElCountIO + fpclass :: fppack - 1 ) /
			fpclass :: fppack;

		const int NewWidthE = NewWidth * ElCount;

		if( SrcScanlineSize < 1 )
		{
			SrcScanlineSize = SrcWidth * ElCountIO;
		}

		Vars.ElCount = ElCount;
		Vars.ElCountIO = ElCountIO;
		Vars.fppack = fpclass :: fppack;
		Vars.fpalign = fpclass :: fpalign;
		Vars.elalign = fpclass :: elalign;
		Vars.IntFltOrder = FilterBank.getOrder();
		Vars.IntFltLen = FilterBank.getFilterLen();

		// Horizontal scanline filtering and resizing.

		CFilterSteps FltSteps;
		buildFilterSteps( FltSteps, Vars, SrcWidth, NewWidth, kx, ox,
			OutMul );

		const int ThreadCount = ThreadPool.getSuggestedWorkloadCount();
			// Includes the current thread.

		CThreadData< Tin, Tout > td[ ThreadCount ];
		int i;

		for( i = 0; i < ThreadCount; i++ )
		{
			if( i > 0 )
			{
				ThreadPool.addWorkload( &td[ i ]);
			}

			td[ i ].init( i, ThreadCount, FltSteps, Vars );

			td[ i ].initScanlineQueue( td[ i ].sopResizeH, SrcHeight,
				SrcWidth );
		}

		CBuffer< fptype > FltBuf( NewWidthE * SrcHeight, fpclass :: fpalign );
			// Temporary buffer that receives horizontally-filtered and
			// resized image.

		for( i = 0; i < SrcHeight; i++ )
		{
			td[ i % ThreadCount ].addScanlineToQueue(
				(void*) &SrcBuf[ i * SrcScanlineSize ],
				&FltBuf[ i * NewWidthE ]);
		}

		ThreadPool.startAllWorkloads();
		td[ 0 ].processScanlineQueue();
		ThreadPool.waitAllWorkloadsToFinish();

		// Horizontal scanline filtering and resizing.

		buildFilterSteps( FltSteps, Vars, SrcHeight, NewHeight, ky, oy );

		if( IsOutFloat && sizeof( FltBuf[ 0 ]) == sizeof( Tout ) &&
			fpclass :: elalign == 1 )
		{
			// In-place output.

			for( i = 0; i < ThreadCount; i++ )
			{
				td[ i ].initScanlineQueue( td[ i ].sopResizeV, NewWidth,
					SrcHeight, NewWidthE, NewWidthE );
			}

			for( i = 0; i < NewWidth; i++ )
			{
				td[ i % ThreadCount ].addScanlineToQueue(
					&FltBuf[ i * ElCount ], (fptype*) &NewBuf[ i * ElCount ]);
			}

			ThreadPool.startAllWorkloads();
			td[ 0 ].processScanlineQueue();
			ThreadPool.waitAllWorkloadsToFinish();
			ThreadPool.removeAllWorkloads();

			return;
		}

		CBuffer< fptype > ResBuf( NewWidthE * NewHeight, fpclass :: fpalign );

		for( i = 0; i < ThreadCount; i++ )
		{
			td[ i ].initScanlineQueue( td[ i ].sopResizeV, NewWidth,
				SrcHeight, NewWidthE, NewWidthE );
		}

		for( i = 0; i < NewWidth; i++ )
		{
			td[ i % ThreadCount ].addScanlineToQueue(
				&FltBuf[ i * ElCount ], &ResBuf[ i * ElCount ]);
		}

		ThreadPool.startAllWorkloads();
		td[ 0 ].processScanlineQueue();
		ThreadPool.waitAllWorkloadsToFinish();

		if( IsOutFloat )
		{
			// Perform output, but skip dithering.

			for( i = 0; i < ThreadCount; i++ )
			{
				td[ i ].initScanlineQueue( td[ i ].sopUnpackH,
					NewHeight, NewWidth );
			}

			for( i = 0; i < NewHeight; i++ )
			{
				td[ i % ThreadCount ].addScanlineToQueue(
					&ResBuf[ i * NewWidthE ],
					&NewBuf[ i * NewWidth * ElCountIO ]);
			}

			ThreadPool.startAllWorkloads();
			td[ 0 ].processScanlineQueue();
			ThreadPool.waitAllWorkloadsToFinish();
			ThreadPool.removeAllWorkloads();
			return;
		}

		// Perform output with dithering (for integer output only).

		int TruncBits; // The number of lower bits to truncate and dither.
		fptype PkOut;

		if( sizeof( Tout ) == 1 )
		{
			TruncBits = 8 - ResBitDepth;
			PkOut = 255.0;
		}
		else
		{
			TruncBits = 16 - ResBitDepth;
			PkOut = 65535.0;
		}

		const fptype TrMul = (fptype) ( TruncBits > 0 ? 1 << TruncBits : 1 );

		if( CDitherer :: isRecursive() )
		{
			td[ 0 ].getDitherer().init( NewWidth, Vars, TrMul, PkOut );

			for( i = 0; i < NewHeight; i++ )
			{
				fptype* const ResScanline = &ResBuf[ i * NewWidthE ];

				td[ 0 ].getDitherer().dither( ResScanline );

				CFilterStep :: unpackScanline( ResScanline,
					&NewBuf[ i * NewWidth * ElCountIO ], NewWidth, Vars );
			}
		}
		else
		{
			for( i = 0; i < ThreadCount; i++ )
			{
				td[ i ].initScanlineQueue( td[ i ].sopDitherAndUnpackH,
					NewHeight, NewWidth );

				td[ i ].getDitherer().init( NewWidth, Vars, TrMul, PkOut );
			}

			for( i = 0; i < NewHeight; i++ )
			{
				td[ i % ThreadCount ].addScanlineToQueue(
					&ResBuf[ i * NewWidthE ],
					&NewBuf[ i * NewWidth * ElCountIO ]);
			}

			ThreadPool.startAllWorkloads();
			td[ 0 ].processScanlineQueue();
			ThreadPool.waitAllWorkloadsToFinish();
		}

		ThreadPool.removeAllWorkloads();
	}

private:
	typedef typename fpclass :: fptype fptype; ///< Floating-point type to use
		///< during processing.
		///<
	typedef typename fpclass :: CFilterStep CFilterStep; ///< Filter step
		///< class to use during processing.
		///<
	typedef typename fpclass :: CDitherer CDitherer; ///< Ditherer class to
		///< use during processing.
		///<
	CImageResizerParams Params; ///< Algorithm's parameters currently in use.
		///<
	int ResBitDepth; ///< Bit resolution of the resulting image.
		///<
	int FilterFracs; ///< The number of fractional delay filters sampled by
		///< the filter bank. This variable affects the signal-to-noise ratio
		///< at interpolation stage. Theoretically, 8-bit image resizing
		///< requires 66.2 dB SNR or 10. 16-bit resizing requires 114.4 dB SNR
		///< or 150.
		///<
	CDSPFracFilterBankLin< fptype > FilterBank; ///< Fractional delay filter
		///< bank.
		///<

	/**
	 * @brief Filtering steps array.
	 *
	 * The object of this class stores filtering steps together.
	 */

	typedef CStructArray< CFilterStep > CFilterSteps;

	/**
	 * Function allocates filter buffer taking "fpclass" alignments into
	 * account. The allocated buffer may be higher than the requested size: in
	 * this case the additional elements will be zeroed by this function.
	 *
	 * @param Flt Filter buffer.
	 * @param ReqCapacity The required filter buffer's capacity.
	 * @param FltExt If non-NULL this variable will receive the number of
	 * elements the filter was extended by.
	 */

	static void allocFilter( CBuffer< fptype >& Flt, const int ReqCapacity,
		int* const FltExt = NULL )
	{
		int UseCapacity = ( ReqCapacity + fpclass :: elalign - 1 ) &
			~( fpclass :: elalign - 1 );

		Flt.alloc( UseCapacity, fpclass :: fpalign );
		int Ext = UseCapacity - ReqCapacity;

		if( FltExt != NULL )
		{
			*FltExt = Ext;
		}

		while( Ext > 0 )
		{
			Ext--;
			Flt[ ReqCapacity + Ext ] = (fptype) 0.0;
		}
	}

	/**
	 * Function assigns filter parameters to the specified filter step object.
	 *
	 * @param fs Filter step to assign parameter to. This step cannot be the
	 * last step if ResampleFactor greater than 1 was specified.
	 * @param IsUpsample "True" if upsampling step. Should be set to "false"
	 * if FltCutoff is negative.
	 * @param ResampleFactor Resampling factor of this filter (>=1).
	 * @param FltCutoff Filter cutoff point. If zero value was specified,
	 * the "half-band" predefined filter will be created. If a negative value
	 * was specified, filter's kernel will not be calculated, and a
	 * pre-calculated kernel is assumed to be available. This value will be
	 * divided by the ResampleFactor if IsUpsample equals "true".
	 * @param k Resizing coefficient at the previous processing step, may be
	 * adjusted on return.
	 * @param o Starting pixel offset inside the source image, may be adjusted
	 * on return.
	 * @param prevfs Previous parameters to derive buffer index and length
	 * from. If NULL was specified, then SrcLen and Vars should be
	 * provided explicitly.
	 * @param SrcLen Source scanline length in pixels.
	 * @param Vars Image resizing-related variables. If prevfs is NULL and
	 * Vars is NULL it means that Vars, InLen, InBuf, OutLen and OutBuf
	 * variables were pre-assigned.
	 * @param ElCount The number of elements (channels) in each pixel.
	 * @param DCGain DC gain to apply to the filter. Assigned to filter step's
	 * DCGain variable.
	 * @return Pointer to the "fs" object.
	 */

	CFilterStep* assignFilterParams( CFilterStep& fs, const bool IsUpsample,
		const int ResampleFactor, const double FltCutoff, double& k,
		double& o, const CFilterStep* const prevfs, const int SrcLen,
		const CImageResizerVars* const Vars, const double DCGain ) const
	{
		int FltExt; // Filter's extension due to fpclass :: elalign.

		if( FltCutoff >= 0.0 )
		{
			double FltAlpha;
			double Len2;
			double Freq;

			if( FltCutoff == 0.0 )
			{
				FltAlpha = Params.HBFltAlpha;
				Len2 = 0.5 * Params.HBFltLen;
				Freq = M_PI * Params.HBFltCutoff;
			}
			else
			{
				FltAlpha = Params.LPFltAlpha;
				Len2 = 0.25 * Params.LPFltBaseLen / FltCutoff;
				Freq = M_PI * Params.LPFltCutoffMult * FltCutoff;
			}

			if( IsUpsample )
			{
				Len2 *= ResampleFactor;
				Freq /= ResampleFactor;
				fs.DCGain = DCGain * ResampleFactor;
			}
			else
			{
				fs.DCGain = DCGain;
			}

			CDSPPeakedCosineLPF w( Len2, Freq, FltAlpha );

			CBuffer< double > Kernel( w.KernelLen );
			fs.FltLatency = w.fl2;
			w.generateLPF( &Kernel[ 0 ], 1.0 );
			optimizeFIRFilter( Kernel, fs.FltLatency );
			normalizeFIRFilter( &Kernel[ 0 ], Kernel.getCapacity(),
				fs.DCGain );

			allocFilter( fs.Flt, Kernel.getCapacity(), &FltExt );
			copyArray( &Kernel[ 0 ], &fs.Flt[ 0 ], Kernel.getCapacity() );
		}
		else
		{
			fs.DCGain = DCGain;
			FltExt = 0;
		}

		fs.IsUpsample = IsUpsample;
		fs.ResampleFactor = ResampleFactor;

		if( prevfs == NULL )
		{
			if( Vars != NULL )
			{
				fs.Vars = Vars;
				fs.InLen = SrcLen;
				fs.InBuf = 0;
			}
		}
		else
		{
			fs.Vars = prevfs -> Vars;
			fs.InLen = prevfs -> OutLen;
			fs.InBuf = prevfs -> OutBuf;
		}

		if( IsUpsample )
		{
			fs.InPrefix = 0;
			fs.InSuffix = 0;
			fs.OutLen = fs.InLen * ResampleFactor;
			fs.OutPrefix = fs.FltLatency;
			fs.OutSuffix = fs.Flt.getCapacity() - fs.FltLatency -
				ResampleFactor;

			k *= ResampleFactor;
			o *= ResampleFactor;

			int l = fs.OutSuffix - FltExt;
			const fptype* ip = &fs.Flt[ fs.FltLatency + ResampleFactor ];
			allocFilter( fs.PrefixDC, l );
			copyArray( ip, &fs.PrefixDC[ 0 ], l );

			while( true )
			{
				ip += ResampleFactor;
				l -= ResampleFactor;

				if( l <= 0 )
				{
					break;
				}

				addArray( ip, &fs.PrefixDC[ 0 ], l );
			}

			l = fs.FltLatency;
			allocFilter( fs.SuffixDC, l );
			fptype* op = &fs.SuffixDC[ 0 ];
			copyArray( &fs.Flt[ 0 ], op, l );

			while( true )
			{
				op += ResampleFactor;
				l -= ResampleFactor;

				if( l <= 0 )
				{
					break;
				}

				addArray( &fs.Flt[ 0 ], op, l );
			}

			int l0 = fs.OutPrefix + fs.OutLen + fs.OutSuffix;
			l = fs.InLen * ResampleFactor + fs.SuffixDC.getCapacity();

			if( l > l0 )
			{
				fs.OutSuffix += l - l0;
			}

			l0 = fs.OutLen + fs.OutSuffix;

			if( fs.PrefixDC.getCapacity() > l0 )
			{
				fs.OutSuffix += fs.PrefixDC.getCapacity() - l0;
			}
		}
		else
		{
			fs.InPrefix = fs.FltLatency;
			fs.InSuffix = fs.Flt.getCapacity() - fs.FltLatency - 1;

			if( ResampleFactor > 1 )
			{
				fs.InSuffix += ResampleFactor - fs.InLen % ResampleFactor;
				fs.OutLen = fs.InLen / ResampleFactor + 1;
				k /= ResampleFactor;
				o /= ResampleFactor;
			}
			else
			{
				fs.OutLen = fs.InLen;
			}
		}

		fs.OutBuf = ( fs.InBuf + 1 ) & 1;

		return( &fs );
	}

	/**
	 * Function adds a correction filter that tries to achieve a linear
	 * frequency response at all frequencies. The actual resulting response
	 * may feature a slight damping of the highest frequencies since a
	 * suitably short correction filter cannot fix steep high-frequency
	 * damping.
	 *
	 * This function assumes that the resizing step is currently the last
	 * step, even if it was not inserted yet: this allows placement of the
	 * correction filter both before and after the resizing step.
	 *
	 * @param Steps Filtering steps.
	 * @param bw Resulting bandwidth relative to the original bandwidth (which
	 * is 1.0), usually 1/k. Should be <= 1.0.
	 * @param k Resizing coefficient at the last processing step, may be
	 * adjusted on return.
	 * @param o Starting pixel offset inside the source image, may be
	 * adjusted on return.
	 * @param IsPreCorrection "True" if the filter step was already created
	 * and it is first in the Steps array.
	 */

	void addCorrectionFilter( CFilterSteps& Steps, const double bw, double& k,
		double& o, const bool IsPreCorrection ) const
	{
		const int BinCount = 65; // Frequency response bins to control.
		const int BinCount1 = BinCount - 1;
		double curbw = 1.0; // Bandwidth of the filter at the current step.
		int i;
		int j;
		double re;
		double im;

		CBuffer< double > Bins( BinCount ); // Adjustment introduced by all
			// steps at all frequencies of interest.

		for( j = 0; j < BinCount; j++ )
		{
			Bins[ j ] = 1.0;
		}

		for( i = ( IsPreCorrection ? 1 : 0 ); i < Steps.getItemCount(); i++ )
		{
			const CFilterStep& fs = Steps[ i ];

			if( fs.ResampleFactor == 0 )
			{
				break;
			}

			const double dcg = 1.0 / fs.DCGain; // DC gain correction.

			if( fs.IsUpsample )
			{
				curbw *= fs.ResampleFactor;
			}

			// Calculate frequency response adjustment introduced by the
			// filter at this step, within the bounds of bandwidth of
			// interest.

			for( j = 0; j < BinCount; j++ )
			{
				const double th = M_PI * bw / curbw * j / BinCount1;

				calcFIRFilterResponse( &fs.Flt[ 0 ], fs.Flt.getCapacity(), th,
					re, im );

				Bins[ j ] /= sqrt( re * re + im * im ) * dcg;
			}

			if( !fs.IsUpsample && fs.ResampleFactor > 1 )
			{
				curbw /= fs.ResampleFactor;
			}
		}

		// Add correction for the resizing step, even if this step was not
		// added yet. Resizing step is accompanied by a filter stored in a
		// filter bank.

		for( j = 0; j < BinCount; j++ )
		{
			const double th = M_PI * bw / curbw * j / BinCount1;

			calcFIRFilterResponse( &FilterBank[ 0 ],
				FilterBank.getFilterLen(), th, re, im );

			Bins[ j ] /= sqrt( re * re + im * im );
		}

		// Calculate filter kernel.

		CDSPFIREQ EQ;
		EQ.init( bw * 2.0, Params.CorrFltLen, BinCount, 0.0, bw, false,
			Params.CorrFltAlpha );

		const CFilterStep* prevfs;
		CFilterStep* fs;

		if( IsPreCorrection )
		{
			prevfs = NULL;
			fs = &Steps[ 0 ];
		}
		else
		{
			prevfs = &Steps[ Steps.getItemCount() - 1 ];
			fs = &Steps.add();
		}

		fs -> FltLatency = EQ.getKernelLatency();

		CBuffer< double > Kernel( EQ.getKernelLength() );
		EQ.buildKernel( Bins, &Kernel[ 0 ]);
		normalizeFIRFilter( &Kernel[ 0 ], Kernel.getCapacity(), 1.0 );
		optimizeFIRFilter( Kernel, fs -> FltLatency );
		normalizeFIRFilter( &Kernel[ 0 ], Kernel.getCapacity(), 1.0 );

		allocFilter( fs -> Flt, Kernel.getCapacity() );
		copyArray( &Kernel[ 0 ], &fs -> Flt[ 0 ], Kernel.getCapacity() );

		// Print a theoretically achieved final frequency response at various
		// feature sizes (from DC to 1 pixel). Values above 255 means features
		// become brighter, values below 255 means features become dimmer.

/*		const double sbw = ( bw > 1.0 ? 1.0 / bw : 1.0 );

		for( j = 0; j < BinCount; j++ )
		{
			const double th = M_PI * sbw * j / BinCount1;

			calcFIRFilterResponse( &fs -> Flt[ 0 ], fs -> Flt.getCapacity(),
				th, re, im );

			printf( "%f\n", sqrt( re * re + im * im ) / Bins[ j ] * 255 );
		}

		printf( "***\n" );
*/
		assignFilterParams( *fs, false, 1, -1.0, k, o, prevfs, 0, NULL, 1.0 );
	}

	/**
	 * Function adds a sharpening filter if image is being upsized. Such
	 * sharpening allows to spot interpolation filter's stop-band attenuation:
	 * if attenuation is too weak, a "dark grid" artifact may become visible.
	 *
	 * It is assumed that 40 decibel stop-band attenuation should be
	 * considered a required minimum: this allows application of 64X
	 * sharpening without spotting any artifacts.
	 *
	 * @param Steps Filtering steps.
	 * @param bw Resulting bandwidth relative to the original bandwidth (which
	 * is 1.0), usually 1/k.
	 * @param k Resizing coefficient at the last processing step, may be
	 * adjusted on return.
	 * @param o Starting pixel offset inside the source image, may be
	 * adjusted on return.
	 */

	void addSharpenTest( CFilterSteps& Steps, const double bw, double& k,
		double& o ) const
	{
		if( bw < 1.0 )
		{
			return;
		}

		const int BinCount = 200;
		CBuffer< double > Bins( BinCount );
		int Thresh = (int) round( BinCount / bw * 1.75 );

		if( Thresh > BinCount )
		{
			Thresh = BinCount;
		}

		int j;

		for( j = 0; j < Thresh; j++ )
		{
			Bins[ j ] = 1.0;
		}

		for( j = Thresh; j < BinCount; j++ )
		{
			Bins[ j ] = 256.0;
		}

		CDSPFIREQ EQ;
		EQ.init( bw * 2.0, 10 * bw, BinCount, 0.0, bw, false, 1.7 );

		const CFilterStep* prevfs = &Steps[ Steps.getItemCount() - 1 ];
		CFilterStep* fs = &Steps.add();

		fs -> FltLatency = EQ.getKernelLatency();

		CBuffer< double > Kernel( EQ.getKernelLength() );
		EQ.buildKernel( Bins, &Kernel[ 0 ]);
		normalizeFIRFilter( &Kernel[ 0 ], Kernel.getCapacity(), 1.0 );
		optimizeFIRFilter( Kernel, fs -> FltLatency );
		normalizeFIRFilter( &Kernel[ 0 ], Kernel.getCapacity(), 1.0 );

		allocFilter( fs -> Flt, Kernel.getCapacity() );
		copyArray( &Kernel[ 0 ], &fs -> Flt[ 0 ], Kernel.getCapacity() );

/*		for( j = 0; j < BinCount; j++ )
		{
			const double th = M_PI * j / ( BinCount - 1 );
			double re;
			double im;

			calcFIRFilterResponse( &fs -> Flt[ 0 ], fs -> Flt.getCapacity(),
				th, re, im );

			printf( "%f\n", sqrt( re * re + im * im ));
		}

		printf( "***\n" );
*/
		assignFilterParams( *fs, false, 1, -1.0, k, o, prevfs, 0, NULL, 1.0 );
	}

	/**
	 * Function calculates an optimal intermediate buffer length that will
	 * cover all needs of the specified filtering steps.
	 *
	 * @param Steps Filtering steps.
	 * @param[out] Vars Variables object, will receive buffer size and length.
	 */

	static void calcBufLen( CFilterSteps& Steps, CImageResizerVars& Vars )
	{
		int MaxPrefix = 0;
		int MaxLen = 0;
		int i;

		for( i = 0; i < Steps.getItemCount(); i++ )
		{
			const CFilterStep& fs = Steps[ i ];

			int l = fs.InPrefix;

			if( l > MaxPrefix )
			{
				MaxPrefix = l;
			}

			l = fs.InLen + fs.InSuffix;

			if( l > MaxLen )
			{
				MaxLen = l;
			}

			if( fs.IsUpsample )
			{
				l = fs.OutPrefix;

				if( l > MaxPrefix )
				{
					MaxPrefix = l;
				}

				l = fs.OutLen + fs.OutSuffix;

				if( l > MaxLen )
				{
					MaxLen = l;
				}
			}
			else
			{
				l = fs.OutLen;

				if( l > MaxLen )
				{
					MaxLen = l;
				}
			}
		}

		Vars.BufLen = MaxPrefix + MaxLen;
		Vars.BufOffset = MaxPrefix;

		if( Vars.elalign == 1 )
		{
			Vars.BufOffset *= Vars.ElCount;
		}
		else
		{
			Vars.BufIncr = Vars.BufLen;
		}

		Vars.BufLen *= Vars.ElCount;
	}

	/**
	 * Function builds sequence of filtering steps depending on the specified
	 * resizing coefficient. The last steps included are always the resizing
	 * step then the correction step.
	 *
	 * @param Steps Array that receives filtering steps.
	 * @param[out] Vars Variables object, will receive buffer size and length.
	 * @param SrcLen Source scanline's length in pixels.
	 * @param NewLen New scanline's length in pixels.
	 * @param k Resizing coefficient.
	 * @param o Starting pixel offset inside the source image.
	 * @param DCGain The overall DC gain to apply. This DC gain is applied to
	 * the first filter step only (upsampling or filtering step).
	 */

	void buildFilterSteps( CFilterSteps& Steps, CImageResizerVars& Vars,
		const int SrcLen, const int NewLen, double k, double o,
		const double DCGain = 1.0 ) const
	{
		Steps.clear();

		const double bw = 1.0 / k; // Resulting bandwidth.
		bool IsPreCorrection; // "True" if the correction filter is applied
			// first.
		double FltCutoff = ( k <= 1.0 ? 1.0 : 1.0 / k ); // Cutoff frequency
			// of the first filter step.
		const int UpsampleFactor = ( (int) floor( k ) < 2 ? 2 : 1 );
		CFilterStep* prevfs;

		// Add 1 upsampling or several downsampling filters.

		if( UpsampleFactor > 1 )
		{
			IsPreCorrection = true;
			prevfs = &Steps.add();
			prevfs -> Vars = &Vars;
			prevfs -> InLen = SrcLen;
			prevfs -> InBuf = 0;
			prevfs -> OutLen = SrcLen;
			prevfs -> OutBuf = 1;

			prevfs = assignFilterParams( Steps.add(), true, UpsampleFactor,
				FltCutoff, k, o, prevfs, 0, NULL, DCGain );
		}
		else
		{
			IsPreCorrection = false;
			prevfs = NULL;
			double UseDCGain = DCGain;

			while( true )
			{
				if( FltCutoff > 0.125 )
				{
					prevfs = assignFilterParams( Steps.add(), false,
						( FltCutoff > 0.25 ? 1 : 2 ), FltCutoff, k, o,
						prevfs, SrcLen, &Vars, UseDCGain );

					break;
				}

				prevfs = assignFilterParams( Steps.add(), false, 2, 0.0, k, o,
					prevfs, SrcLen, &Vars, UseDCGain );

				FltCutoff *= 2.0;
				UseDCGain = 1.0;
			}
		}

		// Insert resizing and correction steps.

		CFilterStep& fs = Steps.add();
		fs.IsUpsample = false;
		fs.ResampleFactor = 0;
		fs.Vars = &Vars;
		fs.InLen = prevfs -> OutLen;
		fs.InBuf = prevfs -> OutBuf;

		const int FilterLenD2 = FilterBank.getFilterLen() / 2;
		const int FilterLenD21 = FilterLenD2 - 1;

		const int ResizeLPix = (int) floor( o ) - FilterLenD21;
		fs.InPrefix = ( ResizeLPix < 0 ? -ResizeLPix : 0 );
		const int ResizeRPix = (int) floor( o + ( NewLen - 1 ) * k ) +
			FilterLenD2 + 1;

		fs.InSuffix = ( ResizeRPix > fs.InLen ? ResizeRPix - fs.InLen : 0 );
		fs.OutLen = NewLen;
		fs.OutBuf = ( fs.InBuf + 1 ) & 1;

		// Fill resizing positions buffer.

		fs.RPosBuf.alloc( fs.OutLen );
		typename CFilterStep :: CResizePos* rpos = &fs.RPosBuf[ 0 ];
		int i;
		const int em = ( fpclass :: elalign == 1 ? Vars.ElCount : 1 );

		for( i = 0; i < fs.OutLen; i++ )
		{
			const double SrcPos = o + k * i;
			const int SrcPosInt = (int) floor( SrcPos );
			double x = ( SrcPos - SrcPosInt ) * FilterFracs;
			const int fti = (int) x;
			rpos -> x = (typename fpclass :: fptypeatom) ( x - fti );
			rpos -> ftp = &FilterBank[ fti ];
			rpos -> SrcOffs = ( SrcPosInt - FilterLenD21 ) * em;
			rpos++;
		}

		addCorrectionFilter( Steps, ( bw >= 1.0 ? 1.0 : bw ), k, o,
			IsPreCorrection );

		//addSharpenTest( Steps, bw, k, o );

		CFilterStep& lastfs = Steps[ Steps.getItemCount() - 1 ];
		lastfs.OutBuf = 2;

		if( lastfs.OutLen > NewLen )
		{
			lastfs.OutLen = NewLen;
		}

		if( UpsampleFactor > 1 )
		{
			Steps[ 1 ].extendUpsample( Steps[ 2 ]);
		}

		calcBufLen( Steps, Vars );
	}

	/**
	 * @brief Thread-isolated data used for scanline processing.
	 *
	 * This structure holds data necessary for image's horizontal or vertical
	 * scanline processing, including scanline processing queue.
	 *
	 * @tparam Tin Source element data type. Intermediate buffers store data
	 * in floating point format.
	 * @tparam Tout Destination element data type. Intermediate buffers store
	 * data in floating point format.
	 */

	template< class Tin, class Tout >
	class CThreadData : public CImageResizerThreadPool :: CWorkload
	{
	public:
		virtual void process()
		{
			processScanlineQueue();
		}

		/**
		 * This enumeration lists possible scanline operations.
		 */

		enum EScanlineOperation
		{
			sopResizeH, ///< Resize horizontal scanline.
			sopResizeV, ///< Resize vertical scanline.
			sopDitherAndUnpackH, ///< Dither and unpack scanline.
			sopUnpackH ///< Unpack scanline.
		};

		/**
		 * Function initializes *this thread data object and assigns certain
		 * variables provided by the higher level code.
		 *
		 * @param aThreadIndex Index of this thread data (0-based).
		 * @param aThreadCount Total number of threads used during processing.
		 * @param aSteps Filtering steps.
		 * @param aVars Image resizer variables.
		 */

		void init( const int aThreadIndex, const int aThreadCount,
			const CFilterSteps& aSteps, const CImageResizerVars& aVars )
		{
			ThreadIndex = aThreadIndex;
			ThreadCount = aThreadCount;
			Steps = &aSteps;
			Vars = &aVars;
		}

		/**
		 * Function initializes scanline processing queue, and updates
		 * capacities of intermediate buffers.
		 *
		 * @param aOp Operation to perform over scanline.
		 * @param TotalLines The total number of scanlines that will be
		 * processed by all threads.
		 * @param aSrcLen Source scanline length in pixels.
		 * @param aSrcIncr Source scanline buffer increment. Ignored in
		 * horizontal scanline processing.
		 * @param aResIncr Resulting scanline buffer increment. Ignored in
		 * horizontal scanline processing.
		 */

		void initScanlineQueue( const EScanlineOperation aOp,
			const int TotalLines, const int aSrcLen, const int aSrcIncr = 0,
			const int aResIncr = 0 )
		{
			Bufs.alloc( Vars -> BufLen * 2, fpclass :: fpalign );
			BufPtrs[ 0 ] = Bufs + Vars -> BufOffset;
			BufPtrs[ 1 ] = BufPtrs[ 0 ] + Vars -> BufLen;

			ScanlineOp = aOp;
			SrcLen = aSrcLen;
			SrcIncr = aSrcIncr;
			ResIncr = aResIncr;
			Queue.alloc(( TotalLines + ThreadCount - 1 ) / ThreadCount );
			QueueLen = 0;
		}

		/**
		 * Function adds a scanline to the queue buffer. The
		 * initScanlineQueue() function should be called before calling this
		 * function. The number of calls to this add function should not
		 * exceed the TotalLines spread over all threads.
		 *
		 * @param SrcBuf Source scanline buffer.
		 * @param ResBuf Resulting scanline buffer.
		 */

		void addScanlineToQueue( void* const SrcBuf, void* const ResBuf )
		{
			Queue[ QueueLen ].SrcBuf = SrcBuf;
			Queue[ QueueLen ].ResBuf = ResBuf;
			QueueLen++;
		}

		/**
		 * Function processes all queued scanlines.
		 */

		void processScanlineQueue()
		{
			int i;

			switch( ScanlineOp )
			{
				case sopResizeH:
				{
					for( i = 0; i < QueueLen; i++ )
					{
						resizeScanlineH( (Tin*) Queue[ i ].SrcBuf,
							(fptype*) Queue[ i ].ResBuf );
					}

					break;
				}

				case sopResizeV:
				{
					for( i = 0; i < QueueLen; i++ )
					{
						resizeScanlineV( (fptype*) Queue[ i ].SrcBuf,
							(fptype*) Queue[ i ].ResBuf );
					}

					break;
				}

				case sopDitherAndUnpackH:
				{
					for( i = 0; i < QueueLen; i++ )
					{
						Ditherer.dither( (fptype*) Queue[ i ].SrcBuf );

						CFilterStep :: unpackScanline(
							(fptype*) Queue[ i ].SrcBuf,
							(Tout*) Queue[ i ].ResBuf, SrcLen, *Vars );
					}

					break;
				}

				case sopUnpackH:
				{
					for( i = 0; i < QueueLen; i++ )
					{
						CFilterStep :: unpackScanline(
							(fptype*) Queue[ i ].SrcBuf,
							(Tout*) Queue[ i ].ResBuf, SrcLen, *Vars );
					}

					break;
				}
			}
		}

		/**
		 * Function returns ditherer object associated with *this thread data
		 * object.
		 */

		CDitherer& getDitherer()
		{
			return( Ditherer );
		}

	private:
		int ThreadIndex; ///< Thread index.
			///<
		int ThreadCount; ///< Thread count.
			///<
		const CFilterSteps* Steps; ///< Filtering steps.
			///<
		const CImageResizerVars* Vars; ///< Image resizer variables.
			///<
		CBuffer< fptype > Bufs; ///< Flip-flop intermediate buffers.
			///<
		fptype* BufPtrs[ 3 ]; ///< Flip-flop buffer pointers (referenced by
			///< filtering step's InBuf and OutBuf indices).
			///<
		EScanlineOperation ScanlineOp; ///< Operation to perform over
			///< scanline.
			///<
		int SrcLen; ///< Source scanline length in the last queue.
			///<
		int SrcIncr; ///< Source scanline buffer increment in the last queue.
			///<
		int ResIncr; ///< Resulting scanline buffer increment in the last
			///< queue.
			///<
		CDitherer Ditherer; ///< Ditherer object to use.
			///<

		/**
		 * @brief Scanline processing queue item.
		 *
		 * Scanline processing queue item.
		 */

		struct CQueueItem
		{
			void* SrcBuf; ///< Source scanline buffer, will by typecasted to
				///< Tin or fptype*.
				///<
			void* ResBuf; ///< Resulting scanline buffer.
				///<
		};

		CBuffer< CQueueItem > Queue; ///< Scanline processing queue.
			///<
		int QueueLen; ///< Queue length.
			///<

		/**
		 * Function resizes a single horizontal scanline.
		 *
		 * @param SrcBuf Source scanline buffer. Can be either horizontal or
		 * vertical.
		 * @param ResBuf Resulting scanline buffer.
		 */

		void resizeScanlineH( const Tin* const SrcBuf, fptype* const ResBuf )
		{
			CFilterStep :: packScanline( SrcBuf, BufPtrs[ 0 ], SrcLen,
				*Vars );

			BufPtrs[ 2 ] = ResBuf;
			int j;

			for( j = 0; j < Steps -> getItemCount(); j++ )
			{
				const CFilterStep& fs = (*Steps)[ j ];
				fs.prepareInBuf( BufPtrs[ fs.InBuf ]);

				if( fs.ResampleFactor != 0 )
				{
					if( fs.IsUpsample )
					{
						fs.doUpsample( BufPtrs[ fs.InBuf ],
							BufPtrs[ fs.OutBuf ]);
					}
					else
					{
						fs.doFilter( BufPtrs[ fs.InBuf ],
							BufPtrs[ fs.OutBuf ], Vars -> ElCount );
					}
				}
				else
				{
					fs.doResize( BufPtrs[ fs.InBuf ], BufPtrs[ fs.OutBuf ],
						Vars -> ElCount );
				}
			}
		}

		/**
		 * Function resizes a single vertical scanline.
		 *
		 * @param SrcBuf Source scanline buffer. Can be either horizontal or
		 * vertical.
		 * @param ResBuf Resulting scanline buffer.
		 */

		void resizeScanlineV( const fptype* const SrcBuf,
			fptype* const ResBuf )
		{
			CFilterStep :: convertVtoH( SrcBuf, BufPtrs[ 0 ], SrcLen, SrcIncr,
				*Vars );

			BufPtrs[ 2 ] = ResBuf;
			int j;

			for( j = 0; j < Steps -> getItemCount(); j++ )
			{
				const CFilterStep& fs = (*Steps)[ j ];
				fs.prepareInBuf( BufPtrs[ fs.InBuf ]);
				const int DstIncr = ( fs.OutBuf == 2 ?
					ResIncr : Vars -> ElCount );

				if( fs.ResampleFactor != 0 )
				{
					if( fs.IsUpsample )
					{
						fs.doUpsample( BufPtrs[ fs.InBuf ],
							BufPtrs[ fs.OutBuf ]);
					}
					else
					{
						fs.doFilter( BufPtrs[ fs.InBuf ],
							BufPtrs[ fs.OutBuf ], DstIncr );
					}
				}
				else
				{
					fs.doResize( BufPtrs[ fs.InBuf ], BufPtrs[ fs.OutBuf ],
						DstIncr );
				}
			}
		}
	};
};

} // namespace avir

#endif // AVIR_CIMAGERESIZER_INCLUDED
