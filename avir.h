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
 * @section license License
 *
 * ! This License Agreement is applicable to non-commercial and low-volume
 * commercial use only. You (individual person or a company) have to sign a
 * special license agreement with AVIR copyright holder if you or your company
 * sells and/or licenses (on a commercial basis, directly or through
 * resellers) software products where AVIR or its parts are included, and if
 * combined sales of such software products reach more than 1000 (one
 * thousand) end user units per year (be it individual or employee users) on
 * desktop and laptop computers and any other devices (e.g. TVs, video
 * projectors, video players, etc), or 5000 (five thousand) end user units per
 * year on smartphones and tablets.
 *
 * If you sublicense to a non-end user third party, or transfer rights to your
 * software where AVIR is used, to a third party, this third party should
 * agree to this License Agreement or sign a special license agreement with
 * AVIR copyright holder, or remove AVIR dependence from the received
 * software.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * Please credit the creator of this library in your documentation in the
 * following way: "AVIR image resizing algorithm designed by Aleksey Vaneev"
 *
 * @version 1.0
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
 * @param fltlat Filter's latency in samples.
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
 * construction (initialization) operations.
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
		, Capacity( 0 )
	{
	}

	/**
	 * Constructor creates the buffer with the specified capacity.
	 *
	 * @param aCapacity Buffer's capacity.
	 */

	CBuffer( const int aCapacity )
	{
		allocinit( aCapacity );
	}

	CBuffer( const CBuffer& Source )
	{
		allocinit( Source.Capacity );
		memcpy( Data, Source.Data, Capacity * sizeof( T ));
	}

	~CBuffer()
	{
		freeData();
	}

	CBuffer& operator = ( const CBuffer& Source )
	{
		alloc( Source.Capacity );
		memcpy( Data, Source.Data, Capacity * sizeof( T ));
		return( *this );
	}

	/**
	 * Function allocates memory so that the specified number of elements
	 * can be stored in *this buffer object.
	 *
	 * @param aCapacity Storage for this number of elements to allocate.
	 */

	void alloc( const int aCapacity )
	{
		freeData();
		allocinit( aCapacity );
	}

	/**
	 * Function deallocates any previously allocated buffer.
	 */

	void free()
	{
		freeData();
		Data = NULL;
		Capacity = 0;
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
	 * performed.
	 *
	 * @param NewCapacity New (increased) capacity.
	 */

	void increaseCapacity( const int NewCapacity )
	{
		if( NewCapacity < Capacity )
		{
			return;
		}

		Data = (T*) :: realloc( Data, NewCapacity * sizeof( T ));
		Capacity = NewCapacity;
	}

	/**
	 * Function "truncates" (reduces) capacity of the buffer without
	 * reallocating it.
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
	 * number of elements is available.
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
		return( Data );
	}

private:
	T* Data; ///< Element buffer pointer.
		///<
	int Capacity; ///< Element buffer capacity.
		///<

	/**
	 * Internal element buffer allocation function used during object
	 * construction.
	 *
	 * @param aCapacity Storage for this number of elements to allocate.
	 */

	void allocinit( const int aCapacity )
	{
		Data = (T*) :: malloc( aCapacity * sizeof( T ));
		Capacity = aCapacity;
	}

	/**
	 * Function frees a previously allocated Data buffer.
	 */

	void freeData()
	{
		:: free( Data );
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
	 * @param aLen2 Half filter length (non-truncated).
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
	double Len2; ///< Half length of the window.
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
	 * @param aKernelLength Required kernel length in samples. The actual
	 * kernel length is truncated to the nearest integer.
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
	 * @return Kernel's length, in samples.
	 */

	int getKernelLength() const
	{
		return( z2 - 1 );
	}

	/**
	 * @return Kernel's output latency, in samples.
	 */

	int getKernelLatency() const
	{
		return( z - 1 );
	}

	/**
	 * Function creates FIR kernel with the specified gain levels at band
	 * crossover points.
	 *
	 * @param BandGains Array of gain levels, count=BandCount specified in
	 * the init() function.
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
 * windowed by the Peaked Cosine window function, for image processing
 * applications.
 */

class CDSPPeakedCosineLPF
{
public:
	int fl2; ///< Half filter's length, excluding the peak value. This value
		///< can be also used as filter's latency in samples (taps).
		///<
	int KernelLen; ///< Filter kernel length in samples (taps).
		///<

	/**
	 * Constructor initalizes *this object.
	 *
	 * @param aLen2 Half-length (non-truncated) of low-pass filter, in
	 * samples.
	 * @param aFreq2 Low-pass filter's corner frequency [0; pi].
	 * @param aAlpha Peaked Cosine window function Alpha parameter.
	 */

	CDSPPeakedCosineLPF( const double aLen2, const double aFreq2,
		const double aAlpha )
		: fl2( (int) floor( aLen2 ))
		, KernelLen( fl2 + fl2 + 1 )
		, Len2( aLen2 )
		, Freq2( aFreq2 )
		, Alpha( aAlpha )
	{
	}

	/**
	 * Function generates a linear-phase low-pass filter windowed by Peaked
	 * Cosine window function.
	 *
	 * @param[out] op Output buffer, length = fl2 * 2 + 1.
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
		///< samples.
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
 * low-pass filter.
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
	 * @param BaseLen Low-pass filter's base length, in samples.
	 * @param Cutoff Low-pass filter's normalized cutoff frequency.
	 * @param WFAlpha Peaked Cosine window function's Alpha parameter.
	 */

	void init( const int ReqFracCount, const double BaseLen,
		const double Cutoff, const double WFAlpha )
	{
		CDSPPeakedCosineLPF p( 0.5 * BaseLen * ReqFracCount,
			M_PI * Cutoff / ReqFracCount, WFAlpha );

		const int ElementSize = 1 + 1;
		FilterLen = ( p.fl2 / ReqFracCount + 1 ) * 2;
		FilterSize = FilterLen * ElementSize;

		const int BufLen = FilterLen * ReqFracCount + InterpPoints - 1;
		const int BufOffs = InterpPoints / 2 - 1;
		const int BufCenter = FilterLen * ReqFracCount / 2 + BufOffs;

		CBuffer< double > Buf( BufLen );
		memset( Buf, 0, ( BufCenter - p.fl2 ) * sizeof( double ));
		int i = BufLen - BufCenter - p.fl2 - 1;
		memset( &Buf[ BufLen - i ], 0, i * sizeof( double ));

		p.generateLPF( &Buf[ BufCenter - p.fl2 ], ReqFracCount );

		Table.alloc( ReqFracCount * FilterSize );
		fptype* op = Table;
		int j;

		for( i = ReqFracCount; i > 0; i-- )
		{
			double* p = Buf + BufOffs + i;

			for( j = 0; j < FilterLen; j++ )
			{
				op[ 0 ] = (fptype) p[ 0 ];
				op[ 1 ] = (fptype) ( p[ -1 ] - p[ 0 ]);
				op += ElementSize;
				p += ReqFracCount;
			}
		}
	}

	/**
	 * Operator returns the pointer to the specified interpolation table
	 * filter.
	 *
	 * @param i Filter (fractional delay) index, in the range 0 to
	 * ReqFracCount - 1, inclusive.
	 * @return Reference to the filter elements.
	 */

	const fptype& operator []( const int i ) const
	{
		return( Table[ i * FilterSize ]);
	}

	/**
	 * @return The length of each filter, in samples.
	 */

	int getFilterLen() const
	{
		return( FilterLen );
	}

private:
	static const int InterpPoints = 2; ///< The maximal number of points the
		///< interpolation is based on.
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
 * The image resizer algorithm makes calls to functions of this class.
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
 */

struct CImageResizerParams
{
	double CorrFltAlpha; ///< Alpha parameter of the Peaked Cosine window
		///< function used on the correction filter. The "usable" values are
		///< above 1.0.
		///<
	double CorrFltBaseLen; ///< Correction filter's base length in samples
		///< (taps). When upsizing is performed, this value is multiplied by
		///< the upsizing factor. The "usable" range is narrow, 6 to 8, as to
		///< minimize the "overcorrection" which is mathematically precise,
		///< but visually unacceptable.
		///<
	double IntFltAlpha; ///< Alpha parameter of the Peaked Cosine window
		///< function used on the interpolation low-pass filter. The "usable"
		///< values are above 1.0.
		///<
	double IntFltBaseLen; ///< Interpolation low-pass filter's base length in
		///< samples (taps). Since resizing algorithm makes sure the image is
		///< at least 2X upsampled before interpolation, this base length
		///< should be at least 14 or otherwise a "grid" artifact will be
		///< introduced due to discontinuity in interpolation function.
		///< IntFltBaseLen together with other parameters should be tuned in a
		///< way that produces the flattest frequency response in 0-0.5
		///< normalized frequency range (this range is due to 2X upsampling).
		///<
	double IntFltCutoff; ///< Interpolation low-pass filter's cutoff frequency
		///< (normalized). The "usable" range is 0.5 to 1.0.
		///<
	double LPFltAlpha; ///< Alpha parameter of the Peaked Cosine window
		///< function used on the low-pass filter. The "usable" values are
		///< above 1.0.
		///<
	double LPFltBaseLen; ///< Base length of the low-pass (aka anti-aliasing
		///< or reconstruction) filter, in samples, further adjusted by the
		///< actual cutoff frequency, upsampling and downsampling factors. The
		///< "usable" range is between 6 and 12.
		///<
	double LPFltCutoffMult; ///< Low-pass filter's cutoff frequency
		///< multiplier. This value can be both below and above 1.0 as
		///< low-pass filters are inserted on downsampling and upsampling
		///< steps and always have corner frequency equal to or below 0.5pi.
		///< This multiplier shifts low-pass filter's corner frequency towards
		///< lower (if below 1.0) or higher (if above 1.0) frequencies. This
		///< multiplier can be way below 1.0 since any additional
		///< high-frequency damping will be partially corrected by the
		///< correction filter. The "usable" range is 0.7 to 1.0.
		///<

	/**
	 * Constructor assigns default parameter values.
	 */

	CImageResizerParams()
		: CorrFltAlpha( 1.0 )
		, CorrFltBaseLen( 7.0 )
		, IntFltAlpha( 4.5 )
		, IntFltBaseLen( 14.0 )
		, IntFltCutoff( 0.75 )
		, LPFltAlpha( 4.5 )
		, LPFltBaseLen( 8.0 )
		, LPFltCutoffMult( 0.80 )
	{
	}
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
 * @tparam fptype Floating point type to use for storing intermediate data and
 * variables. For variables that are not used in intensive calculations the
 * "double" type is always used. On the latest Intel processors (like
 * i7-4770K) there is almost no performance difference between "double" and
 * "float". Image quality differences between "double" and "float" are not
 * apparent on 8-bit images. At the same time the "float" uses half amount of
 * working memory the "double" type uses.
 */

template< class fptype = float >
class CImageResizer
{
public:
	/**
	 * Constructor initializes the resizer.
	 *
	 * @param aBitDepth Required bit depth of resulting image (4-16). Also
	 * affects resolution of the interpolation filter. If integer value output
	 * is used (e.g. uint8_t), the bit depth also affects rounding: for
	 * example, if aBitDepth=6 and "Tout" is uint8_t, the result will be
	 * rounded to 6 most significant bits (2 least significant bits
	 * truncated, with dithering applied).
	 * @param aParams Resizing algorithm's parameters to use. Leave out for
	 * default values. Can be useful when performing automatic optimization of
	 * parameters.
	 */

	CImageResizer( const int aBitDepth = 8,
		const CImageResizerParams& aParams = CImageResizerParams() )
		: Params( aParams )
		, BitDepth( aBitDepth )
	{
		const double SNR = -6.02 * ( BitDepth + 3 );
		FilterFracs = (int) ceil( 0.2077685412 * exp( -0.05744712479 * SNR ));

		if( FilterFracs < 2 )
		{
			FilterFracs = 2;
		}

		FilterBank.init( FilterFracs, Params.IntFltBaseLen,
			Params.IntFltCutoff, Params.IntFltAlpha );
	}

	/**
	 * Function resizes image.
	 *
	 * @param SrcBuf Source image buffer.
	 * @param SrcWidth Source image width.
	 * @param SrcHeight Source image height.
	 * @param SrcScanlineSize Physical size of source scanline in elements. If
	 * this value is below 1, SrcWidth * ElCount0 will be used as the physical
	 * source scanline size.
	 * @param NewBuf Buffer to accept the resized image.
	 * @param NewWidth New image width.
	 * @param NewHeight New image height.
	 * @param ElCount Number of elements (channels) in each pixel (1-4).
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
	 * @param aThreadPool Thread pool to use to perform resizing using
	 * multiple threads. Leave at default NULL value for no multi-threaded
	 * operation.
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
		const int NewWidth, const int NewHeight, const int ElCount,
		const double k, double ox = 0.0, double oy = 0.0,
		CImageResizerThreadPool* const aThreadPool = NULL ) const
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

		CImageResizerThreadPool DefThreadPool;
		CImageResizerThreadPool& ThreadPool = ( aThreadPool == NULL ?
			DefThreadPool : *aThreadPool );

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

		const int SrcWidthE = SrcWidth * ElCount;
		const int NewWidthE = NewWidth * ElCount;

		if( SrcScanlineSize < 1 )
		{
			SrcScanlineSize = SrcWidthE;
		}

		// Horizontal scanline filtering and resizing.

		CStructArray< CFilterStep > FltSteps;
		buildFilterSteps( FltSteps, SrcWidth, NewWidth, ElCount, kx, ox );

		const int ThreadCount = ThreadPool.getSuggestedWorkloadCount();
			// Includes the current thread.

		CThreadData< Tin > td[ ThreadCount ];
		int i;

		for( i = 0; i < ThreadCount; i++ )
		{
			if( i > 0 ) ThreadPool.addWorkload( &td[ i ]);
			td[ i ].init( i, ThreadCount, ElCount, &FltSteps,
				FilterBank.getFilterLen() );

			td[ i ].initScanlineQueue( false, SrcHeight, SrcWidth );
		}

		CBuffer< fptype > FltBuf( NewWidthE * SrcHeight ); // Temporary buffer
			// that receives horizontally-filtered and resized image.

		for( i = 0; i < SrcHeight; i++ )
		{
			td[ i % ThreadCount ].addScanlineToQueue(
				(void*) &SrcBuf[ i * SrcScanlineSize ],
				&FltBuf[ i * NewWidthE ]);
		}

		ThreadPool.startAllWorkloads();
		td[ 0 ].processScanlineQueue();
		ThreadPool.waitAllWorkloadsToFinish();

		// Evaluate pre-multipliers used on output stage.

		const bool IsInFloat = ( (Tin) 0.4 != 0 );
		const bool IsOutFloat = ( (Tout) 0.4 != 0 );
		fptype OutMul; // Output multiplier.
		int j;

		if( IsOutFloat )
		{
			OutMul = 1.0;
		}
		else
		if( sizeof( Tout ) == 1 )
		{
			OutMul = 255.0;
		}
		else
		{
			OutMul = 65535.0;
		}

		if( !IsInFloat )
		{
			if( sizeof( Tin ) == 1 )
			{
				OutMul /= 255.0;
			}
			else
			{
				OutMul /= 65535.0;
			}
		}

		// Horizontal scanline filtering and resizing.

		buildFilterSteps( FltSteps, SrcHeight, NewHeight, ElCount, ky, oy );

		if( IsOutFloat && sizeof( FltBuf[ 0 ]) == sizeof( Tout ))
		{
			// In-place output.

			for( i = 0; i < ThreadCount; i++ )
			{
				td[ i ].initScanlineQueue( true, NewWidth, SrcHeight,
					NewWidthE, NewWidthE );
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

			if( OutMul != 1.0 )
			{
				for( i = 0; i < NewWidthE * NewHeight; i++ )
				{
					NewBuf[ i ] = (Tout) ( NewBuf[ i ] * OutMul );
				}
			}

			return;
		}

		CBuffer< fptype > ResBuf(( NewWidthE + ElCount ) * NewHeight );
			// Resulting resized image, includes 1 more pixel for dithering.

		for( i = 0; i < ThreadCount; i++ )
		{
			td[ i ].initScanlineQueue( true, NewWidth, SrcHeight,
				NewWidthE, NewWidthE + ElCount );
		}

		for( i = 0; i < NewWidth; i++ )
		{
			td[ i % ThreadCount ].addScanlineToQueue(
				&FltBuf[ i * ElCount ], &ResBuf[ i * ElCount ]);
		}

		ThreadPool.startAllWorkloads();
		td[ 0 ].processScanlineQueue();
		ThreadPool.waitAllWorkloadsToFinish();
		ThreadPool.removeAllWorkloads();

		// Perform output with dithering (for integer output only).

		if( IsOutFloat )
		{
			for( i = 0; i < NewHeight; i++ )
			{
				fptype* const ResScanline =
					&ResBuf[ i * ( NewWidthE + ElCount )];

				Tout* op = &NewBuf[ i * NewWidthE ];

				if( OutMul == 1.0 )
				{
					for( j = 0; j < NewWidthE; j++ )
					{
						op[ j ] = (Tout) ResScanline[ j ];
					}
				}
				else
				{
					for( j = 0; j < NewWidthE; j++ )
					{
						op[ j ] = (Tout) ( ResScanline[ j ] * OutMul );
					}
				}
			}

			return;
		}

		int TruncBits; // The number of lower bits to truncate and dither.
		int PeakOutVal;

		if( sizeof( Tout ) == 1 )
		{
			TruncBits = 8 - BitDepth;
			PeakOutVal = 255;
		}
		else
		{
			TruncBits = 16 - BitDepth;
			PeakOutVal = 65535;
		}

		const int TrMul = ( TruncBits > 0 ? 1 << TruncBits : 1 );

		CBuffer< fptype > ResScanlineDith0( NewWidthE + ElCount * 2 );
			// Error propagation array for dithering, first and last pixels
			// unused.
		fptype* const ResScanlineDith = ResScanlineDith0 + ElCount;
			// Offsetted pointer to support -ElCount offset addressing.

		for( i = -ElCount; i < NewWidthE + ElCount; i++ )
		{
			ResScanlineDith[ i ] = 0.0;
		}

		for( i = 0; i < NewHeight; i++ )
		{
			fptype* const ResScanline = &ResBuf[ i * ( NewWidthE + ElCount )];

			for( j = 0; j < ElCount; j++ )
			{
				ResScanline[ NewWidthE + j ] = 0.0;
			}

			if( OutMul == 1.0 )
			{
				for( j = 0; j < NewWidthE; j++ )
				{
					ResScanline[ j ] += ResScanlineDith[ j ];
					ResScanlineDith[ j ] = 0.0;
				}
			}
			else
			{
				for( j = 0; j < NewWidthE; j++ )
				{
					ResScanline[ j ] = ResScanline[ j ] * OutMul +
						ResScanlineDith[ j ];

					ResScanlineDith[ j ] = 0.0;
				}
			}

			Tout* op = &NewBuf[ i * NewWidthE ];

			for( j = 0; j < NewWidthE; j++ )
			{
				const int z = (int) round( ResScanline[ j ] / TrMul ) * TrMul;
				const fptype Noise = ResScanline[ j ] - z;

				// Classic quasi-random dithering. Compressed by PNG well,
				// looks OK.

				ResScanline[ j + ElCount ] += Noise * (fptype) 0.4375;
				ResScanlineDith[ j + ElCount ] += Noise * (fptype) 0.0625;
				ResScanlineDith[ j ] += Noise * (fptype) 0.3125;
				ResScanlineDith[ j - ElCount ] += Noise * (fptype) 0.1875;

				if( z < 0 )
				{
					*op = 0;
				}
				else
				if( z > PeakOutVal )
				{
					*op = (Tout) PeakOutVal;
				}
				else
				{
					*op = (Tout) z;
				}

				op++;
			}
		}
	}

private:
	CImageResizerParams Params; ///< Algorithm's parameters currently in use.
		///<
	int BitDepth; ///< Bit resolution of the resulting image.
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
	 * @brief Resizing position structure.
	 *
	 * Structure holds resizing position and pointer to fractional delay
	 * filter.
	 */

	struct CResizePos
	{
		fptype x; ///< "X" interpolation coefficient.
			///<
		const fptype* ftp; ///< Fractional delay filter pointer.
			///<
		int SrcOffs; ///< Source scanline offset.
			///<
	};

	/**
	 * @brief Filtering step structure.
	 *
	 * Structure defines a single filtering step. 1 or more steps may be
	 * performed before the actual resizing takes place. Each step ensures
	 * that scanline data contains enough pixels to perform the next step
	 * (which may be resizing) without exceeding scanline's bounds.
	 */

	struct CFilterStep
	{
		bool IsUpsample; ///< "True" if this step is an upsampling step,
			///< "false" if downsampling step. Should be set to "false" if
			///< ResampleFactor equals 0.
			///<
		int ResampleFactor; ///< Resample factor (>=1). If 0, this is a
			///< resizing step. This value should be >1 if IsUpsample equals
			///< "true".
			///<
		CBuffer< fptype > Flt; ///< Filter kernel to use on this step.
			///<
		int FltLatency; ///< Filter's latency (shift) in pixels.
			///<
		int ElCount; ///< The number of channels in each source and
			///< destination pixel.
			///<
		int InLen; ///< Input scanline's length in pixels.
			///<
		int InBuf; ///< Input buffer index, 0 or 1.
			///<
		int InPrefix; ///< Required input prefix pixels. These prefix pixels
			///< will be filled with source scanline's first pixel value.
			///<
		int InSuffix; ///< Required input suffix pixels. These suffix pixels
			///< will be filled with source scanline's last pixel value.
			///<
		int OutLen; ///< Length of the resulting scanline.
			///<
		int OutBuf; ///< Output buffer index. 0 or 1; 2 for the last step.
			///<
		int OutPrefix; ///< Required output prefix pixels. These prefix pixels
			///< will not be pre-filled with any values. Value is valid only
			///< if IsUpsample equals "true".
			///<
		int OutSuffix; ///< Required input suffix pixels. These suffix pixels
			///< will not be pre-filled with any values. Value is valid only
			///< if IsUpsample equals "true".
			///<
		CBuffer< fptype > PrefixDC; ///< DC component fluctuations added at
			///< the start of the resulting scanline, used when IsUpsample
			///< equals "true".
			///<
		CBuffer< fptype > SuffixDC; ///< DC component fluctuations added at
			///< the end of the resulting scanline, used when IsUpsample
			///< equals "true".
			///<
		CBuffer< CResizePos > RPosBuf; ///< Resizing positions buffer.
			///<

		/**
		 * Function that replicates a set of adjacent elements several times
		 * in a row. This operation is usually used to replicate pixels at the
		 * start or end of image's scanline.
		 *
		 * @param ip Source array.
		 * @param ipl Source array length (usually 1..4, but can be any
		 * number).
		 * @param[out] op Destination buffer.
		 * @param l Number of times the source array should be replicated (the
		 * destination buffer should be able to hold ipl * l number of
		 * elements).
		 * @param opinc Destination buffer position increment after
		 * replicating the source array. This value should be equal to at
		 * least ipl.
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
			if( InPrefix + InSuffix == 0 )
			{
				return;
			}

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
			fptype* op0 = &Dst[ -OutPrefix * ElCount ];
			memset( op0, 0, ( OutPrefix + OutLen + OutSuffix ) * ElCount *
				sizeof( fptype ));

			const fptype* const f = Flt;
			const int flen = Flt.getCapacity();
			const fptype* ip = Src;
			fptype* op;
			const int opstep = ElCount * ResampleFactor;
			int l = InLen;
			int i;

			if( ElCount == 1 )
			{
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
			}
			else
			if( ElCount == 4 )
			{
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
			}
			else
			if( ElCount == 3 )
			{
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
			}
			else
			if( ElCount == 2 )
			{
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
			}

			ip = &Src[( InLen - 1 ) * ElCount ];
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
			op = Dst;
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
	};

	/**
	 * Function assigns filter parameters to the specified filter step object.
	 *
	 * @param fs Filter step to assign parameter to. This step cannot be the
	 * last step if ResampleFactor greater than 1 was specified.
	 * @param IsUpsample "True" if upsampling step. Should be set to "false"
	 * if FltCutoff is negative.
	 * @param ResampleFactor Resampling factor of this filter (>=1).
	 * @param FltCutoff Filter cutoff point (<= 1.0). If a negative value is
	 * specified, filter's kernel will not be calculated, and a pre-calculated
	 * kernel is assumed to be available. This value will be divided by the
	 * ResampleFactor if IsUpsample equals "true".
	 * @param k Resizing coefficient at the previous processing step, may be
	 * adjusted on return.
	 * @param o Starting pixel offset inside the source image, may be adjusted
	 * on return.
	 * @param prevfs Previous parameters to derive buffer index and length
	 * from. If NULL was specified, then SrcLen and ElCount should be
	 * provided explicitly.
	 * @param SrcLen Source scanline length in pixels.
	 * @param ElCount The number of elements (channels) in each pixel.
	 */

	void assignFilterParams( CFilterStep& fs, const bool IsUpsample,
		const int ResampleFactor, const double FltCutoff, double& k,
		double& o, const CFilterStep* const prevfs, const int SrcLen = 0,
		const int ElCount = 0 ) const
	{
		if( FltCutoff >= 0.0 )
		{
			double Len2 = 0.25 * Params.LPFltBaseLen / FltCutoff;
			double Freq = M_PI * Params.LPFltCutoffMult * FltCutoff;

			if( IsUpsample )
			{
				Len2 *= ResampleFactor;
				Freq /= ResampleFactor;
			}

			CDSPPeakedCosineLPF w( Len2, Freq, Params.LPFltAlpha );

			fs.Flt.alloc( w.KernelLen );
			fs.FltLatency = w.fl2;
			w.generateLPF( &fs.Flt[ 0 ], 1.0 );
			optimizeFIRFilter( fs.Flt, fs.FltLatency );
			normalizeFIRFilter( &fs.Flt[ 0 ], fs.Flt.getCapacity(),
				( IsUpsample ? ResampleFactor : 1.0 ));
		}

		fs.IsUpsample = IsUpsample;
		fs.ResampleFactor = ResampleFactor;

		if( prevfs == NULL )
		{
			fs.ElCount = ElCount;
			fs.InLen = SrcLen;
			fs.InBuf = 0;
		}
		else
		{
			fs.ElCount = prevfs -> ElCount;
			fs.InLen = prevfs -> OutLen;
			fs.InBuf = prevfs -> OutBuf;
		}

		if( IsUpsample )
		{
			fs.InPrefix = 0;
			fs.InSuffix = 0;
			fs.OutLen = SrcLen * ResampleFactor;
			fs.OutPrefix = fs.FltLatency;
			fs.OutSuffix = fs.Flt.getCapacity() - fs.FltLatency -
				ResampleFactor;

			k *= ResampleFactor;
			o *= ResampleFactor;

			fs.PrefixDC.alloc( fs.OutSuffix );
			const fptype* ip = &fs.Flt[ fs.FltLatency + ResampleFactor ];
			int l = fs.OutSuffix;
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

			fs.SuffixDC.alloc( fs.FltLatency );
			fptype* op = &fs.SuffixDC[ 0 ];
			l = fs.FltLatency;
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
	}

	/**
	 * Function adds a correction filter that tries to achieve a linear
	 * frequency response at all frequencies. The actual resulting response
	 * may feature a slight damping of the highest frequencies since a
	 * suitably short correction filter cannot fix steep high-frequency
	 * damping.
	 *
	 * @param Steps Filtering steps.
	 * @param bw Resulting bandwidth relative to the original bandwidth (which
	 * is 1.0), usually 1/k. If image's size was increased this value should
	 * be >1.0, otherwise <1.0.
	 * @param k Resizing coefficient at the last processing step, may be
	 * adjusted on return.
	 * @param o Starting pixel offset inside the source image, may be
	 * adjusted on return.
	 */

	void addCorrectionFilter( CStructArray< CFilterStep >& Steps,
		const double bw, double& k, double& o ) const
	{
		const int BinCount = 65; // Frequency response bins to control.
		const int BinCount1 = BinCount - 1;
		const double bwi = ( bw >= 1.0 ? 1.0 : bw ); // Spectrum bandwidth of
			// interest which is covered by spectral bins.
		const double bwiPerBin = bwi / BinCount1; // Bandwidth controlled by a
			// single spectral bin.
		const int AddedBin = ( bwi + bwiPerBin <= bw ? 1 : 0 ); // If the
			// total bandwidth is much larger than bandwidth of interest, add
			// an additional bin which will not be corrected: it works to
			// bypass correction of frequencies beyond interest.
		double curbw = 1.0; // Bandwidth of the filter at the current step.
		int i;
		int j;

		CBuffer< double > Bins( BinCount + AddedBin ); // Adjustment
			// introduced by all steps at all frequencies of interest.

		for( j = 0; j < BinCount + AddedBin; j++ )
		{
			Bins[ j ] = 1.0;
		}

		for( i = 0; i < Steps.getItemCount(); i++ )
		{
			const CFilterStep& fs = Steps[ i ];

			if( fs.IsUpsample )
			{
				curbw *= fs.ResampleFactor;
			}

			CBuffer< fptype > TmpBuf;
			const fptype* UseFlt;
			int UseFltLen;

			if( fs.ResampleFactor == 0 )
			{
				// Resizing step is accompanied by a filter stored in a filter
				// bank - it must be extracted first. This is also the last
				// resizing stage, so "curbw" will not be adjusted.

				UseFltLen = FilterBank.getFilterLen();
				TmpBuf.alloc( UseFltLen );
				copyArray( &FilterBank[ 0 ], &TmpBuf[ 0 ], UseFltLen, 2, 1 );
				UseFlt = &TmpBuf[ 0 ];
			}
			else
			{
				UseFlt = &fs.Flt[ 0 ];
				UseFltLen = fs.Flt.getCapacity();
			}

			// Calculate frequency response adjustment introduced by the
			// filter on this step, within the bounds of bandwidth of
			// interest.

			for( j = 0; j < BinCount; j++ )
			{
				const double th = M_PI * bwi / curbw * j / BinCount1;
				double re;
				double im;

				calcFIRFilterResponse( UseFlt, UseFltLen, th, re, im );

				Bins[ j ] /= sqrt( re * re + im * im );
			}

			if( !fs.IsUpsample && fs.ResampleFactor > 1 )
			{
				curbw /= fs.ResampleFactor;
			}
		}

		// Calculate filter kernel.

		CDSPFIREQ EQ;
		EQ.init( bw * 2.0, Params.CorrFltBaseLen * ( bw >= 1.0 ? bw : 1.0 ),
			BinCount + AddedBin, 0.0, bwi + bwiPerBin * AddedBin, false,
			Params.CorrFltAlpha );

		const CFilterStep* const prevfs = &Steps[ Steps.getItemCount() - 1 ];
		CFilterStep& fs = Steps.add();
		fs.FltLatency = EQ.getKernelLatency();

		CBuffer< double > Kernel( EQ.getKernelLength() );
		EQ.buildKernel( Bins, &Kernel[ 0 ]);
		normalizeFIRFilter( &Kernel[ 0 ], Kernel.getCapacity(), 1.0 );
		optimizeFIRFilter( Kernel, fs.FltLatency );

		fs.Flt.alloc( Kernel.getCapacity() );
		copyArray( &Kernel[ 0 ], &fs.Flt[ 0 ], Kernel.getCapacity() );
		normalizeFIRFilter( &fs.Flt[ 0 ], fs.Flt.getCapacity(), 1.0 );

/*		// Print a theoretically achieved final frequency response at various
		// feature sizes (from DC to 1 pixel). Values above 255 means features
		// become brighter, values below 255 means features become dimmer.

		for( j = 0; j < BinCount; j++ )
		{
			const double th = M_PI * j / ( BinCount - 1 );
			double re;
			double im;

			calcFIRFilterResponse( &fs.Flt[ 0 ], fs.Flt.getCapacity(), th,
				re, im );

			printf( "%f\n", sqrt( re * re + im * im ) / Bins[ j ] * 255 );
		}

		printf( "***\n" );
*/
		assignFilterParams( fs, false, 1, -1.0, k, o, prevfs );
	}

	/**
	 * Function builds sequence of filtering steps depending on the specified
	 * resizing coefficient. The last steps included are always the resizing
	 * step then the correction step.
	 *
	 * @param Steps Array that receives filtering steps.
	 * @param SrcLen Source scanline's length in pixels.
	 * @param NewLen New scanline's length in pixels.
	 * @param ElCount Number of elements in each pixel (1..4), corresponds to
	 * the number of channels in both source and destination images.
	 * @param k Resizing coefficient, will be adjusted on return.
	 * @param o Starting pixel offset inside the source image, will be
	 * adjusted on return.
	 */

	void buildFilterSteps( CStructArray< CFilterStep >& Steps,
		const int SrcLen, const int NewLen, const int ElCount, double& k,
		double& o ) const
	{
		Steps.clear();

		const double bw = 1.0 / k; // Resulting bandwidth.
		double FltCutoff = ( k <= 1.0 ? 1.0 : 1.0 / k ); // Cutoff frequency
			// of the first filter step.
		const int UpsampleFactor = ( (int) floor( k ) >= 2 ? 1 : 2 );
		CFilterStep* prevfs;

		// Add 1 upsampling or several downsampling filters.

		if( UpsampleFactor > 1 )
		{
			assignFilterParams( Steps.add(), true, UpsampleFactor, FltCutoff,
				k, o, NULL, SrcLen, ElCount );

			prevfs = &Steps[ 0 ];
		}
		else
		{
			int DownsampleSteps = 0; // Number of additional 2x downsampling
				// steps to perform.
			double tk = k;

			while( (int) floor( tk ) >= 8 )
			{
				tk *= 0.5;
				DownsampleSteps++;
				FltCutoff *= 2.0;
			}

			assignFilterParams( Steps.add(), false,
				( (int) floor( k ) >= 4 ? 2 : 1 ), FltCutoff, k, o, NULL,
				SrcLen, ElCount );

			// Insert downsampling steps (if any).

			prevfs = &Steps[ 0 ];

			while( DownsampleSteps > 0 )
			{
				DownsampleSteps--;
				CFilterStep& fs = Steps.add();
				assignFilterParams( fs, false, 2, 0.25, k, o, prevfs );
				prevfs = &fs;
			}
		}

		// Insert resizing and correction steps.

		CFilterStep& fs = Steps.add();
		fs.IsUpsample = false;
		fs.ResampleFactor = 0;
		fs.ElCount = ElCount;
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
		CResizePos* rpos = &fs.RPosBuf[ 0 ];
		int i;

		for( i = 0; i < fs.OutLen; i++ )
		{
			const double SrcPos = o + k * i;
			const int SrcPosInt = (int) floor( SrcPos );
			double x = ( SrcPos - SrcPosInt ) * FilterFracs;
			const int fti = (int) x;
			rpos -> x = (fptype) ( x - fti );
			rpos -> ftp = &FilterBank[ fti ];
			rpos -> SrcOffs = ( SrcPosInt - FilterLenD21 ) * ElCount;
			rpos++;
		}

		addCorrectionFilter( Steps, bw, k, o );

		CFilterStep& lastfs = Steps[ Steps.getItemCount() - 1 ];
		lastfs.OutBuf = 2;

		if( lastfs.OutLen > NewLen )
		{
			lastfs.OutLen = NewLen;
		}
	}

	/**
	 * Function calculates an optimal buffer size that will cover all needs of
	 * the specified filtering steps.
	 *
	 * @param Steps Filtering steps.
	 * @param BufSize This variable receives the required buffer size.
	 * @param BufOffset This variable receives buffer offset that should be
	 * used when passing the buffer pointer to the filtering functions.
	 */

	static void calcBufSize( const CStructArray< CFilterStep >& Steps,
		int& BufSize, int& BufOffset )
	{
		int MaxPrefix = 0;
		int MaxLen = 0;
		int i;

		for( i = 0; i < Steps.getItemCount(); i++ )
		{
			const CFilterStep& fs = Steps[ i ];

			int l = fs.InPrefix * fs.ElCount;

			if( l > MaxPrefix )
			{
				MaxPrefix = l;
			}

			l = ( fs.InLen + fs.InSuffix ) * fs.ElCount;

			if( l > MaxLen )
			{
				MaxLen = l;
			}

			if( fs.IsUpsample )
			{
				l = fs.OutPrefix * fs.ElCount;

				if( l > MaxPrefix )
				{
					MaxPrefix = l;
				}

				l = ( fs.OutLen + fs.OutSuffix ) * fs.ElCount;

				if( l > MaxLen )
				{
					MaxLen = l;
				}
			}
			else
			{
				l = fs.OutLen * fs.ElCount;

				if( l > MaxLen )
				{
					MaxLen = l;
				}
			}
		}

		BufSize = MaxPrefix + MaxLen;
		BufOffset = MaxPrefix;
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
	 * @param DstLineLen Required length of destination scanline, in pixels.
	 * @param DstLineInc Destination scanline position increment, should be
	 * divisible by ElCount.
	 * @param ElCount Number of elements in each pixel (1..4), corresponds to
	 * the number of channels in the image.
	 * @param rpos Resizing positions buffer (length = DstLineLen).
	 * @param aFilterLen Interpolation filter length.
	 */

	static void resizeLine( const fptype* SrcLine, fptype* DstLine,
		int DstLineLen, const int DstLineInc, const int ElCount,
		const CResizePos* rpos, const int aFilterLen )
	{
#define AVIR_RESIZE_PART1 \
			while( DstLineLen > 0 ) \
			{ \
				const fptype x = rpos -> x; \
				const fptype* ftp = rpos -> ftp; \
				const fptype* Src = SrcLine + rpos -> SrcOffs; \
				int l = aFilterLen;

#define AVIR_RESIZE_PART2 \
				DstLineLen--; \
				DstLine += DstLineInc; \
				rpos++; \
			}

		if( ElCount == 1 )
		{
			AVIR_RESIZE_PART1

			fptype sum = 0.0;

			while( l > 0 )
			{
				sum += ( ftp[ 0 ] + ftp[ 1 ] * x ) * Src[ 0 ];
				ftp += 2;
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
				const fptype xx = ftp[ 0 ] + ftp[ 1 ] * x;
				ftp += 2;
				sum[ 0 ] += xx * Src[ 0 ];
				sum[ 1 ] += xx * Src[ 1 ];
				sum[ 2 ] += xx * Src[ 2 ];
				sum[ 3 ] += xx * Src[ 3 ];
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
				const fptype xx = ftp[ 0 ] + ftp[ 1 ] * x;
				ftp += 2;
				sum[ 0 ] += xx * Src[ 0 ];
				sum[ 1 ] += xx * Src[ 1 ];
				sum[ 2 ] += xx * Src[ 2 ];
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
				const fptype xx = ftp[ 0 ] + ftp[ 1 ] * x;
				ftp += 2;
				sum[ 0 ] += xx * Src[ 0 ];
				sum[ 1 ] += xx * Src[ 1 ];
				Src += 2;
				l--;
			}

			DstLine[ 0 ] = sum[ 0 ];
			DstLine[ 1 ] = sum[ 1 ];

			AVIR_RESIZE_PART2
		}
	}
#undef AVIR_RESIZE_PART1
#undef AVIR_RESIZE_PART2

	/**
	 * @brief Thread-isolated data used for scanline processing.
	 *
	 * This structure holds data necessary for image's horizontal or vertical
	 * scanline processing, including scanline processing queue.
	 *
	 * @tparam Tin Source element data type. Intermediate buffers store data
	 * in floating point format.
	 */

	template< class Tin >
	class CThreadData : public CImageResizerThreadPool :: CWorkload
	{
	public:
		virtual void process()
		{
			processScanlineQueue();
		}

		/**
		 * Function initializes *this thread data object and assigns certain
		 * variables provided by the higher level code.
		 *
		 * @param aThreadIndex Index of this thread data (0-based).
		 * @param aThreadCount Total number of threads used during processing.
		 * @param aElCount Pixel element count.
		 * @param aSteps Filtering steps.
		 * @param aFilterLen Interpolation filter length.
		 */

		void init( const int aThreadIndex, const int aThreadCount,
			const int aElCount, const CStructArray< CFilterStep >* aSteps,
			const int aFilterLen )
		{
			ThreadIndex = aThreadIndex;
			ThreadCount = aThreadCount;
			ElCount = aElCount;
			Steps = aSteps;
			FilterLen = aFilterLen;
		}

		/**
		 * Function initializes scanline processing queue, and updates
		 * capacities of intermediate buffers.
		 *
		 * @param aIsVertical "True" if this queue uses vertical scanline
		 * processing.
		 * @param TotalLines The total number of scanlines that will be
		 * processed by all threads.
		 * @param aSrcLen Source scanline length in pixels.
		 * @param aSrcIncr Source scanline buffer increment. Ignored in
		 * horizontal scanline processing.
		 * @param aResIncr Resulting scanline buffer increment. Ignored in
		 * horizontal scanline processing.
		 */

		void initScanlineQueue( const bool aIsVertical, const int TotalLines,
			const int aSrcLen, const int aSrcIncr = 0,
			const int aResIncr = 0 )
		{
			int BufSize;
			int BufOffset;
			calcBufSize( *Steps, BufSize, BufOffset );

			Bufs.alloc( BufSize * 2 );
			BufPtrs[ 0 ] = Bufs + BufOffset;
			BufPtrs[ 1 ] = BufPtrs[ 0 ] + BufSize;

			IsVertical = aIsVertical;
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

		void addScanlineToQueue( void* const SrcBuf, fptype* const ResBuf )
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

			if( IsVertical )
			{
				for( i = 0; i < QueueLen; i++ )
				{
					processScanlineV( (fptype*) Queue[ i ].SrcBuf,
						Queue[ i ].ResBuf );
				}
			}
			else
			{
				for( i = 0; i < QueueLen; i++ )
				{
					processScanlineH( (Tin*) Queue[ i ].SrcBuf,
						Queue[ i ].ResBuf );
				}
			}
		}

	private:
		int ThreadIndex; ///< Thread index.
			///<
		int ThreadCount; ///< Thread count.
			///<
		int ElCount; ///< Pixel element count.
			///<
		const CStructArray< CFilterStep >* Steps; ///< Filtering steps.
			///<
		int FilterLen; ///< Interpolation filter len.
			///<
		CBuffer< fptype > Bufs; ///< Flip-flop intermediate buffers.
			///<
		fptype* BufPtrs[ 3 ]; ///< Flip-flop buffer pointers (referenced by
			///< filtering step's InBuf and OutBuf indices).
			///<
		bool IsVertical; ///< "True" if the last queue uses vertical scanline
			///< processing.
			///<
		int SrcLen; ///< Source scanline length in the last queue.
			///<
		int SrcIncr; ///< Source scanline buffer increment in the last queue.
			///<
		int ResIncr; ///< Resulting scanline buffer increment in the last
			///< queue.
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
			fptype* ResBuf; ///< Resulting scanline buffer.
				///<
		};

		CBuffer< CQueueItem > Queue; ///< Scanline processing queue.
			///<
		int QueueLen; ///< Queue length.
			///<

		/**
		 * Function processes a single vertical scanline.
		 *
		 * @param SrcBuf Source scanline buffer. Can be either horizontal or
		 * vertical.
		 * @param ResBuf Resulting scanline buffer.
		 */

		void processScanlineV( const fptype* const SrcBuf,
			fptype* const ResBuf )
		{
			const fptype* ip = SrcBuf;
			fptype* op = BufPtrs[ 0 ];
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

			BufPtrs[ 2 ] = ResBuf;

			for( j = 0; j < Steps -> getItemCount(); j++ )
			{
				const CFilterStep& fs = (*Steps)[ j ];
				fs.prepareInBuf( BufPtrs[ fs.InBuf ]);
				const int DstIncr = ( fs.OutBuf == 2 ? ResIncr : ElCount );

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
					resizeLine( BufPtrs[ fs.InBuf ], BufPtrs[ fs.OutBuf ],
						fs.OutLen, DstIncr, ElCount, fs.RPosBuf, FilterLen );
				}
			}
		}

		/**
		 * Function processes a single horizontal scanline.
		 *
		 * @param SrcBuf Source scanline buffer. Can be either horizontal or
		 * vertical.
		 * @param ResBuf Resulting scanline buffer.
		 */

		void processScanlineH( const Tin* const SrcBuf,
			fptype* const ResBuf )
		{
			bool IsDirectInput;

			if( (*Steps)[ 0 ].IsUpsample && sizeof( Tin ) == sizeof( fptype ))
			{
				// This is not so reliable way to check for equality between
				// Tin and fptype (typeinfo would be perfect), but should work
				// in most practical cases.

				IsDirectInput = true;
			}
			else
			{
				IsDirectInput = false;
				copyArray( SrcBuf, BufPtrs[ 0 ], SrcLen * ElCount );
			}

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
						fs.doUpsample(( j == 0 && IsDirectInput ?
							(fptype*) SrcBuf : BufPtrs[ fs.InBuf ]),
							BufPtrs[ fs.OutBuf ]);
					}
					else
					{
						fs.doFilter( BufPtrs[ fs.InBuf ],
							BufPtrs[ fs.OutBuf ], ElCount );
					}
				}
				else
				{
					resizeLine( BufPtrs[ fs.InBuf ], BufPtrs[ fs.OutBuf ],
						fs.OutLen, ElCount, ElCount, fs.RPosBuf, FilterLen );
				}
			}
		}
	};
};

} // namespace avir

#endif // AVIR_CIMAGERESIZER_INCLUDED
