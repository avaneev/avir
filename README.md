# AVIR #
## Introduction ##
Me, Aleksey Vaneev, is happy to offer you an open source image resizing
library which has reached a production level of quality, and is ready to be
incorporated into any project. This library features routines for both down-
and upsizing of 8- and 16-bit, 1 to 4-channel images. Image resizing routines
were implemented in multi-platform C++ code, and have a high level of
optimality. Beside resizing, this library offers a sub-pixel shift operation.

The resizing algorithm at first produces 2X upsized image (relative to the
source image size, or relative to the destination image size if downsizing is
performed) and then performs interpolation using a bank of sinc function-based
fractional delay filters. On the last stage a correction filter is applied
which fixes smoothing introduced on previous steps.

An important element utilized by this algorithm is the so called Peaked Cosine
window function, which is applied over sinc function in all filters. Please
consult the documentation for more details.

*AVIR is devoted to women. Your digital photos can look good at any size!*

## Requirements ##
C++ compiler and system with efficient "float" floating point (24-bit
mantissa) type support. This library can also internally use the "double" and
SIMD floating point types during resizing if needed. This library does not
have dependencies beside the standard C library.

## Links ##
* [Documentation](http://avaneev.atspace.cc/avir/Documentation/)

## Usage Information ##
The image resizer is represented by the **avir::CImageResizer<>** class, which
is a single front-end class for the whole library. You do not basically need
to use nor understand any other classes beside this class.

The code of the library resides in the "avir" C++ namespace, effectively
isolating it from all other code. The code is thread-safe. You need just
a single resizer object per running application, at any time, even when
resizing images concurrently.

To resize images in your application, simply add 3 lines of code:
* # include "avir.h"
* avir :: CImageResizer<> ImageResizer( 8 );
* ImageResizer.resizeImage( ... );
(multi-threaded operation requires additional coding, see the documentation)

The library is able to process images of any bit depth: this includes 8-bit,
16-bit, float and double types. Larger integer and signed integer types are
not supported. Supported source and destination image sizes are up to 2.1
gigapixels (46300x46300 or equivalent dimensions, e.g. 53467x40100).

The code of this library was commented in the [Doxygen](http://www.doxygen.org/)
style. To generate the documentation locally you may run the
"doxygen ./other/avirdoxy.txt" command from the library's directory. Note that
the code was suitably documented allowing you to make modifications, and to
gain full understanding of the algorithm.

Preliminary tests show that this library can resize 8-bit RGB 5184x3456
(17.9 Mpixel) image down to 1037x691 (0.7 Mpixel) image in 355 milliseconds,
utilizing a single thread, on a typical Intel Core i7-4770K processor-based
system without overclocking. This scales down to 110 milliseconds if 4 threads
are utilized. This time can be reduced further down to 80 milliseconds by
utilizing SIMD floating point processing. This library's performance has a big
potential to grow together with evolving processor architectures as currently
performance is clearly limited by memory bandwidth, not by algorithm's
mathematical operations and overhead.

Multi-threaded operation is not provided by this library "out of the box".
The multi-threaded infrastructure is fully available, but requires additional
system-specific interfacing code for engagement.

## SIMD Usage Information ##
This library is capable of using SIMD floating point types for internal
variables. This means that up to 4 color channels can be processed in
parallel. For example, this gives 40% performance boost when resizing
3-channel images. Since the processing algorithm itself remains non-SIMD, the
use of SIMD types is not practical for 1-channel image resizing (due to
overhead). SIMD type can be used this way:
* # include "avir_float4_sse.h"
* avir :: CImageResizer< avir :: fpclass_float4 > ImageResizer( 8 );

## Notes ##
This library was tested for compatibility with [GNU C++](http://gcc.gnu.org/),
[Microsoft Visual C++](http://www.microsoft.com/visualstudio/eng/products/visual-studio-express-products)
and [Intel C++](http://software.intel.com/en-us/c-compilers) compilers, on 32-
and 64-bit Windows, Mac OS X and CentOS Linux. The code was also tested with
Dr.Memory/Win32 for the absence of uninitialized or unaddressable memory
accesses.

All code is fully "inline", without the need to compile any source files. The
memory footprint of the library itself is very modest, except that the size of
the temporary image buffers depends on the input and output image sizes, and
is proportionally large.

The "heart" of resizing algorithm's quality resides in the parameters defined
via the **avir::CImageResizerParams** structure. While the default set of
parameters that offers a good quality was already provided, there is still a
place for improvement exists, and the default parameters may change in a
future update. If you need to recall an exact set of parameters, simply save
them locally for a later use.

When the algorithm is run with no resizing applied (k=1), the result of
resizing will not be an exact, but a very close copy of the source image. The
reason for such inexactness is that the image is always low-pass filtered at
first to reduce aliasing during subsequent resizing, and at last filtered by a
correction filter. Such approach allows algorithm to maintain a stable level
of quality regardless of the resizing "k" factor used.

This library includes a binary command line tool "imageresize" for major
desktop platforms. This tool was designed to be used as a demonstration of
library's performance and as a reference, it uses 4 threads and float4 SIMD
internal type during processing. This tool uses the following libraries:
* turbojpeg Copyright (c) 2009-2013 D. R. Commander
* libpng Copyright (c) 1998-2013 Glenn Randers-Pehrson
* zlib Copyright (c) 1995-2013 Jean-loup Gailly and Mark Adler

## Interpolation Discussion ##
The use of certain low-pass filters and 2X upsampling in this library is
hardly debatable, because they are needed to attain a certain anti-aliasing
effect and keep ringing artifacts low. But the use of sinc function-based
interpolation filter that is 18 taps-long (may be higher in practice) can be
questioned, because even in 0th order case such interpolation filter requires
18 multiply-add operations. Comparatively, an optimal Hermit or cubic
interpolation spline requires 8 multiply and 11 add operations.

One of the reasons 18-tap filter is preferred, is because due to memory
bandwidth limitations using a lower-order filter does not provide any
significant performance increase (e.g. 14-tap filter is less than 5% more
efficient). At the same time, in comparison to cubic spline, 18-tap filter
embeds a low-pass filter that rejects signal above 0.5*pi (provides additional
anti-aliasing filtering), and this filter has a consistent shape at all
fractional offsets. Splines have a varying low-pass filter shape at
different fractional offsets (e.g. no low-pass filtering at 0.0 offset,
and maximal low-pass filtering at 0.5 offset). 18-tap filter also offers a
superior stop-band attenuation which almost guarantees absence of artifacts if
the image is considerably sharpened afterwards.

## Why 2X upsizing in AVIR? ##
Classic approaches to image resizing do not perform an additional 2X upsizing.
So, why such upsizing is needed at all in AVIR? Indeed, image resizing can be
implemented using a single interpolation filter which is applied to the source
image directly. However, such approach has limitations:

First of all, interpolation filter (during upsizing) has to be tuned to a
frequency close to pi (Nyquist) in order to reduce high-frequency smoothing:
this reduces the space left for filter optimization. Beside that, filters
with a corner frequency tuned to Nyquist frequency may become distorted in
comparison to lower-frequency tuning (during downsizing). It's usually a good
idea to have filter's stop-band begin below Nyquist so that the transition
band's shape remains stable at any setting. At the same time this complicates
a further correction filtering, because correction filter may become too
steep at the point where the stop-band beings.

Secondly, filter has to be very short (5-7 taps) or otherwise the ringing
artifacts will be very strong: it is a general rule that the steeper the
filter is around signal frequencies being removed the higher the ringing
artifacts. That is why it is preferred to move steep transitions into the
spectral area with a quieter signal. A short filter also means it cannot
provide a strong stop-band attenuation, so an interpolated image will look a
bit edgy or not very clean due to stop-band artifacts.

To sum up, only additional controlled 2X upsizing provides enough spectral
space to design interpolation filter without visible ringing artifacts yet
providing a strong stop-band attenuation and stable spectral characteristics.

## Why Peaked Cosine in AVIR? ##
First of all, AVIR is a general solution to image resizing problem. That is
why it should not be directly compared to "spline interpolation" or "Lanczos
resampling", because the latter two are only means to design interpolation
filters, and they can be implemented in a variety of ways, even in sub-optimal
ways. Secondly, with only a minimal effort AVIR can be changed to use any
existing interpolation formula and any window function, but this is just not
needed.

An effort was made to compare Peaked Cosine to Lanczos window function, and
here is the author's opinion. Peaked Cosine has two degrees of freedom whereas
Lanczos has one degree of freedom. While both functions can be used with
acceptable results, Peaked Cosine window function used in automatic parameter
optimization really pushes the limits of frequency response linearity,
anti-aliasing strength (stop-band attenuation) and low-ringing performance
which Lanczos cannot usually achieve. This is true at least when using a
general-purpose downhill simplex optimization method. Lanczos window has good
(but not better) characteristics in several special cases which makes it of
limited use in a general solution such as AVIR.

Among other window functions (Kaiser, Gaussian, Cauchy, Poisson, generalized
cosine windows) there are no better candidates as well. It looks like Peaked
Cosine function's scalability (it retains stable, almost continously-variable
spectral characteristics at any window parameter values), and its ability to
create "desirable" pass-band ripples in the frequency response near the cutoff
point contribute to its better overall quality. Somehow Peaked Cosine window
function optimization manages to converge to reasonable states in most cases
(that is why AVIR library comes with a set of equally robust, but distinctive
parameter sets) whereas all other window functions tend to produce
unpredictable optimization results.

The only disadvantage of Peaked Cosine window function is that usable filters
windowed by this function tend to be longer than "usual" (with Kaiser window
being the "golden standard" for filter length per decibel of stop-band
attenuation). This is a price that should be paid for stable spectral
characteristics.

## Users ##
This library is used by:

Please drop me a note at aleksey.vaneev@gmail.com and I will include a link to
your software product to the list of users. This list is important at
maintaining confidence in this library among the interested parties.
