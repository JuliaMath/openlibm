// Copyright (C) 1988-1994 Sun Microsystems, Inc. 2550 Garcia Avenue
// Mountain View, California  94043 All rights reserved.
//
// Any person is hereby authorized to download, copy, use, create bug fixes,
// and distribute, subject to the following conditions:
//
// 	1.  the software may not be redistributed for a fee except as
// 	    reasonable to cover media costs;
// 	2.  any copy of the software must include this notice, as well as
// 	    any other embedded copyright notices; and
// 	3.  any distribution of this software or derivative works thereof
// 	    must comply with all applicable U.S. export control laws.
//
// THE SOFTWARE IS MADE AVAILABLE "AS IS" AND WITHOUT EXPRESS OR IMPLIED
// WARRANTY OF ANY KIND, INCLUDING BUT NOT LIMITED TO THE IMPLIED
// WARRANTIES OF DESIGN, MERCHANTIBILITY, FITNESS FOR A PARTICULAR
// PURPOSE, NON-INFRINGEMENT, PERFORMANCE OR CONFORMANCE TO
// SPECIFICATIONS.
//
// BY DOWNLOADING AND/OR USING THIS SOFTWARE, THE USER WAIVES ALL CLAIMS
// AGAINST SUN MICROSYSTEMS, INC. AND ITS AFFILIATED COMPANIES IN ANY
// JURISDICTION, INCLUDING BUT NOT LIMITED TO CLAIMS FOR DAMAGES OR
// EQUITABLE RELIEF BASED ON LOSS OF DATA, AND SPECIFICALLY WAIVES EVEN
// UNKNOWN OR UNANTICIPATED CLAIMS OR LOSSES, PRESENT AND FUTURE.
//
// IN NO EVENT WILL SUN MICROSYSTEMS, INC. OR ANY OF ITS AFFILIATED
// COMPANIES BE LIABLE FOR ANY LOST REVENUE OR PROFITS OR OTHER SPECIAL,
// INDIRECT AND CONSEQUENTIAL DAMAGES, EVEN IF IT HAS BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGES.
//
// This file is provided with no support and without any obligation on the
// part of Sun Microsystems, Inc. ("Sun") or any of its affiliated
// companies to assist in its use, correction, modification or
// enhancement.  Nevertheless, and without creating any obligation on its
// part, Sun welcomes your comments concerning the software and requests
// that they be sent to fdlibm-comments@sunpro.sun.com.
// floord(integer) is itself
T(RN,                  0x0p+0,                  0x0p+0,          0x0p+0, 0)
T(RN,                 -0x0p+0,                 -0x0p+0,          0x0p+0, 0)
T(RN,                  0x1p+0,                  0x1p+0,          0x0p+0, 0)
T(RN,                 -0x1p+0,                 -0x1p+0,          0x0p+0, 0)
T(RN,   0x1.fffffffffffffp+52,   0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RN,  -0x1.fffffffffffffp+52,  -0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RN, 0x1.fffffffffffffp+1023, 0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RN,-0x1.fffffffffffffp+1023,-0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RZ,                  0x0p+0,                  0x0p+0,          0x0p+0, 0)
T(RZ,                 -0x0p+0,                 -0x0p+0,          0x0p+0, 0)
T(RZ,                  0x1p+0,                  0x1p+0,          0x0p+0, 0)
T(RZ,                 -0x1p+0,                 -0x1p+0,          0x0p+0, 0)
T(RZ,   0x1.fffffffffffffp+52,   0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RZ,  -0x1.fffffffffffffp+52,  -0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RZ, 0x1.fffffffffffffp+1023, 0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RZ,-0x1.fffffffffffffp+1023,-0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RU,                  0x0p+0,                  0x0p+0,          0x0p+0, 0)
T(RU,                 -0x0p+0,                 -0x0p+0,          0x0p+0, 0)
T(RU,                  0x1p+0,                  0x1p+0,          0x0p+0, 0)
T(RU,                 -0x1p+0,                 -0x1p+0,          0x0p+0, 0)
T(RU,   0x1.fffffffffffffp+52,   0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RU,  -0x1.fffffffffffffp+52,  -0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RU, 0x1.fffffffffffffp+1023, 0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RU,-0x1.fffffffffffffp+1023,-0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RD,                  0x0p+0,                  0x0p+0,          0x0p+0, 0)
T(RD,                 -0x0p+0,                 -0x0p+0,          0x0p+0, 0)
T(RD,                  0x1p+0,                  0x1p+0,          0x0p+0, 0)
T(RD,                 -0x1p+0,                 -0x1p+0,          0x0p+0, 0)
T(RD,   0x1.fffffffffffffp+52,   0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RD,  -0x1.fffffffffffffp+52,  -0x1.fffffffffffffp+52,          0x0p+0, 0)
T(RD, 0x1.fffffffffffffp+1023, 0x1.fffffffffffffp+1023,          0x0p+0, 0)
T(RD,-0x1.fffffffffffffp+1023,-0x1.fffffffffffffp+1023,          0x0p+0, 0)
// integer - ulp
T(RN,   0x1.eeeeeeeeeeeefp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RN,  -0x1.eeeeeeeeeeeefp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RN,    0x1.fffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RN,   -0x1.fffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RN,    0x1.fffffffffffffp-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,   -0x1.fffffffffffffp-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RZ,   0x1.eeeeeeeeeeeefp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RZ,  -0x1.eeeeeeeeeeeefp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RZ,    0x1.fffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RZ,   -0x1.fffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RZ,    0x1.fffffffffffffp-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,   -0x1.fffffffffffffp-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RU,   0x1.eeeeeeeeeeeefp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RU,  -0x1.eeeeeeeeeeeefp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RU,    0x1.fffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RU,   -0x1.fffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RU,    0x1.fffffffffffffp-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,   -0x1.fffffffffffffp-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,   0x1.eeeeeeeeeeeefp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RD,  -0x1.eeeeeeeeeeeefp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RD,    0x1.fffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RD,   -0x1.fffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RD,    0x1.fffffffffffffp-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,   -0x1.fffffffffffffp-1,                 -0x1p+0,          0x0p+0, INEXACT)
// integer + ulp
T(RN,   0x1.eeeeeeeeeeeedp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RN,  -0x1.eeeeeeeeeeeedp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RN,    0x1.0000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RN,   -0x1.0000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RN,               0x1p-1022,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,              -0x1p-1022,                 -0x1p+0,          0x0p+0, INEXACT)
T(RZ,   0x1.eeeeeeeeeeeedp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RZ,  -0x1.eeeeeeeeeeeedp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RZ,    0x1.0000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RZ,   -0x1.0000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RZ,               0x1p-1022,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,              -0x1p-1022,                 -0x1p+0,          0x0p+0, INEXACT)
T(RU,   0x1.eeeeeeeeeeeedp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RU,  -0x1.eeeeeeeeeeeedp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RU,    0x1.0000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RU,   -0x1.0000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RU,               0x1p-1022,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,              -0x1p-1022,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,   0x1.eeeeeeeeeeeedp+50,   0x1.eeeeeeeeeeeecp+50,          0x0p+0, INEXACT)
T(RD,  -0x1.eeeeeeeeeeeedp+50,   -0x1.eeeeeeeeeeefp+50,          0x0p+0, INEXACT)
T(RD,    0x1.0000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RD,   -0x1.0000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RD,               0x1p-1022,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,              -0x1p-1022,                 -0x1p+0,          0x0p+0, INEXACT)
//  half way case, half way case +- ulp
T(RN,    0x1.fffffffffffffp-2,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,                  0x1p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,    0x1.0000000000001p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,   -0x1.fffffffffffffp-2,                 -0x1p+0,          0x0p+0, INEXACT)
T(RN,                 -0x1p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RN,   -0x1.0000000000001p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RN,    0x1.7ffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RN,                0x1.8p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RN,    0x1.8000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RN,   -0x1.7ffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RN,               -0x1.8p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RN,   -0x1.8000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RN,    0x1.3ffffffffffffp+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RN,                0x1.4p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RN,    0x1.4000000000001p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RN,   -0x1.3ffffffffffffp+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RN,               -0x1.4p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RN,   -0x1.4000000000001p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RN,   0x1.eeeeeeeeeeee7p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RN,   0x1.eeeeeeeeeeee8p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RN,   0x1.eeeeeeeeeeee9p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RN,  -0x1.eeeeeeeeeeee7p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RN,  -0x1.eeeeeeeeeeee8p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RN,  -0x1.eeeeeeeeeeee9p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RZ,    0x1.fffffffffffffp-2,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,                  0x1p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,    0x1.0000000000001p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,   -0x1.fffffffffffffp-2,                 -0x1p+0,          0x0p+0, INEXACT)
T(RZ,                 -0x1p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RZ,   -0x1.0000000000001p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RZ,    0x1.7ffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RZ,                0x1.8p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RZ,    0x1.8000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RZ,   -0x1.7ffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RZ,               -0x1.8p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RZ,   -0x1.8000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RZ,    0x1.3ffffffffffffp+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RZ,                0x1.4p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RZ,    0x1.4000000000001p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RZ,   -0x1.3ffffffffffffp+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RZ,               -0x1.4p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RZ,   -0x1.4000000000001p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RZ,   0x1.eeeeeeeeeeee7p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RZ,   0x1.eeeeeeeeeeee8p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RZ,   0x1.eeeeeeeeeeee9p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RZ,  -0x1.eeeeeeeeeeee7p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RZ,  -0x1.eeeeeeeeeeee8p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RZ,  -0x1.eeeeeeeeeeee9p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RU,    0x1.fffffffffffffp-2,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,                  0x1p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,    0x1.0000000000001p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,   -0x1.fffffffffffffp-2,                 -0x1p+0,          0x0p+0, INEXACT)
T(RU,                 -0x1p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RU,   -0x1.0000000000001p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RU,    0x1.7ffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RU,                0x1.8p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RU,    0x1.8000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RU,   -0x1.7ffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RU,               -0x1.8p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RU,   -0x1.8000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RU,    0x1.3ffffffffffffp+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RU,                0x1.4p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RU,    0x1.4000000000001p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RU,   -0x1.3ffffffffffffp+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RU,               -0x1.4p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RU,   -0x1.4000000000001p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RU,   0x1.eeeeeeeeeeee7p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RU,   0x1.eeeeeeeeeeee8p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RU,   0x1.eeeeeeeeeeee9p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RU,  -0x1.eeeeeeeeeeee7p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RU,  -0x1.eeeeeeeeeeee8p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RU,  -0x1.eeeeeeeeeeee9p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RD,    0x1.fffffffffffffp-2,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,                  0x1p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,    0x1.0000000000001p-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,   -0x1.fffffffffffffp-2,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,                 -0x1p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,   -0x1.0000000000001p-1,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,    0x1.7ffffffffffffp+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RD,                0x1.8p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RD,    0x1.8000000000001p+0,                  0x1p+0,          0x0p+0, INEXACT)
T(RD,   -0x1.7ffffffffffffp+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RD,               -0x1.8p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RD,   -0x1.8000000000001p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RD,    0x1.3ffffffffffffp+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RD,                0x1.4p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RD,    0x1.4000000000001p+1,                  0x1p+1,          0x0p+0, INEXACT)
T(RD,   -0x1.3ffffffffffffp+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RD,               -0x1.4p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RD,   -0x1.4000000000001p+1,               -0x1.8p+1,          0x0p+0, INEXACT)
T(RD,   0x1.eeeeeeeeeeee7p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RD,   0x1.eeeeeeeeeeee8p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RD,   0x1.eeeeeeeeeeee9p+48,    0x1.eeeeeeeeeeeep+48,          0x0p+0, INEXACT)
T(RD,  -0x1.eeeeeeeeeeee7p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RD,  -0x1.eeeeeeeeeeee8p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
T(RD,  -0x1.eeeeeeeeeeee9p+48,   -0x1.eeeeeeeeeeefp+48,          0x0p+0, INEXACT)
// random arguments between -100,100
T(RN,   -0x1.adeefb2b5006dp+3,               -0x1.cp+3,          0x0p+0, INEXACT)
T(RN,    0x1.1ce3efb825911p+5,               0x1.18p+5,          0x0p+0, INEXACT)
T(RN,    0x1.602e109de7505p+5,                0x1.6p+5,          0x0p+0, INEXACT)
T(RN,   -0x1.0b245fba96889p+5,               -0x1.1p+5,          0x0p+0, INEXACT)
T(RN,   -0x1.b171ee27084ddp+3,               -0x1.cp+3,          0x0p+0, INEXACT)
T(RN,   -0x1.f6eff1b093c41p+0,                 -0x1p+1,          0x0p+0, INEXACT)
T(RN,    0x1.ceaa3d18455f5p+4,                0x1.cp+4,          0x0p+0, INEXACT)
T(RN,    0x1.560914a51b239p+5,                0x1.5p+5,          0x0p+0, INEXACT)
T(RN,   -0x1.0ce901079de4dp+3,               -0x1.2p+3,          0x0p+0, INEXACT)
T(RN,   -0x1.7f35b3103b871p+5,               -0x1.8p+5,          0x0p+0, INEXACT)
// inf,nan, and subnormal number
T(RN,               0x1p-1074,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,              -0x1p-1074,                 -0x1p+0,          0x0p+0, INEXACT)
T(RN,                     inf,                     inf,          0x0p+0, 0)
T(RN,                    -inf,                    -inf,          0x0p+0, 0)
T(RN,                     nan,                     nan,          0x0p+0, 0)
T(RZ,               0x1p-1074,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,              -0x1p-1074,                 -0x1p+0,          0x0p+0, INEXACT)
T(RZ,                     inf,                     inf,          0x0p+0, 0)
T(RZ,                    -inf,                    -inf,          0x0p+0, 0)
T(RZ,                     nan,                     nan,          0x0p+0, 0)
T(RU,               0x1p-1074,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,              -0x1p-1074,                 -0x1p+0,          0x0p+0, INEXACT)
T(RU,                     inf,                     inf,          0x0p+0, 0)
T(RU,                    -inf,                    -inf,          0x0p+0, 0)
T(RU,                     nan,                     nan,          0x0p+0, 0)
T(RD,               0x1p-1074,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,              -0x1p-1074,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,                     inf,                     inf,          0x0p+0, 0)
T(RD,                    -inf,                    -inf,          0x0p+0, 0)
T(RD,                     nan,                     nan,          0x0p+0, 0)
T(RD,               0x1.2p+12,               0x1.2p+12,          0x0p+0, 0)
T(RD,                 0x1p+23,                 0x1p+23,          0x0p+0, 0)
T(RD,   0x1.ffffffffffffep+51,   0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RD,                 0x1p+52,                 0x1p+52,          0x0p+0, 0)
T(RD,   0x1.0000000000001p+52,   0x1.0000000000001p+52,          0x0p+0, 0)
T(RD, 0x1.fffffffffffeep+1014, 0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RD, 0x1.ffffffffffff7p+1014, 0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RD, 0x1.fffffffffffffp+1014, 0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RD,               0x1p+1015,               0x1p+1015,          0x0p+0, 0)
T(RD,              -0x1.2p+12,              -0x1.2p+12,          0x0p+0, 0)
T(RD,                -0x1p+23,                -0x1p+23,          0x0p+0, 0)
T(RD,  -0x1.ffffffffffffep+51,  -0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RD,                -0x1p+52,                -0x1p+52,          0x0p+0, 0)
T(RD,  -0x1.0000000000001p+52,  -0x1.0000000000001p+52,          0x0p+0, 0)
T(RD,-0x1.fffffffffffeep+1014,-0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RD,-0x1.ffffffffffff7p+1014,-0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RD,-0x1.fffffffffffffp+1014,-0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RD,              -0x1p+1015,              -0x1p+1015,          0x0p+0, 0)
T(RD, 0x1.ffffffffffffep-1023,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,                0x1.ep-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RD,                0x1.2p+2,                  0x1p+2,          0x0p+0, INEXACT)
T(RD,    0x1.fffffffffffffp+2,                0x1.cp+2,          0x0p+0, INEXACT)
T(RD,    0x1.0000000000001p+3,                  0x1p+3,          0x0p+0, INEXACT)
T(RD,    0x1.0000000000008p+9,                  0x1p+9,          0x0p+0, INEXACT)
T(RD,   0x1.0000000000001p+18,                 0x1p+18,          0x0p+0, INEXACT)
T(RD,   0x1.0000000000001p+23,                 0x1p+23,          0x0p+0, INEXACT)
T(RD,   0x1.ffffffffffffdp+51,   0x1.ffffffffffffcp+51,          0x0p+0, INEXACT)
T(RD,   0x1.fffffffffffffp+51,   0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RD,-0x1.ffffffffffffep-1023,                 -0x1p+0,          0x0p+0, INEXACT)
T(RD,               -0x1.2p+2,               -0x1.4p+2,          0x0p+0, INEXACT)
T(RD,   -0x1.fffffffffffffp+2,                 -0x1p+3,          0x0p+0, INEXACT)
T(RD,   -0x1.0000000000001p+3,               -0x1.2p+3,          0x0p+0, INEXACT)
T(RD,   -0x1.ffffffffffff8p+8,                 -0x1p+9,          0x0p+0, INEXACT)
T(RD,  -0x1.fffffffffffffp+17,                -0x1p+18,          0x0p+0, INEXACT)
T(RD,  -0x1.ffffffffffffdp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RD,  -0x1.ffffffffffffep+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RD,  -0x1.fffffffffffffp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RD,  -0x1.ffffffffffffdp+51,  -0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RD,  -0x1.fffffffffffffp+51,                -0x1p+52,          0x0p+0, INEXACT)
T(RD,                     nan,                     nan,          0x0p+0, 0)
T(RN,               0x1.2p+12,               0x1.2p+12,          0x0p+0, 0)
T(RN,                 0x1p+23,                 0x1p+23,          0x0p+0, 0)
T(RN,   0x1.ffffffffffffep+51,   0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RN,                 0x1p+52,                 0x1p+52,          0x0p+0, 0)
T(RN,   0x1.0000000000001p+52,   0x1.0000000000001p+52,          0x0p+0, 0)
T(RN, 0x1.fffffffffffeep+1014, 0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RN, 0x1.ffffffffffff7p+1014, 0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RN, 0x1.fffffffffffffp+1014, 0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RN,               0x1p+1015,               0x1p+1015,          0x0p+0, 0)
T(RN,              -0x1.2p+12,              -0x1.2p+12,          0x0p+0, 0)
T(RN,                -0x1p+23,                -0x1p+23,          0x0p+0, 0)
T(RN,  -0x1.ffffffffffffep+51,  -0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RN,                -0x1p+52,                -0x1p+52,          0x0p+0, 0)
T(RN,  -0x1.0000000000001p+52,  -0x1.0000000000001p+52,          0x0p+0, 0)
T(RN,-0x1.fffffffffffeep+1014,-0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RN,-0x1.ffffffffffff7p+1014,-0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RN,-0x1.fffffffffffffp+1014,-0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RN,              -0x1p+1015,              -0x1p+1015,          0x0p+0, 0)
T(RN, 0x1.ffffffffffffep-1023,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,                0x1.ep-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RN,                0x1.2p+2,                  0x1p+2,          0x0p+0, INEXACT)
T(RN,    0x1.2000000000001p+2,                  0x1p+2,          0x0p+0, INEXACT)
T(RN,    0x1.0000000000001p+3,                  0x1p+3,          0x0p+0, INEXACT)
T(RN,    0x1.0000000000008p+9,                  0x1p+9,          0x0p+0, INEXACT)
T(RN,   0x1.0000000000001p+18,                 0x1p+18,          0x0p+0, INEXACT)
T(RN,   0x1.0000000000001p+23,                 0x1p+23,          0x0p+0, INEXACT)
T(RN,   0x1.ffffffffffffdp+51,   0x1.ffffffffffffcp+51,          0x0p+0, INEXACT)
T(RN,   0x1.fffffffffffffp+51,   0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RN,   -0x1.2000000000001p+2,               -0x1.4p+2,          0x0p+0, INEXACT)
T(RN,   -0x1.fffffffffffffp+2,                 -0x1p+3,          0x0p+0, INEXACT)
T(RN,   -0x1.ffffffffffff8p+8,                 -0x1p+9,          0x0p+0, INEXACT)
T(RN,  -0x1.fffffffffffffp+17,                -0x1p+18,          0x0p+0, INEXACT)
T(RN,  -0x1.ffffffffffffdp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RN,  -0x1.ffffffffffffep+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RN,  -0x1.fffffffffffffp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RN,  -0x1.ffffffffffffdp+51,  -0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RN,  -0x1.fffffffffffffp+51,                -0x1p+52,          0x0p+0, INEXACT)
T(RN,                     nan,                     nan,          0x0p+0, 0)
T(RU,               0x1.2p+12,               0x1.2p+12,          0x0p+0, 0)
T(RU,                 0x1p+23,                 0x1p+23,          0x0p+0, 0)
T(RU,   0x1.ffffffffffffep+51,   0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RU,                 0x1p+52,                 0x1p+52,          0x0p+0, 0)
T(RU,   0x1.0000000000001p+52,   0x1.0000000000001p+52,          0x0p+0, 0)
T(RU, 0x1.fffffffffffeep+1014, 0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RU, 0x1.ffffffffffff7p+1014, 0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RU, 0x1.fffffffffffffp+1014, 0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RU,               0x1p+1015,               0x1p+1015,          0x0p+0, 0)
T(RU,              -0x1.2p+12,              -0x1.2p+12,          0x0p+0, 0)
T(RU,                -0x1p+23,                -0x1p+23,          0x0p+0, 0)
T(RU,  -0x1.ffffffffffffep+51,  -0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RU,                -0x1p+52,                -0x1p+52,          0x0p+0, 0)
T(RU,  -0x1.0000000000001p+52,  -0x1.0000000000001p+52,          0x0p+0, 0)
T(RU,-0x1.fffffffffffeep+1014,-0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RU,-0x1.ffffffffffff7p+1014,-0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RU,-0x1.fffffffffffffp+1014,-0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RU,              -0x1p+1015,              -0x1p+1015,          0x0p+0, 0)
T(RU,                0x1.ep-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RU,   0x1.0000000000001p+23,                 0x1p+23,          0x0p+0, INEXACT)
T(RU,   0x1.ffffffffffffdp+51,   0x1.ffffffffffffcp+51,          0x0p+0, INEXACT)
T(RU,   0x1.fffffffffffffp+51,   0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RU,  -0x1.ffffffffffffdp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RU,  -0x1.ffffffffffffep+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RU,  -0x1.fffffffffffffp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RU,  -0x1.ffffffffffffdp+51,  -0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RU,  -0x1.fffffffffffffp+51,                -0x1p+52,          0x0p+0, INEXACT)
T(RU,                     nan,                     nan,          0x0p+0, 0)
T(RZ,               0x1.2p+12,               0x1.2p+12,          0x0p+0, 0)
T(RZ,                 0x1p+23,                 0x1p+23,          0x0p+0, 0)
T(RZ,   0x1.ffffffffffffep+51,   0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RZ,                 0x1p+52,                 0x1p+52,          0x0p+0, 0)
T(RZ,   0x1.0000000000001p+52,   0x1.0000000000001p+52,          0x0p+0, 0)
T(RZ, 0x1.fffffffffffeep+1014, 0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RZ, 0x1.ffffffffffff7p+1014, 0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RZ, 0x1.fffffffffffffp+1014, 0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RZ,               0x1p+1015,               0x1p+1015,          0x0p+0, 0)
T(RZ,              -0x1.2p+12,              -0x1.2p+12,          0x0p+0, 0)
T(RZ,                -0x1p+23,                -0x1p+23,          0x0p+0, 0)
T(RZ,  -0x1.ffffffffffffep+51,  -0x1.ffffffffffffep+51,          0x0p+0, 0)
T(RZ,                -0x1p+52,                -0x1p+52,          0x0p+0, 0)
T(RZ,  -0x1.0000000000001p+52,  -0x1.0000000000001p+52,          0x0p+0, 0)
T(RZ,-0x1.fffffffffffeep+1014,-0x1.fffffffffffeep+1014,          0x0p+0, 0)
T(RZ,-0x1.ffffffffffff7p+1014,-0x1.ffffffffffff7p+1014,          0x0p+0, 0)
T(RZ,-0x1.fffffffffffffp+1014,-0x1.fffffffffffffp+1014,          0x0p+0, 0)
T(RZ,              -0x1p+1015,              -0x1p+1015,          0x0p+0, 0)
T(RZ, 0x1.ffffffffffffep-1023,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,                0x1.ep-1,                  0x0p+0,          0x0p+0, INEXACT)
T(RZ,                0x1.2p+2,                  0x1p+2,          0x0p+0, INEXACT)
T(RZ,    0x1.fffffffffffffp+2,                0x1.cp+2,          0x0p+0, INEXACT)
T(RZ,    0x1.0000000000001p+3,                  0x1p+3,          0x0p+0, INEXACT)
T(RZ,    0x1.0000000000008p+9,                  0x1p+9,          0x0p+0, INEXACT)
T(RZ,   0x1.0000000000001p+18,                 0x1p+18,          0x0p+0, INEXACT)
T(RZ,   0x1.0000000000001p+23,                 0x1p+23,          0x0p+0, INEXACT)
T(RZ,   0x1.ffffffffffffdp+51,   0x1.ffffffffffffcp+51,          0x0p+0, INEXACT)
T(RZ,   0x1.fffffffffffffp+51,   0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RZ,  -0x1.ffffffffffffdp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RZ,  -0x1.ffffffffffffep+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RZ,  -0x1.fffffffffffffp+22,                -0x1p+23,          0x0p+0, INEXACT)
T(RZ,  -0x1.ffffffffffffdp+51,  -0x1.ffffffffffffep+51,          0x0p+0, INEXACT)
T(RZ,  -0x1.fffffffffffffp+51,                -0x1p+52,          0x0p+0, INEXACT)
T(RZ,                     nan,                     nan,          0x0p+0, 0)
