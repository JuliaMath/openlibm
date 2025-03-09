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
// 0.7max,0.6max
T(RN, 0x1.666666p+127, 0x1.333334p+127, 0x1.d80a6ap+127,   0x1.b63f02p-7, INEXACT)
T(RZ, 0x1.666666p+127, 0x1.333334p+127, 0x1.d80a68p+127,  -0x1.f92704p-1, INEXACT)
T(RU, 0x1.666666p+127, 0x1.333334p+127, 0x1.d80a6ap+127,   0x1.b63f02p-7, INEXACT)
T(RD, 0x1.666666p+127, 0x1.333334p+127, 0x1.d80a68p+127,  -0x1.f92704p-1, INEXACT)
// tiny,huge = huge,tiny = huge
T(RN,          0x0p+0, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RN,        0x1p-149, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, INEXACT)
T(RN,        0x1p-126, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, INEXACT)
T(RN,          0x1p+0, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, INEXACT)
T(RN, 0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RN, 0x1.fffffep+127,        0x1p-149, 0x1.fffffep+127,          0x0p+0, INEXACT)
T(RN, 0x1.fffffep+127,        0x1p-126, 0x1.fffffep+127,          0x0p+0, INEXACT)
T(RN, 0x1.fffffep+127,          0x1p+0, 0x1.fffffep+127,          0x0p+0, INEXACT)
// pythagoras integers test
T(RN,  0x1.ffe39cp+22,     0x1.c48p+14,  0x1.ffe464p+22,          0x0p+0, 0)
T(RN,  0x1.ffe2bcp+22,     0x1.974p+15,  0x1.ffe544p+22,          0x0p+0, 0)
T(RN,  0x1.ffe15cp+22,     0x1.262p+16,  0x1.ffe6a4p+22,          0x0p+0, 0)
T(RN,  0x1.ffdf7cp+22,     0x1.80ap+16,  0x1.ffe884p+22,          0x0p+0, 0)
T(RN,  0x1.001f3ap+23,     0x1.6a2p+13,  0x1.001f4ap+23,          0x0p+0, 0)
T(RN,  0x1.001efap+23,    0x1.0f98p+15,  0x1.001f8ap+23,          0x0p+0, 0)
T(RN,  0x1.001e7ap+23,    0x1.c4a8p+15,  0x1.00200ap+23,          0x0p+0, 0)
T(RN,  0x1.001dbap+23,    0x1.3cdcp+16,  0x1.0020cap+23,          0x0p+0, 0)
T(RN,  0x1.001cbap+23,    0x1.9764p+16,  0x1.0021cap+23,          0x0p+0, 0)
T(RN,  0x1.004c86p+23,     0x1.6a4p+12,  0x1.004c8ap+23,          0x0p+0, 0)
T(RN,  0x1.004c56p+23,     0x1.c4dp+14,  0x1.004cbap+23,          0x0p+0, 0)
T(RN,  0x1.004be6p+23,    0x1.9788p+15,  0x1.004d2ap+23,          0x0p+0, 0)
T(RN,  0x1.004b36p+23,    0x1.2654p+16,  0x1.004ddap+23,          0x0p+0, 0)
T(RN,  0x1.004a46p+23,    0x1.80e4p+16,  0x1.004ecap+23,          0x0p+0, 0)
T(RN,  0x1.0079cap+23,     0x1.6a6p+13,  0x1.0079dap+23,          0x0p+0, 0)
T(RN,  0x1.00798ap+23,    0x1.0fc8p+15,  0x1.007a1ap+23,          0x0p+0, 0)
T(RN,  0x1.00790ap+23,    0x1.c4f8p+15,  0x1.007a9ap+23,          0x0p+0, 0)
T(RN,  0x1.00784ap+23,    0x1.3d14p+16,  0x1.007b5ap+23,          0x0p+0, 0)
T(RN,  0x1.00774ap+23,    0x1.97acp+16,  0x1.007c5ap+23,          0x0p+0, 0)
T(RN,  0x1.00a71ep+23,     0x1.6a8p+12,  0x1.00a722p+23,          0x0p+0, 0)
// radom argument in (-10,10)
T(RN,  -0x1.57f25cp+1,    0x1.c7d31p+2,   0x1.e72fc4p+2,   -0x1.2ed8fp-2, INEXACT)
T(RN,    0x1.19be7p+3,   -0x1.ab6d7p+2,   0x1.61a0ecp+3,  -0x1.a978c8p-6, INEXACT)
T(RN,  -0x1.5ac18ep+1,  -0x1.925982p-2,   0x1.5e6268p+1,  -0x1.55d62cp-3, INEXACT)
T(RN,   0x1.7221cep+2,   0x1.11a0d4p+3,   0x1.4a5602p+3,   0x1.73f5f4p-2, INEXACT)
T(RN,  -0x1.ae41a2p+0,  -0x1.329154p+3,    0x1.373fep+3,   0x1.1a3a8cp-3, INEXACT)
T(RN,   -0x1.0accfp+2,   0x1.d94512p-2,   0x1.0c6f6ap+2,  -0x1.526f22p-3, INEXACT)
T(RN,    -0x1.e564p+2,   0x1.c7cbf2p+2,   0x1.4ceca6p+3,  -0x1.e8dcc2p-4, INEXACT)
T(RN,  -0x1.3ec60ep+3,  -0x1.3fa3cep+3,   0x1.c36d4cp+3,   0x1.10fbdcp-5, INEXACT)
T(RN,  -0x1.236fd2p+2,   0x1.742432p+2,   0x1.d8ad9ap+2,  -0x1.2fb484p-3, INEXACT)
T(RN,   0x1.6f651ep+1,   0x1.3bfd78p+2,   0x1.6d817ep+2,  -0x1.ba32cep-2, INEXACT)
// nan's resutls
T(RN,             nan,          0x1p+0,             nan,          0x0p+0, 0)
T(RN,             nan,          0x1p+0,             nan,          0x0p+0, 0)
T(RN,             nan,             nan,             nan,          0x0p+0, 0)
T(RN,             nan,             nan,             nan,          0x0p+0, 0)
// inf result
T(RN,            -inf,             nan,             inf,          0x0p+0, 0)
T(RZ,             nan,            -inf,             inf,          0x0p+0, 0)
// inf result with snan argument raise invalid flag
T(RN,             nan,             inf,             inf,          0x0p+0, 0)
// overflow
T(RN,-0x1.fffffep+127, 0x1.fddddcp+127,             inf,          0x0p+0, INEXACT|OVERFLOW)
T(RZ,-0x1.fffffep+127, 0x1.fddddcp+127, 0x1.fffffep+127,         -0x1p+0, INEXACT|OVERFLOW)
T(RU,-0x1.fffffep+127, 0x1.fddddcp+127,             inf,          0x0p+0, INEXACT|OVERFLOW)
T(RD,-0x1.fffffep+127, 0x1.fddddcp+127, 0x1.fffffep+127,         -0x1p+0, INEXACT|OVERFLOW)
// subnormal number
T(RN,          0x0p+0,        0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RN,        0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RN,        0x1p-126,          0x0p+0,        0x1p-126,          0x0p+0, 0)
T(RN,          0x0p+0,       -0x1p-126,        0x1p-126,          0x0p+0, 0)
T(RN,        0x1p-149,        0x1p-149,        0x1p-149,  -0x1.a8279ap-2, INEXACT|UNDERFLOW)
T(RN,       -0x1p-148,       -0x1p-148,      0x1.8p-148,   0x1.5f619ap-3, INEXACT|UNDERFLOW)
T(RD,          0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RD,          0x0p+0,        0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RD,          0x0p+0, 0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RD,          0x0p+0,          0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RD,          0x0p+0,        0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RD,          0x0p+0, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RD,          0x0p+0,             inf,             inf,          0x0p+0, 0)
T(RD,          0x0p+0,         -0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RD,          0x0p+0,       -0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RD,          0x0p+0,-0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RD,          0x0p+0,         -0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RD,          0x0p+0,       -0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RD,          0x0p+0,-0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RD,          0x0p+0,            -inf,             inf,          0x0p+0, 0)
T(RD,        0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RD, 0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RD,          0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RD,        0x1.8p+3,       -0x1.4p+2,        0x1.ap+3,          0x0p+0, 0)
T(RD,        0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RD, 0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RD,             inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RD,             inf,             nan,             inf,          0x0p+0, 0)
T(RD,             inf,             nan,             inf,          0x0p+0, 0)
T(RD,             nan,             inf,             inf,          0x0p+0, 0)
T(RD,             nan,            -inf,             inf,          0x0p+0, 0)
T(RD,         -0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RD,       -0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RD,-0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RD,         -0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RD,       -0x1.8p+1,         -0x1p+2,        0x1.4p+2,          0x0p+0, 0)
T(RD,       -0x1.8p+4,        0x1.cp+2,        0x1.9p+4,          0x0p+0, 0)
T(RD,       -0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RD,-0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RD,            -inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RD,            -inf,             nan,             inf,          0x0p+0, 0)
T(RD,            -inf,             nan,             inf,          0x0p+0, 0)
T(RD,             nan,             inf,             inf,          0x0p+0, 0)
T(RD,             nan,            -inf,             inf,          0x0p+0, 0)
T(RD,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RD,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RD,        0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RD,      0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RD,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RD,             nan,       -0x1p-149,             nan,          0x0p+0, 0)
T(RD,             nan,     -0x1.8p-148,             nan,          0x0p+0, 0)
T(RD,       -0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RD,     -0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RD,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RD,             nan,        0x1p-149,             nan,          0x0p+0, 0)
T(RD,             nan,      0x1.8p-148,             nan,          0x0p+0, 0)
T(RN,          0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RN,          0x0p+0, 0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RN,          0x0p+0,          0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RN,          0x0p+0,        0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RN,          0x0p+0,             inf,             inf,          0x0p+0, 0)
T(RN,          0x0p+0,         -0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RN,          0x0p+0,       -0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RN,          0x0p+0,-0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RN,          0x0p+0,         -0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RN,          0x0p+0,       -0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RN,          0x0p+0,-0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RN,          0x0p+0,            -inf,             inf,          0x0p+0, 0)
T(RN, 0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RN,          0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RN,        0x1.8p+3,       -0x1.4p+2,        0x1.ap+3,          0x0p+0, 0)
T(RN,        0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RN,             inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RN,             inf,             nan,             inf,          0x0p+0, 0)
T(RN,             inf,             nan,             inf,          0x0p+0, 0)
T(RN,             nan,             inf,             inf,          0x0p+0, 0)
T(RN,             nan,            -inf,             inf,          0x0p+0, 0)
T(RN,         -0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RN,       -0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RN,-0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RN,         -0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RN,       -0x1.8p+1,         -0x1p+2,        0x1.4p+2,          0x0p+0, 0)
T(RN,       -0x1.8p+4,        0x1.cp+2,        0x1.9p+4,          0x0p+0, 0)
T(RN,       -0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RN,-0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RN,            -inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RN,            -inf,             nan,             inf,          0x0p+0, 0)
T(RN,             nan,            -inf,             inf,          0x0p+0, 0)
T(RN,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RN,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RN,        0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RN,      0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RN,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RN,             nan,       -0x1p-149,             nan,          0x0p+0, 0)
T(RN,             nan,     -0x1.8p-148,             nan,          0x0p+0, 0)
T(RN,       -0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RN,     -0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RN,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RN,             nan,        0x1p-149,             nan,          0x0p+0, 0)
T(RN,             nan,      0x1.8p-148,             nan,          0x0p+0, 0)
T(RU,          0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RU,          0x0p+0,        0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RU,          0x0p+0, 0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RU,          0x0p+0,          0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RU,          0x0p+0,        0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RU,          0x0p+0, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RU,          0x0p+0,             inf,             inf,          0x0p+0, 0)
T(RU,          0x0p+0,         -0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RU,          0x0p+0,       -0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RU,          0x0p+0,-0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RU,          0x0p+0,         -0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RU,          0x0p+0,       -0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RU,          0x0p+0,-0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RU,          0x0p+0,            -inf,             inf,          0x0p+0, 0)
T(RU,        0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RU, 0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RU,          0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RU,        0x1.8p+3,       -0x1.4p+2,        0x1.ap+3,          0x0p+0, 0)
T(RU,        0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RU, 0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RU,             inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RU,             inf,             nan,             inf,          0x0p+0, 0)
T(RU,             inf,             nan,             inf,          0x0p+0, 0)
T(RU,             nan,             inf,             inf,          0x0p+0, 0)
T(RU,             nan,            -inf,             inf,          0x0p+0, 0)
T(RU,         -0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RU,       -0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RU,-0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RU,         -0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RU,       -0x1.8p+1,         -0x1p+2,        0x1.4p+2,          0x0p+0, 0)
T(RU,       -0x1.8p+4,        0x1.cp+2,        0x1.9p+4,          0x0p+0, 0)
T(RU,       -0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RU,-0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RU,            -inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RU,            -inf,             nan,             inf,          0x0p+0, 0)
T(RU,            -inf,             nan,             inf,          0x0p+0, 0)
T(RU,             nan,             inf,             inf,          0x0p+0, 0)
T(RU,             nan,            -inf,             inf,          0x0p+0, 0)
T(RU,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RU,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RU,        0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RU,      0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RU,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RU,             nan,       -0x1p-149,             nan,          0x0p+0, 0)
T(RU,             nan,     -0x1.8p-148,             nan,          0x0p+0, 0)
T(RU,       -0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RU,     -0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RU,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RU,             nan,        0x1p-149,             nan,          0x0p+0, 0)
T(RU,             nan,      0x1.8p-148,             nan,          0x0p+0, 0)
T(RZ,          0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RZ,          0x0p+0,        0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RZ,          0x0p+0, 0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RZ,          0x0p+0,          0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RZ,          0x0p+0,        0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RZ,          0x0p+0, 0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RZ,          0x0p+0,             inf,             inf,          0x0p+0, 0)
T(RZ,          0x0p+0,         -0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RZ,          0x0p+0,       -0x1p-149,        0x1p-149,          0x0p+0, 0)
T(RZ,          0x0p+0,-0x1.fffffcp-127, 0x1.fffffcp-127,          0x0p+0, 0)
T(RZ,          0x0p+0,         -0x1p+0,          0x1p+0,          0x0p+0, 0)
T(RZ,          0x0p+0,       -0x1p+127,        0x1p+127,          0x0p+0, 0)
T(RZ,          0x0p+0,-0x1.fffffep+127, 0x1.fffffep+127,          0x0p+0, 0)
T(RZ,          0x0p+0,            -inf,             inf,          0x0p+0, 0)
T(RZ,        0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RZ, 0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RZ,          0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RZ,        0x1.8p+3,       -0x1.4p+2,        0x1.ap+3,          0x0p+0, 0)
T(RZ,        0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RZ, 0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RZ,             inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RZ,             inf,             nan,             inf,          0x0p+0, 0)
T(RZ,             inf,             nan,             inf,          0x0p+0, 0)
T(RZ,             nan,             inf,             inf,          0x0p+0, 0)
T(RZ,             nan,            -inf,             inf,          0x0p+0, 0)
T(RZ,         -0x0p+0,          0x0p+0,          0x0p+0,          0x0p+0, 0)
T(RZ,       -0x1p-149,          0x0p+0,        0x1p-149,          0x0p+0, 0)
T(RZ,-0x1.fffffcp-127,          0x0p+0, 0x1.fffffcp-127,          0x0p+0, 0)
T(RZ,         -0x1p+0,          0x0p+0,          0x1p+0,          0x0p+0, 0)
T(RZ,       -0x1.8p+1,         -0x1p+2,        0x1.4p+2,          0x0p+0, 0)
T(RZ,       -0x1.8p+4,        0x1.cp+2,        0x1.9p+4,          0x0p+0, 0)
T(RZ,       -0x1p+127,          0x0p+0,        0x1p+127,          0x0p+0, 0)
T(RZ,-0x1.fffffep+127,          0x0p+0, 0x1.fffffep+127,          0x0p+0, 0)
T(RZ,            -inf,          0x0p+0,             inf,          0x0p+0, 0)
T(RZ,            -inf,             nan,             inf,          0x0p+0, 0)
T(RZ,            -inf,             nan,             inf,          0x0p+0, 0)
T(RZ,             nan,             inf,             inf,          0x0p+0, 0)
T(RZ,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RZ,          0x0p+0,             nan,             nan,          0x0p+0, 0)
T(RZ,        0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RZ,      0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RZ,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RZ,             nan,       -0x1p-149,             nan,          0x0p+0, 0)
T(RZ,             nan,     -0x1.8p-148,             nan,          0x0p+0, 0)
T(RZ,       -0x1p-149,             nan,             nan,          0x0p+0, 0)
T(RZ,     -0x1.8p-148,             nan,             nan,          0x0p+0, 0)
T(RZ,             nan,          0x0p+0,             nan,          0x0p+0, 0)
T(RZ,             nan,        0x1p-149,             nan,          0x0p+0, 0)
T(RZ,             nan,      0x1.8p-148,             nan,          0x0p+0, 0)
