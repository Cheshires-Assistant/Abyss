 /** Author:       Plyashkevich Viatcheslav <plyashkevich@yandex.ru> 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library or "Lesser" General Public License version 3.0 (LGPLv3)
 * All rights reserved. 
 */


#include "VoiceDetector.hpp"
#include <stdlib.h>

enum
{
	D_LOW_THRESHOLD =   16,
	D_HIGH_THRESHOLD = 60
};

enum 
{
	D_NOT_SPEECH = -1,
	D_SPEECH = 1
};


VoiceDetector :: Wu  VoiceDetector :: wu =
{
     0x0000,  0x7FFF,  0x0000,  0x7FFF,  0x0000,  0x7FFF,  0x0C8B,  0x7F61,
     0x18F8,  0x7D89,  0x2527,  0x7A7C,  0x18F8,  0x7D89,  0x30FB,  0x7641,
     0x471C,  0x6A6D,  0x2527,  0x7A7C,  0x471C,  0x6A6D,  0x62F1,  0x5133,
     0x30FB,  0x7641,  0x5A82,  0x5A82,  0x7641,  0x30FB,  0x3C56,  0x70E2,
     0x6A6D,  0x471C,  0x7F61,  0x0C8B,  0x471C,  0x6A6D,  0x7641,  0x30FB,
     0x7D89, -0x18F8,  0x5133,  0x62F1,  0x7D89,  0x18F8,  0x70E2, -0x3C56,
     0x5A82,  0x5A82,  0x7FFF,  0x0000,  0x5A82, -0x5A82,  0x62F1,  0x5133,
     0x7D89, -0x18F8,  0x3C56, -0x70E2,  0x6A6D,  0x471C,  0x7641, -0x30FB,
     0x18F8, -0x7D89,  0x70E2,  0x3C56,  0x6A6D, -0x471C, -0x0C8B, -0x7F61,
     0x7641,  0x30FB,  0x5A82, -0x5A82, -0x30FB, -0x7641,  0x7A7C,  0x2527,
     0x471C, -0x6A6D, -0x5133, -0x62F1,  0x7D89,  0x18F8,  0x30FB, -0x7641,
    -0x6A6D, -0x471C,  0x7F61,  0x0C8B,  0x18F8, -0x7D89, -0x7A7C, -0x2527,
     0x0000,  0x7FFF,  0x0000,  0x7FFF,  0x0000,  0x7FFF,  0x30FB,  0x7641,
     0x5A82,  0x5A82,  0x7641,  0x30FB,  0x5A82,  0x5A82,  0x7FFF,  0x0000,
     0x5A82, -0x5A82,  0x7641,  0x30FB,  0x5A82, -0x5A82, -0x30FB, -0x7641
};

# define DIGIT_PROCESS(i, m, j)                                                   \
    do {                                                                    \
        unsigned _ = (i);                                                   \
        _ = ((_ & 0x33333333) <<  2) | ((_ & ~0x33333333) >>  2);           \
        _ = ((_ & 0x0F0F0F0F) <<  4) | ((_ & ~0x0F0F0F0F) >>  4);           \
        _ = ((_ & 0x00FF00FF) <<  8) | ((_ & ~0x00FF00FF) >>  8);           \
        _ = ((_ & 0x0000FFFF) << 16) | ((_ & ~0x0000FFFF) >> 16);           \
        (j) = _ >> (m);                                                     \
    } while (0)

#define MPY16X32R(x,y) \
    (((int)((short)(x) * (unsigned short)(y) + 0x4000) >> 15) + \
     ((int)((short)(x) * (short)((y) >> 16)) << 1))

void VoiceDetector::fft16x32(const short * ptr_w, int npoints, int * ptr_x, int * ptr_y)
{
    int i, j, l1, l2, h2, predj, tw_offset, stride, fft_jmp;
    int xt0_0, yt0_0, xt1_0, yt1_0, xt2_0, yt2_0;
    int xt0_1, yt0_1, xt1_1, yt1_1, xt2_1, yt2_1;
    int xh0_0, xh1_0, xh20_0, xh21_0, xl0_0, xl1_0, xl20_0, xl21_0;
    int xh0_1, xh1_1, xh20_1, xh21_1, xl0_1, xl1_1, xl20_1, xl21_1;
    int x_0, x_1, x_2, x_3, x_l1_0, x_l1_1, x_l1_2, x_l1_3, x_l2_0, x_l2_1;
    int xh0_2, xh1_2, xl0_2, xl1_2, xh0_3, xh1_3, xl0_3, xl1_3;
    int x_4, x_5, x_6, x_7, x_l2_2, x_l2_3, x_h2_0, x_h2_1, x_h2_2, x_h2_3;
    int x_8, x_9, x_a, x_b, x_c, x_d, x_e, x_f;
    short si10, si20, si30, co10, co20, co30;
    short si11, si21, si31, co11, co21, co31;
    const short *w;
    int *x, *x2, *x0;
    int * y0, * y1, * y2, *y3;
    int n00, n10, n20, n30, n01, n11, n21, n31;
    int n02, n12, n22, n32, n03, n13, n23, n33;
    int n0, j0;
    int radix, m;
    int y0r, y0i,  y4r, y4i;
    int norm;

    for (i = 31, m = 1; (npoints & (1 << i)) == 0; i--, m++) ;
    radix     = m & 1 ? 2 :  4;
    norm      = m - 2;

    stride     =   npoints;
    tw_offset  =   0;
    fft_jmp    =   6 * stride;

    while (stride > 4)
    {
        j         = 0;
        fft_jmp >>= 2;

        h2 = stride >> 1;
        l1 = stride;
        l2 = stride + (stride >> 1);

        x = ptr_x;
        w = ptr_w + tw_offset;
        tw_offset += fft_jmp;

        stride  >>= 2;

        for (i = 0; i < (npoints >> 3); i ++)
        {
            co10 = w[j+1];    si10 = w[j+0];
            co20 = w[j+3];    si20 = w[j+2];
            co30 = w[j+5];    si30 = w[j+4];
            co11 = w[j+7];    si11 = w[j+6];
            co21 = w[j+9];    si21 = w[j+8];
            co31 = w[j+11];   si31 = w[j+10];

            x_0 = x[0];       x_1 = x[1];
            x_2 = x[2];       x_3 = x[3];

            x_l1_0 = x[l1  ]; x_l1_1 = x[l1+1];
            x_l1_2 = x[l1+2]; x_l1_3 = x[l1+3];

            x_l2_0 = x[l2  ]; x_l2_1 = x[l2+1];
            x_l2_2 = x[l2+2]; x_l2_3 = x[l2+3];

            x_h2_0 = x[h2  ]; x_h2_1 = x[h2+1];
            x_h2_2 = x[h2+2]; x_h2_3 = x[h2+3];

            xh0_0  = x_0    + x_l1_0;        xh1_0  = x_1    + x_l1_1;
            xh0_1  = x_2    + x_l1_2;        xh1_1  = x_3    + x_l1_3;

            xl0_0  = x_0    - x_l1_0;        xl1_0  = x_1    - x_l1_1;
            xl0_1  = x_2    - x_l1_2;        xl1_1  = x_3    - x_l1_3;

            xh20_0 = x_h2_0 + x_l2_0;        xh21_0 = x_h2_1 + x_l2_1;
            xh20_1 = x_h2_2 + x_l2_2;        xh21_1 = x_h2_3 + x_l2_3;

            xl20_0 = x_h2_0 - x_l2_0;        xl21_0 = x_h2_1 - x_l2_1;
            xl20_1 = x_h2_2 - x_l2_2;        xl21_1 = x_h2_3 - x_l2_3;

            x0 = x;
            x2 = x0;

            j += 12;
            x += 4;

            predj = (j - fft_jmp);
            if (!predj) x += fft_jmp;
            if (!predj) j = 0;

            y0r   = xh0_0 + xh20_0; y0i   = xh1_0 + xh21_0;
            y4r   = xh0_1 + xh20_1; y4i   = xh1_1 + xh21_1;

            xt0_0 = xh0_0 - xh20_0;  yt0_0 = xh1_0 - xh21_0;
            xt0_1 = xh0_1 - xh20_1;  yt0_1 = xh1_1 - xh21_1;

            xt1_0 = xl0_0 + xl21_0;  yt2_0 = xl1_0 + xl20_0;
            xt2_0 = xl0_0 - xl21_0;  yt1_0 = xl1_0 - xl20_0;

            xt1_1 = xl0_1 + xl21_1;  yt2_1 = xl1_1 + xl20_1;
            xt2_1 = xl0_1 - xl21_1;  yt1_1 = xl1_1 - xl20_1;

            x2[0] = y0r;             x2[1] = y0i;
            x2[2] = y4r;             x2[3] = y4i;

            x2[h2  ] = MPY16X32R(si10 , yt1_0) + MPY16X32R(co10 , xt1_0);
            x2[h2+1] = MPY16X32R(co10 , yt1_0) - MPY16X32R(si10 , xt1_0);

            x2[h2+2] = MPY16X32R(si11 , yt1_1) + MPY16X32R(co11 , xt1_1);
            x2[h2+3] = MPY16X32R(co11 , yt1_1) - MPY16X32R(si11 , xt1_1);

            x2[l1  ] = MPY16X32R(si20 , yt0_0) + MPY16X32R(co20 , xt0_0);
            x2[l1+1] = MPY16X32R(co20 , yt0_0) - MPY16X32R(si20 , xt0_0);

            x2[l1+2] = MPY16X32R(si21 , yt0_1) + MPY16X32R(co21 , xt0_1);
            x2[l1+3] = MPY16X32R(co21 , yt0_1) - MPY16X32R(si21 , xt0_1);

            x2[l2  ] = MPY16X32R(si30 , yt2_0) + MPY16X32R(co30 , xt2_0);
            x2[l2+1] = MPY16X32R(co30 , yt2_0) - MPY16X32R(si30 , xt2_0);

            x2[l2+2] = MPY16X32R(si31 , yt2_1) + MPY16X32R(co31 , xt2_1);
            x2[l2+3] = MPY16X32R(co31 , yt2_1) - MPY16X32R(si31 , xt2_1);
        }
    }

    y0 = ptr_y;
    y2 = ptr_y + (int) npoints;
    x0 = ptr_x;
    x2 = ptr_x + (int) (npoints >> 1);

    if (radix == 2)
    {
        y1 = y0 + (int) (npoints >> 2);
        y3 = y2 + (int) (npoints >> 2);
        l1 = norm + 1;
        j0 = 8;
        n0 = npoints >> 1;
    }
    else
    {
        y1 = y0 + (int) (npoints >> 1);
        y3 = y2 + (int) (npoints >> 1);
        l1 = norm + 2;
        j0 = 4;
        n0 = npoints >> 2;
    }

    j = 0;

    for (i = 0; i < npoints; i += 8)
    {
        DIGIT_PROCESS(j, l1, h2);

        x_0 = x0[0]; x_1 = x0[1];
        x_2 = x0[2]; x_3 = x0[3];
        x_4 = x0[4]; x_5 = x0[5];
        x_6 = x0[6]; x_7 = x0[7];
        x0 += 8;

        xh0_0 = x_0 + x_4; xh1_0 = x_1 + x_5;
        xl0_0 = x_0 - x_4; xl1_0 = x_1 - x_5;
        xh0_1 = x_2 + x_6; xh1_1 = x_3 + x_7;
        xl0_1 = x_2 - x_6; xl1_1 = x_3 - x_7;

        n00 = xh0_0 + xh0_1; n01 = xh1_0 + xh1_1;
        n10 = xl0_0 + xl1_1; n11 = xl1_0 - xl0_1;
        n20 = xh0_0 - xh0_1; n21 = xh1_0 - xh1_1;
        n30 = xl0_0 - xl1_1; n31 = xl1_0 + xl0_1;

        if (radix == 2)
        {
           n00 = x_0 + x_2;     n01 = x_1 + x_3;
           n20 = x_0 - x_2;     n21 = x_1 - x_3;
           n10 = x_4 + x_6;     n11 = x_5 + x_7;
           n30 = x_4 - x_6;     n31 = x_5 - x_7;
        }

        y0[2*h2] = n00;   y0[2*h2 + 1] = n01;
        y1[2*h2] = n10;   y1[2*h2 + 1] = n11;
        y2[2*h2] = n20;   y2[2*h2 + 1] = n21;
        y3[2*h2] = n30;   y3[2*h2 + 1] = n31;

        x_8 = x2[0]; x_9 = x2[1];
        x_a = x2[2]; x_b = x2[3];
        x_c = x2[4]; x_d = x2[5];
        x_e = x2[6]; x_f = x2[7];
        x2 += 8;

        xh0_2 = x_8 + x_c; xh1_2  = x_9 + x_d;
        xl0_2 = x_8 - x_c; xl1_2  = x_9 - x_d;
        xh0_3 = x_a + x_e; xh1_3 = x_b + x_f;
        xl0_3 = x_a - x_e; xl1_3 = x_b - x_f;

        n02 = xh0_2 + xh0_3; n03 = xh1_2 + xh1_3;
        n12 = xl0_2 + xl1_3; n13 = xl1_2 - xl0_3;
        n22 = xh0_2 - xh0_3; n23 = xh1_2 - xh1_3;
        n32 = xl0_2 - xl1_3; n33 = xl1_2 + xl0_3;

        if (radix == 2)
        {
          n02 = x_8 + x_a;     n03 = x_9 + x_b;
          n22 = x_8 - x_a;     n23 = x_9 - x_b;
          n12 = x_c + x_e;     n13 = x_d + x_f;
          n32 = x_c - x_e;     n33 = x_d - x_f;
        }

        y0[2*h2+2] = n02;   y0[2*h2+3] = n03;
        y1[2*h2+2] = n12;   y1[2*h2+3] = n13;
        y2[2*h2+2] = n22;   y2[2*h2+3] = n23;
        y3[2*h2+2] = n32;   y3[2*h2+3] = n33;

        j += j0;
        if (j == n0)
        {
          j  += n0;
          x0 += (int) npoints >> 1;
          x2 += (int) npoints >> 1;
        }
      }
}



int VoiceDetector::voiceProcessing(short *ptr_x)
{
	for(int ii = 0; ii < FOURIE_COEFF << 1; ++ii)
         input.inp[ii] = ptr_x[ii];
             
    fft16x32(VoiceDetector::wu.w, FOURIE_COEFF, input.inp, inner.inr);

    // Find max harmonic's value
    int storeMaxValue = 0;

    for(int ii = 0; ii < FOURIE_COEFF << 1; ++ii)
	   if(storeMaxValue < abs(inner.inr[ii]))
	   {
	       storeMaxValue = abs(inner.inr[ii]);
	   }
    
    #if 1
    int shiftCount = 0;
    while(storeMaxValue > 32767){
	      storeMaxValue >>= 1,
	      ++shiftCount;
    }
    
    if(shiftCount > 0){
	      for(int ii = 0; ii < FOURIE_COEFF << 1; ++ii)	     
		        inner.inr[ii] >>= shiftCount;
	     
    }
    #endif          
    
    for(int ii = 0, jj = FOURIE_COEFF - 1; ii < FOURIE_COEFF << 1; ii+=2, --jj){         
         int temp1 = inner.inr[ii] * inner.inr[ii];
             temp1 >>= 1;
         int temp2 = inner.inr[ii + 1] * inner.inr[ii + 1];
             temp2 >>= 1;             
         input.inp[jj]      = temp1 + temp2;             
    }

    // find average value
    long int AverageValue = input.inp[0];

	for(int ii = 1; ii < FOURIE_COEFF; ++ii){
			AverageValue += input.inp[ii] + 1;
			AverageValue >>= 1;
	}
    ///////////////////////////////////////////

    int Max0Val = 0, indxMax0Val;

    for(int ii = 1; ii < FOURIE_COEFF - (FOURIE_COEFF >> 2); ++ii)   
	    if(input.inp[ii] > Max0Val){
		      Max0Val = input.inp[ii];
		      indxMax0Val = ii;
	    }
   
    int Max1Val = 0, indxMax1Val;

    for(int ii = 1; ii < FOURIE_COEFF - (FOURIE_COEFF >> 2); ++ii)
	    if(input.inp[ii] > Max1Val){
		      Max1Val = input.inp[ii];
		      indxMax1Val = ii;
	    }
       
    if((indxMax0Val > (FOURIE_COEFF >> 2)) ||
       (indxMax1Val > (FOURIE_COEFF >> 2)))
		return D_NOT_SPEECH;

	int Rel = ((Max0Val + Max1Val) >> 1) / AverageValue;

    if((Rel < D_LOW_THRESHOLD) || (Rel > D_HIGH_THRESHOLD))
      return D_NOT_SPEECH;

    return D_SPEECH;
}

bool VoiceDetector::voiceDetection(short *in)
{
	if(voiceProcessing(in) > 0){
		countSpeechDetection = TIME_DURATION;
	}
	else{
		--countSpeechDetection;
	}	

	if(countSpeechDetection > 1){
		for(int ii = 0; ii < FOURIE_COEFF << 1; ++ii)
			storePrevious[ii] = in[ii];
		return true;
	}
	else{
		for(int ii = 0; ii < FOURIE_COEFF << 1; ++ii)
			storePrevious[ii] = 0;
		countSpeechDetection = -1;
		return false;
	}
}
