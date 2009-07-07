#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "arch.h"
#include "dct.h"


#define PI  3.1415926535897932384626433832795
#define IDCTI_AMP 512


uint64_t dctclk = 0;
uint64_t idctclk = 0;


// DCT basis functions coeficients
int dct_tbl_i[8][8];

CACHE_ALIGN short dct_tbl_s[8][8] =
{
	{362, 362, 362, 362, 362, 362, 362, 362},
	{502, 425, 284,  99, -99,-284,-425,-502},
	{473, 195,-195,-473,-473,-195, 195, 473},
	{425, -99,-502,-284, 284, 502,  99,-425},
	{362,-362,-362, 362, 362,-362,-362, 362},
	{284,-502, 99,  425,-425, -99, 502,-284},
	{195,-473, 473,-195,-195, 473,-473, 195},
	{ 99,-284, 425,-502, 502,-425, 284, -99}
};

CACHE_ALIGN short idct_tbl_s[8][8] =
{
	{362, 502, 473, 425, 362, 284, 195,  99},
	{362, 425, 195, -99,-362,-502,-473,-284},
	{362, 284,-195,-502,-362,  99, 473, 425},
	{362,  99,-473,-284, 362, 425,-195,-502},
	{362, -99,-473, 284, 362,-425,-195, 502},
	{362,-284,-195, 502,-362, -99, 473,-425},
	{362,-425, 195,  99,-362, 502,-473, 284},
	{362,-502, 473,-425, 362,-284, 195, -99}
};

/******************************************************************************
**  dct
**  --------------------------------------------------------------------------
**  Fast DCT - Discrete Cosine Transform.
**  This function converts 8x8 pixel block into frequencies.
**  Lowest frequencies are at the upper-left corner.
**  The input and output could point at the same array, in this case the data
**  will be overwritten.
**  
**  ARGUMENTS:
**      pixels  - 8x8 pixel array;
**      data    - 8x8 freq block;
**
**  RETURN: -
******************************************************************************/
void dct(conv pixels[8][8], conv data[8][8])
{
	short rows[8][8];
	unsigned i;

	static const int
				c1 = 1004,  /* cos(pi/16) << 10 */
				s1 = 200,   /* sin(pi/16) */
				c3 = 851,   /* cos(3pi/16) << 10 */
				s3 = 569,   /* sin(3pi/16) << 10 */
				r2c6 = 554, /* sqrt(2)*cos(6pi/16) << 10 */
				r2s6 = 1337,/* sqrt(2)*sin(6pi/16) << 10 */
				r2 = 181;   /* sqrt(2) << 7*/

	uint64_t a = __rdtsc();

	/* transform rows */
	for (i = 0; i < 8; i++)
	{
		int x0,x1,x2,x3,x4,x5,x6,x7,x8;

		x0 = pixels[i][0];
		x1 = pixels[i][1];
		x2 = pixels[i][2];
		x3 = pixels[i][3];
		x4 = pixels[i][4];
		x5 = pixels[i][5];
		x6 = pixels[i][6];
		x7 = pixels[i][7];

		/* Stage 1 */
		x8=x7+x0;
		x0-=x7;
		x7=x1+x6;
		x1-=x6;
		x6=x2+x5;
		x2-=x5;
		x5=x3+x4;
		x3-=x4;

		/* Stage 2 */
		x4=x8+x5;
		x8-=x5;
		x5=x7+x6;
		x7-=x6;
		x6=c1*(x1+x2);
		x2=(-s1-c1)*x2+x6;
		x1=(s1-c1)*x1+x6;
		x6=c3*(x0+x3);
		x3=(-s3-c3)*x3+x6;
		x0=(s3-c3)*x0+x6;

		/* Stage 3 */
		x6=x4+x5;
		x4-=x5;
		x5=r2c6*(x7+x8);
		x7=(-r2s6-r2c6)*x7+x5;
		x8=(r2s6-r2c6)*x8+x5;
		x5=x0+x2;
		x0-=x2;
		x2=x3+x1;
		x3-=x1;

		/* Stage 4 and output */
		rows[i][0]=x6;
		rows[i][4]=x4;
		rows[i][2]=x8>>10;
		rows[i][6]=x7>>10;
		rows[i][7]=(x2-x5)>>10;
		rows[i][1]=(x2+x5)>>10;
		rows[i][3]=(x3*r2)>>17;
		rows[i][5]=(x0*r2)>>17;
	}

	/* transform columns */
	for (i = 0; i < 8; i++)
	{
		int x0,x1,x2,x3,x4,x5,x6,x7,x8;

		x0 = rows[0][i];
		x1 = rows[1][i];
		x2 = rows[2][i];
		x3 = rows[3][i];
		x4 = rows[4][i];
		x5 = rows[5][i];
		x6 = rows[6][i];
		x7 = rows[7][i];

		/* Stage 1 */
		x8=x7+x0;
		x0-=x7;
		x7=x1+x6;
		x1-=x6;
		x6=x2+x5;
		x2-=x5;
		x5=x3+x4;
		x3-=x4;

		/* Stage 2 */
		x4=x8+x5;
		x8-=x5;
		x5=x7+x6;
		x7-=x6;
		x6=c1*(x1+x2);
		x2=(-s1-c1)*x2+x6;
		x1=(s1-c1)*x1+x6;
		x6=c3*(x0+x3);
		x3=(-s3-c3)*x3+x6;
		x0=(s3-c3)*x0+x6;

		/* Stage 3 */
		x6=x4+x5;
		x4-=x5;
		x5=r2c6*(x7+x8);
		x7=(-r2s6-r2c6)*x7+x5;
		x8=(r2s6-r2c6)*x8+x5;
		x5=x0+x2;
		x0-=x2;
		x2=x3+x1;
		x3-=x1;

		/* Stage 4 and output */
		data[0][i]=((x6+16)>>3);
		data[4][i]=((x4+16)>>3);
		data[2][i]=((x8+16384)>>13);
		data[6][i]=((x7+16384)>>13);
		data[7][i]=((x2-x5+16384)>>13);
		data[1][i]=((x2+x5+16384)>>13);
		data[3][i]=(((x3>>8)*r2+8192)>>12);
		data[5][i]=(((x0>>8)*r2+8192)>>12);
	}

	dctclk += __rdtsc() - a;
}
//
void dct_fill_tab()
{
	unsigned u,x;

	for (u = 0; u < 8; u++)
	{
		printf("%d: ", u);
		for (x = 0; x < 8; x++)
		{
			double Cu = (u==0)? 1.0/sqrt(2.0): 1.0;

			double K = Cu * cos((double)(2*x+1) * (double)u * PI/16.0);
			dct_tbl_i[u][x] = K * IDCTI_AMP;
			//dct_tbl_s[u][x] = K * IDCTI_AMP;
			//idct_tbl_s[x][u] = K * IDCTI_AMP; // different order

			printf("%f(%d),", K, idct_tbl_s[u][x]);
		}
		printf("\n");
	}
} 

/* real DCT
void dct2(conv pixel[8][8], conv data[8][8])
{
	unsigned x, y, n;
	float tmp[8][8];

//	uint64_t a = __rdtsc();

	for (y = 0; y < 8; y++)
	for (x = 0; x < 8; x++)
	{
		float q = 0.0f;

		for (n = 0; n < 8; n++)
			q += pixel[y][n] * dct_tbl[x][n];

		tmp[y][x] = q/2;
	}
		
	for (x = 0; x < 8; x++)
	for (y = 0; y < 8; y++)
	{
		float q = 0.0f;

		for (n = 0; n < 8; n++)
			q += tmp[n][x] * dct_tbl[y][n];

		data[y][x] = q/2;
	}

//	idctclk += __rdtsc() - a;
}*/

// integer DCT 
void dct2_i(conv pixel[8][8], conv data[8][8])
{
	unsigned x, y, n;
	conv tmp[8][8];

	//uint64_t a = __rdtsc();

	// process rows
	for (y = 0; y < 8; y++)
	for (x = 0; x < 8; x++)
	{
		int q = 0;

		for (n = 0; n < 8; n++)
			q += pixel[y][n] * dct_tbl_i[x][n];

		tmp[y][x] = (q + ((q<0)? -IDCTI_AMP: IDCTI_AMP))/(IDCTI_AMP*2);
	}
		
	// process columns
	for (x = 0; x < 8; x++)
	for (y = 0; y < 8; y++)
	{
		int q = 0;

		for (n = 0; n < 8; n++)
			q += tmp[n][x] * dct_tbl_i[y][n];

		data[y][x] = (q + ((q<0)? -IDCTI_AMP: IDCTI_AMP))/(IDCTI_AMP*2);
	}

	//idctclk += __rdtsc() - a;
} 

#ifdef _MSC_VER

void dct2_s(conv pixel[8][8], conv data[8][8])
{
	CACHE_ALIGN conv tmp[8][8];
	unsigned x, y;

	//uint64_t t = __rdtsc();

	// process rows
	for (y = 0; y < 8; y++) {

		__m128i a = _mm_loadu_si128 ((__m128i*)pixel[y]);

		for (x = 0; x < 8; x++) {
			__m128i b, c;

			b = _mm_load_si128 ((__m128i*)dct_tbl_s[x]);
			b = _mm_madd_epi16 (a, b);
			c = _mm_shuffle_epi32 (b, 0xB1);
			c = _mm_add_epi32 (b, c);
			b = _mm_shuffle_epi32 (c, 0x27);
			c = _mm_add_epi32 (b, c);
			tmp[x][y] = _mm_cvtsi128_si32(c)/(IDCTI_AMP*2);
		}
	}

	// process columns
	for (y = 0; y < 8; y++) {

		__m128i a = _mm_loadu_si128 ((__m128i*)tmp[y]);

		for (x = 0; x < 8; x++) {
			__m128i b, c;

			b = _mm_load_si128 ((__m128i*)dct_tbl_s[x]);
			b = _mm_madd_epi16 (a, b);
			c = _mm_shuffle_epi32 (b, 0xB1);
			c = _mm_add_epi32 (b, c);
			b = _mm_shuffle_epi32 (c, 0x27);
			c = _mm_add_epi32 (b, c);

			data[x][y] = _mm_cvtsi128_si32(c)/(IDCTI_AMP*2);
		}
	}

	//idctclk += __rdtsc() - t;
}

#endif//_MSC_VER


__inline static int sdiv(const int data, const int quant, const unsigned mag)
{
	return data >> mag;
	//return (data + quant) >> mag;
	//return (data + ((data<0)? -quant: quant))/(1<<mag);
}

// simple but fast DCT
void dct3(conv pixels[8][8], conv data[8][8])
{
	CACHE_ALIGN int rows[8][8];
	unsigned i;

	static const int
		C1 = 1004,// cos(1*Pi/16) = 0.98078528 * 1024
		C2 = 946, // cos(2*Pi/16) = 0.92387953
		C3 = 852, // cos(3*Pi/16) = 0.83146961
		C4 = 724, // cos(4*Pi/16) = 0.70710678
		C5 = 569, // cos(5*Pi/16) = 0.55557023
		C6 = 392, // cos(6*Pi/16) = 0.38268343
		C7 = 200; // cos(7*Pi/16) = 0.19509032

	uint64_t a = __rdtsc();

	/* transform rows */
	for (i = 0; i < 8; i++)
	{
		int s07,s16,s25,s34,s0734,s1625;
		int d07,d16,d25,d34,d0734,d1625;

		s07 = pixels[i][0] + pixels[i][7];
		d07 = pixels[i][0] - pixels[i][7];
		s16 = pixels[i][1] + pixels[i][6];
		d16 = pixels[i][1] - pixels[i][6];
		s25 = pixels[i][2] + pixels[i][5];
		d25 = pixels[i][2] - pixels[i][5];
		s34 = pixels[i][3] + pixels[i][4];
		d34 = pixels[i][3] - pixels[i][4];

		rows[i][1] = sdiv(C1*d07 + C3*d16 + C5*d25 + C7*d34, 512, 10);
		rows[i][3] = sdiv(C3*d07 - C7*d16 - C1*d25 - C5*d34, 512, 10);
		rows[i][5] = sdiv(C5*d07 - C1*d16 + C7*d25 + C3*d34, 512, 10);
		rows[i][7] = sdiv(C7*d07 - C5*d16 + C3*d25 - C1*d34, 512, 10);

		s0734 = s07 + s34;
		d0734 = s07 - s34;
		s1625 = s16 + s25;
		d1625 = s16 - s25;

		rows[i][0] = sdiv(C4*(s0734 + s1625), 512, 10);
		rows[i][4] = sdiv(C4*(s0734 - s1625), 512, 10);

		rows[i][2] = sdiv(C2*d0734 + C6*d1625, 512, 10);
		rows[i][6] = sdiv(C6*d0734 - C2*d1625, 512, 10);
	}

	/* transform columns */
	for (i = 0; i < 8; i++)
	{
		int s07,s16,s25,s34,s0734,s1625;
		int d07,d16,d25,d34,d0734,d1625;

		s07 = rows[0][i] + rows[7][i];
		d07 = rows[0][i] - rows[7][i];
		s16 = rows[1][i] + rows[6][i];
		d16 = rows[1][i] - rows[6][i];
		s25 = rows[2][i] + rows[5][i];
		d25 = rows[2][i] - rows[5][i];
		s34 = rows[3][i] + rows[4][i];
		d34 = rows[3][i] - rows[4][i];

		data[1][i] = sdiv(C1*d07 + C3*d16 + C5*d25 + C7*d34, 2048, 12);
		data[3][i] = sdiv(C3*d07 - C7*d16 - C1*d25 - C5*d34, 2048, 12);
		data[5][i] = sdiv(C5*d07 - C1*d16 + C7*d25 + C3*d34, 2048, 12);
		data[7][i] = sdiv(C7*d07 - C5*d16 + C3*d25 - C1*d34, 2048, 12);

		s0734 = s07 + s34;
		d0734 = s07 - s34;
		s1625 = s16 + s25;
		d1625 = s16 - s25;

		data[0][i] = sdiv(C4*(s0734 + s1625), 2048, 12);
		data[4][i] = sdiv(C4*(s0734 - s1625), 2048, 12);

		data[2][i] = sdiv(C2*d0734 + C6*d1625, 2048, 12);
		data[6][i] = sdiv(C6*d0734 - C2*d1625, 2048, 12);
	}

	dctclk += __rdtsc() - a;
}

// fast DCT, Vetterli & Ligtenberg
void dct4(conv pixels[8][8], conv data[8][8])
{
	int rows[8][8];
	unsigned i;

	static const int
		C1 = 16069,// cos(1*Pi/16) = 0.9808 * 16384
		C2 = 15137,// cos(2*Pi/16) = 0.9239
		C3 = 13623,// cos(3*Pi/16) = 0.8315
		C4 = 11585,// cos(4*Pi/16) = 0.7071
		C5 = 9102, // cos(5*Pi/16) = 0.5556
		C6 = 6270, // cos(6*Pi/16) = 0.3827
		C7 = 3197; // cos(7*Pi/16) = 0.1951

	uint64_t a = __rdtsc();

	/* transform rows */
	for (i = 0; i < 8; i++)
	{
		int s07,s12,s34,s56;
		int d07,d12,d34,d56;
		int x, y;

		s07 = pixels[i][0] + pixels[i][7];
		d07 = pixels[i][0] - pixels[i][7];

		s12 = pixels[i][1] + pixels[i][2];
		d12 = pixels[i][1] - pixels[i][2];
		
		s34 = pixels[i][3] + pixels[i][4];
		d34 = pixels[i][3] - pixels[i][4];

		s56 = pixels[i][5] + pixels[i][6];
		d56 = pixels[i][5] - pixels[i][6];

		x = s07 + s34;
		y = s12 + s56;
		rows[i][0] = C4*(x + y)/32768;
		rows[i][4] = C4*(x - y)/32768;

		x = d12 - d56;
		y = s07 - s34;
		rows[i][2] = (C6*x + C2*y)/32768;
		rows[i][6] = (C6*y - C2*x)/32768;

		x = d07 - (C4*(s12 - s56) >> 14);
		y = d34 - (C4*(d12 + d56) >> 14);
		rows[i][3] = (C3*x - C5*y)/32768;
		rows[i][5] = (C5*x + C3*y)/32768;

		x = d07*2 - x;
		y = d34*2 - y;
		rows[i][1] = (C1*x + C7*y)/32768;
		rows[i][7] = (C7*x - C1*y)/32768;
	}

	/* transform columns */
	for (i = 0; i < 8; i++)
	{
		int s07,s12,s34,s56;
		int d07,d12,d34,d56;
		int x, y;

		s07 = rows[0][i] + rows[7][i];
		d07 = rows[0][i] - rows[7][i];

		s12 = rows[1][i] + rows[2][i];
		d12 = rows[1][i] - rows[2][i];

		s34 = rows[3][i] + rows[4][i];
		d34 = rows[3][i] - rows[4][i];

		s56 = rows[5][i] + rows[6][i];
		d56 = rows[5][i] - rows[6][i];

		x = s07 + s34;
		y = s12 + s56;
		data[0][i] = C4*(x + y)/32768;
		data[4][i] = C4*(x - y)/32768;

		x = d12 - d56;
		y = s07 - s34;
		data[2][i] = (C6*x + C2*y)/32768;
		data[6][i] = (C6*y - C2*x)/32768;

		x = d07 - (C4*(s12 - s56) >> 14);
		y = d34 - (C4*(d12 + d56) >> 14);
		data[3][i] = (C3*x - C5*y)/32768;
		data[5][i] = (C5*x + C3*y)/32768;

		x = d07*2 - x;
		y = d34*2 - y;
		data[1][i] = (C1*x + C7*y)/32768;
		data[7][i] = (C7*x - C1*y)/32768;
	}

	dctclk += __rdtsc() - a;
}

void dct5(conv pixels[8][8], conv data[8][8])
{
	short rows[8][8];
	unsigned i;

	static const int
				c1 = 1004,  /* cos(pi/16) << 10 */
				s1 = 200,   /* sin(pi/16) */
				c3 = 851,   /* cos(3pi/16) << 10 */
				s3 = 569,   /* sin(3pi/16) << 10 */
				r2c6 = 554, /* sqrt(2)*cos(6pi/16) << 10 */
				r2s6 = 1337,/* sqrt(2)*sin(6pi/16) << 10 */
				r2 = 181;   /* sqrt(2) << 7*/

	uint64_t a = __rdtsc();

	/* transform rows */
	for (i = 0; i < 8; i++)
	{
		int x0,x1,x2,x3,x4,x5,x6,x7,x8;

		x0 = pixels[i][0];
		x1 = pixels[i][1];
		x2 = pixels[i][2];
		x3 = pixels[i][3];
		x4 = pixels[i][4];
		x5 = pixels[i][5];
		x6 = pixels[i][6];
		x7 = pixels[i][7];

		/* Stage 1 */
		x8=x7+x0;		// s07 = x0 + x7
		x0-=x7;			// d07 = x0 - x7
		x7=x1+x6;		// s16 = x1 + x6
		x1-=x6;			// d16 = x1 - x6
		x6=x2+x5;		// s25 = x2 + x5
		x2-=x5;			// d25 = x2 - x5
		x5=x3+x4;		// s34 = x3 + x4
		x3-=x4;			// d34 = x3 - x4

		/* Stage 2 */
		x4=x8+x5;		// s0734 = s07 + s34
		x8-=x5;			// d0734 = s07 - s34
		x5=x7+x6;		// s1625 = s16 + s25
		x7-=x6;			// d1625 = s16 - s25
		x6=c1*(x1+x2);		// C1*(d16 + d25)
		x2=(-s1-c1)*x2+x6;	// (-S1 - C1)*d25 + C1*(d16 + d25) = C1*d16 - S1*d25
		x1=(s1-c1)*x1+x6;	//  (S1 - C1)*d16 + C1*(d16 + d25) = C1*d25 + S1*d16
		x6=c3*(x0+x3);		// C3*(d07 + d34)
		x3=(-s3-c3)*x3+x6;	// (-S3 - C3)*d34 + C3*(d07 + d34) = C3*d07 - S3*d34
		x0=(s3-c3)*x0+x6;	//  (S3 - C3)*d07 + C3*(d07 + d34) = C3*d34 + S3*d07

		/* Stage 3 */
		x6=x4+x5;				// s0734 + s1625
		x4-=x5;					// s0734 - s1625
		x5=r2c6*(x7+x8);		// 1024*sqrt(2)*cos(6pi/16)*(d1625 + d0734)
		x7=(-r2s6-r2c6)*x7+x5;	// -1024*sqrt(2)*(sin(6pi/16)+cos(6pi/16))*d1625 + 1024*sqrt(2)*cos(6pi/16)*(d1625 + d0734)
		x8=(r2s6-r2c6)*x8+x5;	// 1024*sqrt(2)*(sin(6pi/16)-cos(6pi/16))*d0734 + 1024*sqrt(2)*cos(6pi/16)*(d1625 + d0734)
		x5=x0+x2;	// (S3 - C3)*d34 + C3*(d07 + d34) + (-S1 - C1)*d25 + C1*(d16 + d25) = C3*d34 + S3*d07 + C1*d16 - S1*d25
		x0-=x2;		// (S3 - C3)*d34 + C3*(d07 + d34) - (-S1 - C1)*d25 - C1*(d16 + d25) = C3*d34 + S3*d07 - C1*d16 + S1*d25
		x2=x3+x1;	// (-S3 - C3)*d34 + C3*(d07 + d34) + (S1 - C1)*d16 + C1*(d16 + d25) = C3*d07 - S3*d34 + C1*d25 + S1*d16
		x3-=x1;		// (-S3 - C3)*d34 + C3*(d07 + d34) - (S1 - C1)*d16 - C1*(d16 + d25) = C3*d07 - S3*d34 - C1*d25 - S1*d16

		/* Stage 4 and output */
		rows[i][0]=x6;
		rows[i][4]=x4;
		rows[i][2]=x8>>10;
		rows[i][6]=x7>>10;
		rows[i][7]=(x2-x5)>>10; // C3*d07 - S3*d34 + C1*d25 + S1*d16 - C3*d34 - S3*d07 - C1*d16 + S1*d25
		rows[i][1]=(x2+x5)>>10;
		rows[i][3]=(x3*r2)>>17;
		rows[i][5]=(x0*r2)>>17;
	}

	/* transform columns */
	for (i = 0; i < 8; i++)
	{
		int x0,x1,x2,x3,x4,x5,x6,x7,x8;

		x0 = rows[0][i];
		x1 = rows[1][i];
		x2 = rows[2][i];
		x3 = rows[3][i];
		x4 = rows[4][i];
		x5 = rows[5][i];
		x6 = rows[6][i];
		x7 = rows[7][i];

		/* Stage 1 */
		x8=x7+x0;
		x0-=x7;
		x7=x1+x6;
		x1-=x6;
		x6=x2+x5;
		x2-=x5;
		x5=x3+x4;
		x3-=x4;

		/* Stage 2 */
		x4=x8+x5;
		x8-=x5;
		x5=x7+x6;
		x7-=x6;
		x6=c1*(x1+x2);
		x2=(-s1-c1)*x2+x6;
		x1=(s1-c1)*x1+x6;
		x6=c3*(x0+x3);
		x3=(-s3-c3)*x3+x6;
		x0=(s3-c3)*x0+x6;

		/* Stage 3 */
		x6=x4+x5;
		x4-=x5;
		x5=r2c6*(x7+x8);
		x7=(-r2s6-r2c6)*x7+x5;
		x8=(r2s6-r2c6)*x8+x5;
		x5=x0+x2;
		x0-=x2;
		x2=x3+x1;
		x3-=x1;

		/* Stage 4 and output */
		data[0][i]=((x6+16)>>3);
		data[4][i]=((x4+16)>>3);
		data[2][i]=((x8+16384)>>13);
		data[6][i]=((x7+16384)>>13);
		data[7][i]=((x2-x5+16384)>>13);
		data[1][i]=((x2+x5+16384)>>13);
		data[3][i]=(((x3>>8)*r2+8192)>>12);
		data[5][i]=(((x0>>8)*r2+8192)>>12);
	}

	dctclk += __rdtsc() - a;
}

/* inverse real DCT
void idct2(conv data[8][8], conv pixel[8][8])
{
	unsigned x, y, n;
	float tmp[8][8];

	uint64_t a = __rdtsc();

	for (y = 0; y < 8; y++)
	{
		for (x = 0; x < 8; x++)
		{
			float q = 0.0f;

			for (n = 0; n < 8; n++)
				q += data[y][n] * dct_tbl[n][x];

			tmp[y][x] = q/2;
		}
	}
		
	for (x = 0; x < 8; x++)
	{
		for (y = 0; y < 8; y++)
		{
			float q = 0.0f;

			for (n = 0; n < 8; n++)
				q += tmp[n][x] * dct_tbl[n][y];

			pixel[y][x] = q/2;
		}
	}

	idctclk += __rdtsc() - a;
}*/


/* inverse integer DCT 
void idct2_i(conv data[8][8], conv pixel[8][8])
{
	unsigned x, y, n;
	conv tmp[8][8];

	uint64_t a = __rdtsc();

	// process rows
	for (y = 0; y < 8; y++)
	for (x = 0; x < 8; x++)
	{
		int q = 0;

		for (n = 0; n < 8; n++)
			q += data[y][n] * dct_tbl_i[n][x];

		tmp[y][x] = q/(IDCTI_AMP*2);
	}
		
	// process columns
	for (x = 0; x < 8; x++)
	for (y = 0; y < 8; y++)
	{
		int q = 0;

		for (n = 0; n < 8; n++)
			q += tmp[n][x] * dct_tbl_i[n][y];

		pixel[y][x] = q/(IDCTI_AMP*2);
	}

	idctclk += __rdtsc() - a;
} 
*/

#ifdef _MSC_VER

/*
void test_idct2_s()
{
	int data[4] = {1,2,3,4};
	int res;

	// process rows
	__m128i b = _mm_loadu_si128 ((__m128i*)data);

	__m128i c = _mm_shuffle_epi32 (b, 0xB1);
	c = _mm_add_epi32 (b, c);

	b = _mm_shuffle_epi32 (c, 0x27);
	c = _mm_add_epi32 (b, c);

	res = _mm_cvtsi128_si32 (c);
}
*/

void idct2_s(conv data[8][8], conv pixel[8][8])
{
	CACHE_ALIGN conv tmp[8][8];
	unsigned x, y;

	uint64_t t = __rdtsc();

	// process rows
	for (y = 0; y < 8; y++) {

		__m128i r0 = _mm_loadu_si128 ((__m128i*)data[y]);

		for (x = 0; x < 8; x++)
		{
			__m128i r1, r2, r3, r4;

			r1 = _mm_load_si128 ((__m128i*)idct_tbl_s[x]);
			r2 = _mm_madd_epi16 (r0, r1);
			r3 = _mm_shuffle_epi32 (r2, 0xB1);
			r4 = _mm_add_epi32 (r2, r3);
			r1 = _mm_shuffle_epi32 (r4, 0x27);
			r2 = _mm_add_epi32 (r1, r4);
			tmp[x][y] = _mm_cvtsi128_si32(r2)/(IDCTI_AMP*2);
		}
	}

	// process columns
	for (y = 0; y < 8; y++) {

		__m128i r0 = _mm_loadu_si128 ((__m128i*)tmp[y]);

		for (x = 0; x < 8; x += 2)
		{
			__m128i r1, r2, r3, r4, r5, r6;

			r1 = _mm_load_si128 ((__m128i*)idct_tbl_s[x]);
			r2 = _mm_madd_epi16 (r0, r1);
			r4 = _mm_load_si128 ((__m128i*)idct_tbl_s[x+1]);
			r1 = _mm_shuffle_epi32 (r2, 0xB1);
			r5 = _mm_madd_epi16 (r0, r4);
			r3 = _mm_add_epi32 (r1, r2);
			r4 = _mm_shuffle_epi32 (r5, 0xB1);
			r2 = _mm_shuffle_epi32 (r3, 0x27);
			r6 = _mm_add_epi32 (r4, r5);
			r1 = _mm_add_epi32 (r2, r3);
			r5 = _mm_shuffle_epi32 (r6, 0x27);
			pixel[x][y] = _mm_cvtsi128_si32(r1)/(IDCTI_AMP*2);
			r4 = _mm_add_epi32 (r5, r6);
			pixel[x+1][y] = _mm_cvtsi128_si32(r4)/(IDCTI_AMP*2);
		}
	}
	idctclk += __rdtsc() - t;
}

#endif//_MSC_VER
