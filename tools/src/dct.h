#ifndef __DCT_H__
#define __DCT_H__

#ifdef __cplusplus
extern "C" {
#endif

void dct_fill_tab();

// integer DCT 
void dct(conv pixel[8][8], conv data[8][8]);
void dct2(conv pixel[8][8], conv data[8][8]);
void dct3(conv pixel[8][8], conv data[8][8]);
void dct4(conv pixel[8][8], conv data[8][8]);
void dct5(conv pixel[8][8], conv data[8][8]);
void dct2_i(conv pixel[8][8], conv data[8][8]);
void dct2_s(conv pixel[8][8], conv data[8][8]);

// inverse real DCT
void idct(conv data[8][8], conv pixel[8][8]);
void idct2(conv data[8][8], conv pixel[8][8]);
// inverse integer DCT 
void idct2_i(conv data[8][8], conv pixel[8][8]);
void idct2_s(conv data[8][8], conv pixel[8][8]);

void test_idct2_s();

#ifdef __cplusplus
}
#endif

#endif//__DCT_H__
