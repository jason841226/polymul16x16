#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


//TO TRY: 256bit shift(alignr and srli)


#define v4591_16 _mm256_set1_epi16(4591)
#define v15631_16 _mm256_set1_epi16(15631)
#define v4158_16 _mm256_set1_epi16(4158)
#define v1_16 _mm256_set1_epi16(1)

void print256_num(__m256i var)
{
    uint16_t *val = (uint16_t*) &var;
    for(int i=0;i<16;i++)
    	if(val[i]>60945)
    		val[i]-=60945;
    printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n", 
           val[15], val[14], val[13], val[12], val[11], val[10], val[9], val[8],
           val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

static inline void print_cmp(__m256i _x,__m256i _y)
{
	uint16_t *x=(uint16_t*) &_x;
	uint16_t *y=(uint16_t*) &_y;
	for(int i=15;i>=0;i--)
	{
		if(x[i]==y[i]||(x[i]-60945)==y[i])
			printf("True ");
		else
			printf("False!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ");
	}
	printf("\n");

}

static inline __m256i montproduct(__m256i x,__m256i y)
{
  __m256i lo,hi,d,e;

  lo = _mm256_mullo_epi16(x,y);
  hi = _mm256_mulhi_epi16(x,y);
  d = _mm256_mullo_epi16(lo,v15631_16);
  e = _mm256_mulhi_epi16(d,v4591_16);
  return _mm256_sub_epi16(hi,e);
}

static inline __m256i polymullo(__m256i x,__m256i y)
{
	__m256i tmp,z;
	__m256i br_reg = _mm256_set1_epi16(0x0100);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	z=montproduct(x,tmp);
	for(int i=1;i<8;++i)
	{
		x=_mm256_slli_si256(x,2);
		y = _mm256_srli_si256(y,2);
		tmp = _mm256_shuffle_epi8( y , br_reg );
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;
}
static inline __m256i polymullo2(__m256i x,__m256i y)//rotate y (other same as polymullo)
{
	__m256i tmp,z;
	__m256i br_reg = _mm256_set1_epi16(0x0100);
	y=_mm256_permute2x128_si256(y,y,0x01);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	z=montproduct(x,tmp);

	for(int i=1;i<8;++i)
	{
		x=_mm256_slli_si256(x,2);
		y = _mm256_srli_si256(y,2);
		tmp = _mm256_shuffle_epi8( y , br_reg );
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;	
}

static inline __m256i polymulhi(__m256i x,__m256i y)
{
	__m256i tmp,z;
	__m256i br_reg = _mm256_set1_epi16(0x0F0E);
	x=_mm256_srli_si256(x,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	z=montproduct(x,tmp);
	for(int i=1;i<7;++i)
	{
		x=_mm256_srli_si256(x,2);
		y=_mm256_slli_si256(y,2);
		tmp=_mm256_shuffle_epi8(y,br_reg);
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;
}
static inline __m256i polymulhi2(__m256i x,__m256i y)
{
	__m256i tmp,z;
	__m256i br_reg = _mm256_set1_epi16(0x0F0E);
	x=_mm256_srli_si256(x,2);
	y=_mm256_permute2x128_si256(y,y,0x01);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	z=montproduct(x,tmp);
	
	for(int i=1;i<7;++i)
	{
		x=_mm256_srli_si256(x,2);
		y=_mm256_slli_si256(y,2);
		tmp=_mm256_shuffle_epi8(y,br_reg);
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;
}


//input(x,y) : polynomial with 16 montgomery domain coef (each coef store in 16bit)
//output (z) : polynomial with 32 montgomery domain coef (first coef=0)
//z=x*y
static void montmut_16x16_readable(__m256i * zl, __m256i *zh , __m256i x , __m256i y , 
		__m256i mask_lo , __m256i mask_hi )
{
	__m256i z0_lo,z0_hi,z1_lo,z1_hi,z_lo,z_hi;

	z0_lo=polymullo(x,y);
	z0_hi=polymulhi(x,y);
	z1_lo=polymullo2(x,y);
	z1_hi=polymulhi2(x,y);

	z_lo=_mm256_permute2x128_si256(z0_lo,z0_hi,0x20);
	z_hi=_mm256_permute2x128_si256(z0_lo,z0_hi,0x31);

	z0_lo=_mm256_permute2x128_si256(z1_lo,z1_hi,0x02);
	z0_hi=_mm256_permute2x128_si256(z1_lo,z1_hi,0x13);

	z0_lo=_mm256_add_epi16(z0_lo,z0_hi);

	z0_hi=_mm256_and_si256(z0_lo,mask_lo);
	z_hi=_mm256_add_epi16(z_hi,z0_hi);

	z0_lo=_mm256_and_si256(z0_lo,mask_hi);
	z_lo=_mm256_add_epi16(z_lo,z0_lo);


_mm256_store_si256(zl , z_lo);
_mm256_store_si256(zh , z_hi);
}

static void montmut_16x16(__m256i * zl, __m256i *zh , __m256i x , __m256i y , 
		__m256i mask_lo , __m256i mask_hi )
{
	__m256i z0_lo,z0_hi,z1_lo,z1_hi,z_lo,z_hi;
	__m256i _y=_mm256_permute2x128_si256(y,y,0x01);
	__m256i tmp,tmp2;
	__m256i br_reg = _mm256_set1_epi16(0x0100);

	//i=0
	tmp = _mm256_shuffle_epi8( y , br_reg );
	z0_lo=montproduct(x,tmp);
	
	//i=1
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,2);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	z0_hi=_mm256_srli_si256(tmp,14);

	//i=2
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,4);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,12);
	z0_hi=_mm256_add_epi16(z0_hi,tmp2);

	//i=3
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,6);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,10);
	z0_hi=_mm256_add_epi16(z0_hi,tmp2);

	//i=4
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,8);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,8);
	z0_hi=_mm256_add_epi16(z0_hi,tmp2);

	//i=5
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,10);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,6);
	z0_hi=_mm256_add_epi16(z0_hi,tmp2);

	//i=6
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,12);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,4);
	z0_hi=_mm256_add_epi16(z0_hi,tmp2);

	//i=7
	y = _mm256_srli_si256(y,2);
	tmp = _mm256_shuffle_epi8( y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,14);
	z0_lo=_mm256_add_epi16(z0_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,2);
	z0_hi=_mm256_add_epi16(z0_hi,tmp2);

	//i=0
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	z1_lo=montproduct(x,tmp);
	
	//i=1
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,2);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	z1_hi=_mm256_srli_si256(tmp,14);


	//i=2
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,4);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,12);
	z1_hi=_mm256_add_epi16(z1_hi,tmp2);

	//i=3
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,6);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,10);
	z1_hi=_mm256_add_epi16(z1_hi,tmp2);

	//i=4
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,8);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,8);
	z1_hi=_mm256_add_epi16(z1_hi,tmp2);

	//i=5
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,10);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,6);
	z1_hi=_mm256_add_epi16(z1_hi,tmp2);

	//i=6
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,12);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,4);
	z1_hi=_mm256_add_epi16(z1_hi,tmp2);

	//i=7
	_y = _mm256_srli_si256(_y,2);
	tmp = _mm256_shuffle_epi8( _y , br_reg );
	tmp=montproduct(x,tmp);
	tmp2=_mm256_slli_si256(tmp,14);
	z1_lo=_mm256_add_epi16(z1_lo,tmp2);
	tmp2=_mm256_srli_si256(tmp,2);
	z1_hi=_mm256_add_epi16(z1_hi,tmp2);

	//combine
	z_lo=_mm256_permute2x128_si256(z0_lo,z0_hi,0x20);
	z_hi=_mm256_permute2x128_si256(z0_lo,z0_hi,0x31);

	tmp=_mm256_permute2x128_si256(z1_lo,z1_hi,0x02);
	tmp2=_mm256_permute2x128_si256(z1_lo,z1_hi,0x13);

	tmp=_mm256_add_epi16(tmp,tmp2);

	z0_hi=_mm256_and_si256(tmp,mask_lo);
	z_hi=_mm256_add_epi16(z_hi,z0_hi);

	z0_lo=_mm256_and_si256(tmp,mask_hi);
	z_lo=_mm256_add_epi16(z_lo,z0_lo);

_mm256_store_si256(zl , z_lo);
_mm256_store_si256(zh , z_hi);
}


#include "benchmark.h"

int main()

{
	//TODO:reduce register
	__m256i x,y,z_lo,z_hi;
	__m256i mask_lo=_mm256_set_epi16(0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1);
	__m256i mask_hi=_mm256_set_epi16(-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0);

	x=_mm256_set_epi16(1650,2855,279,878,477,2971,1138,4437,2730,1293,3360,3586,4330,219,2514,3189);
	y=_mm256_set_epi16(3700,3714,378,2596,1606,3397,1877,2113,2571,1332,1284,3924,3663,1608,3361,1427);

	x=montproduct(x,v4158_16);
	y=montproduct(y,v4158_16);

	// __m256i x2,y2,z_lo2,z_hi2;
	// x2=_mm256_set_epi16(160,255,79,87,47,291,138,447,230,123,330,386,430,29,254,319);
	// y2=_mm256_set_epi16(370,371,37,296,106,397,177,213,251,132,124,324,363,168,331,147);
	struct benchmark bm;
	char msg[256];
	bm_init( & bm);
	for(int i=0;i<10000;i++){
BENCHMARK( bm , {
	montmut_16x16( &z_lo , &z_hi , x , y , mask_lo , mask_hi);
});
	}
	bm_dump( msg , 256 , &bm );
	printf("%s\n", msg );

	z_hi=montproduct(z_hi,v1_16);
	print256_num(z_hi);
	z_lo=montproduct(z_lo,v1_16);
	print256_num(z_lo);

	// print256_num(z_hi2);
	// print256_num(z_lo2);
	
	__m256i z_lo_ans=_mm256_set_epi16(3501,311,4232,3640,830,816,1626,2630,2554,4512,2724,3550,567,2204,151,1022);
	__m256i z_hi_ans=_mm256_set_epi16(0,3561,3315,1500,1705,1115,4287,116,337,3953,2431,782,548,3224,3051,2230);

	print_cmp(z_hi,z_hi_ans);
	print_cmp(z_lo,z_lo_ans);
	return 0;
}