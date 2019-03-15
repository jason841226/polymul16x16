#ifndef POLYMUL16X16
#define POLYMUL16X16

#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define v4591_16 _mm256_set1_epi16(4591)
#define v15631_16 _mm256_set1_epi16(15631)
#define v4158_16 _mm256_set1_epi16(4158)
#define v1_16 _mm256_set1_epi16(1)

void print256_num(__m256i var);
static inline void print_cmp(__m256i _x,__m256i _y);
static inline __m256i montproduct(__m256i x,__m256i y);
static inline __m256i polymullo(__m256i x,__m256i _y);
static inline __m256i polymullo2(__m256i x,__m256i _y);
static inline __m256i polymulhi(__m256i x,__m256i _y);
static inline __m256i polymulhi2(__m256i x,__m256i _y);
static void montmut_16x16(__m256i * zl, __m256i *zh , __m256i x , __m256i y , 
		__m256i mask_lo , __m256i mask_hi );



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

//TODO: 256bit shift

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
//aligr
//TODO:combine polymullo and polymulhi
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
#endif