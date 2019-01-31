#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>



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
static inline __m256i polymullo(__m256i x,__m256i _y)
{
	__m256i tmp,z;
	//TODO:change y (shuffle)
	uint16_t *y = (uint16_t*) &_y;
	tmp=_mm256_set_epi16(y[8],y[8],y[8],y[8],y[8],y[8],y[8],y[8],y[0],y[0],y[0],y[0],y[0],y[0],y[0],y[0]);
	//__m256i br_reg = _mm256_set1_epi16(0x0100);
	//__m256i y = _y;
	//tmp = _mm256_shuffle_epi8( y , br_reg );
	z=montproduct(x,tmp);
	
	//TODO:change j
	int i,j=9;
	for(i=1;i<8;++i)
	{
		x=_mm256_slli_si256(x,2);
		//y = _mm256_srli_si256(y,2);
		//tmp = _mm256_shuffle_epi8( y , br_reg );

		tmp=_mm256_set_epi16(y[j],y[j],y[j],y[j],y[j],y[j],y[j],y[j],y[i],y[i],y[i],y[i],y[i],y[i],y[i],y[i]);
		++j;
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;
}
static inline __m256i polymullo2(__m256i x,__m256i _y)//rotate y (other same as polymullo)
{
	__m256i tmp,z;
	uint16_t *y = (uint16_t*) &_y;
	tmp=_mm256_set_epi16(y[0],y[0],y[0],y[0],y[0],y[0],y[0],y[0],y[8],y[8],y[8],y[8],y[8],y[8],y[8],y[8]);
	z=montproduct(x,tmp);
	
	int i,j=9;
	for(i=1;i<8;++i)
	{
		x=_mm256_slli_si256(x,2);
		tmp=_mm256_set_epi16(y[i],y[i],y[i],y[i],y[i],y[i],y[i],y[i],y[j],y[j],y[j],y[j],y[j],y[j],y[j],y[j]);
		++j;
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;	
}

static inline __m256i polymulhi(__m256i x,__m256i _y)
{
	__m256i tmp,z;
	uint16_t *y = (uint16_t*) &_y;
	x=_mm256_srli_si256(x,2);
	tmp=_mm256_set_epi16(y[15],y[15],y[15],y[15],y[15],y[15],y[15],y[15],y[7],y[7],y[7],y[7],y[7],y[7],y[7],y[7]);
	z=montproduct(x,tmp);
	
	int i,j=14;
	for(i=6;i>0;--i)
	{
		x=_mm256_srli_si256(x,2);
		tmp=_mm256_set_epi16(y[j],y[j],y[j],y[j],y[j],y[j],y[j],y[j],y[i],y[i],y[i],y[i],y[i],y[i],y[i],y[i]);
		--j;
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;
}
static inline __m256i polymulhi2(__m256i x,__m256i _y)
{
	__m256i tmp,z;
	uint16_t *y = (uint16_t*) &_y;

	x=_mm256_srli_si256(x,2);
	tmp=_mm256_set_epi16(y[7],y[7],y[7],y[7],y[7],y[7],y[7],y[7],y[15],y[15],y[15],y[15],y[15],y[15],y[15],y[15]);
	z=montproduct(x,tmp);
	
	int i,j=14;
	for(i=6;i>0;--i)
	{
		x=_mm256_srli_si256(x,2);
		tmp=_mm256_set_epi16(y[i],y[i],y[i],y[i],y[i],y[i],y[i],y[i],y[j],y[j],y[j],y[j],y[j],y[j],y[j],y[j]);
		--j;
		tmp=montproduct(x,tmp);
		z=_mm256_add_epi16(z,tmp);
	}
	return z;
}


static void montmut_16x16(__m256i * zl, __m256i *zh , __m256i x , __m256i y , 
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


#include "benchmark.h"

int main()

{
	//TODO:reduce register
	__m256i x,y,z_lo,z_hi;
	__m256i z0_lo,z0_hi,z1_lo,z1_hi;
	__m256i mask_lo=_mm256_set_epi16(0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1);
	__m256i mask_hi=_mm256_set_epi16(-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0);


	/****/
	int _x[16];
	int _y[16];
	int _z[31]={0};
	FILE *in = fopen("random_input.txt", "r");
	FILE *out = fopen("polymul_avx2.txt", "w");

	for(int j=0;j<10000;j++)
	{
		for(int i=0;i<16;i++)
			fscanf(in,"%d ",&_x[i]);
		for(int i=0;i<16;i++)
			fscanf(in,"%d ",&_y[i]);
		x=_mm256_set_epi16(_x[15],_x[14],_x[13],_x[12],_x[11],_x[10],_x[9],_x[8],_x[7],_x[6],_x[5],_x[4],_x[3],_x[2],_x[1],_x[0]);
		y=_mm256_set_epi16(_y[15],_y[14],_y[13],_y[12],_y[11],_y[10],_y[9],_y[8],_y[7],_y[6],_y[5],_y[4],_y[3],_y[2],_y[1],_y[0]);
		x=montproduct(x,v4158_16);
		y=montproduct(y,v4158_16);
		montmut_16x16( &z_lo , &z_hi , x , y , mask_lo , mask_hi );


		z_hi=montproduct(z_hi,v1_16);
		z_lo=montproduct(z_lo,v1_16);

		uint16_t *lol = (uint16_t*) &z_lo;
    	for(int i=0;i<16;i++)
    	{
    		if(lol[i]>60945)
    			lol[i]-=60945;//-2295~2295 -> 0~4590
    		fprintf(out, "%d," , lol[i]);
    	}
    	uint16_t *hih = (uint16_t*) &z_hi;
    	for(int i=0;i<15;i++)
    	{
    		if(hih[i]>60945)
    			hih[i]-=60945;//-2295~2295 -> 0~4590
    		fprintf(out, "%d," , hih[i]);
    	}
    	fprintf(out, "\n");
	}



	return 0;
}

