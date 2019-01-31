#include <stdlib.h>
#include <stdio.h>

void polymul_ref(int *x, int *y,int *z)
{
	for(int i=0;i<31;i++)
		z[i]=0;
	for(int i=0;i<16;i++)
		for(int j=0;j<16;j++)
			z[i+j]+=x[i]*y[j];
	for(int i=0;i<31;i++)
	{
		z[i]=z[i]%4591;
	}
}

int main()
{
	int x[16];
	int y[16];
	int z[31];
	FILE *in = fopen("random_input.txt", "r");
	FILE *out = fopen("polymul_ref.txt", "w");
	for(int j=0;j<10000;j++)
	{
		for(int i=0;i<16;i++)
			fscanf(in,"%d ",&x[i]);
		for(int i=0;i<16;i++)
			fscanf(in,"%d ",&y[i]);
		polymul_ref(x,y,z);
		for(int i=0;i<31;i++)
			fprintf(out, "%d,", z[i]);
		fprintf(out, "\n" );
	}
	return 0;
}