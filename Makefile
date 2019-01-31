all:
	gcc -mavx2 -O3 polymul_test.c -o polymul_test
	gcc -mavx2 -O3 polymul_avx2.c -o polymul_avx2
	gcc polymul_ref.c -o polymul_ref