import random
f=open('random_input.txt','w');
for i in range(0,16*10000):
	f.write(str(random.randint(0,4591)));
	f.write(" ");
	if i%16==15:
		f.write("\n");
