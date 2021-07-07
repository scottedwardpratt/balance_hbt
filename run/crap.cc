#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	if (argc != 3) {
		printf("Usage: crap a b\n");
		exit(-1);
	}
	double a,b,c;
	a=atof(argv[1]);
	b=atof(argv[2]);
	c=a+b;
	printf("%g\n",c);
	return 0;
}


