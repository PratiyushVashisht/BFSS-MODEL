
#include "unitary.cpp"
extern "C" void zgeev_( char*, char*, int*, double at[], int *,double b[], double dummy[],
                       int *, double dummy2[], int*, double work[], int *, double work2[], int *);
                       
Complex  eigen (Umatrix &x,Complex Eig[NCOLOR])
{	
	
   double aa[2*NCOLOR*NCOLOR],cc[4*NCOLOR*NCOLOR],dd[2],E[2*NCOLOR];
	char nn='N';
	int c1=NCOLOR,c2=2*NCOLOR,c3=1,k;
	
	for(int i=0;i<NCOLOR;i++)
		for(int j=0;j<NCOLOR;j++)
		{
			aa[2*(NCOLOR*i+j)]=x.get(i,j).real();
			aa[2*(NCOLOR*i+j)+1]=x.get(i,j).imag();
                      //  cout<<aa[2*(M*i+j)]<<endl;
		}
	
	zgeev_(&nn,&nn,&c1,aa,&c1,E,dd,&c3,dd,&c3,cc,&c2,cc,&k);   


 for (int j = 0; j < NCOLOR; j++) {
    Eig[j] = Complex (E[2 * j],E[2 * j + 1]);
}

}
