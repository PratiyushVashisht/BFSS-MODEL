
#include "matrix.cpp"


void srand48(long seedval) {
    srand(seedval);
}
long lrand48() {
    return rand();
}
double drand48() {
    return rand() / (double)RAND_MAX;
}


extern "C" void zgesvd_(char *A1, char *A2, int *N1, int *N2, double *store,int *lda, double *junk, double *left, int *Nl,double *right, int *Nr, double *work, int *Nwork,        double *Rwork, int *stat);

Umatrix unitary(Umatrix x_new)
{
  //From SUSY LATTICE code
  register int i;
  char A = 'A';     // Ask LAPACK for all singular values
  int row, col, Npt = NCOLOR, stat = 0, Nwork = 3 * NCOLOR, err = 0;
  double store[2*NCOLOR*NCOLOR], left[2*NCOLOR*NCOLOR], right[2*NCOLOR*NCOLOR], reunit_work[6*NCOLOR], reunit_Rwork[5*NCOLOR], junk[NCOLOR];

  Umatrix  lmat ,rdagmat  ,prod ;
  Umatrix C;
	  for (row = 0; row < NCOLOR; row++)
	    {
	      for (col = 0; col < NCOLOR; col++)
		{
                 i = 2 * (col * NCOLOR+ row);
                 store[i] = (x_new.get(row,col)).real();
                 store[i + 1] = (x_new.get(row,col)).imag();
                 //cout << "store i, i+1 " << store[i] << " " << store[i+1] << endl;
		}
    }

zgesvd_(&A, &A, &Npt, &Npt, store, &Npt, junk, left, &Npt, right, &Npt,          reunit_work, &Nwork, reunit_Rwork, &stat);    
  
 for (row = 0; row < NCOLOR; row++) {
    for (col = 0; col < NCOLOR; col++) {
      i = 2 * (col * NCOLOR + row);
      lmat.set(row,col, Complex (left[i], left[i + 1]));
      rdagmat.set(row,col ,Complex (right[i], right[i + 1]));
    }
  }  
    
  prod = (lmat * rdagmat);
for (row = 0; row < NCOLOR; row++){
	      for (col = 0; col < NCOLOR; col++){
      x_new.set(row,col, Complex( (prod.get(row,col).real()) , (prod.get(row,col).imag()) ));
}
}
   return(x_new) ;
}


