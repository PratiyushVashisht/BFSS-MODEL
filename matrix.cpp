// Defining X matrix :
#include"matrix.h"

Complex::Complex(void) //This operator gives you more freedom in naming your variables 
                      //by letting you distinguish between variables with the same name
{
    re=0.0; im = 0.0;
}

Complex::Complex(double x,double y ){
    re=x;
    im = y;
}

double Complex::real(void) const{
    return(re);
}

double Complex::imag(void) const {
  return(im);
}

double Complex::norm(void){
    return(sqrt(re*re + im*im));
}

void Complex::print(void) const{
    cout<<"("<< re <<","<< im << ")";
}

ostream& operator<<(ostream& out ,Complex c){ // defining output operator
    out <<c.real()<<"\t"<<c.imag()<<endl;
    return out;
}

istream& operator>>(istream& in, Complex & c) {
  double x, y;
  in >> x >> y;
  c = Complex(x, y);
  return in;
}


// **2**Unitary Matrix U******* Traceless,Hermitian,Unitary

////************************************************************************
Matrix::Matrix(void){
    for(int s= 0; s<SITES; s++)
        for(int d=0; d<D ; d++)
            for(int row =0 ; row < NCOLOR; row++)
                for(int col = 0; col < NCOLOR; col++)
                    X[s][d][row][col]= Complex();
                }
                
Matrix::Matrix(int SITES ,int D,int NCOLOR ){
    for(int s= 0; s<SITES; s++){
        for(int d=0; d<D ; d++){
            for(int row =0 ; row < NCOLOR; row++)
                for(int col = 0; col < NCOLOR; col++){
                    X[s][d][row][col]= Complex();
                }
        }
    }
}

void Matrix::print(void) const{
    for(int s= 0; s<SITES; s++){
        for(int d=0; d<D ; d++){
            for(int row =0 ; row < NCOLOR; row++){
                for(int col = 0; col <NCOLOR; col++){
                    X[s][d][row][col].print()  ;
                }
                cout<<"\t";
            }
            cout<<"\t";
        }
        cout<<"\t";
    }
    cout<<"\t";
}


void Matrix::set(int s, int d, int row ,int col, const Complex o){
    X[s][d][row][col] = o;
}

Complex Matrix::get(int s,int d,int row ,int col)  const { // element of particular row,column
    return(X[s][d][row][col]);
}

/*Complex Matrix::loc(const int s) const {
    for (int d=0; d<D; d++)
        for(int row =0 ; row < NCOLOR; row++)
            for(int col = 0; col < NCOLOR; col++){
                return(X[s][d][row][col]);
                }
}
Complex Matrix::scalar(const int s ,const int d ) const {
    for(int row =0 ; row < NCOLOR; row++)
            for(int col = 0; col < NCOLOR; col++){
                return(X[s][d][row][col]);
                }
}

ostream& operator << (ostream& out , Matrix X){
   for(int s= 0; s<SITES; s++)
        for(int d=0; d<D ; d++)
            for(int row =0 ; row < NCOLOR; row++)
                for(int col = 0; col < NCOLOR; col++){
                 out<<X.get(s,d,row,col)<<"\t";  // element wise output
        }
    return out;
}*/

Matrix operator * (const  Matrix &o1 , const Matrix &o2 ){
  Matrix S;  
  Umatrix u ;
  for (int s = 0; s < SITES; s++){
    for (int d = 0; d< D ; d++) {
      Umatrix u1(s,d,o1);
      Umatrix u2(s,d,o2);
      u = u1*u2;
       for(int row =0 ; row < NCOLOR; row++){
        for(int col = 0; col < NCOLOR; col++){
           S.set(s,d,row,col,u.get(row,col));
         }
       }    
      }
  }  
  return(S);  

}

Matrix operator *(const double &o2 , Matrix &o1){

  Matrix S;  
  Umatrix dum;
  for (int s = 0; s < SITES; s++){
    for (int d = 0; d< D ; d++) {
      for(int row =0 ; row < NCOLOR; row++){
        for(int col = 0; col < NCOLOR; col++){
          S.set(s,d,row,col, o1.get(s,d,row,col)* Complex(o2,0.0));
          }
      }
    } 
  }  
  return(S);
}

/*Matrix operator /(const  Matrix &o1 , const double o2){
  Matrix S;  
  Umatrix dum;
  for (int s = 0; s < SITES; s++){
    for (int d = 0; d< D ; d++) {
      for(int row =0 ; row < NCOLOR; row++){
        for(int col = 0; col < NCOLOR; col++){
          S.set(s,d,row,col, o1.get(s,d,row,col)/ Complex(o2,0.0));
          }
      }
    } 
  }  
  return(S);
}*/

Matrix operator +(const Matrix &o1 ,const  Matrix& o2){
  Matrix S;  
    Umatrix dum;
    for (int s = 0; s < SITES; s++){
      for (int d = 0; d< D ; d++) {
        for(int row =0 ; row < NCOLOR; row++){
          for(int col = 0; col < NCOLOR; col++){
            S.set(s,d,row,col, o1.get(s,d,row,col) + o2.get(s,d,row,col));
            }
        }
      } 
    }  
    return(S);
  }     

Matrix operator -(const Matrix &o1 ,const  Matrix& o2){
  Matrix S;  
    Umatrix dum;
    for (int s = 0; s < SITES; s++){
      for (int d = 0; d< D ; d++) {
        for(int row =0 ; row < NCOLOR; row++){
          for(int col = 0; col < NCOLOR; col++){
            S.set(s,d,row,col, o1.get(s,d,row,col) - o2.get(s,d,row,col));
            }
        }
      } 
    }  
    return(S);
  }     



Complex Trace(const Matrix X ){        // Trace as sum over all sites and scalars
 Complex z;
  for(int s =0 ; s< SITES; s++){
    for (int d=0 ; d< D ; d++){
    
    Complex Utrace = Complex();
    for (int i = 0; i < NCOLOR; i++){
      Utrace = Utrace + X.get(s,d,i, i);
    }
    z = z+ Utrace ;
    }  
    }
    return(z);
  }
  

////*********************
Umatrix ::Umatrix(void){

for(int i = 0; i < NCOLOR; i++) 
    for (int j = 0; j < NCOLOR; j++) 
      mat[i][j] = Complex(); //all complex entries 0.0
}

Umatrix::Umatrix(int s,int d, Matrix X){
  
    for(int i = 0; i < NCOLOR; i++) {
      for (int j = 0; j < NCOLOR; j++) {
        mat[i][j] = X.get(s,d,i,j); //all complex entries 0.0
      }
    }   
}

Umatrix::Umatrix(Complex m[NCOLOR][NCOLOR]){
    for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) mat[i][j] = m[i][j];
}

Complex Umatrix::get(int i ,int j)  const { // element of particular row,column
    return(mat[i][j]);
}

void Umatrix::set(int i, int j, const Complex o){
    mat[i][j] = o;
}

void Umatrix::print(void){
     for (int i = 0; i < NCOLOR; i++) {
    for (int j = 0; j < NCOLOR; j++) {
      mat[i][j].print(); // from complex numbers
      cout<<"\t";       // arranging as row
      }
    cout<<"\n"; // arranging as column
    }
}   

/*ostream& operator<<(ostream& out , Umatrix s){
    for (int i = 0; i < NCOLOR; i++)
        for (int j = 0; j < NCOLOR; j++) {
        out<< s.get(i,j)<<"\t";  // element wise output
        }
    return out;
}

istream& operator>>(istream& in , Umatrix& s ){
    Complex v[NCOLOR][NCOLOR];              // v[i][j] = mat[i][j]
    for (int i = 0; i < NCOLOR; i++)
        for (int j = 0; j < NCOLOR; j++) {
             in >> v[i][j]; // v[i][j] = mat[i][j];
        }
    s= Umatrix(v);
    return in;
}*/

Umatrix Adj(const Umatrix &u){
    
    Umatrix res; //defining res as umatrix
   for (int i = 0; i < NCOLOR; i++)
       for (int j = 0; j < NCOLOR; j++) {
            res.set(i, j , conjug(u.get(j, i)));
       }
   return(res);
}

Umatrix operator *(const Umatrix &o1, const Umatrix &o2) {
 
  Umatrix r;
  Complex dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum = Complex();
      for (int k = 0; k < NCOLOR; k++) {
        dum = dum + o1.get(i, k) * o2.get(k, j);
      }
      r.set(i, j, dum);
    }
  return(r);
}

Umatrix operator *(const Umatrix &o1, const Complex &o2) {
   
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum.set(i, j, o1.get(i, j)*o2);
    }
  return(dum);
}

Umatrix operator *(const Complex &o2, const Umatrix &o1) {
  
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum.set(i, j, o1.get(i, j)*o2);
    }
  return(dum);
}

Umatrix operator *(const Umatrix &o1, const double o2) {
    
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum.set(i, j, o1.get(i, j) * o2);
    }
  return(dum);
}

Umatrix operator *(const double o2, const Umatrix &o1) {
    
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum.set(i, j, o1.get(i, j) * o2);
    }
  return(dum);
}
Umatrix comm(const Umatrix &o1, const Umatrix &o2) {
  return(o1*o2 - o2*o1);
}

Umatrix operator /(const Umatrix &o1, const Complex &o2 ){
 
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum.set(i, j, o1.get(i, j) / o2);
    }
  return(dum);
}

Umatrix operator/ (const Umatrix &o1, const double o2)
{
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      dum.set(i, j, o1.get(i, j)/ o2);
    }
  return(dum);
}

Umatrix operator +(const Umatrix &x, const Umatrix &y) {
 
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, x.get(i, j) + y.get(i, j));
  return(dum);
}

Umatrix operator -(const Umatrix &x, const Umatrix &y) {
    
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, x.get(i, j) - y.get(i, j));
  return(dum);
}


/*Umatrix exp(const Umatrix &u) {
    
  Umatrix c, del, prod;
  double fac = 1.0;
  int i = 1;

  static int sum = 0, counter = 0;
  
  do {
    fac = fac * (double)i;
    prod = prod * u;
    del = prod * (1.0/fac);
    c = c + del;
    i++;
  }
  while ( sqrt(Tr(del*Adj(del)).real()) > GAUGETOL );
  sum += i;
  counter++;
  if (counter == 1000) {
    cout << "Mean no. of terms in exp() "
    << (double)sum/counter << "\n" << flush;
    counter = 0;
    sum = 0;
  }
  return(c);
}*/

Complex Tr(const Umatrix &o) {
  Complex dum = Complex();
  for (int i = 0; i < NCOLOR; i++)
    dum = dum + o.get(i, i);
  return(dum);
} 



