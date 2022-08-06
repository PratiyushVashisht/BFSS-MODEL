// BFSS MODEL SIMULATION PART 3
#include"eigen.cpp"


// parameters
const double  m = 9/(SITES*TEMP);        //mass term
const int SWEEPS =1000000;
const int L = 3;                      // no. of steps leap forg
const double EPS =0.0009;              //step size

// ************function: 

Complex action(const Matrix,const Matrix);
Matrix RandM_uni(Matrix); // update Matrix uniform
Matrix RandM_G(Matrix X);
double evolve(Matrix&, Matrix& ,Complex& ,Complex& );
double energy(Matrix X);
double Extent_space(Matrix X);
Complex Poly_loop(Matrix X);

// 
int main()
{   
     //1 seed setting:
     double seed = 41;
     srand48(seed);

    //2 initial set up 
    
    Complex S_i,S_f;
    double dS,exp_dS,r;
    int sweep ,count =0,accept = 0;
    double E , EX;
    Complex pkl;
  
    double acc_rate =0.0;


      //Initial Configuration for x******************
    
    const double v = 0.07;
    Matrix I = RandM_uni(I); 
    Matrix Z =  RandM_G(Z) ;
    for (int s=0; s<SITES ; s++){
        Umatrix F(s,0,Z);
        Umatrix u = unitary(F);
        for(int row =0 ; row < NCOLOR; row++){
            for(int col = 0 ; col <NCOLOR; col++){
                   Z.set(s,0,row, col, u.get(row,col));
            }
        }
    }
    Matrix x = v*I;
    Matrix g = Z;

    Matrix g_old;
    Matrix x_old;

    static int first_time = 1;
    static int second_time =1;
    static ofstream fobs;
    static ofstream Eobs;    
    if (first_time)
    {
        fobs.open("gmc2m9.txt");
        if (fobs.bad()){
            cout << "Failed to open observable file\n"<<flush;
        }
        first_time=0;
    }
       if (second_time)
    {
        Eobs.open("eigen9.txt");
        if (Eobs.bad()){
            cout << "Failed to open observable file\n"<<flush;
        }
        second_time=0;
    }
    ///****************************************************************
    double E_coeff = (1.0*TEMP/(NCOLOR*4.0*LAMBDA));
    double R_coeff = (1.0/(NCOLOR*LAMBDA*2.0*SITES*SITES*TEMP));

    for (int sweep =0; sweep< SWEEPS; sweep++)
    {
      g_old = g;
      x_old = x;
      evolve(x, g, S_i, S_f);
      
      r =drand48();
      dS = (S_f-S_i).real() ;
      exp_dS = exp(-1.0*dS);
      if (exp_dS > r)
      {
        
        accept++ ; //proposed state is accepted
      }
      else
      {
          x = x_old; // proposed states rejected
          g = g_old;
      }
      count++;
      
      if (count%10 ==0)
      {
        acc_rate =double(accept)/count;           //accepted state ./. proposed state
        //cout<<"acceptance rate ="<< acc_rate<<endl;
        E=0.0 , EX =0.0 , pkl = Complex(0.0,0.0) ;
        E = E_coeff* energy(x);
        EX = R_coeff* Extent_space(x);
        pkl =  Poly_loop(g);
        double modp = sqrt(pkl.real()*pkl.real() + pkl.imag()*pkl.imag());
        count= 0;
        accept =0 ;
        Umatrix e(0,0,x);
        double y = action(x,g).real()/(NCOLOR*NCOLOR);
        fobs<< sweep<<"\t"<<y<<"\t"<<acc_rate<<"\t"<<exp_dS<<"\t"<< E <<"\t"<< EX <<"\t"<<modp<<endl;   
	Eobs<<eigen(e)<<endl;
      //cout<<sweep<<"\t"<<y<<"\t"<<acc_rate<<"\t"<<exp_dS<<"\t"<< E <<"\t"<< EX <<"\t"<< pkl<<"\t"<<modp <<endl;  
      }
    } 
    return 0;
}

//     GUAss Matrix ****************************************************

Matrix RandM_uni(Matrix X){ 
  int l = NCOLOR;
    for(int s= 0; s<SITES; s++){
        for(int d=0; d<D ; d++){
            for(int row =0 ; row < NCOLOR; row++){
                for(int col = row +1 ; col <NCOLOR; col++){ 
                   
                    X.set(s, d, row, col , Complex(drand48()-0.5,drand48()-0.5)); 
                    X.set(s,d,col,row, conjug(X.get(s,d,row,col)));         //HErmitian :
                }
                X.set(s,d,row,row, Complex(drand48()-0.5,0));
            }
        Umatrix m(s,d,X);
        Complex z = Tr(m); 
        Complex avg = z/l;
        for(int n=0; n<NCOLOR ; n++){
           X.set(s, d, n, n ,X.get(s,d,n,n)- avg);      //Traceless **************
        }
      }
    }
    return(X);
} 

Matrix RandM_G(Matrix X){ 
  int l = NCOLOR;
    for(int s= 0; s<SITES; s++){
            for(int row =0 ; row < NCOLOR; row++){
                for(int col = 0 ; col <NCOLOR; col++){ 
                    X.set(s, 0, row, col , Complex(drand48()-0.5,drand48()-0.5)); 
               }
            }
   /* Umatrix m(s,0,X);
        Complex z = Tr(m); 
        Complex avg = z/l;
        for(int n=0; n<NCOLOR ; n++){
           X.set(s, 0, n, n ,X.get(s,0,n,n)- avg);      //Traceless **************
        }
	*/	
    //********** all scalars should be same:************//
        for (int d =0; d<D; d++){
          for(int row =0 ; row < NCOLOR; row++){
            for(int col = 0 ; col <NCOLOR; col++){
                X.set(s,d,row,col, X.get(s,0,row,col)); 
            }
          }
        }
    }
   
  return(X);
} 
// Action S -------------------------------

Complex action(const Matrix X, const Matrix G){   

    Umatrix X_1 ;
    Umatrix X_2 ;
    Umatrix X_3 ;
    Complex z   ;
    double k = NCOLOR/(4.0*LAMBDA);

    for(int s=0; s<SITES; s++){
      for(int d=0; d <D; d++){

        Umatrix U1(s,d,X) ; 
        Umatrix U2((s+1)%SITES,d,X);

        Umatrix G1(s,0,G);
        Umatrix G2 = Adj(G1);

        X_1 =  (U1*G1*U2*G2);

        for(int j=d+1; j<D; j++){
          Umatrix X_2;
          Umatrix U4(s,j,X);
          X_2 = X_2 + 0.5*(comm(U1,U4))*(comm(U1,U4));
          }

        X_3 =  (m*m+1)*(U1)*(U1);
        z = z+ Tr(-1.0*X_1 - X_2 + X_3);
        
        X_2 = X_2 -X_2;
      }
    } 
  Complex S = k*z ;
  z=z-z;
  return( S);
}

///**********************
double evolve(Matrix &X, Matrix &G, Complex& S_i ,Complex& S_f){
  int i;
  S_i= action(X,G);
  
   //Step2,3,4,........,l
    for (i=0; i<L; i++)
    {
        Matrix P1 = RandM_uni(X);
        Matrix P2 = RandM_G(G);

        X = X + EPS*P1 ;
        G = G + EPS*P2 ; 
        for (int s=0; s<SITES; s++){
          Umatrix F(s,0,G);
          Umatrix u = unitary(F);
                  
          for(int row =0 ; row < NCOLOR; row++){
              for(int col = 0 ; col <NCOLOR; col++){
                    G.set(s,0,row, col, u.get(row,col));
            }
          }
          Umatrix F_1(s,0,G);
          //cout <<(F_1)*Adj(F_1)<< endl;
        }
    }
    
    // Calculate final action H_f
    S_f = action(X,G);
  return  (0.0) ;  
}

//Internal Energy :
  
double energy(Matrix X){
  Umatrix x_1;
  Umatrix x_2;
  Complex z;

  for(int s=0;s<SITES; s++){
    for(int d=0; d <D; d++){
      Umatrix U1(s,d,X) ;  
        
     x_2 = x_2-x_2;
    for(int j=d+1; j<D; j++)
    {
    Umatrix U4(s,j,X);
    x_2 = x_2 + 0.5*(comm(U1,U4)*comm(U1,U4));
    }
    x_1 =  (m*m)*(U1*U1);
      z = z+ Tr(-1.5*(x_2)+2.0*(x_1));
      x_2 = x_2-x_2;
  } 
  }
  double E = z.real() ;
  z = z-z;
  return(E);
}

double Extent_space(Matrix X){
    Umatrix X_3;
    Complex z ;
    z = z-z;
    for(int s=0;s<SITES; s++){
      for(int d=0; d <D; d++){
        Umatrix U1(s,d,X) ; 
        
        X_3 =  (U1*U1);
        z = z+ Tr(X_3);
        // Umatrix X_3;
      }
    } 
  double ES = z.real() ;
  z = z-z;
  return( ES);
}


Complex  Poly_loop(Matrix X){
  Complex t;
  Umatrix y;
  for(int s=0; s<SITES; s++)
  {
    Umatrix g1(s,0,X);
    if (s == 0) {
     y =g1;
  }
  else {
    y =(g1)*y;
  }
  }
  t= Tr(y)/NCOLOR;
 
  return(t);
}
