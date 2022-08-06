// BFSS MODEL SIMULATION PART 2
#include"matrix.cpp"

void srand48(long seedval) {
    srand(seedval);
}
long lrand48() {
    return rand();
}
double drand48() {
    return rand() / (double)RAND_MAX;
}


// parameters
const int  m = 1.0/(SITES*TEMP);     //mass;     //mass
const int SWEEPS =200000;
const int L = 10;   // no. of steps
const double EPS =0.0003;  //step size

// ************function: 

Complex action(const Matrix);
Matrix RandM_uni(Matrix); // update Matrix uniform
double evolve(Matrix&, Complex& ,Complex& );
double energy(Matrix X);
double Extent_space(Matrix X);
//Complex EV[NCOLOR];

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
    double E= 0.0, EX = 0.0;

    double obs1 ,obs1_e  ,std_err_obs1;
    double obs2 ,obs2_e  ,std_err_obs2;
  
    double acc_rate =0.0, avg_acc_rate = 0.0, total_acc_rate = 0.0;
    double tot_count = 0.0 ,avg_E =0.0, avg_EX = 0.0 ;

      //Initial Configuration for x******************
    
    const double v =  0.1;
    Matrix I = RandM_uni(I); 
    Matrix x = v*I;

    Matrix x_old;
    static int first_time = 1;
    static ofstream fobs;

    if (first_time)
    {
        fobs.open("monted1.txt");
        if (fobs.bad()){
            cout << "Failed to open observable file\n"<<flush;
        }
        first_time=0;
    }
    ///****************************************************************
    double G = (1.0/(NCOLOR*4.0*pow(LAMBDA,1.333)*SITES));
    double R =0.5/(NCOLOR*pow(LAMBDA,0.667)*SITES);

 for (int sweep =0; sweep< SWEEPS; sweep++)
  {
        
      x_old = x;
      evolve(x,S_i,S_f);
      
      r =drand48();
      dS =(S_f-S_i).real();
      exp_dS = exp(-1.0*dS);
      if (exp_dS > r)
      {
        accept++ ; //proposed state is accepted
      }
      else
      {
          x = x_old; // proposed states rejected
      }
      count++;
      
      if (count%10 ==0)
      {
        acc_rate =double(accept)/count;           //accepted state ./. proposed state
        //cout<<"acceptance rate ="<< acc_rate<<endl;
        total_acc_rate =total_acc_rate+acc_rate;        // adding all acceptance rates
      
        E = G*energy(x);
        EX =R* Extent_space(x);

        tot_count++;                     
          // total number of counts, No. of time loop runs
        count= 0;
        accept =0;   

        double y = action(x).real()/(NCOLOR*NCOLOR);
        fobs<< sweep<<"\t"<<y<<"\t"<<acc_rate<<"\t"<<exp_dS<<"\t"<< E <<"\t"<< EX <<"\t" <<endl;   
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
           X.set(s, d, n, n ,X.get(s,d,n,n)- avg); //real and imaginary part
        }
           
      }   
    }
    return(X);
} 
// Action S -------------------------------
    Matrix se;
    Umatrix X_1;
    Umatrix X_2;
    Umatrix X_3;
Complex action(const Matrix X){   
    Matrix se;
    Umatrix X_1;
    Umatrix X_2;
    Umatrix X_3;
    Complex z ;
    double k = NCOLOR/(4.0*LAMBDA);
    for(int s=0; s<SITES; s++){
      for(int d=0; d <D; d++){
        Umatrix U1(s,d,X) ; 
        Umatrix U2((s-1+SITES)%SITES,d,X);
        X_1 =  ((U1 - U2) *(U1 - U2));
        
          for(int j=d+1; j<D; j++){
        
          Umatrix U4(s,j,X);
        X_2 = X_2 + (comm(U1,U4)*comm(U1,U4));
        }
        X_3 =  (m*m*U1*U1);
        z = z+ Tr(X_1 -X_2+X_3);
        Umatrix X_2;
      }
    } 
  Complex S = k*z ;
  return( S);
}

///**********************
double evolve(Matrix &X, Complex& S_i ,Complex& S_f){
  int i;
  
  S_i= action(X);
  
   //Step2,3,4,........,l
    for (i=0; i<L; i++)
    {
        Matrix y;
        y = RandM_uni(y);
        X = X + EPS*y ;
    }

    // Calculate final action S_f
    S_f = action(X);
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

    for(int j=d+1; j<D; j++)
    {
    Umatrix U4(s,j,X);
    x_2 = x_2 + 0.5*(comm(U1,U4)*comm(U1,U4));
    }
    X_1 =  (m*m*U1*U1);
      z = z+ Tr(-1.5*(x_2 +X_1));
      Umatrix x_2;
  } 
  }
  double E = z.real() ;
  return(E);
}

double Extent_space(Matrix X){
    Umatrix X_3;
    Complex z ;
    
    for(int s=0;s<SITES; s++){
      for(int d=0; d <D; d++){
        Umatrix U1(s,d,X) ; 
        
        X_3 =  (U1*U1);
        z = z+ Tr(X_3);
         Umatrix X_3;
      }
    } 
  double ES = z.real() ;
  return( ES);
}