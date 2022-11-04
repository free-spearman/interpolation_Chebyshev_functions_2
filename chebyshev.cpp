#include "chebyshev.hpp"
#include "auto_parallelization.hpp"
#define PI 3.1415926535897932384626433832795

void chebyshevpoints(double *cx, double *cy, int nx, int ny, double a, double b, double c, double d){
    double stepx=PI/nx;
    double anglex=stepx/2.0;
    double semisumx=(a+b)/2;
    double semidifx=(b-a)/2;
    double stepy=PI/ny;
    double angley=stepy/2.0;
    double semisumy=(c+d)/2;
    double semidify=(d-c)/2;
    for(int i=nx+1;i>=0;--i){
        if(i==nx+1)
            cx[nx+1]=b;
        else if(i==0)
            cx[0]=a;
        else{
            cx[i]=cos(anglex);
            cx[i]*=semidifx;
            cx[i]+=semisumx;
            anglex+=stepx;
        }
    }
    for(int j=ny+1;j>=0;--j){
        if(j==ny+1)
            cy[ny+1]=d;
        else if(j==0)
            cy[0]=c;
        else{
            cy[j]=cos(angley);
            cy[j]*=semidify;
            cy[j]+=semisumy;
            angley+=stepy;
        }
    }
}


void Fill_F(double *F, int nx, int ny, double *cx, double *cy, double (*f)(double, double), size_t nthreads){
    if (nthreads > 1){
        Fill_F_arg arg;
        arg.F = F; arg.cx = cx; arg.cy = cy;
        arg.nx = nx; arg.ny = ny;
        arg.f = f; 
        auto_parallelization((void*) &arg,
        nthreads,
        nx+1+1,
        pf_Fill_F
        );
        return ;
    }
        
    for(int i=0;i<=nx+1;++i){
        for(int j=0;j<=ny+1;++j){
            F[i*(ny+2)+j]=f(cx[i], cy[j]);
        }
    }
}

double chebyshevtrans(double a, double b, double x){
    return (2*x-(b+a))/(b-a);
}

void fill_Fx(double *Fx, double *cx, int nx, double a, double b){
    for(int j=0; j<nx; ++j){
        Fx[j]=1;
        Fx[nx+j]=chebyshevtrans(a, b, cx[j+1]);
    }
    for(int j=0; j<nx; ++j){
        double trans=2*chebyshevtrans(a, b, cx[j+1]);
        for(int i=2; i<nx; ++i){
            Fx[i*nx+j]=trans*Fx[(i-1)*nx+j]-Fx[(i-2)*nx+j];
        }
    }
}

void fill_Fy(double *Fy, double *cy, int ny, double c, double d){
    for(int i=0; i<ny; ++i){
        Fy[i*ny]=1;
        Fy[i*ny+1]=chebyshevtrans(c, d, cy[i+1]);
    }
    for(int i=0; i<ny; ++i){
        double trans=2*chebyshevtrans(c, d, cy[i+1]);
        for(int j=2; j<ny; ++j){
            Fy[i*ny+j]=trans*Fy[i*ny+j-1]-Fy[i*ny+j-2];
        }
    }
}

double fscalar_ij(int i, int nx, double *Fx, double *F, double *Fy, int ny, int j){
    double res=0;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l)
            res+=(Fx[i*nx+k]*F[(k+1)*(ny+2)+l+1]*Fy[l*ny+j]);
    }
    res/=nx;
    res/=ny;
    if(i)
        res*=2;
    if(j)
        res*=2;
    return res;
}

double scalar_ij(int i, int nx, double *Fx, double *T, double *Fy, int ny, int j){
    double res=0;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l)
            res+=(Fx[k*nx+i]*T[k*ny+l]*Fy[j*ny+l]);
    }
    return res;
}

double scalar_aj(int nx, double *T, double *Fy, int ny, int j){
    double res=0, u=1;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l)
            res+=(u*T[k*ny+l]*Fy[j*ny+l]);
        u*=(-1);
    }
    return res;
}

double scalar_bj(int nx, double *T, double *Fy, int ny, int j){
    double res=0;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l)
            res+=(T[k*ny+l]*Fy[j*ny+l]);
    }
    return res;
}

double scalar_ic(int i, int nx, double *Fx, double *T, int ny){
    double res=0, u;
    for(int k=0; k<nx; ++k){
        u=1;
        for(int l=0; l<ny; ++l){
            res+=(Fx[k*nx+i]*T[k*ny+l]*u);
            u*=(-1);
        }
    }
    return res;
}

double scalar_id(int i, int nx, double *Fx, double *T, int ny){
    double res=0;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l)
            res+=(Fx[k*nx+i]*T[k*ny+l]);
    }
    return res;
}

double scalar_ac(int nx, double *T, int ny){
    double res=0, u=1, v;
    for(int k=0; k<nx; ++k){
        v=1;
        for(int l=0; l<ny; ++l){
            res+=(u*T[k*ny+l]*v);
            v*=(-1);
        }
        u*=(-1);
    }
    return res;
}

double scalar_ad(int nx, double *T, int ny){
    double res=0, u=1;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l){
            res+=(u*T[k*ny+l]);
        }
        u*=(-1);
    }
    return res;
}

double scalar_bc(int nx, double *T, int ny){
    double res=0, v;
    for(int k=0; k<nx; ++k){
        v=1;
        for(int l=0; l<ny; ++l){
            res+=(T[k*ny+l]*v);
            v*=(-1);
        }
    }
    return res;
}

double scalar_bd(int nx, double *T, int ny){
    double res=0;
    for(int k=0; k<nx; ++k){
        for(int l=0; l<ny; ++l)
            res+=T[k*ny+l];
    }
    return res;
}

void interpolation_tensor(double *T, double *Fx, double *F, double *Fy, int nx, int ny, size_t nthreads){
    if (nthreads > 1){
        itensor_arg arg;
        
        arg.nx = nx; arg.ny = ny;
        arg.T = T; arg.Fx = Fx;
        arg.F = F; arg.Fy = Fy;

        auto_parallelization((void*) &arg,
            nthreads,
            nx,
            pf_interpolation_tensor
        );
        return ;
    }
    for(int i=0; i<nx; ++i){
        for(int j=0; j<ny; ++j){
            T[i*ny+j]=fscalar_ij(i, nx, Fx, F, Fy, ny, j);
        }
    }
}

void Fill_TT(double *Fx, int nx, double *T, double *Fy, int ny, double *TT, size_t nthreads){
    if ( nthreads > 1){
        TT_arg arg;
        arg.Fx = Fx; arg.nx = nx; arg.T = T;
        arg.Fy = Fy; arg.ny = ny; arg.TT = TT;
        auto_parallelization((void*) &arg,
            nthreads,
            nx,
            pf_TT_scalar_ij
        );
        auto_parallelization((void*) &arg,
            nthreads,
            ny,
            pf_TT_scalar_aj_bj
        );
        auto_parallelization((void*) &arg,
        nthreads,
        nx,
        pf_TT_scalar_ic_id
        );
    }
    else {
        
        for(int i=0;i<nx;++i){
            for(int j=0;j<ny;++j){
                TT[(i+1)*(ny+2)+j+1]=scalar_ij(i, nx, Fx, T, Fy, ny, j);
            }
        } 
    
        for(int j=0;j<ny;++j){
            TT[j+1]=scalar_aj(nx, T, Fy, ny, j);
            TT[(nx+1)*(ny+2)+j+1]=scalar_bj(nx, T, Fy, ny, j);
        }
    
        for(int i=0;i<nx;++i){
            TT[(i+1)*(ny+2)]=scalar_ic(i, nx, Fx, T, ny);
            TT[(i+1)*(ny+2)+ny+1]=scalar_id(i, nx, Fx, T, ny);
        }

    }
        
    TT[0]=scalar_ac(nx, T, ny);
    TT[ny+1]=scalar_ad(nx, T, ny);
    TT[(nx+1)*(ny+2)]=scalar_bc(nx, T, ny);
    TT[(nx+2)*(ny+2)-1]=scalar_bd(nx, T, ny);
}

