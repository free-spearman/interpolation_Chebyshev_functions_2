#ifndef auto_parallelization_hpp
#define auto_parallelization_hpp
#include "chebyshev.hpp"
#include <pthread.h>
typedef struct {
  void *arg;
  size_t begin; //начало куска
  size_t end; // конец куска 
} arg_thread;

//void Fill_TT(double *Fx, int nx, double *T, double *Fy, int ny, double *TT)
typedef struct {
  size_t nx;
  size_t ny;
  double *Fx;
  double *Fy;
  double *T;
  double *TT;
} TT_arg; 

typedef struct {
  size_t nx;
  size_t ny;
  double *F; 
  double *cx;
  double *cy;
  double (*f)(double, double);
} Fill_F_arg;

typedef struct {
  double *Fx;
  double *cx;
  double * Fy;
  double *cy;
  double a;
  double b; 
  double c; 
  double d;
  size_t nx; 
  size_t ny;
} F_xy_arg;

typedef struct {
  size_t nx; 
  size_t ny;
  double *T; 
  double *Fx;
  double *F;
  double *Fy;
} itensor_arg;

void* pf_TT_scalar_ic_id (void *ext_arg);
void* pf_TT_scalar_aj_bj (void *ext_arg);
void* pf_TT_scalar_ij (void *ext_arg);
void* pf_Fill_F( void* ext_arg);
void fill_Fx_Fy(double *Fx,
  double *cx,
  size_t nx,
  double a,
  double b,
  double * Fy,
  double *cy, 
  size_t ny, 
  double c, 
  double d
  );
void* pf_interpolation_tensor (void* ext_arg);

int auto_parallelization ( void *arg,
  size_t nthreads, 
  size_t niterations,
  void* (*action) ( void *) 
  );
#endif /* auto_parallelization */



