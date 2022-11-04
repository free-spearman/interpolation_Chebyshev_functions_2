#include <auto_parallelization.hpp>
/*
Функция авто параллеливания
1. Вызываем функцию до цикла 
2. Функции передают:
	1.аргументы цикла как void указатель 
	2.передаю количество нитей и итераций кикла
	3.исполняющая функция
3. Вызывается исполняющая функция в отдельную нить 
 
нужна структура со структурой, которая будет
в себя включать параметры рабоыт и параметры нити как void*
*/

/*обвертки */
void* pf_TT_scalar_ij (void *ext_arg){

	arg_thread arg_thr = *((arg_thread*) ext_arg);
	TT_arg arg = *((TT_arg*) arg_thr.arg);
	
	double *TT = arg.TT, *Fx = arg.Fx, 
		*Fy = arg.Fy, *T = arg.T;
	size_t nx = arg.nx, ny = arg.ny;    

	for(size_t i = arg_thr.begin; i < arg_thr.end; ++i){
        for(size_t j = 0; j < ny; ++j){
            TT[(i+1)*(ny+2)+j+1]=scalar_ij(i, nx, Fx, T, Fy, ny, j);
        }
    }
    return 0; 
}

/*
for(int j=0;j<ny;++j){
        TT[j+1]=scalar_aj(nx, T, Fy, ny, j);
        TT[(nx+1)*(ny+2)+j+1]=scalar_bj(nx, T, Fy, ny, j);
    }
*/
void* pf_TT_scalar_aj_bj (void *ext_arg){

	arg_thread arg_thr = *((arg_thread*) ext_arg);
	TT_arg arg = *((TT_arg*) arg_thr.arg);
	
	double *TT = arg.TT, *Fy = arg.Fy, *T = arg.T;
	size_t nx = arg.nx, ny = arg.ny;    

	for(size_t j = arg_thr.begin; j< arg_thr.end; ++j){
		TT[j+1]=scalar_aj(nx, T, Fy, ny, j);
        TT[(nx+1)*(ny+2)+j+1]=scalar_bj(nx, T, Fy, ny, j);	
	}
	return 0; 
}


/*
for(int i=0;i<nx;++i){
        TT[(i+1)*(ny+2)]=scalar_ic(i, nx, Fx, T, ny);
        TT[(i+1)*(ny+2)+ny+1]=scalar_id(i, nx, Fx, T, ny);
    }
*/ 
void* pf_TT_scalar_ic_id (void *ext_arg){
	arg_thread arg_thr = *((arg_thread*) ext_arg);
	TT_arg arg = *((TT_arg*) arg_thr.arg);
	
	double *TT = arg.TT, *Fx = arg.Fx, *T = arg.T;
	size_t nx = arg.nx, ny = arg.ny;    
	//arg_thr.begin; j< arg_thr.end
	for(size_t i=arg_thr.begin; i<arg_thr.end; ++i){
        TT[(i+1)*(ny+2)]=scalar_ic(i, nx, Fx, T, ny);
        TT[(i+1)*(ny+2)+ny+1]=scalar_id(i, nx, Fx, T, ny);
    }
    return 0;
}
//Fill_F в auto_parallelization
void* pf_Fill_F( void* ext_arg){
	arg_thread arg_thr = *((arg_thread*) ext_arg);
	Fill_F_arg arg = *((Fill_F_arg*) arg_thr.arg);
	double *F = arg.F, *cx = arg.cx, *cy = arg.cy;
	size_t nx = arg.nx, ny = arg.ny;        
	for(size_t i=arg_thr.begin; i < arg_thr.end; ++i){
		for(size_t j=0;j<=ny+1;++j){
			F[i*(ny+2)+j]=arg.f(cx[i], cy[j]);
		}
	}
	return 0;
}

void* fill_Fx( void* ext_arg){
//double *Fx, double *cx, int nx, double a, double b
	F_xy_arg arg = *((F_xy_arg*) ext_arg);
	
	double *Fx = arg.Fx , *cx = arg.cx; 
	double a = arg.a, b = arg.b;
	size_t nx = arg.nx;

	for(size_t j=0; j<nx; ++j){
		Fx[j]=1;
		Fx[nx+j]=chebyshevtrans(a, b, cx[j+1]);
    	}
	for(size_t j=0; j<nx; ++j){
		double trans=2*chebyshevtrans(a, b, cx[j+1]);
		for(size_t i=2; i<nx; ++i){
			Fx[i*nx+j]=trans*Fx[(i-1)*nx+j]-Fx[(i-2)*nx+j];
        	}
    	}
    	return 0;
}

void* fill_Fy(void* ext_arg){
	F_xy_arg arg = *((F_xy_arg*) ext_arg);

	double *Fy = arg.Fy, *cy = arg.cy;
	size_t ny = arg.ny;
	double c = arg.c, d = arg.d;

    for(size_t i=0; i<ny; ++i){
        Fy[i*ny]=1;
        Fy[i*ny+1]=chebyshevtrans(c, d, cy[i+1]);
    }
    for(size_t i=0; i<ny; ++i){
        double trans=2*chebyshevtrans(c, d, cy[i+1]);
        for(size_t j=2; j<ny; ++j){
            Fy[i*ny+j]=trans*Fy[i*ny+j-1]-Fy[i*ny+j-2];
        }
    }
    return 0;
}

void fill_Fx_Fy (double * Fx, double *cx, size_t nx, double a, double b,
 	double *Fy, double *cy, size_t ny, double c, double d){
	F_xy_arg arg;
	pthread_t thread_x, thread_y;
	arg.Fx = Fx; arg.cx = cx; arg.nx = nx; arg.a = a; arg.b = b;
	arg.Fy = Fy; arg.cy = cy; arg.ny = ny; arg.c = c; arg.d = d;
	pthread_create( &thread_x, 
			NULL, 
			fill_Fx, 
			(void*) &arg);
	pthread_create( &thread_y, 
			NULL, 
			fill_Fy, 
			(void*) &arg);
	pthread_join( thread_x, NULL);
	pthread_join( thread_y, NULL); 
}

void* pf_interpolation_tensor (void* ext_arg){
	arg_thread arg_thr = *((arg_thread*) ext_arg);
	itensor_arg arg = *((itensor_arg*) arg_thr.arg);

	size_t nx = arg.nx, ny = arg.ny;
	double *T = arg.T, *Fx = arg.Fx, 
	*F = arg.F, *Fy = arg.Fy;

	for(size_t i = arg_thr.begin; i<arg_thr.end; ++i){
        for(size_t j=0; j<ny; ++j){
            T[i*ny+j]=fscalar_ij(i, nx, Fx, F, Fy, ny, j);
        }
    }
    return 0;
}

int auto_parallelization ( void *arg,
	//const type_info arg_type,
	size_t nthreads, 
	size_t niterations,
	void* (*action) ( void *) 
	){
	if (niterations < nthreads ) 
		nthreads = niterations;
	
	arg_thread arg_list[nthreads];
	pthread_t threads[nthreads-1];   

	size_t size_one_piece =  niterations / nthreads; // 10/4 = 2 
	size_t additive =  niterations % nthreads; //10%4 = 2
	// 2х3 + 2x2
	size_t active_threads = 0, distributed = 0;
	//запуск ниток
	for (; active_threads < nthreads && additive > 0; active_threads++, additive--){
		//заполняем arg
		arg_list[active_threads].begin = distributed; //0 //3
		distributed += size_one_piece + 1; // 3 ) //6
		arg_list[active_threads].end = distributed;
		arg_list[active_threads].arg = arg;
		// вызов action
		pthread_create( &threads[active_threads], 
			NULL, 
			action, 
			(void*) &arg_list[active_threads]);
	}
	// additive = 0 active_threads = 2
	for (; active_threads < nthreads - 1; active_threads++){
		//заполняем arg
		arg_list[active_threads].begin = distributed; //6//8 
		distributed += size_one_piece; // 8 
		arg_list[active_threads].end = distributed;
		arg_list[active_threads].arg = arg;
		// вызов action
		pthread_create( &threads[active_threads], 
			NULL, 
			action, 
			(void*) &arg_list[active_threads]);
	}
	// active_threads = 3 
	// даем работу основной нитке
	arg_list[active_threads].begin = distributed; //8 
	distributed += size_one_piece; // 10 
	arg_list[active_threads].end = distributed;
	arg_list[active_threads].arg = arg;
	//вызов action
	action( (void*) &arg_list[active_threads] );
	//собрать потоки и выйти
	for (size_t i = 0; i < nthreads - 1 ; i++){
		pthread_join( threads[i], NULL);
		
	}
	return 0;
}