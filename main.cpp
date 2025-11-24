#include <iostream>
#include "lib.h"
#include "args.h"
#include <chrono>
#include <fenv.h>
int feenableexcept(int excepts);
int fedisableexcept(int excepts);
int fegetexcept(void);


using namespace std;

int main(int argc, char *argv[])
{

    int n, m , r, s, p, task = 10, thr = 0;
    // int _c=0;
    char *filename = nullptr;
    double r1=0, r2=0;
    double *a, *b, *x, *realx;
    // double el;
    double t1=0, t2=0;
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);


    
    filename=filename;

    if(!((argc == 6 || argc == 7) && sscanf(argv[1], "%d", &n)==1 && sscanf(argv[2], "%d", &m)==1 && sscanf(argv[3], "%d", &p)==1 && sscanf(argv[4], "%d", &r)==1 && sscanf(argv[5], "%d", &s)==1)) 
    {
        cout<<"Usage : "<<argv[0]<<" <n> "<<" <m> "<<" <p> "<<" <r> "<<" <s>"<<endl;
      
        return 0;
    }

    if(argc == 7) filename = argv[6];

    

    if(m<=0 || n<0 || r<0 || s<0 || p < 0)
    {
        printf("<n> or <m> or <r> or <s> or <p> <= 0, usage m,n,r,s,p > 0");
        return 0;
    }

    if(argc == 7 && s!=0)
    {
        printf("Wrong usage! If s!=0 dont use initialization from file or s == 0 and file name not specified \n");
        // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        return 0;
    }else if(argc == 6 && (s<1 || s>4))
    {
        printf("bad argument for s, s=1,2,3,4 \n");
        // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        return 0;
    }

    // printf("n = %d m = %d p = %d r = %d s = %d\n",n,m,p,r,s);

    

    a = new double[n*n]; //create matrix a
    b = new double[n];  // create vector b
    x = new double[n];  // create vector x
    realx = new double[n];  // create vector real x
    int *mainblocks = new int[p];
    double *minnorms = new double[p];
    args *ap = new args[p];
    pthread_t *tid = new pthread_t[p];
    pthread_barrier_t barrier;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    
    for(int i = 0 ; i < p; i++)
    {
        mainblocks[i] = -1;
        minnorms[i] = -1;
    }

    pthread_barrier_init(&barrier,0,p);

    for(thr = 0 ; thr < p ; thr ++)
    {
        ap[thr].n = n;
        ap[thr].m = m;
        ap[thr].a = a;
        ap[thr].b = b;
        ap[thr].x = x;
        ap[thr].r = r;
        ap[thr].s = s;
        ap[thr].name = filename;
        ap[thr].thr = thr;
        ap[thr].p = p;
        ap[thr].err = 0;
        ap[thr].mutex = &mutex;
        ap[thr].barrier = &barrier;
        ap[thr].mainblocks = mainblocks;
        ap[thr].minnorms = minnorms;
    }

    double elapsed = get_full_time();

    auto start_sol= std::chrono::high_resolution_clock::now();

    for(thr = 1; thr < p; thr++)
    {
        if(pthread_create(tid+thr,0,parallelSolve,ap+thr))
        {
            printf("ERROR: Cannot create thread %d\n",thr);
            delete []a;
            delete []b;
            delete []x;
            delete []realx;
            delete []tid;
            delete []ap;
            pthread_mutex_destroy(&mutex);
            pthread_barrier_destroy(&barrier);

            return -1;
        }
    }

    parallelSolve(ap+0);


    for(thr = 1; thr < p; thr++)
    {
        if(pthread_join(tid[thr], 0))
        {
            printf("ERROR: Cannot join thread %d\n",thr);
            delete []a;
            delete []b;
            delete []x;
            delete []realx;
            delete []tid;
            delete []ap;
            pthread_mutex_destroy(&mutex);
            pthread_barrier_destroy(&barrier);

            return -1;
        }

        elapsed = get_full_time() - elapsed;
        printf("CPU Time thread %d = %.2lf\n",thr,elapsed);
    }

    auto end_sol= std::chrono::high_resolution_clock::now();

    

    
    
    

    //print vector x
    cout<<"\nSolution vector x : ";
    printlxn(x,n,1,n,r);
    
    

    auto start_res= std::chrono::high_resolution_clock::now();

    

    //reinit


    // if(argc == 5 && s>=1 && s<=4)
    // {
    //     init(a,&f,n,s);
    // }

    // if(argc == 6 && s == 0 ) 
    // {
    //     filename = argv[5];
    //     FILE *file = fopen(filename,"r");
        

    //     while(fscanf(file,"%lf",&el)==1)
    //     {
    //         if(_c<n*n) 
    //         {
    //             a[_c] = el;
    //             _c++;
    //         }

        
    //     }
    //     if(!feof(file) || _c!=n*n){
    //         printf("Bad file\n");
    //         report(argv[0],task,r1,r2,t1,t2,s,n,m); 
    //         fclose(file);
    //         delete []a;
    //         delete []b;
    //         delete []x;
    //         delete []realx;
    //         return 0;
    //     }
    //     fclose(file);
    // }
    // else if(argc == 6 && s!=0)
    // {
    //     printf("Wrong usage! If s!=0 dont use initialization from file or s == 0 and file name not specified \n");
    //     report(argv[0],task,r1,r2,t1,t2,s,n,m); 
    //     delete []a;
    //     delete []b;
    //     delete []x;
    //     delete []realx;
    //     return 0;
    // }

  

    //reinit vector b

    // for (int i = 0; i < n; i++)
    // {   
    //     double sumbi = 0;
    //     for(int k = 0; k <(n-1)/2+1 ; k++)
    //     {
    //         sumbi+= a[i*n+2*k];
            
            
    //     }
        
    //     b[i] = sumbi;
        
    // }

    double *Ax = new double[n];//mat_x_vector(a,x,n);
    double *Ax_b = new double[n];//vectorsub( Ax , b, n);
    double *x_realx = new double[n];//vectorsub( x , realx, n);

    for(int i = 0 ; i < n; i++)
        realx[i] = (i+1)%2;
    

    residuals(r1,r2,a,b,x,realx,n,Ax,Ax_b,x_realx);

    delete []Ax;
    delete []Ax_b;
    delete []x_realx;
    

    auto end_res= std::chrono::high_resolution_clock::now();

     

    t1 = chrono::duration<double>(end_sol - start_sol ).count();
    t2 = chrono::duration<double>(end_res - start_res).count();



    report(argv[0],task,r1,r2,t1,t2,s,n,m,p); 

    delete []a;
    delete []b;
    delete []x;
    delete []realx;

    pthread_mutex_destroy(&mutex);
    pthread_barrier_destroy(&barrier);

    return 0;
}
