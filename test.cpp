#include <omp.h>
#include<stdio.h>
#include<iostream>
using namespace std;

void add(int a)
{
    cout<<"working "<<a<<"\n";
//    return (a+b);
}

int main()
{
    int nthreads, tid;
    omp_set_num_threads(5);
/* Fork a team of threads with each thread having a private tid variable */
    #pragma omp parallel private(tid)
    {
        add(omp_get_thread_num());            

  /* Only master thread does this */
//      if (tid == 0) 
//        {
//        nthreads = omp_get_num_threads();
//        printf("Number of threads = %d\n", nthreads);
//        }
    } 
    return 0;
}


