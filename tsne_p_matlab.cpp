/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

/*
modified by Yinqing Li
remove dependence on cblas
*/

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

/*
extern "C" {
    #include <cblas.h>
}
*/

//using namespace std;

#define S sizeof(double)

void tsne_run(double* P, double* Y, int N, int no_dims);
void computeExactGradient(double* P, double* Y, int N, int D, double* dC, double* DD, double* Q, double* dataSums);
double evaluateError(double* P, double* Y, int N, int D, double* DD, double* Q, double* dataSums);
void computeSquaredEuclideanDistance(double* X, int N, int D, double* DD, double* dataSums);
void zeroMean(double* X, int N, int D, double* mean);
double randn();

// Perform t-SNE
void tsne_run(double* P, double* Y, int N, int no_dims) {
  // Set learning parameters
  float total_time = .0;
  clock_t start, end;
	int max_iter = 1000, stop_lying_iter = 100, mom_switch_iter = 250;
	double momentum = .5, final_momentum = .8;
	double eta = 200.0;
    
  // Allocate some memory
  double* dY    = (double*) mxMalloc(N * no_dims * sizeof(double));
  double* uY    = (double*) mxMalloc(N * no_dims * sizeof(double));
  double* gains = (double*) mxMalloc(N * no_dims * sizeof(double));
  if(dY == NULL || uY == NULL || gains == NULL) { mexErrMsgTxt("Memory allocation failed! tsne_run\n"); }
    
  //allocate all necessary Memory
  double* DD = (double*) mxMalloc(N * N * sizeof(double));
  if(DD == NULL) { mexErrMsgTxt("Memory allocation failed! DD\n");}
  double* Q = (double*) mxMalloc(N * N * sizeof(double));
  if(Q == NULL) { mexErrMsgTxt("Memory allocation failed! Q\n"); }
  double* dataSums = (double*) mxCalloc(N, sizeof(double));
  if(dataSums == NULL) { mexErrMsgTxt("Memory allocation failed! dataSums\n"); }
  double* mean = (double*) mxCalloc(no_dims, sizeof(double));
  if(mean == NULL) { mexErrMsgTxt("Memory allocation failed! mean\n"); }
  
  for(int i = 0; i < N * no_dims; i++)    uY[i] =  .0;
  for(int i = 0; i < N * no_dims; i++) gains[i] = 1.0;

  // Lie about the P-values
  for(int i = 0; i < N * N; i++)        P[i] *= 12.0;

	// Initialize solution (randomly)
	for(int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;
	
	// Perform main training loop	
  start = clock();
	for(int iter = 0; iter < max_iter; iter++) {
        
    // Compute (approximate) gradient
    computeExactGradient(P, Y, N, no_dims, dY, DD, Q, dataSums);

    // Update gains
    for(int i = 0; i < N * no_dims; i++) gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
    for(int i = 0; i < N * no_dims; i++) if(gains[i] < .01) gains[i] = .01;
            
    // Perform gradient update (with momentum and gains)
    for(int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
		for(int i = 0; i < N * no_dims; i++)  Y[i] = Y[i] + uY[i];
        
    // Make solution zero-mean
    zeroMean(Y, N, no_dims, mean);
    
    // Stop lying about the P-values after a while, and switch momentum
    if(iter == stop_lying_iter) {
      for(int i = 0; i < N * N; i++)        P[i] /= 12.0;
    }
    if(iter == mom_switch_iter) momentum = final_momentum;
        
    // Print out progress
    if(iter > 0 && (iter % 100 == 0 || iter == max_iter - 1)) {
      end = clock();
      double C = .0;
      C = evaluateError(P, Y, N, no_dims, DD, Q, dataSums);
      if(iter == 0) mexPrintf("Iteration %d: error is %f\n", iter + 1, C);
      else {
        total_time += (float) (end - start) / CLOCKS_PER_SEC;
        mexPrintf("Iteration %d: error is %f (100 iterations in %4.2f seconds)\n", iter, C, (float) (end - start) / CLOCKS_PER_SEC);
        mexEvalString("drawnow;"); 
      }
      
			start = clock();
    }
  }
  end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;
  
  // Clean up memory
  mxFree(dY);
  mxFree(uY);
  mxFree(gains);
  mxFree(DD); DD = NULL;
  mxFree(Q);  Q  = NULL;
  mxFree(dataSums); dataSums = NULL;
  mxFree(mean); mean = NULL;   
  mexPrintf("Fitting performed in %4.2f seconds.\n", total_time);
  mexEvalString("drawnow;"); 
}


// Compute gradient of the t-SNE cost function (exact)
void computeExactGradient(double* P, double* Y, int N, int D, double* dC, double* DD, double* Q, double* dataSums) {
	
	// Make sure the current gradient contains zeros
	for(int i = 0; i < N * D; i++) dC[i] = 0.0;
    
  // Compute the squared Euclidean distance matrix
  //double* DD = (double*) mxMalloc(N * N * sizeof(double));
  //if(DD == NULL) { mexErrMsgTxt("Memory allocation failed! DD\n");}
  computeSquaredEuclideanDistance(Y, N, D, DD, dataSums);
  
  // Compute Q-matrix and normalization sum
  //double* Q    = (double*) mxMalloc(N * N * sizeof(double));
  //if(Q == NULL) { mexErrMsgTxt("Memory allocation failed! Q\n"); }
  double sum_Q = .0;
  for(int n = 0; n < N; n++) {
  	for(int m = n+1; m < N; m++){
  	//for(int m = 0; m < N; m++){
          if(n != m) {
              Q[n * N + m] = 1 / (1 + DD[n * N + m]);
              Q[m * N + n] = Q[n * N + m];
              sum_Q += Q[n * N + m];
          }
      }
      Q[n * N + n] = DBL_MIN;
  }
  sum_Q = sum_Q*2;
  
	// Perform the computation of the gradient
	for(int n = 0; n < N; n++) {
    	for(int m = 0; m < N; m++) {
            if(n != m) {
                double mult = (P[n * N + m] - (Q[n * N + m] / sum_Q)) * Q[n * N + m];
                for(int d = 0; d < D; d++) {
                    dC[n * D + d] += (Y[n * D + d] - Y[m * D + d]) * mult;
                }
            }
		}
	}
    
  // Free memory
  //mxFree(DD); DD = NULL;
  //mxFree(Q);  Q  = NULL;
}


// Evaluate t-SNE cost function (exactly)
double evaluateError(double* P, double* Y, int N, int D, double* DD, double* Q, double* dataSums) {
    
    // Compute the squared Euclidean distance matrix
    //double* DD = (double*) mxMalloc(N * N * sizeof(double));
    //double* Q = (double*) mxMalloc(N * N * sizeof(double));
    //if(DD == NULL || Q == NULL) { mexErrMsgTxt("Memory allocation failed! DD or Q\n"); }
    computeSquaredEuclideanDistance(Y, N, D, DD, dataSums);

    // Compute Q-matrix and normalization sum
    double sum_Q = DBL_MIN;
    for(int n = 0; n < N; n++) {
      for(int m = n+1; m < N; m++){
      //for(int m = 0; m < N; m++){
            if(n != m) {
                Q[n * N + m] = 1 / (1 + DD[n * N + m]);
                Q[m * N + n] = Q[n * N + m];
                sum_Q += Q[n * N + m];
            }
            else Q[n * N + m] = DBL_MIN;
        }
        Q[n * N + n] = DBL_MIN;
    }
    sum_Q = sum_Q*2;
    for(int i = 0; i < N * N; i++) Q[i] /= sum_Q;
    
    // Sum t-SNE error
    double C = .0;
	for(int n = 0; n < N; n++) {
		for(int m = 0; m < N; m++) {
			C += P[n * N + m] * log(P[n * N + m] / Q[n * N + m]);
		}
	}
    
  // Clean up memory
  //mxFree(DD);
  //mxFree(Q);
	return C;
}


// Compute squared Euclidean distance matrix (not using BLAS)
void computeSquaredEuclideanDistance(double* X, int N, int D, double* DD, double* dataSums) {
    //double* dataSums = (double*) mxCalloc(N, sizeof(double));
    //if(dataSums == NULL) { mexErrMsgTxt("Memory allocation failed! dataSums\n"); }
    
    for(int n = 0; n < N; n++) {
        //initialize dataSums
        dataSums[n] = 0;
        for(int d = 0; d < D; d++) {
            dataSums[n] += (X[n * D + d] * X[n * D + d]);
        }
    }
    
    double datasum;
    for(int n = 0; n < N; n++) {
        for(int m=n+1; m<N; m++){
            DD[n * N + m] = dataSums[n] + dataSums[m];
            
            datasum = 0;
        		for(int d=0; d<D; d++){
        			datasum += X[n * D + d] * X[m * D + d];
        		}
        		datasum = datasum*2;
        		
        		DD[n * N + m] -= datasum;
        		DD[m * N + n] = DD[n * N + m];
        }
        DD[n * N + n] = 0;
    }
    /*
    double datasum;
    for(int n=0; n<N; n++){
    	for(int m=n+1; m<N; m++){
				datasum = 0;
    		for(int d=0; d<D; d++){
    			datasum += X[n * D + d] * X[m * D + d];
    		}
    		datasum = datasum*2;
    		DD[n * N + m] -= datasum;
    		DD[m * N + n] -= datasum;
    	}
    	DD[n * N + n] = 0;
    }
    */
    //mxFree(dataSums); dataSums = NULL;
}


// Makes data zero-mean
void zeroMean(double* X, int N, int D, double* mean) {
	
	// Compute data mean
	//double* mean = (double*) mxCalloc(D, sizeof(double));
  //if(mean == NULL) { mexErrMsgTxt("Memory allocation failed! mean\n"); }
  
  //initialize mean
  for(int d = 0; d < D; d++) {
    mean[d] = 0;
  }
  
	for(int n = 0; n < N; n++) {
		for(int d = 0; d < D; d++) {
			mean[d] += X[n * D + d];
		}
	}
	for(int d = 0; d < D; d++) {
		mean[d] /= (double) N;
	}
	
	// Subtract data mean
	for(int n = 0; n < N; n++) {
		for(int d = 0; d < D; d++) {
			X[n * D + d] -= mean[d];
		}
	}
  //mxFree(mean); mean = NULL;
}


// Generates a Gaussian random number
double randn() {
	double x, y, radius;
	do {
		x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
		y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
		radius = (x * x) + (y * y);
	} while((radius >= 1.0) || (radius == 0.0));
	radius = sqrt(-2 * log(radius) / radius);
	x *= radius;
	y *= radius;
	return x;
}

//transform indexing from matlab to c
void Mat2C(double * X, int M, int N)
{
  double* Y = (double*) mxMalloc(M * N * sizeof(double));
  int m,n;
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      Y[m * N + n] = X[n * M + m];
    }
  }
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      X[m * N + n] = Y[m * N + n];
    }
  }
  mxFree(Y);
}

//transform indexing from c to matlab
void C2Mat(double * X, int M, int N)
{
  double* Y = (double*) mxMalloc(M * N * sizeof(double));
  int m,n;
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      Y[n * M + m] = X[m * N + n];
    }
  }
  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      X[n * M + m] = Y[n * M + m];
    }
  }
  mxFree(Y);
}

//Matlab stores matrix in columnwise
//C stores matrix in rowwise
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Macros for the ouput and input arguments */
	#define YDATA plhs[0]
	#define P prhs[0]
	#define NO_DIMS prhs[2]
	
	double *ydata, *p;
	double no_dims_in;
	
	int nrow, ncol, origN;
	int no_dims;
	
	nrow = mxGetM(P);
	ncol = mxGetN(P);
	p = mxGetPr(P);
	no_dims_in = (double)mxGetScalar(NO_DIMS);	
	no_dims = (int)floor(no_dims_in);
	
	YDATA = mxCreateDoubleMatrix(nrow, no_dims, mxREAL);
  
	if(YDATA==NULL){
		mexErrMsgTxt("cannot allocate memory for output\n");
	}
	
	ydata = mxGetPr(YDATA);

	//translate to tsne convention
	origN = nrow;
	srand(time(NULL));
	
	Mat2C(p, nrow, ncol);
	tsne_run(p, ydata, origN, no_dims);
	C2Mat(ydata, nrow, no_dims);
}