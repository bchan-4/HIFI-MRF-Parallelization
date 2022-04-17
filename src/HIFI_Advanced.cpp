#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <list>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include "HIFI_options.h"
#include "HIFI_misc.h"
#include "HIFI_MyMatrix.h"
#include "mpi.h"

using namespace std;
//using namespace __gnu_cxx;
 
HIFI_options options;
HIFI_matrixProperties matrixProperties;

UpperDiag<float> biasMatrix;
float *bias;
double *precomputedNormPDF;
double *precomputedNormLogPDF;

#define NORM_RES 10000
int fullRange=100*NORM_RES;
int floatFullRange=2*100*NORM_RES;


// retrieve normpdf value with mean m and sigma s
float normpdf(float x, float m, float s) {
  int i=((x-m)/s)*NORM_RES + fullRange;
  if (i<0 || i>floatFullRange) return 0;
  return precomputedNormPDF[i];
}

void outputSparseBoundaryMatrix(char *fn, UpperDiag<char> &horizontalBoundary );

// for boundary identification. These boundaries set "domains" within the matrix across which elements cannot be part of the same neighborhood.
// Iterates through columns and rows and searches for where distributions of entries between OoverBias(i,j:j') and OoverBias(i+1,j:j') are significantly different, based on Kolmogorov-Smirnov (KS) test.
#define MAX_NBNONMAX 10    
void setBoundaries(UpperDiag<float> &OoverBias, UpperDiag<char> &horizontalBoundary) {
  
  float sqrtLeftOver2[100000];
  for (int i=0;i<100000;i++) sqrtLeftOver2[i]=sqrt(i/2.0);

  // horizontal
  for (int i=OoverBias.startRow();i<OoverBias.endRow()-1;i++) {

    // consider a boundary starting down from (i,j)
    int distLeft[10000];
    int distRight[10000];
    int maxj=-1;
    int maxj2=-1;
    float maxks=-1;
    
    for (int j=OoverBias.startCol(i)+1;j<OoverBias.endCol(i);j++) {
      if (OoverBias.get(i,j)==0 && j+1>= OoverBias.startCol(i) && j+1 < OoverBias.endCol(i) && OoverBias.get(i+1,j)==0) continue;
      
      memset(distLeft,0,10000*sizeof(int));
      memset(distRight,0,10000*sizeof(int));
      int nLeft=0;
      int nRight=0;
      int nbNonMax=0;
      float maxSearch=0;
      int maxEntry=0;

      int j2=0;


      for (j2=j;j2<OoverBias.endCol(i) && nbNonMax<MAX_NBNONMAX;j2++) { // end point of boundary
	      nbNonMax++;
	      // only consider end points where there is data
	      if (horizontalBoundary.get(i,j2)) break;
	      int entryLeft=(int)(OoverBias.get(i,j2)+0.49);
	      int entryRight=(int)(OoverBias.get(i+1,j2)+0.49);

	      distLeft[entryLeft]++; nLeft++;
	      distRight[entryRight]++; nRight++;
	      if (entryLeft>maxEntry) maxEntry=entryLeft;
	      if (entryRight>maxEntry) maxEntry=entryRight;

	      if (entryLeft==0 && entryRight==0) continue;

	      int cumLeft2=0;
	      int cumRight2=0;
	      int maxDiff2=0;
	      for (int a=0;a<maxEntry;a++) {
	        cumLeft2+=distLeft[a];
	        cumRight2+=distRight[a];
	        if (abs(cumLeft2-cumRight2)>maxDiff2) { maxDiff2=abs(cumLeft2-cumRight2);}
	      }

	      float ks=((float)maxDiff2/nLeft)*sqrtLeftOver2[nLeft];

	      if (ks>maxSearch) {nbNonMax=0;maxSearch=ks;}

	      if (ks>maxks) {

	      maxks=ks;
	      maxj=j;
	      maxj2=j2;
	      }
	      else {nbNonMax++;}
        }
      }

    if (maxks>options.boundaryKS) {
      for (int j=maxj;j<=maxj2;j++) {
	      horizontalBoundary.set(i,j,1);
      }
      for (int x=i;x<=maxj;x++) {
	      horizontalBoundary.set(x,maxj2,1);
      }
	    //	verticalBoundary.set(j,i,1); 
      //Note from Brian: the line above was commented out in the original source code. The authors thought the horizontal boundaries and 
      // vertical boundaries should be symmetric. However, based on how they are setting the boundaries,
      // They really should have just set vertical boundaries as being *perpendicular* to horizontal boundaries. Biologically, there is no reason for the boundaries to be at different places.
      // This would allow us to calculate boundaries per task instead of all at the root, since the matrices are held by rows and we do not
      // need an entire column to decide vertical boundaries. Also, this means we do NOT have to hold a verticalBoundary matrix in memory.
      // We can just use horizontalBoundary for BOTH (really should name it boundarymatrix, not horizontalBoundary)
      
      i--;
    }
  }

  // NOTE FROM BRIAN: There is no reason for the horizontal boundaries to be different from the vertical boundaries biologically. Do not
  // need to calculate twice. See comment above. Vertical boundary calculation below is deleted.
  
  if (options.boundaryOutput[0]) {
    outputSparseBoundaryMatrix(options.boundaryOutput, horizontalBoundary);
  }

}

// Retrieve bias value from either bias matrix if supplied or just uses identity matrix.
float getBias(int i, int j) {
  if (options.normalizationMatrix[0]) {
    return biasMatrix.get(i,j);
  }
  return bias[i]*bias[j];
}

// precomputing OoverBias matrix 
void computeOoverBias(UpperDiag<float> &O, UpperDiag<float> &OoverBias) {                                                                                                     

  for (int i=O.startRow();i<O.endRow();i++) {
    for (int j=O.startCol(i);j<O.endCol(i);j++) {
      int x=O.get(i,j);
      if (x>0) OoverBias.set(i,j,x/getBias(i,j));
    }
  }
}
	   
// adaptive kernel density estimation.
void computeTMatrix_kde(UpperDiag<float> &O, UpperDiag<float> &OoverBias, UpperDiag<float> &T) {
  
  int nbStdDev=3;
  UpperDiag<unsigned long> cumO(options.firstRow,options.firstCol, options.lastRow,options.lastCol, options.bandSize+2*nbStdDev*options.kdeMaxBandwidth+10);
  UpperDiag<int> nextNonEmpty(options.firstRow,options.firstCol, options.lastRow, options.lastCol, options.bandSize+2*nbStdDev*options.kdeMaxBandwidth+10);

  // build sparse matrix representation
  for (int i=O.startRow();i<O.endRow();i++) {
    int last=O.startCol(i);
    for (int j=O.startCol(i);j<O.endCol(i);j++) {
      nextNonEmpty.set(i,j,-1);
      if (O.get(i,j)!=0) {
	for (int jj=last;jj<j;jj++) {
	  nextNonEmpty.set(i,jj,j);
	}
	last=j;
      }
    }
  }

  //   initialize T matrix
  int sumObs=0;
  
  for (int i=O.startRow();i<O.endRow();i++) {
    for (int j=O.startCol(i);j<O.endCol(i);j++) {
      T.set(i,j,-1);
      int t=O.get(i,j);
      sumObs+=t;
    }
  }

  // setting up bandwidth range
  float mini_h=0.3;
  float maxi_h=options.kdeMaxBandwidth; //default is 25

  if (options.method_kde) {
    mini_h=options.kdeBandwidth;
    maxi_h=mini_h+0.0001;
  }

  


  bool lastIter=false;

  // consider increasing values of bandwidth h
  double multFactor=1.3;
  for (float h=mini_h;h<maxi_h;h*=multFactor) {

    if (h*multFactor>=maxi_h) lastIter=true;

    int nbStdDevTh=(int)(nbStdDev*h);
    
    // compute cumulative.
    if (h==mini_h) {
      
      for (int i=cumO.startRow();i<cumO.endRow();i++) {
	for (int j=cumO.startCol(i);j<cumO.endCol(i);j++) {
	  unsigned long above=(i>O.startRow() && j<cumO.endCol(i-1))? cumO.get(i-1,j) : 0;
	  
	  unsigned long left=0;
	  if (j>O.startCol(i)) left=cumO.get(i,j-1);
	  else {
	    if (i>O.startRow() && j-1>=O.startCol(i-1)) left=cumO.get(i-1,j-1);
	  }

	  unsigned long aboveleft=0;
	  if (i-1>=O.startRow() && j-1>=O.startCol(i-1)) aboveleft=cumO.get(i-1,j-1);
	  
	  float x=0;
	  if (i>=O.startRow() && i<O.endRow() && j>=O.startCol(i) && j<O.endCol(i))  x= O.get(i,j);
	  cumO.set(i,j,above+left-aboveleft+x);
	}
      }
    }
    
    // 2D gaussian distribution. See Methods section of paper.
    //    fprintf(stderr,"Computing gaussian\n");
    float** gaussian=(float**) malloc(sizeof(float*)*(int)(nbStdDev*2*(maxi_h+1)+1));
    for (int i=0;i<(int)(nbStdDev*2*(maxi_h+1)+1);i++) gaussian[i]=(float*) malloc(sizeof(float)*(int)((nbStdDev*2*(maxi_h+1)+1)));
    
    float s=0;
    for (int i=-nbStdDevTh;i<=nbStdDevTh;i++) {
      for (int j=-nbStdDevTh;j<=nbStdDevTh;j++) {
	gaussian[i+nbStdDevTh][j+nbStdDevTh]=0;
	normpdf(sqrt(i*i+j*j),0,h);
	gaussian[i+nbStdDevTh][j+nbStdDevTh]=normpdf(sqrt(i*i+j*j),0,h);
	s+=gaussian[i+nbStdDevTh][j+nbStdDevTh];
      }
    }
    
    for (int i=-nbStdDevTh;i<=nbStdDevTh;i++) {
      for (int j=-nbStdDevTh;j<=nbStdDevTh;j++) {
	gaussian[i+nbStdDevTh][j+nbStdDevTh]/=s;
      }
    }
    
    int nSet=0;    
    float maxRadius2=nbStdDevTh*nbStdDevTh;
    float maxRadius=nbStdDevTh;
    
    int *minb=(int*)malloc(sizeof(int)*(maxRadius*2+3));
    for (int a=-maxRadius;a<=maxRadius;a++) {
      minb[(int)(a+maxRadius)]=(int)(sqrt(maxRadius2-a*a));
    }

    double gaussianSum=0;
    for (int a=-maxRadius;a<=maxRadius;a++) {
      int range=minb[(int)(a+maxRadius)];
      for (int b=-range;b<=range;b++) {
	gaussianSum+=gaussian[a+nbStdDevTh][b+nbStdDevTh];
      }
    }

    // iterate over matrix to apply gaussian filter
    for (int i=O.startRow();i<O.endRow();i++) {

      for (int j=O.startCol(i);j<O.endCol(i);j++) {
	if (T.get(i,j)<0) {
	  
	  // check if there are enough entries in the square around i,j
	  int lowi= (i-nbStdDevTh)>=O.startRow() ? (i-nbStdDevTh):O.startRow();
	  int highi= (i+nbStdDevTh)<O.endRow()? (i+nbStdDevTh):O.endRow()-1;
	  
	  int lowj= (j-nbStdDevTh)>=O.startCol(lowi)? (j-nbStdDevTh):O.startCol(lowi);
	  int highj= (j+nbStdDevTh)<O.endCol(lowi)? (j+nbStdDevTh):O.endCol(lowi)-1;
	  
	  int nHits=0;
	  if (lowi>O.startRow() && lowj>O.startCol(O.startRow())) {
	    nHits=cumO.get(highi,highj)-cumO.get(lowi-1,highj);
	    if (highi<=lowj-1) {nHits-=cumO.get(highi,lowj-1);}
	    else {nHits-=cumO.get(lowj-1,lowj-1);}
	    
	    if (lowi-1<=lowj-1) {nHits+=cumO.get(lowi-1,lowj-1);}
	    else {nHits+=cumO.get(lowj-1,lowj-1);}
	  }
	  
	  
	  if (lowi>O.startRow() && lowj==O.startCol(O.startRow())) nHits=cumO.get(highi,highj)-cumO.get(lowi-1,highj);
	  if (lowi==O.startRow() && lowj>O.startCol(O.startRow())) {
	    nHits=cumO.get(highi,highj);
	    if (highi<=lowj-1) nHits-=cumO.get(highi,lowj-1);
	    else nHits-=cumO.get(lowj-1,lowj-1);
	  }
	  if (lowi==O.startRow() && lowj==O.startCol(O.startRow())) nHits=cumO.get(highi,highj);
	  
	  if ((lastIter && nHits>0) || nHits>=options.kdeMinCount) {
	    
	      double s=0;
	      int c=0;
	      double sumw=0;
	      double gs=gaussianSum;
	      bool isSafe;
	      // check if entire square around i,j is away from diagonal and borders
	      if (i-maxRadius>=O.startRow() && i+maxRadius<O.endRow() && j+maxRadius<O.endCol(i+maxRadius) && j-maxRadius>=O.startCol(i-maxRadius) && i+maxRadius < j-maxRadius) isSafe=true;
	      else isSafe=false;

	      if (isSafe) {
		int mr=(int) maxRadius;
		for (int ii=i-mr;ii<=i+mr;ii++) {
		  int range=minb[ii-i+mr];
		  int jj=j-range;
		  if (jj>=O.endCol(ii) || O.get(ii,jj)==0) 
		    jj=nextNonEmpty.get(ii,jj);

		  while (jj!=-1 && jj<=j+mr) {
		    if (ii==jj) continue; // skip items on main diagonal
		    if (jj-j<=range) {
		      s+=OoverBias.get(ii,jj)*gaussian[ii-i+nbStdDevTh][jj-j+nbStdDevTh];
		    }
		    jj=nextNonEmpty.get(ii,jj);
		  }
		}
	      }
	      else {
		for (int a=-maxRadius;a<=maxRadius;a++) {
		  int ii=i+a;
		  int range=minb[(int)(a+maxRadius)];
		  for (int b=-range;b<=range;b++) {
		    int jj=j+b;
		    if (ii!=jj && (ii>=O.startRow() && ii<O.endRow() && jj>=ii && jj>=O.startCol(ii) && jj<O.endCol(ii))) {
		      float x=OoverBias.get(ii,jj);
		      if (x>0) {
			s+=x*gaussian[a+nbStdDevTh][b+nbStdDevTh];
		      }
		    }
		    else {gs-=gaussian[a+nbStdDevTh][b+nbStdDevTh];}
		  }
		}
	      }
	      T.set(i,j,s/gs);

	      nSet++;
	    }
	} else {nSet++;}
      }
    }
    
    free(minb);
    fprintf(stderr,"KDE: h = %lf, nSet = %d\n",h,nSet);
    fflush(stderr);    
    for (int i=0;i<(int)((nbStdDev*2*(maxi_h+1)+1));i++) free(gaussian[i]);
    free(gaussian);
    
  }
  

  for (int i=O.startRow();i<O.endRow();i++) {
    for (int j=O.startCol(i);j<O.endCol(i);j++) {
      if (T.get(i,j)==-1) T.set(i,j,0);
    }
  }

  // clean up
  cumO.clear();
  nextNonEmpty.clear();
}


// returns the median of the neighborhood of i,j. Ignore elements that cross the boundaries.
float getMedianNeighbor(UpperDiag<float> &T, UpperDiag<char> &horizontalBoundary, int i, int j) {
 
  int d=9;
  int nNei=0;
  float sumNeighbor=0;
  float allNei[8];
  for (int off_i=-1;off_i<=1;off_i++) {
    for (int off_j=-1;off_j<=1;off_j++) {
      if (off_i==0 && off_j==0) continue;
      
      int  ii=i+off_i;
      int jj=j+off_j;
      if (ii<jj && ii>=T.startRow() && jj>=T.startCol(ii) & ii<T.endRow() && jj<T.endCol(ii)) { 
	
 	if (horizontalBoundary.get(i,j-1) && off_j==-1) continue;
	if (horizontalBoundary.get(i,j) && off_j==+1) continue;
	if ((i==T.startRow() || horizontalBoundary.get(i-1,j)) && off_i==-1) continue;
	if (horizontalBoundary.get(i,j) && off_i==+1) continue;
	
	
	allNei[nNei]=T.get(ii,jj)+(float)rand()/RAND_MAX*0.0001;
	sumNeighbor+=T.get(ii,jj);
	nNei++;
      }
    }
  }
  // calculate median
  float med1=-1;
  float med2=-1;
  float median;
  for (int i=0;i<nNei;i++) {
    int nSmEq=0;
    int nSm=0;
    for (int j=0;j<nNei;j++) {
      if (allNei[j]<=allNei[i]) nSmEq++;
      if (allNei[j]<allNei[i]) nSm++;
    }
    if (nNei%2==1 && nSmEq>=nNei/2+1 && nSm<nNei/2+1) med1=med2=allNei[i];
    if (nNei%2==0 && nSmEq>=nNei/2 && nSm<nNei/2) med1=allNei[i];
    if (nNei%2==0 && nSmEq>=nNei/2+1 && nSm<nNei/2+1) med2=allNei[i];
  }
  median=(med1+med2)/2;
  if (nNei>0) return median;
  else return 0;
}


float *logKFact;

// Evaluates log-likelihood that the observed read count obs occured if the kog of the sparse, optimized matrix entry t is normally distributed with 
// mean=log(median of neighborhood) and std=log(median of neighborhood)*localSigma (given in HIFI_options.h).
// Also assumes that observed read count is poisson distributed with lambda parameter equal to t (the sparse, optimized matrix entry) times to bias factors (standard calculation).
// See Eq (1) and following eq. in "Markov random field estimation" section of the Methods in the paper. 
double evaluateLikelihood(float currentT, int obs, float medNei, float bias) {

  float transformedMedNeighbor;

  transformedMedNeighbor=log(medNei+1);

  // no +1 here?
  float transformedCurrentTbias;
  transformedCurrentTbias=log((currentT+0.00001)*bias);
 
  float localSigma=max(options.minSigma, sqrt(medNei*options.sigmaMultiplier));
  float NORM_RES_over_sigma=NORM_RES/localSigma;
  
  float p=currentT;
  float correctedp=currentT*bias;
  
  // log of Poisson distribution. logKFact is an array prepopulated with the natural logarithm of the k factorial (needed for distribution). Precomputed in "precomputeStuff"
  float loglike1= -correctedp + obs  * (transformedCurrentTbias) - logKFact[obs];
  
  float transformedp;
  transformedp=log(p+1);

  float loglike2=0; 
  int iii=(transformedp-transformedMedNeighbor)*NORM_RES_over_sigma + fullRange;
  if (iii>=0 && iii<=floatFullRange) loglike2= precomputedNormLogPDF[iii];
  return loglike1+loglike2;
}

// Finds log-likelihood of having observed read matrix O given the sparse, optimized matrix T by iterating through each elements and going through the function above.
double evaluateFullLikelihood(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> &horizontalBoundaries) {

  double like=0;
  
  for (int i=O.startRow();i<O.endRow();i++) {
    for (int j=O.startCol(i)+1;j<O.endCol(i);j++) {
       
      float medianNeighbor=getMedianNeighbor(T, horizontalBoundaries, i,j);
      float t=T.get(i,j);
      int o=O.get(i,j);
      double l=evaluateLikelihood(t, o, medianNeighbor,getBias(i,j));
      like+=l;
    }
  }
  return like;
}

// Same as above, but only evaluate between row i1 and i2 inclusive
double evaluateFullLikelihoodSection(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> &horizontalBoundaries,int i1, int i2) {

  double like=0;
  
  for (int i=i1;i<=i2;i++) {
    for (int j=O.startCol(i)+1;j<O.endCol(i);j++) {
       
      float medianNeighbor=getMedianNeighbor(T, horizontalBoundaries, i,j);
      float t=T.get(i,j);
      int o=O.get(i,j);
      double l=evaluateLikelihood(t, o, medianNeighbor,getBias(i,j));
      like+=l;
    }
  }
  return like;
}

// log-likelihood of neighborhood, as above function. neighborhood of ii-1:ii+1 inclusive (unless hitting boundary) and jj-1:jj+1 inclusing (unless hitting boundary).
double evaluateFullLikelihoodPatch(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> &horizontalBoundaries, int ii, int jj) {
  double like=0;
  
  for (int i=max(O.startRow(),ii-1);i<min(O.endRow(),ii+2);i++) {
    for (int j=max(O.startCol(i)+1,jj-1);j<min(O.endCol(i),jj+2);j++) {
       
      float medianNeighbor=getMedianNeighbor(T, horizontalBoundaries, i,j);
      float t=T.get(i,j);
      int o=O.get(i,j);
      double l=evaluateLikelihood(t, o, medianNeighbor,getBias(i,j));
      like+=l;
    }
  }
  return like;
}


#define NFACTORS 13
float factor[NFACTORS]={0.5, 0.7, 0.9, 0.95, 0.97, 0.99, 1.0 ,1.01,1.03,1.05,1.1,1.3,1.5};
float logFactor[NFACTORS];

// returns optimal choice of value for an entry T(i,j), given the number of reads observed (obs) and the median of the neighborhood.
// This section goes through an iterated conditional mode (ICM) algorithm, using the "evaluateFullLikelihoodPatch" as the optimized function and the multiplicative factors above to modify the "guesses" for T(i,j)
// Note that each iteration through the entire matrix goes through ALL of the factors above.
double optimizeTnew(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> horizontalBoundaries,  int i, int j, int obs, float medNei, float localbias) {
   
  // search for the best multiplicative factor to apply to currenT
  float loglike[NFACTORS];
  float old=T.get(i,j);
  for (int f=0;f<NFACTORS;f++) {
    T.set(i,j,old*factor[f]);
    loglike[f]=evaluateFullLikelihoodPatch(O, T, horizontalBoundaries,  i,j);
  }
  T.set(i,j,old);
  // select best factor
  int chosenf=0;
  float bestLogLike=loglike[0];

  for (int f=1;f<NFACTORS;f++) {
    if (loglike[f]>bestLogLike) {
      bestLogLike=loglike[f];
      chosenf=f;
    }
  }

  return old*factor[chosenf];
}

  
// Precomputed distributions and factorials
void precomputeStuff() {

  // log k!
  logKFact=(float*)malloc(1000000*sizeof(float));

  logKFact[0]=0;
  for (int i=1;i<1000000;i++) {
    logKFact[i]=logKFact[i-1]+log(i);
  }


  // normal distribution PDF
  precomputedNormPDF=(double*) malloc((2000000+1)*sizeof(double));
  precomputedNormLogPDF=(double*) malloc((2000000+1)*sizeof(double));
  double s=0;
  float denom=sqrt(6.28318530717959);
  for (int i=-NORM_RES*100;i<=NORM_RES*100;i++) {

    precomputedNormPDF[i+NORM_RES*100] = exp(-((float)i/NORM_RES*(float)i/NORM_RES)/2)/denom;
    precomputedNormLogPDF[i+NORM_RES*100] = -((float)i/NORM_RES*(float)i/NORM_RES)/2-log(denom);
    
    s+=precomputedNormPDF[i+NORM_RES*100];
  }
  
  for (int i=-NORM_RES*100;i<=NORM_RES*100;i++) {
    precomputedNormPDF[i+NORM_RES*100]/=s;
    precomputedNormLogPDF[i+NORM_RES*100]-=log(s);
  }

  // log of multiplicative factors
  for (int i=0;i<NFACTORS;i++) logFactor[i]=log(factor[i]);
}


// Self explanatory.
void getMatrixSize(char *fn) {
  matrixProperties.fullMatrix_firstRow=matrixProperties.fullMatrix_firstCol=999999999;
  matrixProperties.fullMatrix_lastRow=matrixProperties.fullMatrix_lastCol=-1;
  FILE *f=fopen(fn,"r");
  if (!f) {
    fprintf(stderr,"Error: Can't open input file %s\n",fn);
    exit(1);
  }
  int x=fscanf(f,"#%d %d %d %d\n", &matrixProperties.fullMatrix_firstRow, &matrixProperties.fullMatrix_lastRow, &matrixProperties.fullMatrix_firstCol, &matrixProperties.fullMatrix_lastCol);
  matrixProperties.fullMatrix_size=max(matrixProperties.fullMatrix_lastRow,matrixProperties.fullMatrix_lastCol)+1;
  
  if (options.firstRow==-1) {
    options.firstRow=matrixProperties.fullMatrix_firstRow;
    options.firstCol=matrixProperties.fullMatrix_firstCol;
    options.lastRow=matrixProperties.fullMatrix_lastRow;
    options.lastCol=matrixProperties.fullMatrix_lastCol;
  }
  
  fclose(f);
}

// for bias calculation
#define BIAS_PSEUDOCOUNT 5
#define MAX_BIAS 10

// Function to read in .tsv file and populate the observed read count matrix O. Will also set the bias list if setBias=true.
void readFullOMatrix(char *fn, UpperDiag<float> *mat, float *bis, bool setBias) {
  
  FILE *f=fopen(fn,"r");  
  if (!f) {
    fprintf(stderr,"Error: Can't open input file %s\n",fn);
    exit(1);
  }
  
  char line[1000];

  double fullSum=0;
  double *sumRow=(double*)malloc(matrixProperties.fullMatrix_size*sizeof(double));
  for (int i=matrixProperties.fullMatrix_firstRow;i<=matrixProperties.fullMatrix_lastRow;i++) {
    sumRow[i]=BIAS_PSEUDOCOUNT;
    fullSum+=BIAS_PSEUDOCOUNT;
  }
  
  while (fscanf(f,"%[^\n]\n",line)!=EOF) {
    if (line[0]=='#' || line[0]=='0') continue;
    int foo;
    int p1,p2;
    float c;
    char c1[100];
    char c2[100];
    sscanf(line,"%f %s %d %s %d",&c,c1,&p1,c2,&p2);
    
    if (matrixProperties.chr1[0]==0) {strcpy(matrixProperties.chr1,c1);}
    else { if (strcmp(c1,matrixProperties.chr1)) {fprintf(stderr,"ERROR1: input file contains data from multiple chromosomes: %s and %s\n",matrixProperties.chr1, c1);exit(1);}}
    if (matrixProperties.chr2[0]==0) {strcpy(matrixProperties.chr2,c2);}
    else { if (strcmp(c2,matrixProperties.chr2)) {fprintf(stderr,"ERROR2: input file contains data from multiple chromosomes: %s and %s\n",matrixProperties.chr2, c2);exit(1);}}
    
    
    
    sumRow[p1]+=c;
    sumRow[p2]+=c;
    fullSum+=c+c;
    if (abs(p2-p1)<options.bandSize) {
      if (options.lastRow==-1 || (p1>=options.firstRow && p1<=options.lastRow && p2>=options.firstCol && p2<=options.lastCol))
	{
	  if (p1<options.firstRow || p1>options.lastRow) continue;
	  if (p2<options.firstCol || p1>options.lastCol) continue;
	  mat->set(p1,p2,c);
      //mat.set(p1,p2,c);
	}
    }
  }

  fclose(f);
  
  // Calculate bias matrix
  if (setBias) {
    for (int i=matrixProperties.fullMatrix_firstRow;i<matrixProperties.fullMatrix_lastRow;i++) {
      bis[i]=sumRow[i]/((fullSum)/(matrixProperties.fullMatrix_lastRow-matrixProperties.fullMatrix_firstRow+1));
      if (bis[i]>MAX_BIAS) bis[i]=MAX_BIAS;
      if (bis[i]<1.0/MAX_BIAS) bis[i]=1.0/MAX_BIAS;
    }
  }
}


// If a bias matrix is supplied. Will not use this for project.
void readBiasMatrix(char *fn, UpperDiag<float> *mat) {
  
  FILE *f=fopen(fn,"r");  
  if (!f) {
    fprintf(stderr,"Error: Can't open input file %s\n",fn);
    exit(1);
  }
  
  char line[1000];

  // initialize bias to 1
  for (int i=mat->startRow();i<mat->endRow();i++) {
    for (int j=mat->startCol(i)+1;j<mat->endCol(i);j++) {
      mat->set(i,j,1);
    }
  }

  
  while (fscanf(f,"%[^\n]\n",line)!=EOF) {
    if (line[0]=='#') continue;
    int foo;
    int p1,p2;
    float c;
    char c1[100];
    char c2[100];
    sscanf(line,"%f %s %d %s %d",&c,c1,&p1,c2,&p2);
    if (strcmp(c1,matrixProperties.chr1)) {fprintf(stderr,"ERROR1: Chromosome mismatch: %s and %s\n",matrixProperties.chr1, c1);exit(1);}
    if (strcmp(c2,matrixProperties.chr2)) {fprintf(stderr,"ERROR2: Chromosome mismatch: %s and %s\n",matrixProperties.chr2, c2);exit(1);}
    
    if (abs(p2-p1)<options.bandSize) {
      if (options.lastRow==-1 || (p1>=options.firstRow && p1<=options.lastRow && p2>=options.firstCol && p2<=options.lastCol))
	{
	  if (p1<options.firstRow || p1>options.lastRow) continue;
	  if (p2<options.firstCol || p1>options.lastCol) continue;
	  mat->set(p1,p2,c);
	}
    }
  }

  fclose(f);
}

// Fixed binning. Not using this.
void computeTMatrix_fixed(UpperDiag<float> &O, UpperDiag<float> &T, float *bias) {
  float *sumRow=(float*)malloc(O.endRow()*sizeof(float));
  float fullSum=0;
  for (int i=O.startRow();i<O.endRow();i++) {
    sumRow[i]=0;
  }
  
  for (int i=O.startRow();i<O.endRow();i+=options.fragmentsPerBin) {
    for (int j=O.startCol(i);j<O.endCol(i);j+=options.fragmentsPerBin) {
      float sum=0;
      int n=0;
      for (int a=0;a<options.fragmentsPerBin && i+a<O.endRow();a++) {
	for (int b=0;b<options.fragmentsPerBin && j+b<O.endCol(i);b++) {
	  if (i+a<=j+b) {
	    sum+=O.get(i+a,j+b);
	    n++;
	  }
	}
      }
      float x=sum/n;
      for (int a=0;a<options.fragmentsPerBin && i+a<O.endRow();a++) {
        for (int b=0;b<options.fragmentsPerBin && j+b<O.endCol(i);b++) {
	  if (i+a<=j+b) {
	    T.set(i+a,j+b,x);
	    sumRow[i+a]+=x;
	    sumRow[j+b]+=x;
	    fullSum+=x;
	  }
	}
      }
    }
  }
  // deal with biases
  for (int i=0;i<O.endRow();i++) {
    bias[i]=sumRow[i]/((fullSum)/matrixProperties.fullMatrix_size);
    if (bias[i]>MAX_BIAS) bias[i]=MAX_BIAS;
    if (bias[i]<1.0/MAX_BIAS) bias[i]=1.0/MAX_BIAS;
  }
  for (int i=O.startRow();i<O.endRow();i++) {
    for (int j=O.startCol(i);j<O.endCol(i);j++) {
      T.set(i,j,T.get(i,j)/getBias(i,j));
    }
  }
}

// Outputs the final tsv file with optimized IF matrix.
void outputSparseMatrix(char *fn, UpperDiag<float> &T, float *bias) {
  FILE *out=fopen(fn,"w");
  if (!out) {
    fprintf(stderr,"Error: Can't open input file %s\n",fn);
    exit(1);
  }
  if (options.lastRow!=-1 && options.lastCol!=-1) fprintf(out,"# %d %d %d %d\n",options.firstRow,options.lastRow,options.firstCol,	    options.lastCol);
  else fprintf(out,"# %d %d %d %d\n",options.firstRow, options.lastRow, options.firstCol, options.lastCol);

  for (int i=T.startRow();i<T.endRow();i++) {

    for (int j=T.startCol(i);j<T.endCol(i);j++){
      float x=T.get(i,j);
      if (!options.outputNormalized) x*=getBias(i,j);
      if (x>options.minOutput) fprintf(out,"%5.3lf %s %d %s %d\n",x/options.minOutput,matrixProperties.chr1,i,matrixProperties.chr2,j);
    }
  }
  fclose(out);
}


// Outputs the boundaries. See "setBoundaries" fcn
void outputSparseBoundaryMatrix(char *fn, UpperDiag<char> &horizontalBoundary) {
  FILE *out=fopen(fn,"w");
  if (options.lastRow!=-1 && options.lastCol!=-1) fprintf(out,"# %d %d %d %d\n",options.firstRow,options.lastRow,options.firstCol,	    options.lastCol);
  else fprintf(out,"# %d %d %d %d\n",options.firstRow, options.lastRow, options.firstCol, options.lastCol);

  for (int i=horizontalBoundary.startRow();i<horizontalBoundary.endRow();i++) {
    for (int j=horizontalBoundary.startCol(i)+1;j<horizontalBoundary.endCol(i)-1;j++){
      if (horizontalBoundary.get(i,j) || (i>horizontalBoundary.startRow() && horizontalBoundary.get(i-1,j)) ||  horizontalBoundary.get(i,j) || (j>horizontalBoundary.startCol(j) && horizontalBoundary.get(i,j-1))) fprintf(out,"1 %s %d %s %d\n",matrixProperties.chr1,i+options.firstRow,matrixProperties.chr2,j+options.firstCol);

    }
  }
  fclose(out);
}


// Main pipeline
int main(int argc, char *argv[]) {
  if (argc<3) {
    options.help();
    exit(1);
  } 

  options.read(argc, argv);
  options.print();
   
  int rank,size,provided;
  // Initialize MPI
  MPI_Init_thread(&argc, &argv,MPI_THREAD_SINGLE,&provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  fprintf(stderr,"Hello from the beginning, rank %d\n",rank);

  // Start timer
  double StartTime=MPI_Wtime();
  
  //Pointer for time file
  FILE *outfile_fp;

  precomputeStuff(); // Precomputes distributions, etc.

  getMatrixSize(argv[1]); // size of matrix provided. All tasks have this stored.

  // Partition the matrix rows such that number of ELEMENTS per rank is about equal (or as equal as possible s.t. rows are not broken up)
  // Don't want to just divide number of rows evenly because the "matrix" is upper triangular, so the tasks with smaller ranks would
  // have many more elements than the tasks with larger ranks if we just divide by rows.

  // Approximate number of elements in matrix (~N^2 / 2)
  int OverallNum=round((options.lastRow-options.firstRow+1)*(options.lastRow-options.firstRow+1)/2);
  // Approximate number of elements ideally handled by each rank
  int idealNumPerRank=round(OverallNum/size);
  // Number of rows
  int numRows=options.lastRow-options.firstRow+1;
  int rowStart,rowEnd,rowStartNoHalo,rowEndNoHalo;
  // Will hold all row starts/ends at root
  int allrowstart[size];
  int allrowend[size];
  int allrowstartnohalo[size];
  int allrowendnohalo[size];

  // Algorithm to find bounds. rowEnds are inclusive
  int cumulativesum[numRows];
  int tempsum=0;
  int cumdividedbynumperrank[numRows];
  for (int rr=0;rr<numRows;rr++){
      cumulativesum[rr]=tempsum+numRows-rr;   // Add numRows-rr elements for each row index rr
      cumdividedbynumperrank[rr]=round(cumulativesum[rr]/idealNumPerRank);
      tempsum=tempsum+numRows-rr; 
    }
  for (int rr=1;rr<numRows;rr++){
      if (cumdividedbynumperrank[rr]==(rank+1) && cumdividedbynumperrank[rr-1]==(rank)){rowEndNoHalo=rr;} 
      if (cumdividedbynumperrank[rr]==(rank) && cumdividedbynumperrank[rr-1]==(rank-1)){rowStartNoHalo=rr+1;} 
   }
  // Edge cases
  if (rank==0){rowStart=0;rowStartNoHalo=0;}
  if (rank==(size-1)){rowEndNoHalo=numRows-1;rowEnd=numRows-1;}

   if (rank>0){rowStart=rowStartNoHalo-2;}
   if (rank<size-1){rowEnd=rowEndNoHalo+2;};


   // Gather row boundaries
   MPI_Gather(&rowStart,1,MPI_INT,&allrowstart[0],1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&rowStartNoHalo,1,MPI_INT,&allrowstartnohalo[0],1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&rowEnd,1,MPI_INT,&allrowend[0],1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&rowEndNoHalo,1,MPI_INT,&allrowendnohalo[0],1,MPI_INT,0,MPI_COMM_WORLD);
   if (rank==0){fprintf(stderr,"Success gather of row boundaries\n");}
  // define neighbors
  int topneighRank=rank-1;
  int bottomneighRank=rank+1;

  // allocate tables. Only root should have full data.
  // For some reason need to initialize matrices on all tasks or else does not work?
  UpperDiag<float> T;
  UpperDiag<float> O;
  UpperDiag<float> OoverBias;
  UpperDiag<char> horizontalBoundaries;
  //UpperDiag<char> verticalBoundaries;
  bias=(float*) malloc(matrixProperties.fullMatrix_size*sizeof(float));

    if (rank==0){
        T=UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
        O=UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
        OoverBias=UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
        horizontalBoundaries=UpperDiag<char>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);

        for (int i=0;i<matrixProperties.fullMatrix_size;i++) bias[i]=1;
        fprintf(stderr,"Reading matrix\n"); // Reads in .tsv file with read count matrices. Also calculates the bias list (see "Fragment-specific bias calculation" in Methods of paper).
        readFullOMatrix(argv[1], &O, bias, true);
        fprintf(stderr,"After read, rank %d \n",rank);
        // Reads in bias matrix if supplied
        if (options.normalizationMatrix[0]) {
            biasMatrix = UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
            readBiasMatrix(options.normalizationMatrix,&biasMatrix);
        }
        fprintf(stderr,"About to compute OoverBias\n");
        computeOoverBias(O,OoverBias); // calculate bias-corrected read count matrix.
    }   // End read of data and bias computation
    
    MPI_Bcast(bias,matrixProperties.fullMatrix_size,MPI_FLOAT,0,MPI_COMM_WORLD);
    if (rank==0) {  fprintf(stderr,"Success broadcast of bias matrix\n");}

  
  // Allocate T and O for each task
  UpperDiag<float> Ttask = UpperDiag<float>(options.firstRow+rowStart,options.firstCol,options.firstRow+rowEnd,options.lastCol, options.bandSize);
  UpperDiag<float> Otask = UpperDiag<float>(options.firstRow+rowStart,options.firstCol,options.firstRow+rowEnd,options.lastCol, options.bandSize);
 
  // Calculate number of elements held by each task in individual Ttask,Otask, inclduing halos
  int numElements=0;
  int realcolstart=0;
  int realcolend=0;
  for (int whichrow=rowStart;whichrow<=rowEnd;whichrow++){ // for each row handled by the task
      realcolstart=max(whichrow+options.firstRow,options.firstCol); // which column to start at per row
      realcolend=min(options.firstRow+whichrow+options.bandSize,options.lastCol); // which column to end at per row
      numElements=numElements+(realcolend-realcolstart+1);
  }

  // Calculate number of elements held by each task NOT including halos
  int numElNH=0;
  int realcolstartNH=0;
  int realcolendNH=0;
  for (int whichrow=rowStartNoHalo;whichrow<=rowEndNoHalo;whichrow++){ // for each row handled by the task
      realcolstartNH=max(whichrow+options.firstRow,options.firstCol); // which column to start at per row
      realcolendNH=min(options.firstRow+whichrow+options.bandSize,options.lastCol); // which column to end at per row
      numElNH=numElNH+(realcolendNH-realcolstartNH+1);
  }

  // send number of elements per task to root
  int *numElemonEachTask=(int *)malloc(size*sizeof(int)); // each element has number of elements handled per task;
  int *numElEachTaskNH=(int *)malloc(size*sizeof(int));// no halos

  MPI_Gather(&numElements,1,MPI_INT, &numElemonEachTask[0],1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Gather(&numElNH,1,MPI_INT,&numElEachTaskNH[0],1,MPI_INT,0,MPI_COMM_WORLD);
  if (rank==0){
  fprintf(stderr,"Finish gather of number of elements at root rank %d\n",rank);}


  // allocate buffers for receiving data
  float *bufferO;
  bufferO=(float*) malloc(sizeof(float)*numElements);
  float *bufferT;
  bufferT=(float*) malloc(sizeof(float)*numElements);

  // Partition and distribute O data, akde guess for initial T, to all tasks
  MPI_Status status;
  float *bufferOsend;
  float *bufferTsend;
  MPI_Request req1[size];
  MPI_Request req2[size];
  MPI_Request req3;
  MPI_Request req4;
  if (rank==0){
      fprintf(stderr,"Starting KDE initialization\n");
      fflush(stderr);
      computeTMatrix_kde(O,OoverBias,T); // Calculate akde for starting guess at root. 
      fprintf(stderr,"Communicating O and T to other ranks\n");
      fflush(stderr);
      int colstart,colend;
        for (int k=0;k<size;k++) {
            // Assemble the data to send
            bufferOsend=(float*) malloc(sizeof(float)*numElemonEachTask[k]);
            bufferTsend=(float*) malloc(sizeof(float)*numElemonEachTask[k]);
            int buffcounter=0;
            for (int thisrow=allrowstart[k];thisrow<=allrowend[k];thisrow++) { // for the rows on each task
                  colstart=max(thisrow+options.firstRow,options.firstCol); // which column to start at per row
                  colend=min(options.firstRow+thisrow+options.bandSize,options.lastCol); // which column to end at per row
                for (int thiscol=colstart;thiscol<=colend;thiscol++){
                    bufferOsend[buffcounter]=O.get(thisrow+options.firstRow,thiscol); // Access the correct element
                    bufferTsend[buffcounter]=T.get(thisrow+options.firstRow,thiscol);
                    buffcounter++;
                }
            }
            // Send data
            MPI_Isend(&bufferOsend[0],numElemonEachTask[k],MPI_FLOAT,k,0,MPI_COMM_WORLD,&req1[k]);
            MPI_Isend(&bufferTsend[0],numElemonEachTask[k],MPI_FLOAT,k,1,MPI_COMM_WORLD,&req2[k]);
        }
   }
  double KDETime=MPI_Wtime()-StartTime; // Time through the AKDE algorithm
  if (rank==0) {
        // Root has to "receive" partitioned data from itself. The MPI_Wait()s are for receiving the partitioned data.
        // The MPI_Waitall()s are for sending the partitioned data.
        MPI_Irecv(&bufferO[0],numElements,MPI_FLOAT,0,0,MPI_COMM_WORLD,&req3);
        MPI_Irecv(&bufferT[0],numElements,MPI_FLOAT,0,1,MPI_COMM_WORLD,&req4);
        MPI_Wait(&req3,&status);MPI_Wait(&req4,&status);
        MPI_Waitall(size,req1,MPI_STATUSES_IGNORE);MPI_Waitall(size,req2,MPI_STATUSES_IGNORE);
    } else {
        MPI_Irecv(&bufferO[0],numElements,MPI_FLOAT,0,0,MPI_COMM_WORLD,&req3);
        MPI_Irecv(&bufferT[0],numElements,MPI_FLOAT,0,1,MPI_COMM_WORLD,&req4);
        MPI_Wait(&req3,&status);MPI_Wait(&req4,&status);
    }

    // These are the indices over which I actually want to calculate the MRF. endIterRow is NOT inclusive.
    // For some reason I needed to put these variables really early in the code.
    int startIterRow;
    startIterRow=rowStartNoHalo+options.firstRow;
    int endIterRow;
    endIterRow=rowEndNoHalo+options.firstRow+1;

  // Organize received data into the UpperDiag format
  int tempcolstart=0;
  int tempcolend=0;
  int communcounter=0;
  for (int ii=rowStart;ii<=rowEnd;ii++){ // for each row
      tempcolstart=max(ii+options.firstRow,options.firstCol); // which column to start at per row
      tempcolend=min(options.firstRow+ii+options.bandSize,options.lastCol); // which column to end at per row
      for (int jj=tempcolstart;jj<=tempcolend;jj++){
          Otask.set(ii+options.firstRow,jj,bufferO[communcounter]);
          Ttask.set(ii+options.firstRow,jj,bufferT[communcounter]);
          communcounter++; // iterate through elements
      }
    }
   printf("Finish populating initial matrices at rank %d\n",rank);

  // Put MRF into main. Not sure how to pass the objects otherwise...

  // Boundary part of mrf //
    
  // allocate boundary matrices per task
  UpperDiag<char> horzBoundtask=UpperDiag<char>(options.firstRow+rowStart,options.firstCol,options.firstRow+rowEnd,options.lastCol, options.bandSize);
  
  // allocate buffers for receiving boundary data
  char *bufferhorz;
  bufferhorz=(char*) malloc(sizeof(char)*numElements);


  // Calculate, partition and distribute horizontalBoundaries and verticalBoundaries to all tasks
  MPI_Request reqhorz[size];
  MPI_Request req5;
  if (rank==0){
        fprintf(stderr,"Finding boundaries at root\n");
        fflush(stderr);
        setBoundaries(OoverBias, horizontalBoundaries);
        fprintf(stderr,"Communicating boundaries from root\n");
        fflush(stderr);
        int colstart,colend;
        for (int k=0;k<size;k++) {
            // Assemble the data to send
            char *bufferhsend;
            bufferhsend=(char*) malloc(sizeof(char)*numElemonEachTask[k]);
            char *buffervsend;
            buffervsend=(char*) malloc(sizeof(char)*numElemonEachTask[k]);
            int buffcounter=0;
            for (int thisrow=allrowstart[k];thisrow<=allrowend[k];thisrow++) { // for the rows on each task
                  colstart=max(thisrow+options.firstRow,options.firstCol); // which column to start at per row
                  colend=min(options.firstRow+thisrow+options.bandSize,options.lastCol); // which column to end at per row
                for (int thiscol=colstart;thiscol<=colend;thiscol++){
                    bufferhsend[buffcounter]=horizontalBoundaries.get(thisrow+options.firstRow,thiscol); // Access the correct element
                    buffcounter++;
                }
            }
            // Send data
            MPI_Isend(&bufferhsend[0],numElemonEachTask[k],MPI_CHAR,k,3,MPI_COMM_WORLD,&reqhorz[k]);
        }
   }

      if (rank==0) { 
          // Root has to "receive" partitioned data from itself. The MPI_Wait()s are for receiving the partitioned data.
          // The MPI_Waitall()s are for sending the partitioned data.
         MPI_Irecv(&bufferhorz[0],numElements,MPI_CHAR,0,3,MPI_COMM_WORLD,&req5); 
         MPI_Wait(&req3,&status);
         MPI_Waitall(size,reqhorz,MPI_STATUSES_IGNORE);
       } else {
          // Other ranks only have to receive the partitioned data from root. Don't care about req1 and req2.
         MPI_Irecv(&bufferhorz[0],numElements,MPI_CHAR,0,3,MPI_COMM_WORLD,&req5);
         MPI_Wait(&req3,&status);
       }

 
  // Organize received boundary data into the UpperDiag format
  tempcolstart=0;tempcolend=0;communcounter=0;
  for (int ii=rowStart;ii<=rowEnd;ii++){ // for each row
      tempcolstart=max(ii+options.firstRow,options.firstCol); // which column to start at per row
      tempcolend=min(options.firstRow+ii+options.bandSize,options.lastCol); // which column to end at per row
      for (int jj=tempcolstart;jj<=tempcolend;jj++){
          horzBoundtask.set(ii+options.firstRow,jj,bufferhorz[communcounter]);
          communcounter++; // iterate through elements
      }
  }
  MPI_Barrier(MPI_COMM_WORLD); // Make sure everything is set to go

 // Holds temporary optimized T matrices
  UpperDiag<float> Ttemptask = UpperDiag<float>(options.firstRow+rowStart,options.firstCol,options.firstRow+rowEnd,options.lastCol, options.bandSize);

  // Indicates if/where each iteration has updated the T matrix
  UpperDiag<char> hasChangedtask=UpperDiag<char>(options.firstRow+rowStart,options.firstCol,options.firstRow+rowEnd,options.lastCol, options.bandSize);
  
  unsigned long nChanges=1;  
  
  // Buffer arrays for communication

    // Variables for nchanges, sumchanges reductions
    unsigned long nChangesTot;
    double sumChangesTot;
    // Variable for likelihoodtask reduction
    double LikelihoodTot;

    // Define indices for halos to receive
    // First two rows
    int topstart=max(rowStart+options.firstRow,options.firstCol); // which column to start at for first row
    int topend=min(options.firstRow+rowStart+options.bandSize,options.lastCol); // which column to end at for first row
    int toprowsize=topend-topstart+1; // number of elements in top row
    int nextstart=max(rowStart+1+options.firstRow,options.firstCol); // which column to start at for first row
    int nextend=min(options.firstRow+rowStart+1+options.bandSize,options.lastCol); // which column to end at for first row
    int nextrowsize=nextend-nextstart+1; // number of elements in top row
    // Last two rows
    int bottstart=max(rowEnd+options.firstRow,options.firstCol);// which column to start at for bottom row
    int bottend=min(options.firstRow+rowEnd+options.bandSize,options.lastCol); // which column to end at for bottom row
    int bottrowsize=bottend-bottstart+1; //number of elements in bottom row
    int bottnextstart=max(rowEnd-1+options.firstRow,options.firstCol);// which column to start at for bottom row
    int bottnextend=min(options.firstRow+rowEnd-1+options.bandSize,options.lastCol); // which column to end at for bottom row
    int bottnextsize=bottnextend-bottnextstart+1;

    // Indices for halos to send
    // third and fourth rows = last two rows of prev rank
    int topstartSend=max(rowStart+2+options.firstRow,options.firstCol); // which column to start at for first row
    int topendSend=min(options.firstRow+rowStart+2+options.bandSize,options.lastCol); // which column to end at for first row
    int toprowsizeSend=topendSend-topstartSend+1; // number of elements in top row
    int nextstartSend=max(rowStart+3+options.firstRow,options.firstCol); // which column to start at for first row
    int nextendSend=min(options.firstRow+rowStart+3+options.bandSize,options.lastCol); // which column to end at for first row
    int nextrowsizeSend=nextendSend-nextstartSend+1; // number of elements in top row
    // end-2 and end-3 rows = first two rows of next rank
    int bottstartSend=max(rowEnd-2+options.firstRow,options.firstCol);// which column to start at for bottom row
    int bottendSend=min(options.firstRow+rowEnd-2+options.bandSize,options.lastCol); // which column to end at for bottom row
    int bottrowsizeSend=bottendSend-bottstartSend+1; //number of elements in bottom row
    int bottnextstartSend=max(rowEnd-3+options.firstRow,options.firstCol);// which column to start at for bottom row
    int bottnextendSend=min(options.firstRow+rowEnd-3+options.bandSize,options.lastCol); // which column to end at for bottom row
    int bottnextsizeSend=bottnextendSend-bottnextstartSend+1;

    // Arrays for sending halos
    float *bufftoprows;
    bufftoprows=(float*) malloc(sizeof(float)*(toprowsizeSend+nextrowsizeSend));
    float *buffbottrows;
    buffbottrows=(float*) malloc(sizeof(float)*(bottrowsizeSend+bottnextsizeSend));
    // Arrays for recv halos
    float *buffrecvtop;
    buffrecvtop=(float*) malloc(sizeof(float)*(toprowsize+nextrowsize));
    float *buffrecvbott;
    buffrecvbott=(float*) malloc(sizeof(float)*(bottrowsize+bottnextsize));

    // Arrays for entire matrices without halos. 
    // Sending. To send to root after each iteration. All ranks must allocate memory.
    float *buffSendAll;
    buffSendAll=(float*) malloc(sizeof(float)*numElNH); // numElNH is per task
    // Receiving. For receiving all elements after each iteration. Only the root needs to allocate memory.
    float *buffCollectAll;
    int displs[size]; // Count up numElEachTaskNH array and displs array for MPI_Gatherv at root.
    if (rank==0){
        int totnumelementsNH=0; // Total number of elements without halos
        int prevd=0;
        for (int ttt=0;ttt<size;ttt++){
            totnumelementsNH=totnumelementsNH+numElEachTaskNH[ttt];
            displs[ttt]=prevd;
            prevd=prevd+numElEachTaskNH[ttt];
        }
        buffCollectAll=(float*) malloc(sizeof(float)*totnumelementsNH);
        fprintf(stderr,"TotalNumElements without halos is %d\n",totnumelementsNH); // Useful for debugging and troubleshooting.
    }
  // End allocation of buffer arrays

  // Time before main loop. Considers functions for reading data, AKDE algorithm, finding boundaries, distributing data
  double SetupTime=MPI_Wtime()-StartTime;

  // ** START OF MRF **//
  // Iterate through the Markov Random Field to update the optimized Ttask matrix. Default option is 5 iterations.
  // May change this to a tolerance for checking convergance, since I thought that was what was here...
  fprintf(stderr,"Finished allocating halo buffers. Starting MRF iterations, rank %d\n",rank);

  // Clear matrices that do not need to be accessed anymore
  if (rank==0){
  O.clear();
  T.clear();
  horizontalBoundaries.clear();
  OoverBias.clear();
  }

  for (int rep=0;rep<options.mrfMaxIter;rep++) {
      double sumChanges=0;

        /* See line ~1030 for how startIterRow and endIterRow are defined
        int startIterRow;
        startIterRow=rowStartNoHalo+options.firstRow;
        int endIterRow;
        endIterRow=rowEndNoHalo+options.firstRow+1;
        */

        for (int i=startIterRow;i<endIterRow;i++){ // Loop through elements
            for (int j=Otask.startCol(i)+1;j<Otask.endCol(i);j++) {
	
            // check if anything in the neighbourhood has changed. 
            // For some reason the code actually goes to 2 elements out for checking if something changed...
            // This is why halo is increased to two rows.
	
	        int neiLeft=max(startIterRow,i-2);
	        int neiRight=min(endIterRow-1,i+2);
	        int neiBottom=max(Otask.startCol(i),j-2);
	        int neiTop=min(Otask.endCol(i)-1,j+2);

	        bool changed=false;
	        for (int a=neiLeft;a<=neiRight;a++) {
	            for (int b=neiBottom;b<=neiTop;b++) {
	                if (a<=b) changed=changed || hasChangedtask.get(a,b);
	            }
	        }
      
	        float t=Ttask.get(i,j);	
	        Ttemptask.set(i,j,t);
      
            // If I am at the first iteration or something has changed in neighborhood in Ttask, do the optimization again.
	        if (rep==0 || changed) {
	            float medianNeighbor=getMedianNeighbor(Ttask, horzBoundtask, i,j);
                int o=Otask.get(i,j);
	            float newT=optimizeTnew(Otask, Ttask, horzBoundtask, i,j, o, medianNeighbor,getBias(i,j));
                Ttemptask.set(i,j,newT);
	        }
        }
    } // end of for i for j
    
    nChanges=0;
    sumChanges=0;
    hasChangedtask.clear();
    hasChangedtask=UpperDiag<char>(options.firstRow+rowStart,options.firstCol,options.firstRow+rowEnd,options.lastCol, options.bandSize);

    // This is where the code adds up the number of changes (nchanges). I modified so that it is per task and does not consider halo cells. 
    // Also updates the Ttask matrix with what is in Ttemptask.
    // sumChanges is really the sum of squared residuals between iterations.
    for (int i=startIterRow;i<endIterRow;i++){
      for (int j=Otask.startCol(i)+1;j<Otask.endCol(i);j++) {
	        float tt=Ttemptask.get(i,j);
	        float x=Ttask.get(i,j);
	        if (x!=tt) {
	            hasChangedtask.set(i,j,1);
	            sumChanges+=(x-tt)*(x-tt);
	            Ttask.set(i,j,tt);
	            nChanges++;
	        }
        }
    } // end of transfer of Ttemptask to Ttask
    
    // Use new function evalulateFullLikelihoodSection to calculate likelihood per task, without halos
    double LikelihoodTask=evaluateFullLikelihoodSection(Otask, Ttask, horzBoundtask,options.firstRow+rowStartNoHalo,options.firstRow+rowEndNoHalo);

    // Reduce nChanges,sumChanges to root. Post early and do nonblocking. Wait at end of interation.
        MPI_Request reqChanges[3];
        MPI_Ireduce(&nChanges,&nChangesTot,1,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD,&reqChanges[0]);
        MPI_Ireduce(&sumChanges,&sumChangesTot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,&reqChanges[1]);
    // Reduce LikelihoodTask to root.
        MPI_Ireduce(&LikelihoodTask,&LikelihoodTot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,&reqChanges[2]);

    
    // Collect halos to send.
    // Third and fourth rows = last two rows for halo on previous rank
    for (int jjj=topstartSend;jjj<=topendSend;jjj++){
        bufftoprows[jjj-topstartSend]=Ttask.get(rowStart+2+options.firstRow,jjj);
    }
    for (int jjj=nextstartSend;jjj<=nextendSend;jjj++){
        bufftoprows[jjj-nextstartSend+toprowsizeSend]=Ttask.get(rowStart+3+options.firstRow,jjj);
    }
    // end-2 and end-3 rows = first two rows on next rank
    for (int jjj=bottnextstartSend;jjj<=bottnextendSend;jjj++){
        buffbottrows[jjj-bottnextstartSend]=Ttask.get(rowEnd-3+options.firstRow,jjj);
    }
    for (int jjj=bottstartSend;jjj<=bottendSend;jjj++){
        buffbottrows[jjj-bottstartSend+bottnextsizeSend]=Ttask.get(rowEnd-2+options.firstRow,jjj);
    }

    
    MPI_Request reqsRecvTop;
    MPI_Request reqsRecvBott;
    MPI_Request reqsSendBott;
    MPI_Request reqsSendTop;
    // Receive from top neighbor, send to top. Not including rank 0.
    if (rank>0){
        MPI_Irecv(&buffrecvtop[0],toprowsize+nextrowsize,MPI_FLOAT,topneighRank,123,MPI_COMM_WORLD,&reqsRecvTop);
        MPI_Isend(&bufftoprows[0],toprowsizeSend+nextrowsizeSend,MPI_FLOAT,topneighRank,323,MPI_COMM_WORLD,&reqsSendTop);
    }
    // Send to bottom neighbor, receive from bottom. Not including last rank.
    if (rank<size-1) {
        MPI_Isend(&buffbottrows[0],bottrowsizeSend+bottnextsizeSend,MPI_FLOAT,bottomneighRank,123,MPI_COMM_WORLD,&reqsSendBott);
        MPI_Irecv(&buffrecvbott[0],bottrowsize+bottnextsize,MPI_FLOAT,bottomneighRank,323,MPI_COMM_WORLD,&reqsRecvBott);
    }

    // Organize the data received from top neighbor into UpperDiag format. Not including rank 0.
    if (rank>0){
        MPI_Wait(&reqsRecvTop,&status);
        MPI_Wait(&reqsSendTop,&status);
        for (int iii=0;iii<toprowsize;iii++){
            Ttask.set(options.firstRow+rowStart,topstart+iii,buffrecvtop[iii]);
        }
        for (int iii=0;iii<nextrowsize;iii++){
            Ttask.set(options.firstRow+rowStart+1,nextstart+iii,buffrecvtop[iii+toprowsize]);
        }

    }

    // Organize the data received from bottom neighbor into UpperDiag format. Not including last rank.
    if (rank<size-1){
        MPI_Wait(&reqsRecvBott,&status);
        MPI_Wait(&reqsSendBott,&status);
        for (int iii=0;iii<bottnextsize;iii++){
            Ttask.set(options.firstRow+rowEnd-1,bottnextstart+iii,buffrecvbott[iii]);
        }
        for (int iii=0;iii<bottrowsize;iii++){
            Ttask.set(options.firstRow+rowEnd,bottstart+iii,buffrecvbott[iii+bottnextsize]);
        }

    }

    MPI_Waitall(3,reqChanges, MPI_STATUSES_IGNORE); // wait for nChangesTot and sumChangesTot and LikelihoodTot to finish communication
    if (rank==0){     // Print results of reductions for global evaluation of IF matrix
        fprintf(stderr,"****** rep=%d: nChanges=%ld, sumChanges=%lf ******\n",rep, nChangesTot,sumChangesTot);
        fflush(stderr);
        fprintf(stderr,"*******After re-estimation Likelihood: %lf *******\n",LikelihoodTot);
    }
  } // end of MRF iteration loop

    hasChangedtask.clear();
    // In advanced version, do not need to compile at root every iteration. Only at end.
    // Communicate Ttask to root to populate full T 
        //First need to collect elements into a linear array.
        int sendcountall=0;
        for (int xx=rowStartNoHalo;xx<=rowEndNoHalo;xx++){
            for (int yy=Ttask.startCol(xx+options.firstRow);yy<Ttask.endCol(xx+options.firstRow);yy++){
                buffSendAll[sendcountall]=Ttask.get(xx+options.firstRow,yy);
                sendcountall++;
            }
        }
        // Gather at root. Nonblocking gatherv.
        MPI_Request reqFullT;
        MPI_Igatherv(&buffSendAll[0],numElNH,MPI_FLOAT,&buffCollectAll[0],numElEachTaskNH,displs,MPI_FLOAT,0,MPI_COMM_WORLD,&reqFullT);
        
        if (rank==0){ // If root, compile the entire T matrix
            MPI_Wait(&reqFullT,MPI_STATUS_IGNORE);
            // Reallocate T
            T=UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
            int recvcountall=0;
            for (int zz=options.firstRow;zz<=options.lastRow;zz++){ // For each row
                int herestartcol=T.startCol(zz);
                int hereendcol=T.endCol(zz);
                for (int uu=herestartcol;uu<hereendcol;uu++){
                    T.set(zz,uu,buffCollectAll[recvcountall]);
                    recvcountall=recvcountall+1;
                }
            }
        }
    

  Ttemptask.clear();

// ** END OF MRF **//
  if (rank==0){
  outputSparseMatrix(argv[2], T, bias); // Output, only root//
  }
MPI_Wait(&reqFullT,MPI_STATUS_IGNORE);
 // Overall time
 double FinalTime=MPI_Wtime()-StartTime;

 // Main MRF loop time
 double MRFTime=FinalTime-SetupTime;
 
 // Compile times
 double AllTimes[4];
 AllTimes[0]=KDETime;AllTimes[1]=SetupTime;AllTimes[2]=MRFTime;AllTimes[3]=FinalTime;

 // Reduce times and print
 double MaxTimes[4];
 MPI_Reduce(&AllTimes[0],&MaxTimes[0],4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
 if (rank==0){
 fprintf(stderr,"----------------------------------------------------\n");
 fprintf(stderr,"Max AKDE time: %f \n",MaxTimes[0]);
 fprintf(stderr,"Max Setup time: %f \n",MaxTimes[1]);
 fprintf(stderr,"Max MRF time: %f \n",MaxTimes[2]);
 fprintf(stderr,"Max Overall time: %f \n",MaxTimes[3]);
 
 // Print to file
 char fname[15];
 sprintf(fname, "TimeFile_%d.txt",size);
 outfile_fp = fopen(fname, "a");
 fprintf(outfile_fp,"%f %f %f %f\n",MaxTimes[0],MaxTimes[1],MaxTimes[2],MaxTimes[3]);
 fclose(outfile_fp);
 }

  MPI_Finalize();
}
