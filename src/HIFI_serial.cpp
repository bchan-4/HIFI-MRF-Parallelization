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

using namespace std;
//using namespace __gnu_cxx;
 
HIFI_options options;
HIFI_matrixProperties matrixProperties;

UpperDiag<float> biasMatrix;

double *precomputedNormPDF;
double *precomputedNormLogPDF;

#define NORM_RES 10000
int fullRange=100*NORM_RES;
int floatFullRange=2*100*NORM_RES;

float *bias;

// retrieve normpdf value with mean m and sigma s
float normpdf(float x, float m, float s) {
  int i=((x-m)/s)*NORM_RES + fullRange;
  if (i<0 || i>floatFullRange) return 0;
  return precomputedNormPDF[i];
}


void outputSparseBoundaryMatrix(char *fn, UpperDiag<char> &horizontalBoundary,UpperDiag<char> &verticalBoundary );


// for boundary identification. These boundaries set "domains" within the matrix across which elements cannot be part of the same neighborhood.
// Iterates through columns and rows and searches for where distributions of entries between OoverBias(i,j:j') and OoverBias(i+1,j:j') are significantly different, based on Kolmogorov-Smirnov (KS) test.
#define MAX_NBNONMAX 10    
void setBoundaries(UpperDiag<float> &OoverBias, UpperDiag<char> &horizontalBoundary, UpperDiag<char> &verticalBoundary) {

  
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
	//	verticalBoundary.set(j,i,1);
      }
      i--;
    }
  }


  // vertical
  for (int j=OoverBias.startCol(OoverBias.startRow());j<OoverBias.endCol(OoverBias.startRow())-1;j++) {

    // consider a boundary starting down from (i,j)
    int distLeft[10000];
    int distRight[10000];
    int maxi=-1;
    int maxi2=-1;
    float maxks=-1;
    
    for (int i=OoverBias.startRow();i<OoverBias.endRow();i++) {
      if (i>=j) continue;
      if (OoverBias.get(i,j)==0 && OoverBias.get(i,j+1)==0) continue;
      
      memset(distLeft,0,10000*sizeof(int));
      memset(distRight,0,10000*sizeof(int));
      int nLeft=0;
      int nRight=0;
      int nbNonMax=0;
      float maxSearch=0;
      int maxEntry=0;

      int i2=0;


      for (i2=i;i2<=j && i2<OoverBias.endRow() &&  nbNonMax<MAX_NBNONMAX;i2++) { // end point of boundary
	nbNonMax++;
	// only consider end points where there is data
	if (verticalBoundary.get(i2,j)) break;
	int entryLeft=(int)(OoverBias.get(i2,j)+0.49);
	int entryRight=(int)(OoverBias.get(i2,j+1)+0.49);

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
	  maxi=i;
	  maxi2=i2;
	}
	else {nbNonMax++;}
      }
    }

    if (maxks>options.boundaryKS) {
      for (int i=maxi;i<=maxi2;i++) {
	verticalBoundary.set(i,j,1);
      }
      j--;
    }
  }

  if (options.boundaryOutput[0]) {
    outputSparseBoundaryMatrix(options.boundaryOutput, horizontalBoundary, verticalBoundary);
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
	   
// adaptive kernel density estimation. Was not originally planning on parallelizing this, but honestly I may because it probably will not be too much effort. Will test to see if serial version is slow.
// This is just embarrassingly parallel, though.
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
float getMedianNeighbor(UpperDiag<float> &T, UpperDiag<char> &horizontalBoundary, UpperDiag<char> &verticalBoundary, int i, int j) {
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
	
	if (verticalBoundary.get(i,j-1) && off_j==-1) continue;
	if (verticalBoundary.get(i,j) && off_j==+1) continue;
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
double evaluateFullLikelihood(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> &horizontalBoundaries, UpperDiag<char> &verticalBoundaries) {
  double like=0;
  
  for (int i=O.startRow();i<O.endRow();i++) {
    for (int j=O.startCol(i)+1;j<O.endCol(i);j++) {
       
      float medianNeighbor=getMedianNeighbor(T, horizontalBoundaries, verticalBoundaries, i,j);
      float t=T.get(i,j);
      int o=O.get(i,j);
      double l=evaluateLikelihood(t, o, medianNeighbor,getBias(i,j));
      like+=l;
    }
  }
  return like;
}

// log-likelihood of neighborhood, as above function. neighborhood of ii-1:ii+1 inclusive (unless hitting boundary) and jj-1:jj+1 inclusing (unless hitting boundary).
double evaluateFullLikelihoodPatch(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> &horizontalBoundaries, UpperDiag<char> &verticalBoundaries, int ii, int jj) {
  double like=0;
  
  for (int i=max(O.startRow(),ii-1);i<min(O.endRow(),ii+2);i++) {
    for (int j=max(O.startCol(i)+1,jj-1);j<min(O.endCol(i),jj+2);j++) {
       
      float medianNeighbor=getMedianNeighbor(T, horizontalBoundaries, verticalBoundaries, i,j);
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
// This section goes through an iterated conditional mode (ICM) algorithm, using the "evaluateFUllLikelihoodPatch" as the optimized function and the multiplicative factors above to modify the "guesses" for T(i,j)
// Note that each iteration through the entire matrix goes through ALL of the factors above.
double optimizeTnew(UpperDiag<float> &O, UpperDiag<float> &T, UpperDiag<char> horizontalBoundaries, UpperDiag<char> verticalBoundaries,  int i, int j, int obs, float medNei, float localbias) {
  
  // search for the best multiplicative factor to apply to currenT
  float loglike[NFACTORS];
  float old=T.get(i,j);
  for (int f=0;f<NFACTORS;f++) {
    T.set(i,j,old*factor[f]);
    loglike[f]=evaluateFullLikelihoodPatch(O, T, horizontalBoundaries, verticalBoundaries,  i,j);
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
void readFullOMatrix(char *fn, UpperDiag<float> *mat, float *bias, bool setBias) {
  
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
	}
    }
  }

  fclose(f);
  
  // Calculate bias matrix
  if (setBias) {
    for (int i=matrixProperties.fullMatrix_firstRow;i<matrixProperties.fullMatrix_lastRow;i++) {
      bias[i]=sumRow[i]/((fullSum)/(matrixProperties.fullMatrix_lastRow-matrixProperties.fullMatrix_firstRow+1));
      if (bias[i]>MAX_BIAS) bias[i]=MAX_BIAS;
      if (bias[i]<1.0/MAX_BIAS) bias[i]=1.0/MAX_BIAS;
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


// This function goes through the actual Markov Random Field. Inputs are observed read counts (O), read counts/bias (calculated with a function above), 
// buffer UpperDiag to hold the optimized sparse matrix T, and the bias factor array
void computeTMatrix_mrf(UpperDiag<float> &O, UpperDiag<float> &OoverBias, UpperDiag<float> &T, float *bias) {

 // Set up the boundaries in the matrix.
  UpperDiag<char> horizontalBoundaries  = UpperDiag<char>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
  UpperDiag<char> verticalBoundaries =  UpperDiag<char>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize); 

  setBoundaries(OoverBias, horizontalBoundaries, verticalBoundaries);

 // Holds temporary optimized T matrices
  UpperDiag<float> Ttemp = UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);

  // Indicates if/where each iteration has updated the T matrix
  UpperDiag<char> hasChanged=UpperDiag<char>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);

  // initalization of T
  // Start with adaptive KDE as the initial "guess" for the optimized matrix "T". 
  // Initially I was not planning on parallelizing this, but I still might since it seems amenable to parallelization too.
  fprintf(stderr,"Starting KDE initialization\n");
  fflush(stderr);
  computeTMatrix_kde(O, OoverBias, T);
  
  fprintf(stderr,"Starting MRF\n");
  unsigned long nChanges=1;  
  
  // Iterate through the Markov Random Field to update the optimized T matrix. Default uption is 5 iterations.
  // May change this to a tolerance for checking convergance, since I thought that was what was here...
  for (int rep=0;rep<options.mrfMaxIter;rep++) {
    double sumChanges=0;

    // Iterate through all elements of the observed read count matrix O    
    for (int i=O.startRow();i<O.endRow();i++) {
      for (int j=O.startCol(i)+1;j<O.endCol(i);j++) {
	
    // check if anything in the neighbourhood has changed
	
	int neiLeft=max(O.startRow(),i-2);
	int neiRight=min(O.endRow()-1,i+2);
	int neiBottom=max(O.startCol(i),j-2);
	int neiTop=min(O.endCol(i)-1,j+2);

	  bool changed=false;
	  for (int a=neiLeft;a<=neiRight;a++) {
	    for (int b=neiBottom;b<=neiTop;b++) {
	      if (a<=b) changed=changed || hasChanged.get(a,b);
	    }
	  }

        
	  float t=T.get(i,j);	
	  Ttemp.set(i,j,t);
      // If I am at the first iteration or something has changed in T, do the optimization again.
	  if (rep==0 || changed) {

	    float medianNeighbor=getMedianNeighbor(T, horizontalBoundaries, verticalBoundaries, i,j);
	    int o=O.get(i,j);
	    float newT=optimizeTnew(O, T, horizontalBoundaries, verticalBoundaries, i,j, o, medianNeighbor,getBias(i,j));
	    Ttemp.set(i,j,newT);
	  }
      }
    } // end of for i for j
    
    nChanges=0;
    sumChanges=0;
    hasChanged.clear();
    hasChanged=UpperDiag<char>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);

    // This is where the code adds up the number of changes in the entire matrix (nchanges). Also updates the T matrix with what is in Ttemp.
    // sumChanges is really the sum of squared residuals between iterations.
    for (int i=O.startRow();i<O.endRow();i++) {
      for (int j=O.startCol(i)+1;j<O.endCol(i);j++) {
	float tt=Ttemp.get(i,j);
	float x=T.get(i,j);
	if (x!=tt) {
	  hasChanged.set(i,j,1);
	  sumChanges+=(x-tt)*(x-tt);
	  T.set(i,j,tt);
	  nChanges++;
	}
      }
    } // end of transfer of Ttemp to T
    
    fprintf(stderr,"rep=%d: nChanges=%ld, sumChanges=%lf\n",rep, nChanges,sumChanges);
    fflush(stderr);
    fprintf(stderr,"After re-estimation Likelihood: %lf\n",evaluateFullLikelihood(O, T, horizontalBoundaries, verticalBoundaries));

  } // end of MRF iteration loop
  
  Ttemp.clear();

  hasChanged.clear();
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
void outputSparseBoundaryMatrix(char *fn, UpperDiag<char> &horizontalBoundary,UpperDiag<char> &verticalBoundary ) {
  FILE *out=fopen(fn,"w");
  if (options.lastRow!=-1 && options.lastCol!=-1) fprintf(out,"# %d %d %d %d\n",options.firstRow,options.lastRow,options.firstCol,	    options.lastCol);
  else fprintf(out,"# %d %d %d %d\n",options.firstRow, options.lastRow, options.firstCol, options.lastCol);

  for (int i=horizontalBoundary.startRow();i<horizontalBoundary.endRow();i++) {
    for (int j=horizontalBoundary.startCol(i)+1;j<horizontalBoundary.endCol(i)-1;j++){
      if (horizontalBoundary.get(i,j) || (i>horizontalBoundary.startRow() && horizontalBoundary.get(i-1,j)) ||  verticalBoundary.get(i,j) || (j>verticalBoundary.startCol(j) && verticalBoundary.get(i,j-1))) fprintf(out,"1 %s %d %s %d\n",matrixProperties.chr1,i+options.firstRow,matrixProperties.chr2,j+options.firstCol);
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

  // declares matrices for observed read counts O, observed read counts divided by the bias factor OoverBias, and the optimized interaction matrix T.
  UpperDiag<float> O; 
  UpperDiag<float> OoverBias; 
  UpperDiag<float> T; 


  options.read(argc, argv);
  options.print();

  precomputeStuff(); // Precomputes distributions, etc.

  getMatrixSize(argv[1]); // size of matrix provided.

  // allocate tables
  T = UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
  O = UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
  OoverBias = UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
  
  bias=(float*) malloc(matrixProperties.fullMatrix_size*sizeof(float));

  // Initializes bias to 1.
  for (int i=0;i<matrixProperties.fullMatrix_size;i++) bias[i]=1;

  fprintf(stderr,"Reading matrix\n"); // Reads in .tsv file with read count matrices. Also calculates the bias list (see "Fragment-specific bias calculation" in Methods of paper).
  readFullOMatrix(argv[1], &O, bias, true);

  // Reads in bias matrix if supplied
  if (options.normalizationMatrix[0]) {
    biasMatrix = UpperDiag<float>(options.firstRow,options.firstCol,options.lastRow,options.lastCol, options.bandSize);
    readBiasMatrix(options.normalizationMatrix,&biasMatrix);
  }
  
  computeOoverBias(O,OoverBias); // calculate bias-corrected read count matrix.
  if (options.method_fixed) computeTMatrix_fixed(O, T, bias);
  if (options.method_kde || options.method_akde)  computeTMatrix_kde(O, OoverBias, T);
  if (options.method_mrf) computeTMatrix_mrf(O, OoverBias, T, bias); /* This is the method I am focusing on */
  outputSparseMatrix(argv[2], T, bias); // Output //
}
  
