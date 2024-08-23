# BME 590L Final Project: Parallelization of HIFI-MRF    
>**Author**: Brian Chan  
----------
## Summary  
This project parallelizes the HIFI-MRF pipeline described in Cameron, Dostie, and Blanchette, Genome Biology, 2020. The GitHub repository for the original HIFI package can be found at https://github.com/BlanchetteLab/HIFI.  Paralllelization was achieved using MPI. Here, we describe the contents of this repository as well as how to run the parallelized version of HIFI-MRF.   

This project was for the "Computational Foundations of Biomedical Simulations" class at Duke University taught by Dr. Amanda Randles.

-----------

Make sure that within your current working directory you have the following files and directories in addition to this `README.md`:  

- `Makefile`  
- `src`  
- `Data`  
- `Analysis`      

----------   
Within ``` src ``` there should be the following files:    

- `BAMtoSparseMatrix.py`  
- `callPeaks.cpp`  
- `HIFI.cpp`  
- `HIFI_serial.cpp`  
- `HIFI_misc.h`  
- `HIFI_MyMatrix.h`  
- `HIFI_options.h`  
- `HIFI_options.cpp`  
- `HIFI_Advanced.cpp`  
- `parseHIFIoutput.py`  
- `plotHIFIoutput.py`  
- `SparseToFixed.py`  

You will not used any of the Python scripts or `callPeaks.cpp`, but they are part of the original pipeline and can be used for pre-/post-processing.
See original documentation for details.  

``` HIFI_serial.cpp ``` is the serial code. It has some comments in it from the original source, but many comments are added by yours truly to help guide the user.  
``` HIFI.cpp ``` is the parallelized code with comments.    
``` HIFI_Advanced.cpp ``` further optimizes and parallelizes HIFI-MRF.
```HIFI_misc.h```, ```HIFI_MyMatrix.h```, ```HIFI_options.h```, and ```HIFI_options.cpp``` define class templates and HIFI parameters used during runtime.

-------------------  
Within ``` Data ``` there should be the following files:  
- `HindIII.hg19.chr9_chr9.RF.tsv`  
- `MboI.hg38.chr4_chr4.RF.Truncate_1000.tsv`  
- `MboI.hg38.chr4_chr4.RF.Truncate_2000.tsv`  
- `MboI.hg38.NEWBAM.chr4_chr4.RF.Truncate_4000.tsv`  
- `MboI.hg38.NEWBAM.chr4_chr4.RF.Truncate_8000.tsv`  
- `MboI.hg38.chr4_chr4.RF.Truncate_16000.tsv`  


`HindIII.hg19.chr9_chr9.RF.tsv` was used for strong scaling analysis. The `MboI...` files were used for weak scaling analysis The `MboI...` files have restriction fragment read count matrices for either 1000, 2000, 4000, 8000, or 16000 restriction fragments.  The `HindIII` data is the supplied example from the original GitHub repository (see above). The `MboI...` files are processed from the BAM file found in the 4DNucleome project (https://data.4dnucleome.org/files-processed/4DNFIP9ADMXB/#details).

Note that the files in this directory are already pre-processed from BAM files using the `BAMtoSparseMatrix.py` script in the `src` directory on the Duke Compute CLuster. With reference to the original HIFI package documentation, this project starts optimization from step 2 in the Quick Start section of the original `README.md`.

-------------------  
Within ``` Analysis ``` there should be the following files:  
- `CompSCC.m`  
- `PlotHeatMap.m`  
- `plotIFTriangles.m`  
- `CompareIFs.m`  

`CompSCC.m` calculates the stratum-adjusted correlation coefficient between two interaction frequency (IF) matrices (see Yang, et al, Genome Research, 2017).  
`PlotHeatMap.m` visualizes IF matrices as heat maps.  
`plotIFTriangles.m` simultaneously plots two IF matrix heat maps: one in the upper triangular portion of the plot and the second in the lower triangular portion.  
`CompareIFs.m` uses the three previous files to run a full comparison between two IF matrices.  

To use, enter the `Analysis` directory, load MATLAB and run:  

```CompareIFs(<SERIAL_IF_TSV>, <PARALLEL_IF_TSV>, <OUTPUT_FILE_PREFIXES>);```    

where `<SERIAL_IF_TSV>` is the path and filename of the IF matrix `.tsv` file from the serial HIFI-MRF output, `<PARALLEL_IF_TSV>` is the path and filename of the IF matrix `.tsv` file from the pHIFI-MRF output, and `<OUTPUT_FILE_PREFIXES>` is the desired prefix to name the output files. All three arguments must be strings. The function will plot a heat map where the upper triangular portion is the parallel IF matrix and the lower triangular portion is the serial IF matrix. This will be saved in the file `<OUTPUT_FILE_PREFIXES>.pdf`. Another output file called `<OUTPUT_FILE_PREFIXES>.txt` will contain the SSE, SCC, and correlation per stratum (diagonal) for the two IF matrices. See comments within `CompareIFs.m` for returning function values in a MATLAB session.

-------------------  

## Compile and run  

Compilation uses the g++ and mpicc. To compile the original serial HIFI code:

``` make HIFIserial ```  

To run the original serial HIFI code:

``` src/HIFIserial <PATH/TO/INPUT> <OUTPUTNAME>.tsv -method=mrf```  

To compile the parallel version (pHIFI-MRF) used for the majority of the project:  

``` make HIFI ```    

To run:

```ibrun -np <NUMTASKS> src/HIFI <PATH/TO/INPUT> <OUTPUTNAME>.tsv -method=mrf ```  
or  
```mpirun -np <NUMTASKS> src/HIFI <PATH/TO/INPUT> <OUTPUTNAME>.tsv -method=mrf ```  

To compile the advanced parallel version (advanced pHIFI-MRF) used for preliminary testing of continued optimizations:  

``` make HIFIAdvanced ```  

To run:

```ibrun -np <NUMTASKS> src/HIFIAdvanced <PATH/TO/INPUT> <OUTPUTNAME>.tsv -method=mrf ```  
or  
```mpirun -np <NUMTASKS> src/HIFIAdvanced <PATH/TO/INPUT> <OUTPUTNAME>.tsv -method=mrf ```  

For the parallel versions, ``` <NUMTASKS> ``` is the number of MPI tasks you want to use. For all versions, `<PATH/TO/INPUT>` is the path and filename of your desired input file, for example `Data/HindIII.hg19.chr9_chr9.RF.tsv`. For all versions, ``` <OUTPUTNAME>.tsv ``` will be the name of the `.tsv` file containing the optimized interaction frequency matrix. The parallel versions will also create a file called `TimeFile_<NUMTASKS>.txt` with the runtimes for the adaptive kernel density estimation algorithm, the full setup of the code before the Markiv random field component, the Markov random field component, and the entire code. The program will also print periodic status updates to standard error. Note that the parallel code can also be run with a single task without issue.  

