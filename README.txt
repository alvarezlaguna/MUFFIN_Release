Instructions to install easily MUFFIN:

1) Download the code with 

git clone https://github.com/alvarezlaguna/MUFFIN_Release.git ./MUFFIN

2) enter inside the directory and update the submodules with 

git submodule update --init --recursive

IF YOU HAVE ALREADY CONDA IN YOUR SYSTEM, YOU CAN SKIP TO STEP 4

3) Download mamba-forge in order to manage the packages. In the following website download the Mambaforge for your specific distribution:

https://github.com/conda-forge/miniforge

  3.1) Change the permissions to run the script:
	chmod +x Mambaforge-[PUT HERE YOUR DISTRIBUTION]-x86_64.sh
  3.2) Run the script:
	./Mambaforge-[PUT HERE YOUR DISTRIBUTION]-x86_64.sh

4) Open a new terminal and go inside MUFFIN, create the conda environment by running:
conda env create -f conda/environment.yml

5) Change to the created environment:

conda activate muffin-env

(IMPORTANT: YOU NEED TO DO STEP 5 EVERYTIME YOU COMPILE OR RUN THE CODE IN A NEW TERMINAL)

5) Create the build directory where we will compile the code:
mkdir build

6) Enter in the build directory and compile by running:
cd build; CC=mpicc CXX=mpicxx cmake ..; make -j 16; cmake --install .; 

7) In the code we need to use mpi4py. We need to install them:
env MPICC=mpicc python3 -m pip install --no-cache-dir mpi4py

8) We install the muffin library in python
cd ..; python3 -m pip install -e .

9) Do the regression tests:
cd Regression_tests_Euler; python3 Regression_tests.py 

(You can check the plot Results.pdf)

10) To test the code, go to Results/1_Test_Euler:
cd ../testcases/1_Test_Euler

11) Run the script by doing:
python3 Euler_FirstOrder.py

12) You can visualize the results by doing:
jupyter-notebook 
  12.1) Open Tutorial.ipynb and run it

WARNING: You need install the library in your python environment everytime you compile
`cmake --install build`



