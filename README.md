*******
Table of contents:

- [Installing](#install)

- [Launch the main code](#launch)

*******

<div id='install'/>  

# Installing

## 0) Prerequisites

You need python3 and a fortran compiler, as gfortran. 

## 1) Create a container

It is recommended that you create a virtual environment. It is an isolated execution environment. It allows the packages used for a project to be isolated, so that there are different versions for each project.

To do this, use the following command: 

`python3 -m venv <environment_path>`

You can now activate the environment: 

`source <environment_path>/bin/activate`

The `python3` command will now refer to the python3 installed in your environment `<environment_path>/bin/python3` and not the one on the system. And the libraries are now installed in the `<environment_path>/lib/` directory.


If you do not choose this option, simply create an usual directory:

`mkdir <environment_path>`


## 2) Download the github directory  

Now that the container is ready, you can download the gihub directory inside:

`cd <environment_path>`

`git clone https://github.com/lacquema/Oracle`


## 3) Install python packages

All the necessary python packages are listed in the `<environment_path>/Oracle/requirements.txt` file. You can install them all:

`python3 -m pip install -r <environment_path>/Oracle/requirements.txt`

**On computing servers**, the use of modern package managers such as Guix is highly recommended and often mandatory. These tools enable isolated and reproducible management of software dependencies, thus avoiding unintentional changes to the server environment. In addition, Guix avoids unnecessary duplication of packages by storing shared dependencies only once in the server's `/gnu/store`. This considerably reduces the use of storage space. Instead of installing packages directly, you configure them to reference the Guix-managed library. All required Python packages are also explicitly listed in the `<environment_path>/Oracle/manifest.scm` file. You can reference them all:

`guix package -m <environment_path>/Oracle/manifest.scm`

Ensure to consult the Guix documentation provided by your server administrators, as specific configurations or permissions might apply.

Note that software dependencies may includes the python packages we need, as well as python3 itself and the fortran compiler `gfortran-toolchain`.


## 4) Compile fortran code

Open the `<environment_path>/Oracle/Makefile` file and update the following parameters: COMPILF. This corresponds to the paths or simply the command for the executables of the fortran compiler installed on your computer. 

**On compute servers**, you can set it as `gfortran-toolchain`. 

You can now compile all the fortran files:

`make compile`




<div id='launch'/>  

# Launch the main code



`python3 <environment_path>/Oracle/Interface/Main.py`


