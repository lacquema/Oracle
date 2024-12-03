# Installing on your computer

## 0) Create a container

It is recommended that you create a virtual environment. It is an isolated execution environment. It allows the packages used for a project to be isolated, so that there are different versions for each project.

To do this, use the following command: 

`python3 -m venv <environment_path>`

Note that the path to the python3 executable is in this environment: `<environment_path>/bin/python3`, and that the libraries are installed in the `<environment_path>/lib/` directory.

If you do not choose this option, simply create a directory:

`mkdir <environment_path>`


## 1) Download the github directory

Now that the container is ready, you can download the gihub directory inside:

`git clone https://github.com/lacquema/Oracle`


## 2) Compile and install packages

First, you need to open the newly installed github directory in your new environment:

`cd <environment_path>/Oracle/`

Now, you need to open the `<environment_path>/Oracle/Makefile` file and update the following parameters: COMPILF and PYTHON3. They correspond respectively to the paths to the executables of the fortran compiler installed on your computer and of python3 in your environment.

And finally:

`make all`

This command compiles all the code and installs all the necessary packages.




# Installing on a computing server

## 0) Create a container

Simply create a directory:

`mkdir <environment_path>`

## 1) Download the github directory

Now that the container is ready, you can download the gihub directory inside:

`git clone https://github.com/lacquema/Oracle`

## 2) Compile and install packages

The use of modern package managers such as Nix or Guix is recommended. These tools enable dependencies to be managed in an isolated and reproducible way, without polluting the overall system environment.






# Launch the main code

`python3 <environment_path>/Oracle/Interface/Main.py`


