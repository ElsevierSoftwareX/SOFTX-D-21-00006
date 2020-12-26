<div style="text-align: justify">

# **Optimized Microstructure Generator**

|<img src="./images/tessellations.png" title="Overview of generating tessellations"> |
| :---: |
| **Overview of generating tessellations** |

|  |
| :---: |
|<img src="./images/meshing.png" title="Overview of meshing process"> |
| **Overview of meshing process** |

In this package an attempt is made to create a tool to optimize these Tessellations based on user requirements. Basically user can provide the desired distribution along with the total number of seeds required, as input and the tool will optimize the tessellations generated to suite the user requirements. The final configuration obtained can then be meshed using hexahedral or tetrahedral elements.

The tool is developed at TU Bergakademie Freiberg, Germany under the supervision and guidance of Dr. Arun Prakash.

# **Quick Start**

The tool would be compiled into a package in the near future. Till then you can directly clone the git repository to your desired directory and then navigate in the terminal to the directory `ppp2019_OptimizedMicrostructureGeneration`.

Before starting, please ensure that you have statisfied all the [**requirements**](https://gitlab.com/arun.prakash.mimm/ppp2019_optimizedmicrostructuregeneration/-/wikis/Home/Requirements#requirements) and most importantly, ensure that you have installed all the [**required libraries**](https://gitlab.com/arun.prakash.mimm/ppp2019_optimizedmicrostructuregeneration/-/wikis/Home/Requirements/installing_libraries#installing-libraries).

## **Basic Usage**

For detailed reference usage manual please [**click here**](https://gitlab.com/arun.prakash.mimm/ppp2019_optimizedmicrostructuregeneration/-/wikis/home#optimized-micro-structure-generator).

### **Navigate to the directory and changing file permissions**

First ensure that you are in the correct directory. `pwd` command can be used to find the current working directory.

```bash
$ pwd
...../ppp2019_optimizedmicrostructuregeneration
```
Once you are in the correct directory, we need to change the permissions of the file 'execute' and this can be achieved by

```bash
$ chmod 777 execute
```
In order to check that proper permissions are granted, use the following:

```bash
$ ls -al
```
The above command lists all the files in the current directory with its repective permission details in the first column. The first column of the file 'execute' should look like `-rwxrwxrwx`. This means that all the three blocks ie; user, group and others has been granted with read, write and execute permissions for this file.

In order to execute the program, you can now use:

```bash
$ ./execute input-arguments
```


### **Undertanding available options**

It is better to have a grasp of the options available

```bash
$ ./execute main --help
```
The following would be displayed:
```bash
Usage: main.py main [OPTIONS]

Options:
  --s FLOAT...     Size of simulation box along X, Y & Z direction in the
                   format n n n  [required]
  --d INTEGER      Dimension of study ie; 2D or 3D in the format n where n is
                   an integer
  --n INTEGER      Number of seeds in the format n
  --t TEXT         Target distribution file name as a string stored in the
                   same directory where the package is executed from
  --c INTEGER      The characteristic that has to be optimized in the format n
                   where n corresponds to the integer number corresponding to
                   the characteristic
  --m TEXT         The name of the material as a string
  --st INTEGER...  The direction of stress for computing the Schmid Factors
  --o FLOAT...     Required texture common to each grain in the format n n n
                   as provided in the documentation
  --r              Flag to indicate if Random orientations is to be generated
  --noopti         Flag to indicate if optimization is not to be performed
  --f              This flag is to be used to indicate if a closed surface is
                   to be used for visualization files for all grain in one
                   file
  --ss TEXT        Option to indicate if randomly placed seeds are required or
                   regularly spaced seeds ie; cubic_2d, cubic_3d, etc
  --sl FLOAT       Option to specify the spacing length such that it is
                   exactly a multiple of size of simulation box along all
                   three directions
  --method TEXT    Method to be used for optimization
  --sk             Flag to indicate skewed boundary requirement
  --ucf TEXT       Specify the user defined cost function file name. Refer
                   documentation for file specifications
  --msh TEXT       Flag to indicate type of meshing required of the simulation
                   box. Eg. hex (for Hexahedral), tet (for Tetrahedral), vis
                   (for Visualization)
  --gms FLOAT      Provide global mesh size
  --maxf INTEGER   Provide maximum number of function evaluation during
                   optimization
  --rseed INTEGER  Enter the seed value for Numpy random function
  --nbins INTEGER  Specify the number of bins
  --help           Show this message and exit.
```

### **Generate microstructure without optimization**

The microstructure with Cubic spacing in 2D without optimization can be obtained by using:

```bash
$ ./execute main --s 10 10 10 --d 2 --m steel --r --st 0 1 0 --ss cubic_2d --noopti
```
Here `--noopti` indicates that optimization is not required. The `cubic_2d` can also be replaced with `cubic_3d`, `fcc_2d`, `fcc_3d`, `bcc_3d`, or `random_3d` along with appropriate value for `--d` option representing the dimension.

### **Generate microstructure with optimization**

The microstructure with Random spacing in 3D with optimization can be obtained by using:

```bash
$ ./execute main --s 10 10 3 --d 3 --m steel --r --st 0 1 0 --ss random_3d --t user_grain_size_distribution.txt --c 0 --n 10
```
Here `--n` represents the number of seeds required. The results are stored in the directory `visualization_files`.

### **Generate microstructure with meshing but without optimization**

Here we would be again using the case of `cubic_2d` seed spacing.

```bash
$ ./execute main --s 10 10 10 --d 2 --m steel --r --st 0 1 0 --ss cubic_2d --noopti --msh hex
```

Here `--msh hex` represents the meshing of obtained configuration using hexahedral meshes. The `hex` keyword can also be replaced with `tet` for tetrahedral elements or `vis` to generate tetrahedral elements for the meshing of raw configuration obtained directly without preprocessing into a cuboid, just for visualization purposes.

**For more details regarding the possibilities of microstructure generation using this tool, please refer the [reference usage manual](){Update link here}.**

</div>
