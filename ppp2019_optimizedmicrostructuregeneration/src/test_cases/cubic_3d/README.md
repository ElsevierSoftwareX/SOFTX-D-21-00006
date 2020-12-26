<div style="text-align: justify">

# **Test case 2: Cubic 3D**

## **Aim**

To test the following structural characteristics:

1. Grain size of each grain
2. Number of neighbors of each grain
3. Grain boundary area of each grain
4. Junction lengths
5. Junction angles (in degrees)
6. Spacing lengths

## **Expected result**

### **Expected results in brief:**

| Characteristic feature | Expected Result | Spacing length = 1 | Spacing length = 2 | Spacing length = 2.5 | Spacing length = 5 |
|:----------------------------:|:---------------------------------------------------------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Grain size | 1.2407 * spacing length | 1.2407 | 2.4814 | 3.1017 | 6.2035 |
| Grain size | $d = \sqrt[3]{\frac{total \: volume * 6 }{\pi * number\_of\_grains}}$ | 1.2407 | 2.4814 | 3.1017 | 6.2035 |
| Number of Neighbors | 6 | 6 | 6 | 6 | 6 |
| Grain boundary area | Each GB area = $(spacing\: length)^2$ | 1 | 4 | 6.25 | 25 |
| Junction length | Equal to spacing length | 1 | 2 | 2.5 | 5 |
| Junction angles (in degrees) | 90° | 90° | 90° | 90° | 90° |


## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name cubic_3d
```
### **Options:**
1. `--name` refers to the test case name (cubic_2d in this case).
2. `--f` refers to flag if opaque surface is required for VTK and OBJ files of entire configuration.
3. `--rseed` refers to the seed of Numpy random function (for eg: `--rseed 1`).

Please refer section 'Basic usage' section for more details.

### **Parameters used**

All the parameters used for this test case are enlisted below:

* store_folder = "output_of_tests"    
* skewed_boundary_flag = False
* face_flag = None
* number_of_bins = 10
* size_of_simulation_box = 10
* spacing_lengths = [1, 2, 2.5, 5]
* required_texture = np.array([1, 1, 1])
* rand_quat_flag = True
* mesh_flag = 'TET'
* global_mesh_size = 0.5
* dimension = 3
* orientation_data = None

### **Files used**

The test case uses all modules of the package.

## **Obtained result**

The obtained results are summarized below:

| Characteristic feature | Spacing length = 1 | Spacing length = 2 | Spacing length = 2.5 | Spacing length = 5 |
|:----------------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Grain size | 1.2407 | 2.4814 | 3.1017 | 6.2035 |
| Grain size | 1.2407 | 2.4814 | 3.1017 | 6.2035 |
| Number of Neighbors | 6 | 6 | 6 | 6 |
| Grain boundary area | 1 | 4 | 6.25 | 25 |
| Junction length | 1 | 2 | 2.5 | 5 |
| Junction angles (in degrees) | 90° | 90° | 90° | 90° |

The results obtained were matching with the expected results.

</div>