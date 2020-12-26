<div style="text-align: justify">

# **Test case 6: One seed**

## **Aim**

To test the following:

1. Total volume of all grains
2. Periodicity of tessellations generated

## **Expected result**

### **Expected results in brief:**

| Characteristic feature | Expected results | Spacing length = 1 | Spacing length = 2 | Spacing length = 2.5 | Spacing length = 5 |
|:--------------------------------------:|---------------------------------------------|:------------------:|:------------------:|:--------------------:|:------------------:|
| Total volume of all grains | $(size\: of \: simulation \: box)^2 * length\_z$ | 100 | 100 | 100 | 100 |
| Neighbors due to periodicity of tessellations generated | All neighbors must be same (cell number: 0) | 0 | 0 | 0 | 0 |

## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name one_seed
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
* length_z = 1
* dimension = 2
* cell_number = 0
* orientation_data = None

### **Files used**

The test case uses all modules of the package.

## **Obtained result**

The obtained results are summarized below:

| Characteristic feature | Spacing length = 1 | Spacing length = 2 | Spacing length = 2.5 | Spacing length = 5 |
|:--------------------------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Total volume of all grains | 100 | 100 | 100 | 100 |
| Neighbors due to periodicity of tessellations generated | 0 | 0 | 0 | 0 |


The results obtained were matching with the expected results.

</div>