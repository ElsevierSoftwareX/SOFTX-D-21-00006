<div style="text-align: justify">

# **Test case 9: Random 3D**

## **Aim**

To test the following:

1. Read seeds data without orientation data from file
2. Read seeds data with orientation data from file
3. Total number of seeds
4. Total volume of all grains

## **Expected result**

### **Expected results in brief:**

| Characteristic feature | Expected Result | Expected output |
|:--------------------------------------------------:|:----------------------------------:|:-------------------:|
| Read seeds data without orientation data from file | Original seeds data | Original seeds data |
| Read seeds data with orientation data from file | Original seeds data | Original seeds data |
| Total number of seeds | Total input seeds | 100 |
| Total volume of all grains | $(size\: of\: simulation\: box)^3$ | 1000 |

## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name random_3d
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
* number_of_seeds = 100

### **Files used**

The test case uses all modules of the package.

## **Obtained result**

The obtained results are summarized below:

| Characteristic feature | Obtained output |
|:--------------------------------------------------:|:-------------------:|
| Read seeds data without orientation data from file | Original seeds data |
| Read seeds data with orientation data from file | Original seeds data |
| Total number of seeds | 100 |
| Total volume of all grains | 1000 |

The results obtained were matching with the expected results.

</div>


