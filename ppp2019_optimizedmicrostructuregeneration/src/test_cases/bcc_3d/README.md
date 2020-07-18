<div style="text-align: justify">

# **Test case 5: BCC 3D**

## **Aim**

To test the following structural characteristics:

1. Grain size of each grain
2. Number of neighbors of each grain
3. Junction lengths
4. Spacing lengths

## **Expected result**

### **Expected results in brief:**

| Characteristic feature | Expected Result | Spacing length = 1 | Spacing length = 2 | Spacing length = 2.5 | Spacing length = 5 |
|:----------------------:|:---------------------------------------------------------------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Grain size | 0.9847 * spacing length | 0.9847 | 1.9695 | 2.4619 | 4.9235 |
| Grain size | $d = \sqrt[3]{\frac{total \: volume * 6 }{\pi * number\_of\_grains}}$ | 0.9847 | 1.9695 | 2.4619 | 4.9235 |
| Number of Neighbors | 14 | 14 | 14 | 14 | 14 |
| Junction length | All junctions have length and are exactly equal to $0.3536*spacing\_length$ | 0.3536 | 0.7071 | 0.884 | 1.768 |


## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name bcc_3d
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
|:----------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Grain size | 0.9847 | 1.9695 | 2.4619 | 4.9235 |
| Grain size | 0.9847 | 1.9695 | 2.4619 | 4.9235 |
| Number of Neighbors | 14 | 14 | 14 | 14 |
| Junction length | 0.3536 | 0.7071 | 0.884 | 1.768 |

The results obtained were matching with the expected results.