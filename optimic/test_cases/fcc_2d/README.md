<div style="text-align: justify">

# **Test case 3: FCC 2D**

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
|:----------------------------:|:---------------------------------------------------------------------------------------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Grain size | 0.7979 * spacing length | 0.7979 | 1.5958 | 1.9947 | 3.9895 |
| Grain size | $d = \sqrt{\frac{\frac{total \: volume}{length\: along \: Z-axis} * 4 }{\pi * number\_of\_grains}}$ | 0.7979 | 1.5958 | 1.9947 | 3.9895 |
| Number of Neighbors | 4 | 4 | 4 | 4 | 4 |
| Grain boundary area | All GB areas must be same and must be equal to $0.7071*spacing\_length*length\_z$ | 0.7071 | 1.4142 | 1.7677 | 3.5355 |
| Junction length | All junction lengths must be  equal to $0.7071*spacing\_length$ or length_z | 0.7071 or 1 | 1.4142 or 1 | 1.7677 or 1 | 3.5355 or 1 |
| Junction angles (in degrees) | 90° | 90° | 90° | 90° | 90° |


## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name fcc_2d
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
* orientation_data = None

### **Files used**

The test case uses all modules of the package.

## **Obtained result**

The obtained results are summarized below:

| Characteristic feature | Expected Result | Spacing length = 1 | Spacing length = 2 | Spacing length = 2.5 | Spacing length = 5 |
|:----------------------------:|:---------------------------------------------------------------------------------------------------:|:------------------:|:------------------:|:--------------------:|:------------------:|
| Grain size | 0.7979 * spacing length | 0.7979 | 1.5958 | 1.9947 | 3.9895 |
| Grain size | $d = \sqrt{\frac{\frac{total \: volume}{length\: along \: Z-axis} * 4 }{\pi * number\_of\_grains}}$ | 0.7979 | 1.5958 | 1.9947 | 3.9895 |
| Number of Neighbors | 4 | 4 | 4 | 4 | 4 |
| Grain boundary area | All GB areas must be same and must be equal to $0.7071*spacing\_length*length\_z$ | 0.7071 | 1.4142 | 1.7677 | 3.5355 |
| Junction length | All junction lengths must be  equal to $0.7071*spacing\_length$ or length_z | 0.7071 or 1 | 1.4142 or 1 | 1.7677 or 1 | 3.5355 or 1 |
| Junction angles (in degrees) | 90° | 90° | 90° | 90° | 90° |

The results obtained were matching with the expected results.

</div>
