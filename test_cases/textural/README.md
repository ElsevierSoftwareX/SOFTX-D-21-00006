<div style="text-align: justify">

# **Test case 8: Textural**

## **Aim**

To test the following textural characteristics:

* 'sharp_texture_quaternions' function returns unit quaternions
* Disorientation Angles
    * Maximum value of disorientation angles
    * Grains with same orientations
* Schmid Factor
    * Maximum value of Schmid factor
    * Schmid Factor for known stress directions
* Known misorientation data 

## **Expected result**

### **Expected results in brief:**

| Characteristic feature | Expected Result |
|:-------------------------------------------------------:|:-------------------------:|
| Norm of sharp texture quaternions | 1 |
| Maximum value of disorientation angles | <= 62.8 |
| Disorientation angles for grains with same orientations | All angles are equal to 0 |
| Maximum value of Schmid factor | <= 0.5 |

The following table shows the heighest schmid factor for the known stress directions:

| **Stress Direction [X, Y, Z]** | **Highest Schmid factor** |
|:------------------------------:|:-------------------------:|
| [1, -1, 0] | 0.408 |
| [1, 0, 0] | 0.408 |
| [1, 1, 0] | 0.408 |

**Known misorientation data**

The misorientation data (angle, axis and GB type) is based on the data given in *'The measurement of Grain Boundary geometry', V Randle, Institute of Physics (IOP) Publishing, 1993,* Page No. 40-43, Table 3.3.

The misorientation axis data is then tested to ensure that all the GB types along with its corresponding misorientation angles and axes has been captured.

## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name textural
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

**Norm of sharp texture quaternions**

* number_of_quaternions = 1

**For disorientation angles**

* orientation_data = None

**For all grains having same crystal orientation**
    
* orientation_quaternion = np.array([1, 0, 0, 0])

**Schmid factors with random orientations**
    
* orientation_data = None
* stress_direction = np.array([1, 0, 0])

**Schmid factors with same orientations but with different stress directions**

* orientation_quaternion = np.array([1, 0, 0, 0])
* stress_direction = np.array([1, -1, 0])
* stress_direction = np.array([1, 0, 0])
* stress_direction = np.array([1, 1, 0])

### **Files used**

The test case uses all modules of the package.

## **Obtained result**

The obtained results are summarized below:

| Characteristic feature | Obtained Results |
|:-------------------------------------------------------:|:-------------------------:|
| Norm of sharp texture quaternions | 1 |
| Maximum value of disorientation angles | <= 62.8 |
| Disorientation angles for grains with same orientations | All angles are equal to 0 |
| Maximum value of Schmid factor | <= 0.5 |

The following table shows the heighest schmid factor for the known stress directions:

| **Stress Direction [X, Y, Z]** | **Obtained  Schmid factor** |
|:------------------------------:|:-------------------------:|
| [1, -1, 0] | 0.408 |
| [1, 0, 0] | 0.408 |
| [1, 1, 0] | 0.408 |

The results obtained were matching with the expected results.

</div>