<div style="text-align: justify">

# **Test case 7: Two seed**

## **Aim**

To test the following:

1. Volume of each grain
2. Total volume of both grains
3. Distance between grains
4. Type of grain boundary
5. Disorientation angle between grains

## **Expected result**

### **Expected results in brief:**

**Volume of each grain, total volume of all grains and distance between grains**

| Characteristic feature | Expected Result | Expected output |
|:--------------------------:|:------------------------------------------------:|:---------------:|
| Volume of each grain | 50% of total volume | 50 |
| Total volume of all grains | $(size\: of \: simulation \: box)^2 * length\_z$ | 100 |
| Distance between grains | 0.5 * size of simulation box | 5 |

**Type of grain boundary**

| **Misorientation in the form of Quaternions (real first format)** | Expected **GB type (Σ)** |
|:-----------------------------------------------------------------:|:------------------------:|
| (0.8166416, 0.4081033, 0.4081033, 0) | 3 |
| (0.6710739, 0.6706129, 0.2235376, 0.2235376) | 5 |
| (0.8017756, 0.5345322, 0.2672661, 0) | 7 |
| (0.7745967, 0.5163978, 0.2581989, 0.2581989) | 15 |

**Disorientation angle between grains**

| **GB type** | **Misorientation in the form of Quaternions (real first format)**| Expected **Disorientation angles** |
| :---: | :---: | :---: |
| Σ3 | (0.8166416, 0.4081033, 0.4081033, 0) | 60°  |
| Σ5 | (0.6710739, 0.6706129, 0.2235376, 0.2235376) | 36.9° |
| Σ7 | (0.8017756, 0.5345322, 0.2672661, 0) | 38.2° |
| Σ15 | (0.7745967, 0.5163978, 0.2581989, 0.2581989) | 48.2° |

## **Command used to run the program**

Please refer section 'Execute test cases' of the documentation for more details.

Navigate to the **‘ppp2019_optimizedmicrostructuregeneration‘** directory using terminal. Once
you are in the appropriate directory, you can execute the test case using:

```bash
$ ./execute test --name two_seed
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

**# Σ = 15 grain boundary orientations** 

* orientation_data = np.array([[1, 0, 0, 0], [0.7745967, 0.5163978, 0.2581989, 0.2581989]])

**# Σ = 5 grain boundary orientations**

* orientation_data = np.array([[1, 0, 0, 0], [0.6710739, 0.6706129, 0.2235376, 0.2235376]])

**# Σ = 3 grain boundary orientations**

* orientation_data = np.array([[1, 0, 0, 0], [0.8166416, 0.4081033, 0.4081033, 0]])

**# Σ = 7 grain boundary orientations**

* orientation_data = np.array([[1, 0, 0, 0], [0.8017756, 0.5345322, 0.2672661, 0]])

### **Files used**

The test case uses all modules of the package.

## **Obtained result**

The obtained results are summarized below:

**Volume of each grain, total volume of all grains and distance between grains**

| Characteristic feature | Obtained output |
|:--------------------------:|:---------------:|
| Volume of each grain | 50 |
| Total volume of all grains | 100 |
| Distance between grains | 5 |

**Type of grain boundary**

| **Misorientation in the form of Quaternions (real first format)** | Obtained **GB type (Σ)** |
|:-----------------------------------------------------------------:|:------------------------:|
| (0.8166416, 0.4081033, 0.4081033, 0) | 3 |
| (0.6710739, 0.6706129, 0.2235376, 0.2235376) | 5 |
| (0.8017756, 0.5345322, 0.2672661, 0) | 7 |
| (0.7745967, 0.5163978, 0.2581989, 0.2581989) | 15 |

**Disorientation angle between grains**

| **GB type** | **Misorientation in the form of Quaternions (real first format)** | Obtained **Disorientation angles** |
|:-----------:|:-----------------------------------------------------------------:|:----------------------------------:|
| Σ3 | (0.8166416, 0.4081033, 0.4081033, 0) | 60° |
| Σ5 | (0.6710739, 0.6706129, 0.2235376, 0.2235376) | 36.9° |
| Σ7 | (0.8017756, 0.5345322, 0.2672661, 0) | 38.2° |
| Σ15 | (0.7745967, 0.5163978, 0.2581989, 0.2581989) | 48.2° |


The results obtained were matching with the expected results.

</div>
