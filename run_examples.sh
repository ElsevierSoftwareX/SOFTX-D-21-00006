pip install .

optimic main --size 15. 15. 15. --dimension 3 --number_seed 500 --target ./optimic/user_grain_size_distribution.txt --characteristic 0 --material example_1_3D --stress_direction 1 0 0 --seed_spacing random_3d --rand_seed 1 --optimization_method COBYLA --max_iter 2 --number_bins 20 --mesh hex

optimic main --size 180 120 15 --dimension 2 --number_seed 122 --target ./optimic/user_grain_size_distribution.txt --characteristic 2 --characteristic 5 --material example_2_2D_10k_so --stress_direction 1 0 0 --sharp_orientation 1 1 1 --seed_spacing random_3d --rand_seed 2 --optimization_method COBYLA --max_iter 10 --number_bins 10 --user_cost_func ./optimic/util/user_cost_function_example2.py --mesh hex --mesh_size 1.0

optimic main --size 10. 10. 1. --spacing_length 1 --dimension 2 --target ./optimic/junc_angle_80160_sc_17.txt --characteristic 4 --material example_3_2D_hcp_longer --stress_direction 1 0 0 --sharp_orientation 1 1 1 --seed_spacing hcp_2d --rand_seed 3 --optimization_method COBYLA --max_iter 7 --number_bins 15 --user_cost_func ./optimic/util/user_cost_function_example3.py --mesh vis --mesh_size 0.2