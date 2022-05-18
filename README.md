# dtu_25314_marine_models_project2


Scripts
single_copepod:
	Is the model running a steady state solution with no seasons. This script generates
	Figure 3.
multi_copepods:
	Is the script to run the full model with seasons and multiple copepods with adaptive
	behaviour. You can change the amount of copepods yourself.
multi_copepods_random:
	Is the script to run the full model with seasons and multiple copepods with
	random movement. You can change the amount of copepods yourself.

Functions
call_param:
	Is the function to call all parameters in a structure.
func_diff:
	Is the function for the NPD model.
func_diff_season:
	Is the function for the NPD model with seasonal changes.
func_ligth:
	The function calculate light.
func_ligth_s:
	The function calculation light with seasonal changes.
