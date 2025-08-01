# gcd_nangate45 IO placement
source "helpers.tcl"
read_lef Nangate45/Nangate45.lef
read_def gcd.def

place_pins -hor_layers {metal3 metal5} -ver_layers {metal2 metal4} \
  -corner_avoidance 0 -min_distance 0.12

set def_file [make_result_file multi_layers.def]

write_def $def_file

diff_file multi_layers.defok $def_file
