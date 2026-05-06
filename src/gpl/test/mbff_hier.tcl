source helpers.tcl

read_lef ./asap7/asap7_tech_1x_201209.lef
read_lef ./SingleBit/asap7sc7p5t_28_L_1x_220121a.lef
read_lib ./SingleBit/asap7sc7p5t_SEQ_LVT_TT_nldm_220123.lib
read_lef ./2BitTrayH2/asap7sc7p5t_DFFHQNV2X.lef
read_lib ./2BitTrayH2/asap7sc7p5t_DFFHQNV2X_LVT_TT_nldm_FAKE.lib
read_lef ./4BitTrayH4/asap7sc7p5t_DFFHQNV4X.lef
read_lib ./4BitTrayH4/asap7sc7p5t_DFFHQNV4X_LVT_TT_nldm_FAKE.lib

read_verilog ./mbff_hier.v
link_design -hier mbff_hier

set block [ord::get_db_block]
set locs {{a/ff_a1 6000 6000} {a/ff_a2 4000 6000}
          {b/ff_b1 4000 4000} {b/ff_b2 6000 4000}}
foreach e $locs {
  set inst [$block findInst [lindex $e 0]]
  $inst setLocation [lindex $e 1] [lindex $e 2]
  $inst setPlacementStatus PLACED
}

create_clock -name clk -period 1000 [get_ports clk1]

proc report_modules { tag } {
  set block [ord::get_db_block]
  puts "===== $tag ====="
  foreach mod [list [$block getTopModule] \
    [$block findModule sub_a] \
    [$block findModule sub_b]] {
    set names [list]
    foreach inst [$mod getInsts] { lappend names [$inst getName] }
    puts "[$mod getName]: [llength $names] inst(s): [lsort $names]"
  }
}

report_modules "BEFORE cluster_flops"

cluster_flops -tray_weight 40.0 -timing_weight 0.0 \
  -max_split_size -1 -num_paths 0

report_modules "AFTER cluster_flops"
