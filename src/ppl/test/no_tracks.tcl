source "helpers.tcl"
# design without routing tracks
read_lef Nangate45/Nangate45.lef
read_lib Nangate45/Nangate45_typ.lib
read_verilog ../../../test/gcd_nangate45.v
link_design gcd

initialize_floorplan -site FreePDK45_38x28_10R_NP_162NW_34O \
  -utilization 30 -core_space 0.0

catch { place_pins -hor_layers metal3 -ver_layers metal2 } error
puts $error
