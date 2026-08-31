# -place_ios with -timing_driven. repair_design creates and destroys
# instances mid-solve, which permutes gCellStor_ and reallocates gPinStor_,
# while the IO pin GCells live in their own storage: every pin must still own
# its GPin afterwards and end up on the die perimeter.

source helpers.tcl
read_liberty ./library/nangate45/NangateOpenCellLibrary_typical.lib

read_lef ./nangate45.lef
read_def ./simple01-td-repair.def

create_clock -name core_clock -period 0.5 clk

set_wire_rc -signal -layer metal3
set_wire_rc -clock -layer metal5

# The def places the ports FIXED, which -place_ios keeps as anchors.
set block [ord::get_db_block]
foreach bterm [$block getBTerms] {
  foreach bpin [$bterm getBPins] {
    $bpin setPlacementStatus PLACED
  }
}
set insts_before [llength [$block getInsts]]

# Configuring the pin placer is what lets the solve legalize its own pins
# mid-run and before each timing-driven pass, so rsz and STA read shapes on
# real tracks.
set_io_pin_placement -hor_layers metal5 -ver_layers metal6
global_placement -timing_driven -timing_driven_repair_timing -place_ios \
  -place_ios_legalize_every 20

puts "instance count: $insts_before -> [llength [$block getInsts]]"

# A pin that lost its GPin in the churn would stop moving, so check that every
# one of them ends on the perimeter. The solve now hands the pins to the pin
# placer before it returns, so the shape sits just inside the die edge rather
# than centred on it: allow the shape's own extent.
set die [$block getDieArea]
foreach bterm [$block getBTerms] {
  set name [$bterm getName]
  set bpins [$bterm getBPins]
  if { [llength $bpins] != 1 } {
    error "$name has [llength $bpins] bpins, expected exactly 1"
  }
  set bpin [lindex $bpins 0]
  if { [$bpin getPlacementStatus] ne "PLACED" } {
    error "$name is [$bpin getPlacementStatus], expected PLACED"
  }
  set box [lindex [$bpin getBoxes] 0]
  set slack [expr { max([$box xMax] - [$box xMin], [$box yMax] - [$box yMin]) }]
  set to_edge [expr {
    min([$box xMin] - [$die xMin], [$die xMax] - [$box xMax],
      [$box yMin] - [$die yMin], [$die yMax] - [$box yMax])
  }]
  if { $to_edge > $slack } {
    error "$name is $to_edge from the die edge, more than its own $slack"
  }
}

estimate_parasitics -placement
report_worst_slack

puts pass
