# -place_ios with nothing breaking the perimeter, so the slots form one cyclic
# arc. The fit rounds against the point the cycle is cut at, and that point is
# a pin position rather than a slot: rounding to it once held every pin on the
# arc a fraction of a slot off the track grid. place_ios01 excludes an edge, so
# its arcs are open and it never reaches this.

source helpers.tcl
set test_name place_ios_cyclic
read_lef ./nangate45.lef
read_def ./simple01.def

set block [ord::get_db_block]
foreach bterm [$block getBTerms] {
  foreach bpin [$bterm getBPins] {
    odb::dbBPin_destroy $bpin
  }
}

set_io_pin_placement -hor_layers metal5 -ver_layers metal6
global_placement -place_ios -init_density_penalty 0.01

foreach bterm [$block getBTerms] {
  set box [lindex [[lindex [$bterm getBPins] 0] getBoxes] 0]
  set layer [$box getTechLayer]
  set grid [$block findTrackGrid $layer]
  if { [$layer getDirection] eq "HORIZONTAL" } {
    set coord [expr { ([$box yMin] + [$box yMax]) / 2 }]
    set tracks [$grid getGridY]
  } else {
    set coord [expr { ([$box xMin] + [$box xMax]) / 2 }]
    set tracks [$grid getGridX]
  }
  if { [lsearch -exact -integer $tracks $coord] < 0 } {
    error "[$bterm getName] is at $coord on [$layer getName], not on a track"
  }
}

puts pass
