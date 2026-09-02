# Ranking slots by wirelength discards a pin placement the caller has already
# optimized, which is what global_placement -place_ios leaves behind.
# -minimize_displacement ranks them by distance instead, so the assignment
# legalizes that placement rather than replacing it.
source "helpers.tcl"
read_lef Nangate45/Nangate45.lef
read_def gcd.def

proc pin_centers { } {
  set centers {}
  foreach bterm [[ord::get_db_block] getBTerms] {
    set box [$bterm getBBox]
    dict set centers [$bterm getName] \
      [list [expr { ([$box xMin] + [$box xMax]) / 2 }] \
        [expr { ([$box yMin] + [$box yMax]) / 2 }]]
  }
  return $centers
}

proc report_move { what from } {
  set to [pin_centers]
  set moved 0
  set max_d 0
  dict for { name pos } $from {
    lassign $pos x0 y0
    lassign [dict get $to $name] x1 y1
    set d [expr { abs($x1 - $x0) + abs($y1 - $y0) }]
    if { $d > 0 } {
      incr moved
      set max_d [expr { max($max_d, $d) }]
    }
  }
  puts "$what: $moved of [dict size $from] pins moved, at most $max_d dbu"
}

place_pins -hor_layers metal3 -ver_layers metal2
set placed [pin_centers]

# Mirror the cell placement across the core. Every pin's wirelength-optimal
# edge is now the opposite one, while the pins themselves have not moved.
set block [ord::get_db_block]
set core [$block getCoreArea]
foreach inst [$block getInsts] {
  lassign [$inst getLocation] x y
  $inst setLocation [expr {
    [$core xMin] + [$core xMax] - $x - [[$inst getMaster] getWidth]
  }] $y
}

place_pins -hor_layers metal3 -ver_layers metal2 -minimize_displacement
report_move "minimize_displacement" $placed

# The switch is per call, so this one ranks by wirelength again and follows
# the cells across the die.
place_pins -hor_layers metal3 -ver_layers metal2
report_move "wirelength" $placed

# Pins with no shape have no location to keep, which has to be said rather
# than measured from the die origin in silence.
foreach bterm [$block getBTerms] {
  foreach bpin [$bterm getBPins] {
    odb::dbBPin_destroy $bpin
  }
}
place_pins -hor_layers metal3 -ver_layers metal2 -minimize_displacement
