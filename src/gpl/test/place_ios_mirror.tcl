# -place_ios with mirrored pins. A follower's position comes from reflecting
# its master, so the spacing projection never sees it: whatever edge the
# reflection lands on already holds pins the projection did space, and the two
# groups know nothing about each other. Constrain a dozen pins to the right
# edge and let the mirrored followers reflect onto it, then check that every
# edge still holds its pins a slot pitch apart.

source helpers.tcl
read_lef ./nangate45.lef
read_def ./simple01.def

set block [ord::get_db_block]
foreach bterm [$block getBTerms] {
  foreach bpin [$bterm getBPins] {
    odb::dbBPin_destroy $bpin
  }
}

set pairs {}
for { set i 0 } { $i < 8 } { incr i } {
  lappend pairs "req_msg\[$i\]" "resp_msg\[$i\]"
}
set_io_pin_constraint -mirrored_pins $pairs

set direct {}
for { set i 10 } { $i < 22 } { incr i } {
  lappend direct "req_msg\[$i\]"
}
set_io_pin_constraint -region right:* -pin_names $direct

set_io_pin_placement -hor_layers metal5 -ver_layers metal6
global_placement -place_ios -init_density_penalty 0.01

# metal5/metal6 both step 560, and the default spacing is two tracks.
set pitch 1120
set die [$block getDieArea]
array set pos {}
foreach bterm [$block getBTerms] {
  set bpins [$bterm getBPins]
  if { [llength $bpins] != 1 } {
    error "[$bterm getName] has [llength $bpins] bpins, expected exactly 1"
  }
  set box [lindex [[lindex $bpins 0] getBoxes] 0]
  set cx [expr { ([$box xMin] + [$box xMax]) / 2 }]
  set cy [expr { ([$box yMin] + [$box yMax]) / 2 }]
  if { $cx == [$die xMin] } {
    lappend pos(left) $cy
  } elseif { $cx == [$die xMax] } {
    lappend pos(right) $cy
  } elseif { $cy == [$die yMin] } {
    lappend pos(bottom) $cx
  } elseif { $cy == [$die yMax] } {
    lappend pos(top) $cx
  } else {
    error "[$bterm getName] is at ($cx $cy), off the die perimeter"
  }
}

foreach edge { left right bottom top } {
  if { ![info exists pos($edge)] } {
    continue
  }
  set v [lsort -integer $pos($edge)]
  for { set i 1 } { $i < [llength $v] } { incr i } {
    set gap [expr { [lindex $v $i] - [lindex $v [expr { $i - 1 }]] }]
    if { $gap < $pitch } {
      error "two pins on the $edge edge are $gap apart, less than the\
             $pitch slot pitch"
    }
  }
}

# The mirrored pins must still be exact reflections of one another.
foreach { master follower } $pairs {
  set mb [lindex [[lindex [[$block findBTerm $master] getBPins] 0] getBoxes] 0]
  set fb [lindex [[lindex [[$block findBTerm $follower] getBPins] 0] getBoxes] 0]
  set mx [expr { ([$mb xMin] + [$mb xMax]) / 2 }]
  set my [expr { ([$mb yMin] + [$mb yMax]) / 2 }]
  set fx [expr { ([$fb xMin] + [$fb xMax]) / 2 }]
  set fy [expr { ([$fb yMin] + [$fb yMax]) / 2 }]
  if { $mx == [$die xMin] || $mx == [$die xMax] } {
    set want [expr { $mx == [$die xMin] ? [$die xMax] : [$die xMin] }]
    if { $fx != $want || $fy != $my } {
      error "$follower is at ($fx $fy), not the reflection of $master\
             at ($mx $my)"
    }
  } else {
    set want [expr { $my == [$die yMin] ? [$die yMax] : [$die yMin] }]
    if { $fy != $want || $fx != $mx } {
      error "$follower is at ($fx $fy), not the reflection of $master\
             at ($mx $my)"
    }
  }
}

puts pass
