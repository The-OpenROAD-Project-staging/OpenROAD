# check route guides for gcd_nangate45. def file from the openroad-flow
source "helpers.tcl"
read_lef "Nangate45/Nangate45.lef"
read_def "gcd.def"

set_global_routing_layer_adjustment metal2 0.7
set_global_routing_layer_adjustment metal3-metal10 1.0

set_routing_layers -signal metal2-metal10

catch { global_route -verbose } error
puts $error
