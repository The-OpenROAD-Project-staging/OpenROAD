# resize to target_slew
source "helpers.tcl"
read_lef sky130hd/sky130hd.tlef
read_lef sky130hd/sky130hd_std_cell.lef
read_liberty sky130hd/sky130hd_tt.lib

read_verilog gain_buffering1_hier.v
link_design top -hier

create_clock -name clk -period 1.0 [get_ports clk]

set_dont_use {sky130_fd_sc_hd__probe_*
    sky130_fd_sc_hd__lpflow_*
    sky130_fd_sc_hd__clkdly*
	sky130_fd_sc_hd__dlygate*
	sky130_fd_sc_hd__dlymetal*
	sky130_fd_sc_hd__clkbuf_*
	sky130_fd_sc_hd__bufbuf_*}

set_dont_touch _53_

report_checks -fields {fanout}
#set_debug_level RSZ early_sizing 2
repair_design -pre_placement
report_checks -fields {fanout}
set verilog_file [make_result_file gain_buffering1_hier_out.v]
write_verilog $verilog_file
diff_files $verilog_file gain_buffering1_hier_out.vok
