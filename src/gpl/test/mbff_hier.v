// Hierarchical netlist for MBFF clustering regression.
// Two child modules each containing two flip-flops.

module sub_a (clk, d1, d2, q1, q2);
  input  clk;
  input  d1, d2;
  output q1, q2;
  DFFHQNx1_ASAP7_75t_L ff_a1 (.D(d1), .CLK(clk), .QN(q1));
  DFFHQNx1_ASAP7_75t_L ff_a2 (.D(d2), .CLK(clk), .QN(q2));
endmodule

module sub_b (clk, d1, d2, q1, q2);
  input  clk;
  input  d1, d2;
  output q1, q2;
  DFFHQNx1_ASAP7_75t_L ff_b1 (.D(d1), .CLK(clk), .QN(q1));
  DFFHQNx1_ASAP7_75t_L ff_b2 (.D(d2), .CLK(clk), .QN(q2));
endmodule

module mbff_hier (clk1, d1, d2, d3, d4, o1, o2, o3, o4);
  input  clk1;
  input  d1, d2, d3, d4;
  output o1, o2, o3, o4;

  sub_a a (.clk(clk1), .d1(d1), .d2(d2), .q1(o1), .q2(o2));
  sub_b b (.clk(clk1), .d1(d3), .d2(d4), .q1(o3), .q2(o4));
endmodule
