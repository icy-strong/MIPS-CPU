`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: PSU CMPEN 331 SP24 Section 02
// Engineer: Jenicy Strong
// 
// Create Date: 03/11/2024 06:25:19 PM
// Design Name: cpu testbench
// Project Name: Lab 3
// 
//////////////////////////////////////////////////////////////////////////////////


module testbench();
    reg clk;
    wire [31:0] pc, eqa, eqb, eimm32, distOut, wr, wdo, mr, mqb;
    wire [4:0] edestReg, wdestReg, mdestReg;
    wire [3:0] ealuc;
    wire ewreg, em2reg, ealuimm, ewmem, wwreg, wm2reg, mwreg, mm2reg, mwmem;
    
    parameter PERIOD = 10;
    
    datapath DATAPATH (clk, pc, distOut, ewreg, em2reg, ewmem, ealuc, ealuimm, edestReg, eqa, eqb, eimm32, wwreg, wm2reg, wdestReg, wr, wdo,
                    mwreg, mm2reg, mwmem, mdestReg, mr, mqb);
    
    initial begin
        clk = 0;
    end
    
    always
        #(PERIOD / 2) clk = ~clk; 
    
    
    
endmodule
