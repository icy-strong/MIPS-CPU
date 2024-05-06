`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 03/11/2024 05:20:17 PM
// Design Name: 
// Module Name: progam_counter
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module progam_counter(nextPC, clk, pc);
    input wire clk;
    input wire [31:0] nextPC;
    output reg [31:0] pc;
    
    initial begin
        pc = 32'd100; //initalize pc to 100dec
    end
    
    
    always @(posedge clk) begin
        pc <= nextPC; //pc = nextPC
    end 

endmodule
