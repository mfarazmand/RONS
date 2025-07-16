# RONS Package:

This computational package implements RONS with examples from the following papers: <br/>
[1] W. Anderson, M. Farazmand. Fisher information and shape-morphing modes for solving the Fokker-Planck equation in higher dimensions. J. Appl. Math. Comput., vol. 467, pp. 128489, 2024  <br/> 
[2] W. Anderson, M. Farazmand. Fast and scalable computation of shape-morphing nonlinear solutions with application to evolutional neural networks. J. Comput. Phys., vol. 498, pp. 112649, 2024 <br/>
[3] W. Anderson, M. Farazmand. Evolution of nonlinear reduced-order solutions for PDEs with conserved quantities. SIAM J. on Scientific Computing, vol. 44, pp. A176-A197, 2022 <br/>


# License
You may use these codes as long as the above papers are appropriately cited.

# SRONS (RONS with symbolic computing)
To use the SRONS package, the user only needs to change: 
1) Input parameters in “Input_File.nb”
2) The PDE parameters and right-hand side in “General_Symbolic_EXP_f.m” or “General_Symbolic_f.m”.

Tip: If your solution uses Gaussian modes, edit General_Symbolic_EXP_f.m which is much faster. Otherwise use General_Symbolic_f.m.

Once the files are edited, running Input_File.nb will create a folder for the corresponding PDE which contains files to build the metric tensor, M, and RHS vector, f. The ODE RHS file is also output as a matlab file. Symbolic expressions for the inner products are also stored in a .mx file which can be opened in Mathematica.
Warning: Enforcement of conserved quantities is not yet implemented in the package, although they are manually added to the examples in “Examples_matlab” directory.

In the “Examples” folder we provide input files for the 1D bistable potential and Duffing oscillator described in [1] to serve as a template for users. Additionally, we provide Mathematica notebooks for the advection-diffusion example in [3], and a heat equation example to help familiarize users with Mathematica and RONS. 

# Example: 
The default Input_File.nb is set up for the Fokker-Planck (FP) equation corresponding to a 1D bistable potential. In order to change this to FP for the Duffing oscillator, one must take the following steps:
1. In Input_File.nb change “dim” to 2
2. In Input_File.nb change “pdename” to “Duffing” or other desired file name
3. In exec_mathematica/General_Symbolic_EXP_f.m change “params” to include all parameters in your PDE
4. In exec_mathematica/General_Symbolic_EXP_f.m change “Assumptions” to specify the type of parameters
5. In exec_mathematica/General_Symbolic_EXP_f.m change “ordernonlin” to the highest order of nonlinearity in the PDE. For instance, ordernonlin=1 for Fokker-Planck and ordernonlin =2 for Burgers equation.
6. In exec_mathematica/General_Symbolic_EXP_f.m under (*1D Bistable Potential*), change F0 to the RHS of Duffing.

# CRONS (RONS with collocation points)
Two examples are provided, the Kuramoto-Sivashinsky example from [2] and the 1D bistable potential from [1]. 
