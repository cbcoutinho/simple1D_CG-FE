# simple1D_CG-FE
A 1D implementation of a continuous Galerkin Model

This model is used to solve a 1D ODE:
  $ a u_x = cos(x) $
  
Which has the known solution of:
  $ u = sin(x)/a $

Testing the model with an odd number of elements produces the correct result, as shown here with 7 elements:

![Sample image](http://i.stack.imgur.com/Qo0ab.png)
