# Spp
Welcome to S++. This is a port of the R/numpy matrix notation (and a few datatype macros) to C++. You can index arrays much more naturally and develoopment speeds up 
rather quickly. This whole system is built on pointers so its fairly fast. It needs some improvement (which I will get to) and the documentation will be uploaded but 
that's going to take a long time so you'll just need to read the header file until I get there.

There are a few files. I need to merge them and clear up the code base so I would use S++2.h but really I would read over them and take whichever you like best. 
  There are some differences like row vs column major order matricies in each file for speed for different projects.
  
Extra Functions:
Tic and Toc. This replicates the tic and toc functions from tictoc R package. Call tic and pass string name as an argument. Call toc and the time since tic will 
      be printed along with the string name passed. Very helpful for opimization testing. May place into its own repo but feels a bit like a waste.
 
Data Types:
num - shorthand for double 
pos - shorthand for unsigned int 
Vector - an array with attribute len
Matrix - a 2D array with attribute shape which returns a 2 element array of (row, col)


Vector Details:

Methods and attributes:
.len will return an unsigned int with the length of the array 


Operators:

Vector - num : does an elementwise subtraction and returns a Vector as the result
num - Vector : again, elementwise subtraction
-Vector : returns the vector with all elements signed flipped 
Vector - Vector : returns error if incompatible legnths. Otherwise, Vector with elementwise subtraction 
num < Vector & num > Vector : returns array with same length as the Vector with a 1 in places where equality is met and 0 elsewhere 
num * Vector & Vector * num & Vector * Vector: returns an Vector with elementwise multiplication as the result
num / Vector & Vector / num : returns a Vector with the elementwise division as the result 
num == Vector & Vector == num : returns a Vector with 1 and 0 in spots where equality is met
Vector == Vector : returns true or false if the vectors are equal after elementwise comparison


