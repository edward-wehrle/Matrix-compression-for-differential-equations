# Matrix-compression-for-differential-equations
This repository shows off the (value, column, row pointer) matrix compression algorithm used for reducing the size of sparse matrices.

When modeling a 2-dimensional Laplacian differential equation, a weighted average of points is used to iterate through snapshots of a model.
That weighted average is computed via a matrix multiplication, one that scales exponentially as the size of the model increases. Due to the 
nature of how the weighted average is laid out, the matrix representing it is, "sparse", meaning it contains a large amount of zeros. During 
the typical multiplication process, many operations are wasted through trying to do a multiplication with these zeros. By retaining location
and value data of the non-zero entries, and implementing the compressed matrix into a solving algorithm, those wasted operations are 
removed. This saves a large amount of time, scaling up with the size of the model.

This code was done in C++ and consists of a Kroenecker multiplication function to establish the initial sparse matrix, the matrix compression
algorithm itself, and 4 separate matrix solving algorithms that utilize the compressed matrix. It is a student project done for a computational
science class, and represents my problem-solving abilities. It has a number of flaws that I did not notice or fix when I was coding it, but
it does serve as a display of my effort to pursue a project that pushed my boundaries.
