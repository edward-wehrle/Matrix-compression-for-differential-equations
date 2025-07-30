#include<iostream>
#include<vector>
#include<cmath>
#include<chrono>
#include<algorithm>
#include<fstream>
using namespace std;


double f(double x, double y){
    double val = 2.0 * pow(M_PI,2) * sin(M_PI * x) * sin(M_PI * y);
    return val;
};

double f_exact(double x, double y){
    double val = sin(M_PI * x) * sin(M_PI * y);
    return val;
};

class SparseMatrix{

private:
    int nRows, nCols;
    vector<double> values;
    vector<int> colIndices;
    vector<int> rowPtr;

public:
    SparseMatrix(int rows, int cols):
    nRows(rows), nCols(cols) {
        rowPtr.resize(rows+1, 0); 
    }

    int rows() const {return nRows;}
    int cols() const {return nCols;}

    void addValue(int r, int c, double val){
        // How to add values into a sparse matrix
        values.push_back(val);
        colIndices.push_back(c);
        rowPtr[r+1]++; // rowPtr indexes the count of elements in
                       // a row and needs to be adjusted later.
    }

    void finalize(){
        for(int i = 1; i < (int)rowPtr.size(); i++){
            rowPtr[i] += rowPtr[i-1]; // Final adjustment for rowPtr
        }
    }

    vector<double> operator*(const vector<double> &x) const{
        // Overloading the * operator to return a vector
        vector<double> y(nRows, 0.0);

        for(int i = 0; i < nRows; i++){
            int start = rowPtr[i];
            int end = rowPtr[i+1];
            double sum = 0.0;

            for(int idx = start; idx < end; idx++){
                // Calculates the dot product of A and x per row
                sum+= values[idx] * x[colIndices[idx]];
            }
            y[i] = sum;
            }
            return y;
        }

    #pragma region Friend functions

    friend void JacobiSolver(const SparseMatrix &A, const vector<double> &b,
        vector<double> &u, int maxIters, double tol){

        int n = A.rows();
        vector<double> diag(n,0.0);

        for(int i = 0; i < n; ++i){
            int start = A.rowPtr[i];
            int end = A.rowPtr[i+1];

            for(int idx = start; idx < end; idx++){
            // Iterate for the number of elements in a row,
            // this statement comes up numerous times in the code

                if(A.colIndices[idx] == i){ // Checks the diagonal by seeing if the column = row for the index
                    diag[i] = A.values[idx]; // Sets the diagonal
                    break;
                }
            }
        }

        vector<double> uOld(n,0.0);

        for(int iter = 0; iter < maxIters; iter++){
            uOld = u;
            double maxDiff = 0.0;

            for(int i = 0; i < n; i++){
                double rowSum = 0.0;
                int start = A.rowPtr[i];
                int end = A.rowPtr[i+1];

                for(int idx = start; idx < end; idx++){
                    int j = A.colIndices[idx];
                    if(j != i){
                        rowSum += A.values[idx] * uOld[j];
                        // Calculates the dot product of A and uOld without
                        // using any values where j = i
                    }
                }
                u[i] = (b[i] - rowSum) / diag[i]; // Adjust u[i]

                maxDiff = max(maxDiff, fabs(u[i] - uOld[i])); // Find the max change in any
                                                              // element of u and record it,
                                                              // if it is under tolerance then
                                                              // break the program and return u
            }

            if(maxDiff < tol){
                cout << "[Jacobi] Converged at iteration " << iter+1 << "\n";
                return;
            }
        }
        cout << "[Jacobi] did not converge. \n";
    };



    friend void GaussSeidelSolver(const SparseMatrix &A, const vector<double> &b,
                                  vector<double> &u, int maxIters, double tol){
        // Beginning section of code is identical to Jacobi,
        // it establishes the diagonal vector of A 
        int n  = A.rows();
        vector<double> diag(n,0.0);
                                    
        for(int i = 0; i < n; ++i){
            int start = A.rowPtr[i];
            int end = A.rowPtr[i+1];
                        
            for(int idx = start; idx < end; idx++){
                if(A.colIndices[idx] == i){
                    diag[i] = A.values[idx];
                    break;
                }
            }
        }

        vector<double> uOld(n,0.0);

        for (int iter = 0; iter < maxIters; ++iter){
            uOld = u;
            double maxDiff = 0.0;

            for(int i = 0; i < n; i++){
                double rowSumOld = 0.0;
                double rowSumNew = 0.0;
                
                int start = A.rowPtr[i];
                int end = A.rowPtr[i+1];

                for(int idx = start; idx < end; idx++){
                    int j = A.colIndices[idx];
                    
                    // Calculate the two needed dot products,
                    // one of which uses updated values of u[j]
                    // to increase accuracy
                    if(j < i){
                        rowSumNew += A.values[idx] * u[j];
                    }
                    if(j > i){
                        rowSumOld += A.values[idx] * uOld[j];
                    }
                }

                // Update u[i] and the maxDiff
                u[i] = (b[i] - rowSumOld - rowSumNew) / diag[i];
                maxDiff = max(maxDiff, fabs(u[i] - uOld[i]));
            }
            if(maxDiff < tol){
                cout << "[Gauss-Seidel] Converged at iteration " << iter+1 << "\n";
                return;
            }
        }
        cout << "[Gauss-Seidel] did not converge. \n";
    };

    friend void SorSolver(const SparseMatrix &A, const vector<double> &b,
        vector<double> &u, int maxIters, double tol, double omega){
        // Reference variable for iter so I can make a vector
        // listing # of iterations with the value of omega.
        //
        // This code is identical to the Gauss-Seidel code except for
        // including omega and using a slightly adjusted formula to update u[i]

        int n = A.rows();
        vector<double> diag(n, 0.0);
                                                        
        for(int i = 0; i < n; ++i){
            int start = A.rowPtr[i];
            int end = A.rowPtr[i+1];
                                            
            for(int idx = start; idx < end; idx++){
                if(A.colIndices[idx] == i){
                    diag[i] = A.values[idx];
                        break;
                }
            }
        }
                    
        vector<double> uOld(n,0.0);
                    
        for (int iter = 0; iter < maxIters; ++iter){
            uOld = u;
            double maxDiff = 0.0;
    
            for(int i = 0; i < n; i++){
                double rowSumOld = 0.0;
                double rowSumNew = 0.0;
                    
                int start = A.rowPtr[i];
                int end = A.rowPtr[i+1];
    
                for(int idx = start; idx < end; idx++){
                    int j = A.colIndices[idx];
                    if(j < i){
                        rowSumNew += A.values[idx] * u[j];
                    }
                    if(j > i){
                        rowSumOld += A.values[idx] * uOld[j];
                    }
                }
                // Below, u[i] is adjusted by a factor of omega. In testing,
                // multiple omegas were used to find which converged in the
                // fewest iterations, the time for that one was reported and 
                // the number of iterations needed as a factor of omega was
                // also reported

                u[i] = (b[i] - rowSumOld - rowSumNew) * (omega/diag[i]) + (1-omega)*uOld[i];
                maxDiff = max(maxDiff, fabs(u[i] - uOld[i]));
            }
            if(maxDiff < tol){
                cout << "[SOR] Converged at iteration " << iter+1 << "\n";
                return;
            }
        }
        cout << "[SOR] did not converge. \n";
    };

    friend void ConjugateGradientSolver(const SparseMatrix &A, const vector<double> &b,
                                        vector<double> &u, int maxIters, double tol){
        // The first thing I did was establish all the necessary vectors and variables
        // needed for the conjugate gradient method.

        int n = A.rows();
        double alpha;
        double beta;
        double maxDiff = 0.0;
        vector<double> uOld(n, 0.0);
        vector<double> r(n, 0.0);
        vector<double> r_new(n, 0.0);
        vector<double> p(n, 0.0);
        vector<double> Ap(n, 0.0);
        
        
        r = (A * u);

        // As shown below, when adding and subtracting vectors,
        // I implemented a loop that ran over each element and 
        // performed the necessary operations.
        // 
        // Wherever possible I iterated over multiple vectors 
        // in one loop to save time.

        for(int i = 0; i < n; i++){
            r[i] = r[i] - b[i];
        }
        p = r;
        // p and r vectors are initialized

        for(int iter = 0; iter < maxIters; iter++){
            maxDiff = 0.0;
            uOld = u; 
            Ap = A * p; // First use of the overloaded * operator

            double sumrtr = 0.0;
            double sumptAp = 0.0;
            // Both sums need to be set to 0 before calculating alpha.
            // I forgot this initially and had them outside the for loop.
            // As a result, alpha (and therefore the rest of the code)
            // were all off, and the program would not converge

            for(int i = 0; i < n; i++){
                // I needed to find r.T @ r and p.T @ Ap, which are both single values
                // and not vectors. I didn't need to transpose any vectors because I
                // performed all operations element-wise in the same way they would be computed
                // using a dot product operation in python.
                sumrtr += r[i] * r[i];
                sumptAp += p[i] * Ap[i];
            }
            alpha = sumrtr / sumptAp; // Establish alpha, the step size

            for(int i = 0; i < n; i++){
                u[i] = u[i] - (alpha*p[i]);
                r_new[i] = r[i] - (Ap[i] * alpha);
                // Here I updated the u and r_new vectors in the same loop
                // to save time
                maxDiff = max(maxDiff, fabs(u[i] - uOld[i]));
            }

            if(maxDiff < tol){
                cout << "[Conjugate Gradient] converged at iteration " << iter+1 << endl;
                // This is the one tolerance check that doesn't occur at the end of the loop.
                // It can be performed before updating beta, p, and r,
                // which does save some time.
                return;
            }

            double sum1 = 0;
            double sum2 = 0;
            for(int i = 0; i < n; i++){
                sum1 += r_new[i] * r_new[i];
                sum2 += r[i] * r[i];
            }
            beta = sum1 / sum2;
            // beta has to be updated before the next loop begins,
            // otherwise I would have combined them.

            for(int i = 0; i < n; i++){
                p[i] = r_new[i] + (beta * p[i]);
                r[i] = r_new[i];
            }
        }
        // This function was the most fun to code out of all 4
        // B)
        cout << "Conjugate Gradient did not converge." << endl;
    };

    #pragma endregion // I added some regions to make editing code easier,
                      // because I'm not very familiar with coding, one of
                      // the things I needed to learn was how to make scrolling
                      // through and editing code as easy as possible.
    };


vector<int> Kroneckerproduct(vector<vector<int>> A, vector<vector<int>> B, int N)
{
    // A function to find the Kronecker product, I used this to implement the dense
    // A matrix

    // Square matrices ONLY!
    vector<int> x(N*N*N*N,0); // setup a flat vector with length N^4
    int p = 0; // total number of iterations

    // i loops
    for (int i = 0; i < N; i++) {
 
        // k loops
        for (int k = 0; k < N; k++) {
 
            // j loops
            for (int j = 0; j < N; j++) {
 
                // l loops
                for (int l = 0; l < N; l++) {
 
                    // Each element of matrix A is
                    // multiplied by whole Matrix B
                    // and stored in the flat vector x.
                    //
                    // Base Kroenecker code was borrowed from online
                    // and only printed the Kroenecker matrix 
                    // by printing each individual element in the terminal.
                    //
                    // The resulting C matrix updated over the same elements
                    // multiple times, so trying to return it to use in the main code
                    // gave an incorrect matrix.
                    //
                    // I spent a lot of time trying to adjust the code to return
                    // a useful matrix before I realized that each individual element
                    // printed out in the correct order. I had a eureka moment and
                    // put the values into a flat vector, which I returned and adjusted
                    // to fit my needs.
                    // Figuring this out took much more time than I'm comfortable to admit haha
                    x.at(p) = A[i][j] * B[k][l];
                    // This is the only time I remembered I could use the .at function for the 
                    // vectors. I was using x[p] but only established N^2 spaces for memory
                    // and kept getting errors until I realized.
                    ++p;
                }
            }
        }
    }
    return x;
}

int main(){
    
    int N = 100;
    double h = 1.0 / N;
    int idx = 0;
    int max_its = 100000;
    double tol = .0000001; // 1e-7 tolerance, chosen arbitrarily
    double omega = 1.9;
    int iter = 0;
    int N2 = N*N;
    int N4 = N*N*N*N;

    // I allocate memory for many vectors and the A matrix at once
    SparseMatrix Ainit(N2,N2);
    vector<double> x_grid(N+1,0); // N+1 points from 0 to N
    vector<double> y_grid(N+1,0);
    vector<double> u(pow(N+1,2),0);
    vector<double> u_exact(pow(N+1,2),0);
    vector<double> b(pow(N+1,2),0);

    for (int i = 0; i < N+1; ++i){
        double temp = h * i;
        x_grid[i] = temp;
        y_grid[i] = temp;
    }
    

    int iters = 0;
    for(int i = 0; i < N+1; ++i){
        for(int j = 0; j < N+1; ++j){ 
                if(!(i == 0 && j == 0 && i == N+1 && j == N+1)){ // Checking for interior points
                    u[iters] = f(x_grid[i], y_grid[j]); // Establish the initial u vector
                    u_exact[iters] = f_exact(x_grid[i], y_grid[j]); // Establish the u_exact vector
                }
                ++iters;
        }
    }

    
    // T and I are initialized as NxN matrices with all zeros
    vector<vector<int>> T(N, vector<int>(N));
    vector<vector<int>> I(N, vector<int>(N));

    //vector<vector<int>> A(N2, vector<int>(N2));  This is an old line I cut

    vector<int> c1(N4,0); // c1 and c2 are flat vectors that will be used in intermediate kroenecker calculations
    vector<int> c2(N4,0);
    vector<int> cflat(N4,0);
 
    for (int i = 0; i < N; ++i){
        I[i][i] = 1; // Establish the identity matrix
        T[i][i] = 2; // Establish the T matrix diagonal
        if(i>=1){
            T[i-1][i] = -1; 
            T[i][i-1] = -1; // sub-diagonals above and below main diagonal
        }
    }
    
    c1 = Kroneckerproduct(I,T,N); // Here's where I implement the Kronecker function
    c2 = Kroneckerproduct(T,I,N); // I mentioned earlier.
    
    for(int i = 0; i < N4; ++i){
        cflat[i] = (c1[i]+c2[i]) /(h*h); // The two flat vectors are combined to get A flattened

        // Now I finally can get my compressed matrix
        if(cflat[i] != 0){ 
            Ainit.addValue(i / N2, i % N2, cflat[i]); // i / N2 truncates down to the row #, i % N2 gives the col
        }
    }

    Ainit.finalize(); // finish the rowptr vector
    
    
    /* Code cut for efficiency
    for(int i = 0; i < N2; ++i){
        for(int j = 0; j < (N2); ++j){
            A[i][j] = cflat[idx]; // Turn A into a full N^2 x N^2 matrix, but I realized
                                  // that I could use the flat vector to get Ainit
                                  // and this step was unnecessary. But I kept this
                                  // in to show my work and what I learned.
            if(A[i][j] != 0){
            Ainit.addValue(i, j, A[i][j]);
            }
        }
    }
    End of the cut code*/ 
    
    // The u vector in each function is a non-constant reference, so I'm assigning space
    // for each function to have their own u vector, which will also make it easier to export 
    // and categorize the data to make my graphs

    int length = pow(N+1,2); // At first I did this calculation in each line, but then
                             // I just added this variable instead.

    vector<double> u_Jacobi(length,0);
    vector<double> u_Gauss_Seidel(length,0);
    vector<double> u_SOR(length,0);
    vector<double> u_Conjugate_Gradient(length,0);
    vector<double> b_exact(length,0);
    b_exact = Ainit*u_exact; // establish the solution b
    u_Jacobi = u; // Copy u over to each vector
    u_Gauss_Seidel = u;
    u_SOR = u;
    u_Conjugate_Gradient = u;

    // Vector starts at omega = 0.1 and iterates by 0.1. Ideal is 1.9
    vector<int> SOR_iters = {-1, 89048, 58778, 42931, 33090, 26343, 21405, 17620, 14616, 
        12165, 10123, 8388, 6893, 5586, 4429, 3392, 2451, 1583, 745, -1, 4222};
    

    /* Code I used to find the ideal omega for SOR

    for(int i = 0; i < 22; i++){
    omega = (i+1) * 0.1;
    u_SOR = u; // Don't forget to reset this between loops
    SorSolver(Ainit, b_exact, u_SOR, max_its, tol, omega);
    } 
    */


    ofstream myfile; 

    
    /*myfile.open("u_exact.txt");
    for(int i = 0; i < u_exact.size(); i++){
        myfile << u_exact[i] << " ";
    }
    myfile.close();
    */
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    JacobiSolver(Ainit, b_exact, u_Jacobi, max_its, tol);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    
    myfile.open("Jacobi.txt");
    for(int i = 0; i < u_exact.size(); i++){
        myfile << u_Jacobi[i] << " ";
    }
    myfile.close();
    
    
    chrono::steady_clock::time_point begin2 = chrono::steady_clock::now();
    GaussSeidelSolver(Ainit, b_exact, u_Gauss_Seidel, max_its, tol);
    chrono::steady_clock::time_point end2 = chrono::steady_clock::now();
    /*
    myfile.open("GaussSeidel.txt");
    for(int i = 0; i < u_exact.size(); i++){
        myfile << u_Gauss_Seidel[i] << " ";
    }
    myfile.close();
    */
    
    chrono::steady_clock::time_point begin3 = chrono::steady_clock::now();
    SorSolver(Ainit, b_exact, u_SOR, max_its, tol, omega);
    chrono::steady_clock::time_point end3 = chrono::steady_clock::now();
     /*
    myfile.open("SOR_data.txt");
    for(int i = 0; i < u_exact.size(); i++){
        myfile << u_SOR[i] << " ";
    }
    myfile.close();
    */
    
    chrono::steady_clock::time_point begin4 = chrono::steady_clock::now();
    ConjugateGradientSolver(Ainit, b_exact, u_Conjugate_Gradient, max_its, tol);
    chrono::steady_clock::time_point end4 = chrono::steady_clock::now();
        /*
    myfile.open("Conjugate_Gradient.txt");
    for(int i = 0; i < u_exact.size(); i++){
        myfile << u_Conjugate_Gradient[i] << " ";
    }
    myfile.close();
    */

    /*myfile.open("Iterations.txt");
    myfile << 22895 << " " << 12165 << " " << 745 << " " << 289;
    myfile.close();
    */

    /*
    myfile.open("Times.txt");
    myfile << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " ";
    myfile << chrono::duration_cast<chrono::milliseconds>(end2 - begin2).count() << " ";
    myfile << chrono::duration_cast<chrono::milliseconds>(end3 - begin3).count() << " ";
    myfile << chrono::duration_cast<chrono::milliseconds>(end4 - begin4).count();
    myfile.close();
    */
   
    return 0;
}