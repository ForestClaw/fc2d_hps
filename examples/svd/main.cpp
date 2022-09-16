#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

#include <HPS/fc2d_hps.hpp>
#include <Structures/fc2d_hps_matrix.hpp>
#include <Structures/fc2d_hps_patchgrid.hpp>
#include <Structures/fc2d_hps_patchsolver.hpp>
#include <petsc.h>
#include <petscmat.h>
#include <slepc.h>
#include <slepcsvd.h>

extern "C" {
    void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);
    void dgesdd_(char* JOBZ, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO);
}

void SVDFromMat(Mat& mat) {

    MatType matType; MatGetType(mat, &matType);
    // Mat matDense;
    if (matType != MATDENSE) {
        std::cout << "Converting..." << std::endl;
        MatConvert(mat, MATDENSE, MAT_INPLACE_MATRIX, &mat);
    }

    char JOBU = 'A';
    char JOBVT = 'A';
    int M, N; MatGetSize(mat, &M, &N);
    int NSIGMA = std::min(M,N);
    int LDA = M;
    int LDU = M;
    int LDVT = N;
    int LWORK = -1;
    double* AMemory = (double*) malloc(M*N*sizeof(double)); MatDenseGetArray(mat, &AMemory);
    std::vector<double> A; A.assign(AMemory, AMemory + M*N);
    std::vector<double> S(NSIGMA);
    std::vector<double> U(LDU*M);
    std::vector<double> VT(LDVT*N);
    std::vector<double> WORK(1);
    int INFO;
    std::cout << "Calling dgesvd_..." << std::endl;
    // Call to dgesvd to get optimal size for WORK
    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, VT.data(), &LDVT, WORK.data(), &LWORK, &INFO);
    LWORK = WORK[0];
    WORK.resize(WORK[0]);
    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, VT.data(), &LDVT, WORK.data(), &LWORK, &INFO);
    if (INFO) {
        std::cerr << "Call to dgesvd_ returned non-zero value for INFO of: " << INFO << std::endl;
    }
    std::cout << "Done with dgesvd_" << std::endl;

    for (auto i = 0; i < NSIGMA; i++) {
        printf("%8.4f,  ", S[i]);
    }
    std::cout << std::endl;

    MatDenseRestoreArray(mat, &AMemory);
    WORK.clear();

}

int main(int argc, char** argv) {

    std::cout << "Hello from SVD!" << std::endl;

    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    SlepcInitialize(&argc, &argv, NULL, NULL);

    int nCells = 64;
    fc2d_hps_patchgrid grid(nCells, nCells, -1, 1, -1, 1);
    fc2d_hps_FISHPACK_solver solver;
    
    fc2d_hps_matrix<double> T = solver.build_dtn(grid);

    // int rows = 4;
    // int cols = 5;
    // int lda = cols;
    // std::vector<double> TData = {
    //     1, 0, 0, 0, 2,
    //     0, 0, 3, 0, 0,
    //     0, 0, 0, 0, 0,
    //     0, 2, 0, 0, 0
    // };
    // fc2d_hps_matrix<double> T(rows, cols, TData);

    // SingularValueDecomposition(rows, cols, lda, T.data());

    // for (auto i = 0; i < T.rows; i++) {
    //     for (auto j = 0; j < T.cols; j++) {
    //         printf("%8.4f,  ", T(i,j));
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<int> I(T.rows);
    std::vector<int> J(T.cols);
    std::iota(I.begin(), I.end(), 0);
    std::iota(J.begin(), J.end(), 0);

    // Create PETSc matrix
    Mat TPetsc;
    MatCreate(MPI_COMM_WORLD, &TPetsc);
    MatSetSizes(TPetsc, PETSC_DECIDE, PETSC_DECIDE, (PetscInt) T.rows, (PetscInt) T.cols);
    MatSetFromOptions(TPetsc);
    MatSetUp(TPetsc);
    // MatCreateSeqAIJ(MPI_COMM_WORLD, T.rows, T.cols, T.rows, NULL, &TPetsc);
    // MatSetUp(TPetsc);
    for (auto i = 0; i < T.rows; i++) {
        for (auto j = 0; j < T.cols; j++) {
            MatSetValue(TPetsc, i, j, T(i,j), INSERT_VALUES);
        }
    }
    // MatSetValues(TPetsc, T.rows, I.data(), T.cols, J.data(), T.data(), INSERT_VALUES);
    MatAssemblyBegin(TPetsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(TPetsc, MAT_FINAL_ASSEMBLY);

    SVDFromMat(TPetsc);

    // Create SVD object
    int nConv = -1;
    int iters = -1;
    double sigma;
    Mat A;
    Vec u, v;
    MatCreateVecs(TPetsc, &v, &u);
    // VecCreate(MPI_COMM_WORLD, &u);
    // VecCreate(MPI_COMM_WORLD, &v);

    SVD svd;
    SVDType svdType;
    int nsv = -1;
    int maxIter = -1;
    double tolerance = -1;
    SVDCreate(MPI_COMM_WORLD, &svd);
    SVDSetFromOptions(svd);
    // SVDSetOperators(svd, TPetsc, NULL);
    // SVDSetDimensions(svd, 1, PETSC_DEFAULT, PETSC_DEFAULT);
    // SVDSolve(svd);
    // SVDGetIterationNumber(svd, &iters);             std::cout << "iters = " << iters << std::endl;
    // SVDGetType(svd, &svdType);                      std::cout << "type = " << svdType << std::endl;
    // SVDGetDimensions(svd, &nsv, NULL, NULL);        std::cout << "nsv = " << nsv << std::endl;
    // SVDGetTolerances(svd, &tolerance, &maxIter);    std::cout << "tolerance = " << tolerance << std::endl; std::cout << "maxIter = " << maxIter << std::endl;
    // SVDGetConverged(svd, &nConv);                   std::cout << "nConv = " << nConv << std::endl;
    // SVDGetSingularTriplet(svd, nConv, &sigma, u, v);
    
    // Clean up
    MatDestroy(&TPetsc);
    VecDestroy(&u);
    VecDestroy(&v);
    SVDDestroy(&svd);

    SlepcFinalize();
    PetscFinalize();
    MPI_Finalize();

    return EXIT_SUCCESS;

}