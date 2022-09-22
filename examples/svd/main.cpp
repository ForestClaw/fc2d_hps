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
// #include <slepc.h>
// #include <slepcsvd.h>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

extern "C" {
    void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);
    void dgesdd_(char* JOBZ, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO);
}

std::vector<double> SVDValuesFromMat(Mat& mat) {

    MatType matType; MatGetType(mat, &matType);
    // Mat matDense;
    if (matType != MATDENSE) {
        // std::cout << "Converting..." << std::endl;
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
    // std::cout << "Calling dgesvd_..." << std::endl;
    // Call to dgesvd to get optimal size for WORK
    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, VT.data(), &LDVT, WORK.data(), &LWORK, &INFO);
    LWORK = WORK[0];
    WORK.resize(WORK[0]);
    dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, VT.data(), &LDVT, WORK.data(), &LWORK, &INFO);
    if (INFO) {
        std::cerr << "Call to dgesvd_ returned non-zero value for INFO of: " << INFO << std::endl;
    }
    // std::cout << "Done with dgesvd_" << std::endl;

    // for (auto i = 0; i < NSIGMA; i++) {
    //     printf("%8.4f,  ", S[i]);
    // }
    // std::cout << std::endl;

    MatDenseRestoreArray(mat, &AMemory);
    WORK.clear();

    return S;

}

int main(int argc, char** argv) {

    std::cout << "Hello from SVD!" << std::endl;

    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    // SlepcInitialize(&argc, &argv, NULL, NULL);

    fc2d_hps_FISHPACK_solver solver;
    std::vector<int> nCellsVector({2, 4, 8, 16, 32, 64, 128});
    std::vector<double> conditionNumbers(nCellsVector.size());
    for (auto i = 0; i < nCellsVector.size(); i++) {
        int nCells = nCellsVector[i];
        fc2d_hps_patchgrid grid(nCells, nCells, -1, 1, -1, 1);
        
        fc2d_hps_matrix<double> T = solver.build_dtn(grid);
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
        MatSetValues(TPetsc, T.rows, I.data(), T.cols, J.data(), T.data(), INSERT_VALUES);
        MatAssemblyBegin(TPetsc, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(TPetsc, MAT_FINAL_ASSEMBLY);

        std::vector<double> sigmas = SVDValuesFromMat(TPetsc);

        double conditionNumber = sigmas[0] / sigmas[sigmas.size()-1];
        conditionNumbers[i] = conditionNumber;

        printf("N = %4i, k = %16.8E\n", nCells, conditionNumber);

        MatDestroy(&TPetsc);
    }

    plt::figure();
    plt::loglog(nCellsVector, conditionNumbers);
    plt::title("Condition Number of DtN Matrix");
    plt::xlabel("Cells per Leaf Side");
    plt::ylabel("Condition Number");
    plt::save("condition_number_dtn.pdf");
    plt::show();
    
    // Clean up

    // SlepcFinalize();
    PetscFinalize();
    MPI_Finalize();

    return EXIT_SUCCESS;

}