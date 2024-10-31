#include "csr_matrix.hpp"

// Convert a dense matrix with only main diagonal and one-off diagonals to CSR format
void sparse_csr(float dense[N][N], CSRMatrix& csr) {
    int value_index = 0;            // Index to track position in values and col_indices
    csr.row_pointers[0] = 0;        // Start the first row pointer at 0

    for (int i = 0; i < N; ++i) {
        // Main diagonal element
        if (dense[i][i] != 0.0) {
            csr.values[value_index] = dense[i][i];
            csr.col_indices[value_index] = i;
            value_index++;
        }

        // Subdiagonal element if within bounds and non-zero
        if (i > 0 && dense[i][i - 1] != 0.0) {
            csr.values[value_index] = dense[i][i - 1];
            csr.col_indices[value_index] = i - 1;
            value_index++;
        }

        // Superdiagonal element if within bounds and non-zero
        if (i < N - 1 && dense[i][i + 1] != 0.0) {
            csr.values[value_index] = dense[i][i + 1];
            csr.col_indices[value_index] = i + 1;
            value_index++;
        }

        // Update row pointers to the current size of values
        csr.row_pointers[i + 1] = value_index;
    }
}