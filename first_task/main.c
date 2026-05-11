#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON __DBL_EPSILON__

typedef struct Matrix {
  size_t size;
  double** table;
} Matrix;

void InitMatrix(Matrix* mat) {
  mat->table = (double**) malloc(mat->size * sizeof(double*));
  scanf("%zu", &mat->size);

  for (size_t i = 0; i < mat->size; ++i) {
    mat->table[i] = (double*) malloc(mat->size * sizeof(double));
    for (size_t j = 0; j < mat->size; ++j) {
      scanf("%lf", &mat->table[i][j]);
    }
  }

  return;
}

void PrintMatrix(Matrix* mat) {
  for (size_t i = 0; i < mat->size; ++i) {
    for (size_t j = 0; j < mat->size; ++j) {
      printf("%10lf ", mat->table[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  return;
}

void DeleteMatrix(Matrix* mat) {
  for (size_t i = 0; i < mat->size; ++i) {
    free(mat->table[i]);
  }
  free(mat->table);
  return;
}

int CheckZeroLines(Matrix* mat) {
  for (size_t i = 0; i < mat->size; ++i) {
    int is_zero = 1;
    for (size_t j = 0; j < mat->size; ++j) {
      if (mat->table[i][j] < -EPSILON || mat->table[i][j] > EPSILON) {
        is_zero = 0;
        break;
      }
    }
    if (is_zero) {
      return 1;
    }
  }

  return 0;
}

void ChangeLines(Matrix* mat, size_t line1, size_t line2) {
  if (line1 == line2) {
    return;
  }

  for (size_t i = 0; i < mat->size; ++i) {
    double tmp = mat->table[line1][i];
    mat->table[line1][i] = mat->table[line2][i];
    mat->table[line2][i] = tmp;
  }
  return;
}

double GaussianElimination(Matrix* mat) {
  double determinant = 1;

  for (size_t step = 0; step < mat->size; ++step) {
    if (CheckZeroLines(mat)) {
      return 0;
    }

    int divide_line = -1;
    for (size_t i = step; i < mat->size; ++i) {
      if (mat->table[i][step] < -EPSILON || mat->table[i][step] > EPSILON) {
        divide_line = i;
        break;
      }
    }

    if (divide_line == -1) {
      return 0;
    }

    ChangeLines(mat, divide_line, step);

    for (size_t i = step + 1; i < mat->size; ++i) {
      double lambda = -mat->table[i][step] / mat->table[step][step];
      mat->table[i][step] = 0;

      if (lambda < -EPSILON || lambda > +EPSILON) {
        for (size_t j = step + 1; j < mat->size; ++j) {
          mat->table[i][j] += lambda * mat->table[step][j];
        }
      }
    }
  }

  for (size_t i = 0; i < mat->size; ++i) {
    determinant *= mat->table[i][i];
  }
  return determinant;
}


int main() {
  Matrix matrix;
  InitMatrix(&matrix);
  printf("\nInput matrix:\n\n");
  PrintMatrix(&matrix);

  double determinant = GaussianElimination(&matrix);
  printf("Matrix after calculations:\n\n");
  PrintMatrix(&matrix);
  printf("Determinant = %lf\n", determinant);

  DeleteMatrix(&matrix);
  return 0;
}
