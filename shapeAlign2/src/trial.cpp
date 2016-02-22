#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <string>

using namespace std;

int main(){
	  int i, j;
  gsl_matrix * m = gsl_matrix_alloc (10, 3);

  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);


  	gsl_matrix * c= gsl_matrix_alloc(10,3);
  	gsl_matrix_memcpy(c,m);

  	gsl_matrix_set(c,0,0,GSL_NAN);

  	gsl_matrix_sub(c,m);


  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      printf ("c(%d,%d) = %g\n", i, j,
              gsl_matrix_get (m, i, j));

  for (i = 0; i < 10 ;i++){
    gsl_vector_view row = gsl_matrix_row(m,i);
    gsl_vector_add_constant(&row.vector,100);
  }

  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      printf ("c(%d,%d) = %g\n", i, j,
              gsl_matrix_get (m, i, j));


  gsl_matrix_free (m);
	gsl_matrix_free(c);
  return 0;
}
