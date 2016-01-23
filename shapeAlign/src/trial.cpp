#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <string>

using namespace std;

int main(){
	//   int i, j;
 //  gsl_matrix * m = gsl_matrix_alloc (10, 3);

 //  for (i = 0; i < 10; i++)
 //    for (j = 0; j < 3; j++)
 //      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);


 //  	gsl_matrix * c= gsl_matrix_alloc(10,3);
 //  	gsl_matrix_memcpy(c,m);

 //  	gsl_matrix_set(c,0,0,GSL_NAN);

 //  	gsl_matrix_sub(c,m);


 //  for (i = 0; i < 10; i++)
 //    for (j = 0; j < 3; j++)
 //      printf ("c(%d,%d) = %g\n", i, j,
 //              gsl_matrix_get (m, i, j));

 //  for (i = 0; i < 10 ;i++){
 //    gsl_vector_view row = gsl_matrix_row(m,i);
 //    gsl_vector_add_constant(&row.vector,100);
 //  }

 //  for (i = 0; i < 10; i++)
 //    for (j = 0; j < 3; j++)
 //      printf ("c(%d,%d) = %g\n", i, j,
 //              gsl_matrix_get (m, i, j));


 //  gsl_matrix_free (m);
	// gsl_matrix_free(c);


  gsl_vector* v = gsl_vector_alloc(5);

  for (int i = 0; i < 5 ; i++)
    gsl_vector_set(v,i,i);

  gsl_vector_set(v,0,GSL_NAN);

  double mu = gsl_stats_mean(v->data,v->stride,v->size);

  cerr << "With Nan: " << mu << endl;

  gsl_vector* x = gsl_vector_alloc(4);
  for (int i = 0; i < 4; i++)
    gsl_vector_set(x,i,i+1);

  mu =gsl_stats_mean(x->data,x->stride,x->size);

  cerr << "Without Nan: " << mu << endl;

  return 0;
}
