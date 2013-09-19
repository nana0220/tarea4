#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h> 

/*metodo para calcular los An------------------------------------*/

void main(int argc,char **argv){

int f=atof(argv[2]);
int n=atoi(argv[2]);

// n es el grado que se recibe por entrada

if(n < 0) {
     printf("\n El numero es negativo");
exit(1);
   }
if(f%1) {
     printf("\n El numero debe ser entero");
exit(1);
   }

 FILE *in;
 in=fopen(argv[1],"r");

 double *x;
 double *y;

int lineas=0;
int letra=0;

  //contar lineas archivo

 do{
    letra=fgetc(in) ;
    if(letra=='\n'){
      lineas++;
    }
  }while(letra!=EOF);
  rewind(in);
  printf("el numero de filas es %d \n", lineas);


//creacion matriz

  double *M;
  double *MT;
  int columnas=n+1;
  int pos;
  int i,j;
  int s; //signo permutaciones

x = malloc(lineas * sizeof(float));
y = malloc(lineas * sizeof(float));

 printf("se crea aloc memoria con los datos");


 //guardar datos en matriz

 for(i=0; i<lineas; i++){
 fscanf(in,"%lf %lf", &x[i],&y[i]);
    }
 printf("se crearon las matrices de datos \n");


  /*creando matriz*/

 //dandole memoria a matriz y transpuesta

M=malloc(sizeof(float)*lineas*(columnas));
MT=malloc(sizeof(float)*lineas*(columnas));

//creando matriz normal

  for(i=0;i<lineas;i++){

    for(j=0;j<columnas;j++){

      pos=j+(lineas*i);
      M[pos]=pow(x[i],j);
   }
  }

  //crenaod matriz transpuesta
   for(i=0;i<lineas;i++){

    for(j=0;j<columnas;j++){

      pos=i+(lineas*j);
      M[pos]=pow(x[i],j);
   }
  }
   printf("las matrices han sido creadas");

   
gsl_matrix_view m= gsl_matrix_view_array (M, lineas, columnas);
gsl_matrix_view mt= gsl_matrix_view_array (MT, columnas, lineas);

 printf("matrices en terminos de gsl");

 gsl_matrix_view v=gsl_matrix_view_array(y,lineas,1);

  printf("se creo v\n");

gsl_matrix * multiplicacion = gsl_matrix_alloc(columnas,columnas);
  gsl_matrix * prueba = gsl_matrix_alloc(columnas,columnas);

  printf("se crearon matrices de multiplicacion y de prueba");

  //inversa
gsl_matrix * inversa = gsl_matrix_alloc(columnas,columnas);
//total
gsl_matrix * total = gsl_matrix_alloc(columnas,columnas);

 printf("se crea matriz total e inversa");

 gsl_matrix *vector;
 vector=&v.matrix;
 gsl_matrix*m=gsl_matrix_alloc(lineas,columnas);
 m=&m.matrix;
 gsl_matrix *mt=gsl_matrix_alloc(columnas,lineas);
 mt=&MT.matrix;

gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m, mt, 0.0, multiplicacion);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mt, m, 0.0, prueba);

  // matriz de permutacion para el algoritmo LU
  gsl_permutation * perm = gsl_permutation_alloc (columnas);
  
  printf("SE CREO la matriz de permutacion \n");

gsl_linalg_LU_decomp(prueba,perm,&s);
gsl_linalg_LU_invert(prueba, perm, inversa);

 printf("se creo descomp lu");

gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inversa, mt, 0.0, total);

gsl_matrix * result = gsl_matrix_alloc(columnas,1);

gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, total, vector, 0.0, result);

gsl_matrix * chi= gsl_matrix_alloc(lineas,1);

gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m, result, 0.0, chi);

printf("se han sacado los coeficientes An");


/*sacando chi cuadrado*/

  double xi=0;	

  for(i=0;i<lineas;i++){
    printf("%d\n",i);
    printf("%lf\n",gsl_matrix_get(chi, i, 0));
    xi=xi+ pow((y[i]- gsl_matrix_get(chi, i, 0)),2)/lineas;
  }

// imprime resultado

   for(i=0;i<lineas;i++){

     printf("m%d = %lf\n",i,gsl_matrix_get(result,i,0));

   }

   printf("X²=%lf\n",xi);





}


