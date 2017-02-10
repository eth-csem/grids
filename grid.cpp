#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "point_clouds.h"

/*! \file
 \brief Main executable.
 */

/*==========================================================================*/
/* Main function to compute various grids. ---------------------------------*/
/*==========================================================================*/

int main(int argc, char *argv[])
{
    /* Local variables. ----------------------------------------------------*/    
    PointCloud C;
    void print_help();

    /* Check which grid to build. ------------------------------------------*/
    
    /* Cubed sphere. -------------------------------------------------------*/
    if (!strcmp(argv[1],"cubed_sphere"))
    {
        printf("\ncubed sphere\n");
        printf("name of refinement list: %s\n",argv[2]);
        printf("output written to: %s\n",argv[3]);
        printf("radius of the sphere: %f km\n",atof(argv[4]));
        printf("minimum point distance: %f km\n",atof(argv[5]));

        C.cubed_sphere(argv[2],atof(argv[4]));
        C.write(argv[3],atof(argv[5]));
    }
    
    /* Cubed ball. ---------------------------------------------------------*/
    else if (!strcmp(argv[1],"cubed_ball"))
    {
        printf("\ncubed ball\n");
        printf("name of refinement list: %s\n",argv[2]);
        printf("output written to: %s\n",argv[3]);
        printf("minimum point distance: %f km\n",atof(argv[4]));
        
        C.cubed_ball(argv[2]);
        C.write(argv[3], atof(argv[4]));
    }

    
    /* Fibonacci sphere. ---------------------------------------------------*/
    else if (!strcmp(argv[1],"fibonacci_sphere"))
    {
        printf("\nFibonacci sphere\n");
        printf("name of refinement list: %s\n",argv[2]);
        printf("output written to: %s\n",argv[3]);
        printf("radius of the sphere: %f km\n",atof(argv[4]));
        printf("minimum point distance: %f km\n",atof(argv[5]));
        
        C.fibonacci_sphere(argv[2],atof(argv[4]));
        C.write(argv[3], atof(argv[5]));
    }
    
    /* Fibonacci ball. -----------------------------------------------------*/
    else if (!strcmp(argv[1],"fibonacci_ball"))
    {
        printf("\nFibonacci ball\n");
        printf("name of refinement list: %s\n",argv[2]);
        printf("output written to: %s\n",argv[3]);
        printf("minimum point distance: %f km\n",atof(argv[4]));
        
        C.fibonacci_ball(argv[2]);
        C.write(argv[3], atof(argv[4]));
    }
    
    /* Regular spherical grid. ---------------------------------------------*/
    else if (!strcmp(argv[1],"regular"))
    {
        printf("\nregular spherical grid\n");
        printf("name of refinement list: %s\n",argv[2]);
        printf("output written to: %s\n",argv[3]);
        printf("minimum point distance: %f km\n",atof(argv[4]));
        
        C.regular(argv[2]);
        C.write(argv[3], atof(argv[4]));
    }
    
    /* Get help. -----------------------------------------------------------*/
    else if (!strcmp(argv[1],"-help"))
    {
        print_help();
    }
    
    /* Error. --------------------------------------------------------------*/
    else
    {
        printf("ERROR! No valid option!\n");
        print_help();
    }
    
    return 0;
}


/*==========================================================================*/
/* Print help message. -----------------------------------------------------*/
/*==========================================================================*/

void print_help()
{
    printf("\nUsage of grid:\n");
    printf("--------------\n");
    printf("grid cubed_sphere [name of refinement list] [output file name] [radius of the sphere in km] [minimum point distance in km]\n");
    printf("grid cubed_ball [name of refinement list] [output file name] [minimum point distance in km]\n");
    printf("grid fibonacci_sphere [name of refinement list] [output file name] [radius of the sphere in km] [minimum point distance in km]\n");
    printf("grid fibonacci_ball [name of refinement list] [output file name] [minimum point distance in km]\n");
    printf("grid regular [name of refinement list] [output file name] [minimum point distance in km]\n");
}