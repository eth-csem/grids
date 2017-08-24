#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "point_clouds.h"

/*! \file
 \brief Main executable. Read input, calls functions to generate specific types of
 point clouds, and write the point clouds into a file..
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
    
    /* Vertical profile. ---------------------------------------------------*/
    else if (!strcmp(argv[1],"profile"))
    {
        printf("\nvertical profile\n");
        printf("latitude=%lg deg, longitude=%lg deg\n",atof(argv[2]),atof(argv[3]));
        printf("radius=%lg:%lg:%lg km\n",atof(argv[4]),atof(argv[6]),atof(argv[5]));
        
        C.profile(atof(argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]));
        C.write(argv[7], 0.0);
    }
    
    /* Vertical slice. -----------------------------------------------------*/
    else if (!strcmp(argv[1],"slice"))
    {
        printf("\nvertical slice\n");
        printf("min. latitude=%lg deg, min. longitude=%lg deg, min. radius=%lg km\n",atof(argv[2]),atof(argv[3]),atof(argv[4]));
        printf("max. latitude=%lg deg, max. longitude=%lg deg, max. radius=%lg km\n",atof(argv[5]),atof(argv[6]),atof(argv[7]));
        printf("radius increment=%lg km, angular increment=%lg deg\n",atof(argv[8]),atof(argv[9]));
        
        C.slice(atof(argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8]),atof(argv[9]));
        C.write(argv[10], 0.0);
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
    printf("grid cubed_sphere [name of refinement list] [output file name] [radius of the sphere (km)] [min. point distance (km)]\n");
    printf("grid cubed_ball [name of refinement list] [output file name] [min. point distance (km)]\n");
    printf("grid fibonacci_sphere [name of refinement list] [output file name] [radius of the sphere (km)] [min. point distance (km)]\n");
    printf("grid fibonacci_ball [name of refinement list] [output file name] [min. point distance (km)]\n");
    printf("grid regular [name of refinement list] [output file name] [minimum point distance (km)]\n");
    printf("grid profile [latitude (deg)] [longitude (deg)] [min. radius (km)] [max. radius (km)] [radius increment (km)]\n");
    printf("grid slice [min. latitude (deg)] [min. longitude (deg)] [min. radius (km)] [max. latitude (deg)] [max. longitude (deg)] [max. radius (km)] [radius increment (km)] [angular increment (deg)]\n");
}