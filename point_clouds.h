#ifndef point_clouds_h
#define point_clouds_h

#include <stdio.h>
#include "auxiliary.h"

/*! \file
 \brief Definition of the point cloud classes.
 */

/********************************************************************************//**
 * Point cloud in 3D.
***********************************************************************************/
class PointCloud
{
    private:
        Point* p;               /**< Array of 3D points to be filled by one of the member functions. */
        unsigned int n_points;  /**< Number of points. */
    
    public:
        PointCloud();   /**< Constructor. */
        ~PointCloud();  /**< Destructor. */
    
        /**
         * Write points to file. Duplicate points are omitted provided that the array of points
         has been sorted before.
         */
        void write(
                    const char* filename    /**< Pointer to filename. */
                    ) const;
    
        /** 
         * Compute a sorted point cloud for a refinable cubed sphere.
         */
        void cubed_sphere(
                        const char* filename,   /**< Name of the file containing the list of refinement regions. */
                        double r                /**< Radius of the sphere [km]. */
                            );
    
        /**
         * Compute a sorted point cloud for a refinable cubed ball.
         */
        void cubed_ball(
                        const char* filename    /**< Name of the file containing the list of refinement regions. */
                        );
    
        /**
         * Compute a point cloud for a refinable Fibonacci sphere.
         */
        void fibonacci_sphere(
                        const char* filename,   /**< Name of the file containing the list of refinement regions. */
                        double r                /**< Radius of the sphere [km]. */
                            );
    
        /**
         * Compute a point cloud for a refinable Fibonacci ball.
         */
        void fibonacci_ball(
                        const char* filename    /**< Name of the file containing the list of refinement regions. */
                            );
    
        /**
         * Compute a point cloud for a regular spherical grid.
         */
        void regular(
                        const char* filename    /**< Name of the file containing the list of refinement regions. */
                            );
};

#endif /* point_clouds_h */
