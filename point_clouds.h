#ifndef point_clouds_h
#define point_clouds_h

#include <stdio.h>
#include "auxiliary.h"

/*! \file
 \brief The source file point_clouds.h defines the PointCloud class. A member of this class is an array of 3D points in space that can be filled by the member functions of the class, which generate a specific type of point cloud. The point cloud can be written to a simple output text file, and in this process, points that are too close to another one can be removed.
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
                    const char* filename,    /**< Pointer to filename. */
                    double d_min             /**< Minimum allowable distance to the nearest grid point. */
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
    
        /** 
         * Compute a point cloud for a vertical profile.
         */
        void profile(
                        double lat,             /**< Latitude [deg]. */
                        double lon,             /**< Longitude [deg]. */
                        double r_min,           /**< Minimum radius [km]. */
                        double r_max,           /**< Maximum radius [km]. */
                        double dr               /**< Radius increment [km]. */
                     );
    
        /** 
         * Compute a point cloud for a vertical slice.
         */
        void slice(
                        double lat_min,         /**< Minimum latitude [deg]. */
                        double lon_min,         /**< Minimum longitude [deg]. */
                        double r_min,           /**< Minimum radius [km]. */
                        double lat_max,         /**< Maximum latitude [deg]. */
                        double lon_max,         /**< Maximum longitude [deg]. */
                        double r_max,           /**< Maximum radius [km]. */
                        double dr,              /**< Radius increment [km]. */
                        double dgamma           /**< Angular increment [deg]. */
                    );
    
};

#endif /* point_clouds_h */
