#ifndef auxiliary_h
#define auxiliary_h

/*! \file
 \brief Definitions for auxiliary classes and functions to compute point clouds.
 */

/********************************************************************************//**
 * Point in 3D.
 ***********************************************************************************/
class Point
{
    public:
        double x;   /**< x coordinate. */
        double y;   /**< y coordinate. */
        double z;   /**< z coordinate. */
    
        Point *next;
    
        Point(double x, double y, double z);    /**< Constructor to initialise with x, y, z values. */
        Point();                                /**< Default constructor initialises all values to 0. */
    
        Point operator=(const Point a)          /**< Assignent operator. */
        {
            x=a.x;
            y=a.y;
            z=a.z;
            return *this;
        }
    
        void print() const;                     /**< Print point to screen. */
        void set(double x, double y, double z); /**< Set x, y, z values of the point. */
};

bool operator==(const Point a, const Point b);
bool operator!=(const Point a, const Point b);
bool operator>(const Point a, const Point b);
bool operator<(const Point a, const Point b);

/********************************************************************************//**
 * List of points.
 ***********************************************************************************/
class Pointlist
{
    public:
        Point *start;               /**< List anchor. */
    
    public:
        int n;                      /**< Number of elements in the list. */
    
        Pointlist();                /**< Constructor. */
        ~Pointlist();               /**< Destructor. Removes all list entries. */
    
        void append(Point p);                       /**< Append a point. */
        void append(double x, double y, double z);  /**< Append a point. */
        void print();               /**< Print the list to the screen. */
    
        void list2array(Point *p);  /**< Write list into an array p of size n. */
};


/********************************************************************************//**
 * Refinement region. Defines the Cartesian extent of a refinement region,
 * as well es the horizontal and vertical refinements. This is done for two types
 * of points clouds: cubed sphere/ball and Fibonacci sphere/ball.
 ***********************************************************************************/
class Refinement_Regions
{
    public:
    
        /* General properties. ----------------------------------------------------*/
    
        int n;           /**< Number of refinement regions. */
    
        int* n_points;   /**<   Cubed sphere:       not used.
                                Fibonacci sphere:   Number of points in the Fibonacci sphere. */
    
        /* Spatial extent of refinement region. -----------------------------------*/
    
        double* xmin;    /**<   Cubed sphere:       Minimum x of the refinement region within the unit cube [-1 1]x[-1 1]x[-1 1].
                                Fibonacci sphere:   Minimum colatitude [deg].
                                Regular grid:       Minimum colatitude [deg]. */
        double* xmax;    /**<   Cubed sphere:       Maximum x of the refinement region within the unit cube [-1 1]x[-1 1]x[-1 1].
                                Fibonacci sphere:   Maximum colatitude [deg].
                                Regular grid:       Maximum colatitude [deg]. */
        double* ymin;    /**<   Cubed sphere:       Minimum y of the refinement region within the unit cube [-1 1]x[-1 1]x[-1 1]. 
                                Fibonacci sphere:   Minimum longitude [deg].
                                Regular grid:       Minimum longitude [deg]. */
        double* ymax;    /**<   Cubed sphere:       Maximum y of the refinement region within the unit cube [-1 1]x[-1 1]x[-1 1]. 
                                Fibonacci sphere:   Maximum longitude [deg].
                                Regular grid:       Maximum longitude [deg]. */
        double* zmin;    /**<   Cubed sphere:       Minimum z of the refinement region within the unit cube [-1 1]x[-1 1]x[-1 1]. 
                                Fibonacci sphere:   Not used.
                                Regular grid:       Minimum radius [km]. */
        double* zmax;    /**<   Cubed sphere:       Maximum z of the refinement region within the unit cube [-1 1]x[-1 1]x[-1 1]. 
                                Fibonacci sphere:   Not used. 
                                Regular grid:       Maximum radius [km]. */
        double* rmin;    /**<   Cubed sphere:       Normalised minimum radius [1].
                                Fibonacci sphere:   Minimum radius [km].
                                Regular grid:       Not used. */
        double* rmax;    /**<   Cubed sphere:       Normalised maximum radius [1].
                                Fibonacci sphere:   Maximum radius [km].
                                Regular grid:       Not used. */
        double* d_horizontal;    /**<   Cubed sphere:       Horizontal refinement within the unit cube.
                                        Fibonacci sphere:   Not used. 
                                        Regular grid:       Horizontal spacing [deg]. */
        double* d_vertical;      /**<   Cubed sphere:       Vertical refinement within the unit cube. 
                                        Fibonacci sphere:   Radial increment. 
                                        Regular grid:       Vertical spacing [km]. */
    
        /* Rotation properties. ---------------------------------------------------*/
    
        double* phi;     /**<   Cubed sphere:       Not used.
                                Fibonacci sphere:   Rotation angle [radians]. 
                                Regular grid:       Rotation angle [radians]. */
        double* b[3];    /**<   Cubed sphere:       Not used.
                                Fibonacci sphere:   Normalised rotation vector. 
                                Regular grid:       Normalised rotation vector. */
    
        /* Functions. -------------------------------------------------------------*/
    
        Refinement_Regions();
        ~Refinement_Regions();
    
        void read_cubed(
                        const char* refinement_list  /**< Name of file containing list of refinement regions. */
                    );
    
        void read_fibonacci(
                        const char* refinement_list  /**< Name of file containing list of refinement regions. */
                    );
    
        void read_regular(
                        const char* refinement_list  /**< Name of file containing list of refinement regions. */
                    );
    
        void print_cubed();     /**< Print refinement regions for cubed sphere on screen. */
        void print_fibonacci(); /**< Print refinement regions for Fibonacci sphere on screen. */
        void print_regular();   /**< Print refinement regions for regular spherical grid on screen. */
};

/********************************************************************************//**
 * Projection from the unit cube onto the sphere.
 ***********************************************************************************/
void project2sphere(
                    Point* p,       /**< Pointer to an array of points. */
                    int n_points,   /**< Number of points in the array. */
                    double R        /**< Radius in km. */
);

Point project2sphere(
                    double x,       /**< x-coordinate in the unit cube. */
                    double y,       /**< y-coordinate in the unit cube. */
                    double z,       /**< z-coordinate in the unit cube. */
                    double R        /**< Radius in km. */
);

/********************************************************************************//**
 * Rotation.
***********************************************************************************/

/**
 * Rotation matrix for points in Cartesian coordinates.
 */
void rotation_matrix(
                     double R[][3],    /**< Rotation matrix, to be computed. */
                     double b[],     /**< Rotation unit vector. Vector about which to rotate. */
                     double phi     /**< Rotation angle in degrees. */
);

/********************************************************************************//**
 * Sorting algorithms for 3D points.
***********************************************************************************/

/**
 * Quicksort algorithm, recursive part.
 */
void quicksort(
               Point* array,        /**< Pointer to the input array. */
               unsigned int start,  /**< Starting index. */
               unsigned int end     /**< Final index (must be less than length of array). */
);

/**
 * Quicksort algorithm, partitioning part.
 */
unsigned int partition(
                       Point* array,        /**< Pointer to the input array. */
                       unsigned int start,  /**< Starting index. */
                       unsigned int end     /**< Final index (must be less than length of array). */
);


/********************************************************************************//**
 * Remove newline character from a string.
***********************************************************************************/

/**
 * Removing newline from a string, read for instance with fgets.
 */
void remove_newline(
                    char *line  /**< Pointer to string. */
);

/********************************************************************************//**
 * Compute number of grid points in a refinement region.
***********************************************************************************/

/**
 * Compute number of grid points in a refinement region.
 */
int number_of_grid_points(
                                Refinement_Regions &r,  /**< List of refinement regions. */
                                int index               /**< Index of a specific refinement region. */
);


#endif /* auxiliary_h */
