#include <stdio.h>
#include <math.h>
#include "point_clouds.h"

/*! \file
 \brief Implement PointCould class..
 */

const double PI = 3.14159265358979323846264338327;
const double RE = 6371.0;

/*=================================================================================*/
/* Point cloud in 3D. -------------------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/

PointCloud::PointCloud()
{
    p=0;
    n_points=0;
}

/* Destructor. --------------------------------------------------------------------*/

PointCloud::~PointCloud()
{
    if (p!=0) delete[] p;
}

/* Write points to a file. --------------------------------------------------------*/

void PointCloud::write(const char* filename, double d_min) const
{
    /* Local variables. */
    FILE* pfile;
    int j, k, i=0;
    int n_accepted=0, n_rejected=0;
    bool* is_far_enough;
    
    /* Initialise the vector indicating if the point is far enough from the others. */
    is_far_enough=new bool[n_points];
    for (int k=0; k<n_points; k++) is_far_enough[k]=true;
    
    pfile=fopen(filename,"w");
    
    /* Write grid points without eliminating points that are too close. */
    if (p && pfile && d_min==0.0)
    {
        for (int i=0; i<n_points; i++)
        {
            fprintf(pfile, "%lg %lg %lg\n",p[i].x,p[i].y,p[i].z);
        }
        n_accepted=n_points;
    }
    
    /* Write grid points, except those that are too close to each other. */
    else if (p && pfile && d_min>0.0)
    {
        /* Sort with increasing x coordinate. */
        if (n_points>1) quicksort(p,0,n_points-1);

        /* Walk through the grid points. */
        while (i<n_points)
        {
            /* Print the current grid point. */
            fprintf(pfile, "%lg %lg %lg\n",p[i].x,p[i].y,p[i].z);
            n_accepted++;
            
            /* Walk from the current grid point until the x coordinate is too far away. */
            for (k=i+1; fabs(p[k].x-p[i].x)<d_min && k<n_points; k++)
            {
                /* Check if the candidate grid points are far enough from the current grid point. */
                if ((fabs(p[k].y-p[i].y)<d_min && fabs(p[k].z-p[i].z)<d_min) && is_far_enough[k]==true)
                {
                    is_far_enough[k]=false;
                    n_rejected++;
                }
            }
            /* Set the first accepted grid point, counted from the current one, as the new current one. */
            for (j=i+1; j<=k; j++)
            {
                if (is_far_enough[j]==true) break;
            }
            i=j;
        }
    }
    
    /* Output about accepted and rejected points. */
    printf("total number of grid points: %d\n",n_points);
    printf("number of grid points written: %d\n",n_accepted);
    printf("number of grid points rejected: %d\n",n_rejected);

    /* Clean up. */
    delete[] is_far_enough;
    if (pfile) fclose(pfile);
}

/* Cubed sphere. ------------------------------------------------------------------*/

void PointCloud::cubed_sphere(const char* filename, double R)
{
    /* Local variables. */
    Pointlist pl;
    Refinement_Regions r;
    
    /* Read refinement regions. */
    
    r.read_cubed(filename);
    
    /* March through refinement regions and make list. */
    for (int idx=0; idx<r.n; idx++)
    {
        /* Check if radius is within radii of the refinement region. */
        if (R/RE>=r.rmin[idx] && R/RE<=r.rmax[idx])
        {
            for (double x=r.xmin[idx]; x<=r.xmax[idx]; x+=r.d_horizontal[idx])
            {
                for (double y=r.ymin[idx]; y<=r.ymax[idx]; y+=r.d_horizontal[idx])
                {
                    for (double z=r.zmin[idx]; z<=r.zmax[idx]; z+=r.d_horizontal[idx])
                    {
                        pl.append(x,y,z);
                    }}}}}
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);

    /* Project onto sphere. */
    project2sphere(p,pl.n,R);
}

/* Cubed ball. --------------------------------------------------------------------*/

void PointCloud::cubed_ball(const char* filename)
{
    /* Local variables. */
    Point p_temp;
    Pointlist pl;
    Refinement_Regions r;
    double R;
    
    /* Read refinement regions. */
    r.read_cubed(filename);
    //r.print_cubed();
    
    /* March through refinement regions and make list. */
    for (int idx=0; idx<r.n; idx++)
    {
        /* Loop over normalised radii in the refinement region. */
        for (double R=r.rmin[idx]; R<=r.rmax[idx]; R+=r.d_vertical[idx])
        {
            printf("%d %lg\n",idx, R);
            /* Loops over the 3 Cartesian coordinates. */
            for (double x=r.xmin[idx]; x<=r.xmax[idx]; x+=r.d_horizontal[idx])
            {
                for (double y=r.ymin[idx]; y<=r.ymax[idx]; y+=r.d_horizontal[idx])
                {
                    for (double z=r.zmin[idx]; z<=r.zmax[idx]; z+=r.d_horizontal[idx])
                    {
                        p_temp=project2sphere(x,y,z,R*RE);
                        pl.append(p_temp.x,p_temp.y,p_temp.z);
                    }}}}}
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);
}


/* Fibonacci sphere. --------------------------------------------------------------*/

void PointCloud::fibonacci_sphere(const char* filename, double radius)
{
    /* Local variables. */
    Pointlist pl;
    Point *temp;
    Refinement_Regions ref;
    double R[3][3], n[3];
    double x, y, z, x_temp, y_temp, z_temp;
    double r, theta, phi;
    double golden_angle, part_of_circle, offset;
    double phi_min, phi_max, theta_min, theta_max, z_start, z_end, z_len;
    
    /* Read refinement regions. */
    ref.read_fibonacci(filename);

    /* March through refinement regions and make a list. */
    for (int idx=0; idx<ref.n; idx++)
    {
        /* Check if this refinement region overlaps with the radius of the sphere. */
        if (radius>=ref.rmin[idx] && radius<=ref.rmax[idx])
        {
            /* Setup for the Fibonacci sphere. */
            phi_min=ref.ymin[idx]*PI/180.0;
            phi_max=ref.ymax[idx]*PI/180.0;
            theta_min=ref.xmin[idx]*PI/180.0;
            theta_max=ref.xmax[idx]*PI/180.0;
            z_start=cos(theta_min);
            z_end=cos(theta_max);
            z_len=z_end-z_start;
            part_of_circle=(phi_max-phi_min)/(2.0*PI);
            golden_angle=PI*(3.0-sqrt(5.0))*part_of_circle;
            offset=z_len/ref.n_points[idx];
            
            /* Compute rotation matrix. */
            if (ref.phi[idx]!=0.0)
            {
                n[0]=ref.b[0][idx]; n[1]=ref.b[1][idx]; n[2]=ref.b[2][idx];
                rotation_matrix(R,n,ref.phi[idx]*PI/180.0);
            }
        
            /* Remove already existing points in the list that fall within this refinement region. */
            if (idx>0)
            {
                /* March through all the points in the list. */
                /* Consider successor of current point, pt->next, to be removed. */
                
                for (Point *pt=pl.start; pt->next; pt=pt->next)
                {
                    /* Rotate the point in reverse direction using transpose of rotation matrix. */
                    if (ref.phi[idx]!=0.0)
                    {
                        x=R[0][0]*pt->next->x+R[1][0]*pt->next->y+R[2][0]*pt->next->z;
                        y=R[0][1]*pt->next->x+R[1][1]*pt->next->y+R[2][1]*pt->next->z;
                        z=R[0][2]*pt->next->x+R[1][2]*pt->next->y+R[2][2]*pt->next->z;
                    }
                    else
                    {
                        x=pt->next->x; y=pt->next->y; z=pt->next->z;
                    }
               
                    /* Compute spherical coordinates (of backward rotated point). */
                    r=sqrt(x*x+y*y+z*z);
                    theta=acos(z/r);
                    phi=atan2(y,x);
            
                    /* Check if point falls into current refinement region. If so, remove from list. */
                    if (theta < theta_max && theta > theta_min && phi < phi_max && phi > phi_min)
                    {
                        temp=pt->next;
                        if (pt->next->next) // Avoid removing last point in the list, which would cause problem in the loop statement.
                        {
                            pt->next=pt->next->next;
                            delete temp;
                            pl.n--;
                        }
                    }
                }
            } /* END: Remove already existing points in the list that fall within this refinement region. */
            
            /* Compute Cartesian coordinates for all the points on the sphere segment. */
            for (int i=0; i<ref.n_points[idx]; i++)
            {
                /* Unrotated coordinates. */
                z=((i*offset)+z_start)+offset/2.0;
                r=sqrt(1.0-z*z);
                phi=fmod(i*golden_angle,part_of_circle*2.0*PI)+phi_min;
                
                x=cos(phi)*r; y=sin(phi)*r;
                x=radius*x; y=radius*y; z=radius*z;
            
                /* Rotate if needed. */
                if (ref.phi[idx]!=0.0)
                {
                    x_temp=R[0][0]*x+R[0][1]*y+R[0][2]*z;
                    y_temp=R[1][0]*x+R[1][1]*y+R[1][2]*z;
                    z_temp=R[2][0]*x+R[2][1]*y+R[2][2]*z;
                    x=x_temp; y=y_temp; z=z_temp;
                }
                
                /* Append to the list of points. */
                pl.append(x,y,z);
            }
        } /* END: Check if this refinement region overlaps with the radius of the sphere. */
    } /* END: March through refinement regions and make a list. */
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);
}

/* Fibonacci sphere. --------------------------------------------------------------*/

void PointCloud::fibonacci_ball(const char* filename)
{
    /* Local variables. */
    Pointlist pl;
    Point *temp;
    Refinement_Regions ref;
    double R[3][3], n[3];
    double x, y, z, x_temp, y_temp, z_temp;
    double r, theta, phi;
    double golden_angle, part_of_circle, offset;
    double phi_min, phi_max, theta_min, theta_max, z_start, z_end, z_len;
    
    /* Read refinement regions. */
    ref.read_fibonacci(filename);
    
    /* March through refinement regions and make a list. */
    for (int idx=0; idx<ref.n; idx++)
    {
        /* Setup for the Fibonacci sphere. */
        phi_min=ref.ymin[idx]*PI/180.0;
        phi_max=ref.ymax[idx]*PI/180.0;
        theta_min=ref.xmin[idx]*PI/180.0;
        theta_max=ref.xmax[idx]*PI/180.0;
        z_start=cos(theta_min);
        z_end=cos(theta_max);
        z_len=z_end-z_start;
        part_of_circle=(phi_max-phi_min)/(2.0*PI);
        golden_angle=PI*(3.0-sqrt(5.0))*part_of_circle;
        offset=z_len/ref.n_points[idx];
        
        /* Compute rotation matrix. */
        if (ref.phi[idx]!=0.0)
        {
            n[0]=ref.b[0][idx]; n[1]=ref.b[1][idx]; n[2]=ref.b[2][idx];
            rotation_matrix(R,n,ref.phi[idx]*PI/180.0);
        }
        
        /* Remove already existing points in the list that fall within this refinement region. */
        if (idx>0)
        {
            /* March through all the points in the list. */
            /* Consider successor of current point, pt->next, to be removed. */
            for (Point *pt=pl.start; pt->next; pt=pt->next)
            {
                /* Rotate the point in reverse direction using transpose of rotation matrix. */
                if (ref.phi[idx]!=0.0)
                {
                    x=R[0][0]*pt->next->x+R[1][0]*pt->next->y+R[2][0]*pt->next->z;
                    y=R[0][1]*pt->next->x+R[1][1]*pt->next->y+R[2][1]*pt->next->z;
                    z=R[0][2]*pt->next->x+R[1][2]*pt->next->y+R[2][2]*pt->next->z;
                }
                else
                {
                    x=pt->next->x; y=pt->next->y; z=pt->next->z;
                }
                
                /* Compute spherical coordinates (of backward rotated point). */
                r=sqrt(x*x+y*y+z*z);
                theta=acos(z/r);
                phi=atan2(y,x);
                
                /* Check if point falls into current refinement region. If so, remove from list. */
                if (theta <= theta_max && theta >= theta_min && phi <= phi_max && phi >= phi_min && r <= ref.rmax[idx] && r >= ref.rmin[idx])
                {
                    temp=pt->next;
                    if (pt->next->next) // Avoid removing last point in the list.
                    {
                        pt->next=pt->next->next;
                        delete temp;
                        pl.n--;
                    }
                }
            }
        } /* END: Remove already existing points in the list that fall within this refinement region. */

        /* March through radii of the refinement region. */
        for (double radius=ref.rmin[idx]; radius<=ref.rmax[idx]; radius+=ref.d_vertical[idx])
        {
            /* Compute Cartesian coordinates for all the points on the sphere segment. */
            for (int i=0; i<ref.n_points[idx]; i++)
            {
                /* Unrotated coordinates. */
                z=((i*offset)+z_start)+offset/2.0;
                r=sqrt(1.0-z*z);
                phi=fmod(i*golden_angle,part_of_circle*2.0*PI)+phi_min;
                
                x=cos(phi)*r; y=sin(phi)*r;
                x=radius*x; y=radius*y; z=radius*z;
                
                /* Rotate if needed. */
                if (ref.phi[idx]!=0.0)
                {
                    x_temp=R[0][0]*x+R[0][1]*y+R[0][2]*z;
                    y_temp=R[1][0]*x+R[1][1]*y+R[1][2]*z;
                    z_temp=R[2][0]*x+R[2][1]*y+R[2][2]*z;
                    x=x_temp; y=y_temp; z=z_temp;
                }
                
                /* Append to the list of points. */
                pl.append(x,y,z);
            }
        } /* END: March through radii in this refinement region. */
    } /* END: March through refinement regions and make a list. */
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);
}

/* Regular spherical grid. --------------------------------------------------------*/

void PointCloud::regular(const char* filename)
{
    /* Local variables. */
    Pointlist pl;
    Refinement_Regions ref;
    double R[3][3], n[3];
    double x, y, z, x_temp, y_temp, z_temp;
    
    /* Read refinement regions. */
    ref.read_regular(filename);
    
    /* March through refinement regions and make list. */
    for (int idx=0; idx<ref.n; idx++)
    {
        /* Compute rotation matrix. */
        if (ref.phi[idx]!=0.0)
        {
            n[0]=ref.b[0][idx]; n[1]=ref.b[1][idx]; n[2]=ref.b[2][idx];
            rotation_matrix(R,n,ref.phi[idx]);
        }

        /* Loop over normalised radii in the refinement region. */
        for (double r=ref.zmin[idx]; r<=ref.zmax[idx]; r+=ref.d_vertical[idx])
        {
            /* Loops over the 3 Cartesian coordinates. */
            for (double theta=ref.xmin[idx]*PI/180.0; theta<=ref.xmax[idx]*PI/180.0; theta+=ref.d_horizontal[idx]*PI/180.0)
            {
                for (double phi=ref.ymin[idx]*PI/180.0; phi<=ref.ymax[idx]*PI/180.0; phi+=ref.d_horizontal[idx]*PI/180.0)
                {
                    x=r*cos(phi)*sin(theta);
                    y=r*sin(phi)*sin(theta);
                    z=r*cos(theta);
                    
                    /* Rotate if needed. */
                    if (ref.phi[idx]!=0.0)
                    {
                        x_temp=R[0][0]*x+R[0][1]*y+R[0][2]*z;
                        y_temp=R[1][0]*x+R[1][1]*y+R[1][2]*z;
                        z_temp=R[2][0]*x+R[2][1]*y+R[2][2]*z;
                        x=x_temp; y=y_temp; z=z_temp;
                    }

                    pl.append(x,y,z);
                }}}}
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);
}

/* Vertical profile. --------------------------------------------------------------*/

void PointCloud::profile(double lat, double lon, double r_min, double r_max, double dr)
{
    /* Local variables. */
    Pointlist pl;
    double x, y, z;
    
    double theta=(90.0-lat)*PI/180.0;
    double phi=lon*PI/180;
    
    /* Make grid points on the profile. */
    for (double r=r_min; r<=r_max; r+=dr)
    {
        x=r*cos(phi)*sin(theta);
        y=r*sin(phi)*sin(theta);
        z=r*cos(theta);
        pl.append(x,y,z);
    }
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);
}

/* Vertical slice. ----------------------------------------------------------------*/

void PointCloud::slice(double lat_min, double lon_min, double r_min, double lat_max, double lon_max, double r_max, double dr, double dgamma)
{
    /* Local variables. */
    Pointlist pl;
    double x, y, z, phi, theta;
    
    double theta_min=(90.0-lat_min)*PI/180.0;
    double theta_max=(90.0-lat_max)*PI/180.0;
    double phi_min=lon_min*PI/180.0;
    double phi_max=lon_max*PI/180.0;
    
    double x_min=cos(phi_min)*sin(theta_min);
    double y_min=sin(phi_min)*sin(theta_min);
    double z_min=cos(theta_min);
    
    double x_max=cos(phi_max)*sin(theta_max);
    double y_max=sin(phi_max)*sin(theta_max);
    double z_max=cos(theta_max);
    
    double gamma=acos(x_min*x_max+y_min*y_max+z_min*z_max); // Angle between the end points.
    double dt=dgamma*PI/(180.0*gamma);  // Angular increment normalised to 1.
    
    printf("dt=%lg\n",dt);
    
    /* March through the grid points. */
    for (double r=r_min; r<=r_max; r+=dr)
    {
        for (double t=0.0; t<=1.0; t+=dt)
        {
            phi=phi_min+t*(phi_max-phi_min);
            theta=theta_min+t*(theta_max-theta_min);
            
            x=r*cos(phi)*sin(theta);
            y=r*sin(phi)*sin(theta);
            z=r*cos(theta);
            pl.append(x,y,z);
        }
    }
    
    /* Write list into an array. */
    p=new Point[pl.n];
    n_points=pl.n;
    pl.list2array(p);
}