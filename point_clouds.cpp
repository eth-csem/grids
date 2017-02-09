#include <stdio.h>
#include <math.h>
#include "point_clouds.h"

/*! \file
 \brief Definition of the point cloud classes.
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

void PointCloud::write(const char* filename) const
{
    FILE* pfile;
    
    pfile=fopen(filename,"w");
    
    if (p!=0)
    {
        fprintf(pfile, "%lg %lg %lg\n",p[0].x,p[0].y,p[0].z);
        for (int i=1; i<n_points; i++)
        {
            if (p[i]!=p[i-1]) fprintf(pfile, "%lg %lg %lg\n",p[i].x,p[i].y,p[i].z);
        }
    }

    fclose(pfile);
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

    /* Sort to later avoid duplications. */
    if (pl.n>1)quicksort(p,0,pl.n-1);
    
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
    
    /* Sort to later avoid duplications. */
    if (pl.n>1) quicksort(p,0,pl.n-1);
    
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
                rotation_matrix(R,n,ref.phi[idx]);
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
                        pt->next=pt->next->next;
                        delete temp;
                        pl.n--;
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
            rotation_matrix(R,n,ref.phi[idx]);
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
                    pt->next=pt->next->next;
                    delete temp;
                    pl.n--;
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
    
    /* Sort to later avoid duplications. */
    if (pl.n>1) quicksort(p,0,pl.n-1);
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
    
    /* Sort to later avoid duplications. */
    if (pl.n>1) quicksort(p,0,pl.n-1);
}
