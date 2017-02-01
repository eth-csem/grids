#include <stdio.h>
#include <string.h>
#include <math.h>
#include "auxiliary.h"

/*! \file
 \brief Definitions for auxiliary classes and functions to compute point clouds.
 */

const double PI = 3.14159265358979323846264338327;

/*=================================================================================*/
/* Point in 3D. -------------------------------------------------------------------*/
/*=================================================================================*/

Point::Point(double x_in, double y_in, double z_in)
{
    x=x_in;
    y=y_in;
    z=z_in;
    
    next=0;
}

Point::Point()
{
    x=0.0;
    y=0.0;
    z=0.0;
    
    next=0;
}


void Point::print() const
{
    printf("x=%lg, y=%lg, z=%lg\n",x,y,z);
}

void Point::set(double x_in, double y_in, double z_in)
{
    x=x_in;
    y=y_in;
    z=z_in;
}

bool operator==(const Point a, const Point b)
{
    return (a.x==b.x && a.y==b.y && a.z==b.z);
}

bool operator!=(const Point a, const Point b)
{
    return !(a==b);
}

bool operator>(const Point a, const Point b)
{
    return (a.x>b.x || (a.x==b.x && a.y>b.y) || (a.x==b.x && a.y==b.y && a.z>b.z));
}

bool operator<(const Point a, const Point b)
{
    return !(a>b || a==b);
}

/*=================================================================================*/
/* List of points. ----------------------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/

Pointlist::Pointlist()
{
    n=0;
    start=0;
}

/* Destructor. Clean up from the list anchor onwards. -----------------------------*/

Pointlist::~Pointlist()
{
    Point *temp;
    
    while (start)
    {
        temp=start->next;
        delete start;
        start=temp;
    }
}

/* Append a new list element. -----------------------------------------------------*/

void Pointlist::append(Point p)
{
    Point *pt;
    pt=new Point;

    n++;
    
    pt->x=p.x;
    pt->y=p.y;
    pt->z=p.z;

    pt->next=start;
    start=pt;
}

void Pointlist::append(double x, double y, double z)
{
    Point *pt;
    pt=new Point;
    
    n++;
    
    pt->x=x;
    pt->y=y;
    pt->z=z;
    
    pt->next=start;
    start=pt;
}

/* Print the list to screen. ------------------------------------------------------*/

void Pointlist::print()
{
    printf("number of list elements: %d\n",n);
    for (Point *p=start; p; p=p->next) p->print();
}

/* Write list to an array. --------------------------------------------------------*/

void Pointlist::list2array(Point *p)
{
    for (Point *pt=start; pt; pt=pt->next)
    {
        *p=*pt;
        p++;
    }
}

/*=================================================================================*/
/* Refinement regions. ------------------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/

Refinement_Regions::Refinement_Regions()
{
    n=0;
    n_points=0;
    xmin=0; xmax=0; ymin=0; ymax=0; zmin=0; zmax=0;
    d_horizontal=0; d_vertical=0;
    
    phi=0;
    b[0]=0; b[1]=0; b[2]=0;
}

/* Destructor. --------------------------------------------------------------------*/

Refinement_Regions::~Refinement_Regions()
{
    if (n_points!=0) delete[] n_points;
    
    if (xmin!=0) delete[] xmin;
    if (xmax!=0) delete[] xmax;
    if (ymin!=0) delete[] ymin;
    if (ymax!=0) delete[] ymax;
    if (zmin!=0) delete[] zmin;
    if (zmax!=0) delete[] zmax;
    if (d_horizontal!=0) delete[] d_horizontal;
    if (d_vertical!=0) delete[] d_vertical;
    
    if (phi!=0) delete[] phi;
    if (b[0]!=0) delete[] b[0];
    if (b[1]!=0) delete[] b[1];
    if (b[2]!=0) delete[] b[2];
}

/* Read refinement regions for cubed sphere/ball from file. -----------------------*/

void Refinement_Regions::read_cubed(const char* refinement_list)
{
    /* Local variables. -----------------------------------------------------------*/
    
    FILE *fid_refine, *fid_list;
    char str[1000], filename[1000];
    int number_of_regions, number_of_subregions, count, face;

    /* Read to get number of subregions. ------------------------------------------*/
    
    fid_list=fopen(refinement_list,"r");
    
    if (fid_list)
    {
        /* Read the list of refinement regions. */
        
        fgets(str,1000,fid_list);
        fscanf(fid_list,"%d",&number_of_regions);
        fgets(str,1000,fid_list);
        
        /* March through refinement regions to read number of subregions. */
        for (int k=0; k<number_of_regions; k++)
        {
            fgets(filename,1000,fid_list);
            remove_newline(filename);
            
            /* Open and read refinement region file. */
            fid_refine=fopen(filename,"r");
            if (fid_refine)
            {
                fgets(str,1000,fid_refine);
                fscanf(fid_refine,"%d",&number_of_subregions);
                n+=number_of_subregions;
                fclose(fid_refine);
            }
            else
            {
                printf("Could not open file %s.\n",str);
            }
        }
        fclose(fid_list);
    }
    else
    {
        printf("file could not be opened.\n");
    }
    
    /* Allocate memory. -----------------------------------------------------------*/

    xmin=new double[n];
    xmax=new double[n];
    ymin=new double[n];
    ymax=new double[n];
    zmin=new double[n];
    zmax=new double[n];
    rmin=new double[n];
    rmax=new double[n];
    d_vertical=new double[n];
    d_horizontal=new double[n];
    
    /* Read to get subregion properties. ------------------------------------------*/

    count=0;
    fid_list=fopen(refinement_list,"r");
    
    if (fid_list)
    {
        /* Read the list of refinement regions. */
        
        fgets(str,1000,fid_list);
        fscanf(fid_list,"%d",&number_of_regions);
        fgets(str,1000,fid_list);
        
        /* March through refinement regions to read subregion data. */
        for (int k=0; k<number_of_regions; k++)
        {
            fgets(filename,1000,fid_list);
            remove_newline(filename);
            
            /* Open and read refinement region file. */
            fid_refine=fopen(filename,"r");
            if (fid_refine)
            {
                fgets(str,1000,fid_refine);
                fscanf(fid_refine,"%d",&number_of_subregions);
                fgets(str,1000,fid_refine);
                
                for (int i=0; i<number_of_subregions; i++)
                {
                    fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%d",&face); fgets(str,1000,fid_refine);
                    if (face==1)
                    {
                        fscanf(fid_refine,"%lg",&ymin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&ymax[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmax[count]); fgets(str,1000,fid_refine);
                        xmin[count]=1.0;
                        xmax[count]=1.0;
                    }
                    else if (face==2)
                    {
                        fscanf(fid_refine,"%lg",&ymin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&ymax[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmax[count]); fgets(str,1000,fid_refine);
                        xmin[count]=-1.0;
                        xmax[count]=-1.0;
                    }
                    else if (face==3)
                    {
                        fscanf(fid_refine,"%lg",&xmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&xmax[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmax[count]); fgets(str,1000,fid_refine);
                        ymin[count]=1.0;
                        ymax[count]=1.0;
                    }
                    else if (face==4)
                    {
                        fscanf(fid_refine,"%lg",&xmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&xmax[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&zmax[count]); fgets(str,1000,fid_refine);
                        ymin[count]=-1.0;
                        ymax[count]=-1.0;
                    }
                    else if (face==5)
                    {
                        fscanf(fid_refine,"%lg",&xmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&xmax[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&ymin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&ymax[count]); fgets(str,1000,fid_refine);
                        zmin[count]=1.0;
                        zmax[count]=1.0;
                    }
                    else if (face==6)
                    {
                        fscanf(fid_refine,"%lg",&xmin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&xmax[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&ymin[count]); fgets(str,1000,fid_refine);
                        fscanf(fid_refine,"%lg",&ymax[count]); fgets(str,1000,fid_refine);
                        zmin[count]=-1.0;
                        zmax[count]=-1.0;
                    }
                    
                    fscanf(fid_refine,"%lg",&rmin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&rmax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&d_horizontal[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&d_vertical[count]); fgets(str,1000,fid_refine);
                    count++;
                }
                
                fclose(fid_refine);
            }
            else
            {
                printf("Could not open file %s.\n",str);
            }
        }
        fclose(fid_list);
    }
    else
    {
        printf("File could not be opened.\n");
    }
}

/* Read refinement regions for Fibonacci sphere from file. ------------------------*/

void Refinement_Regions::read_fibonacci(const char* refinement_list)
{
    /* Local variables. -----------------------------------------------------------*/
    
    FILE *fid_refine, *fid_list;
    char str[1000], filename[1000];
    int number_of_regions, number_of_subregions, count;
    
    /* Read to get number of subregions. ------------------------------------------*/
    
    fid_list=fopen(refinement_list,"r");
    
    if (fid_list)
    {
        /* Read the list of refinement regions. */
        
        fgets(str,1000,fid_list);
        fscanf(fid_list,"%d",&number_of_regions);
        fgets(str,1000,fid_list);
        
        /* March through refinement regions to read number of subregions. */
        for (int k=0; k<number_of_regions; k++)
        {
            fgets(filename,1000,fid_list);
            remove_newline(filename);
            
            /* Open and read refinement region file. */
            fid_refine=fopen(filename,"r");
            if (fid_refine)
            {
                fgets(str,1000,fid_refine);
                fscanf(fid_refine,"%d",&number_of_subregions);
                n+=number_of_subregions;
                fclose(fid_refine);
            }
            else
            {
                printf("Could not open file %s.\n",str);
            }
        }
        fclose(fid_list);
    }
    else
    {
        printf("file could not be opened.\n");
    }
    
    /* Allocate memory. -----------------------------------------------------------*/
    
    n_points=new int[n];
    
    xmin=new double[n];
    xmax=new double[n];
    ymin=new double[n];
    ymax=new double[n];
    rmin=new double[n];
    rmax=new double[n];
    d_vertical=new double[n];
    
    phi=new double[n];
    b[0]=new double[n];
    b[1]=new double[n];
    b[2]=new double[n];
    
    /* Read to get subregion properties. ------------------------------------------*/
    
    count=0;
    fid_list=fopen(refinement_list,"r");
    
    if (fid_list)
    {
        /* Read the list of refinement regions. */
        
        fgets(str,1000,fid_list);
        fscanf(fid_list,"%d",&number_of_regions);
        fgets(str,1000,fid_list);
        
        /* March through refinement regions to read subregion data. */
        for (int k=0; k<number_of_regions; k++)
        {
            fgets(filename,1000,fid_list);
            remove_newline(filename);
            
            /* Open and read refinement region file. */
            fid_refine=fopen(filename,"r");
            if (fid_refine)
            {
                fgets(str,1000,fid_refine);
                fscanf(fid_refine,"%d",&number_of_subregions);
                fgets(str,1000,fid_refine);
                
                for (int i=0; i<number_of_subregions; i++)
                {
                    fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&ymin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&ymax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&xmin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&xmax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%d",&n_points[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&rmin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&rmax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&d_vertical[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg %lg %lg",&b[0][count],&b[1][count],&b[2][count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&phi[count]); fgets(str,1000,fid_refine);
                    
                    count++;
                }

                fclose(fid_refine);
            }
            else
            {
                printf("Could not open file %s.\n",str);
            }
        }
        fclose(fid_list);
    }
    else
    {
        printf("File could not be opened.\n");
    }
}

/* Read refinement regions for regular spherical coordinates. ---------------------*/

void Refinement_Regions::read_regular(const char* refinement_list)
{
    /* Local variables. -----------------------------------------------------------*/
    
    FILE *fid_refine, *fid_list;
    char str[1000], filename[1000];
    int number_of_regions, number_of_subregions, count;
    
    /* Read to get number of subregions. ------------------------------------------*/
    
    fid_list=fopen(refinement_list,"r");
    
    if (fid_list)
    {
        /* Read the list of refinement regions. */
        
        fgets(str,1000,fid_list);
        fscanf(fid_list,"%d",&number_of_regions);
        fgets(str,1000,fid_list);
        
        /* March through refinement regions to read number of subregions. */
        for (int k=0; k<number_of_regions; k++)
        {
            fgets(filename,1000,fid_list);
            remove_newline(filename);
            
            /* Open and read refinement region file. */
            fid_refine=fopen(filename,"r");
            if (fid_refine)
            {
                fgets(str,1000,fid_refine);
                fscanf(fid_refine,"%d",&number_of_subregions);
                n+=number_of_subregions;
                fclose(fid_refine);
            }
            else
            {
                printf("Could not open file %s.\n",str);
            }
        }
        fclose(fid_list);
    }
    else
    {
        printf("file could not be opened.\n");
    }
    
    /* Allocate memory. -----------------------------------------------------------*/
    
    n_points=new int[n];
    
    xmin=new double[n];
    xmax=new double[n];
    ymin=new double[n];
    ymax=new double[n];
    zmin=new double[n];
    zmax=new double[n];
    d_vertical=new double[n];
    d_horizontal=new double[n];
    
    phi=new double[n];
    b[0]=new double[n];
    b[1]=new double[n];
    b[2]=new double[n];
    
    /* Read to get subregion properties. ------------------------------------------*/
    
    count=0;
    fid_list=fopen(refinement_list,"r");
    
    if (fid_list)
    {
        /* Read the list of refinement regions. */
        
        fgets(str,1000,fid_list);
        fscanf(fid_list,"%d",&number_of_regions);
        fgets(str,1000,fid_list);
        
        /* March through refinement regions to read subregion data. */
        for (int k=0; k<number_of_regions; k++)
        {
            fgets(filename,1000,fid_list);
            remove_newline(filename);
            
            /* Open and read refinement region file. */
            fid_refine=fopen(filename,"r");
            if (fid_refine)
            {
                fgets(str,1000,fid_refine);
                fscanf(fid_refine,"%d",&number_of_subregions);
                fgets(str,1000,fid_refine);
                
                for (int i=0; i<number_of_subregions; i++)
                {
                    fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&ymin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&ymax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&xmin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&xmax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&zmin[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&zmax[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&d_horizontal[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&d_vertical[count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg %lg %lg",&b[0][count],&b[1][count],&b[2][count]); fgets(str,1000,fid_refine);
                    fscanf(fid_refine,"%lg",&phi[count]); fgets(str,1000,fid_refine);
                    
                    count++;
                }
                
                fclose(fid_refine);
            }
            else
            {
                printf("Could not open file %s.\n",str);
            }
        }
        fclose(fid_list);
    }
    else
    {
        printf("File could not be opened.\n");
    }
}


/* Print refinement regions for cubed sphere on screen. ---------------------------*/

void Refinement_Regions::print_cubed()
{
    printf("Total number of refinement regions: %d\n",n);
    printf("--------------------------------------\n");
    
    for (int i=0; i<n; i++)
    {
        printf("xmin=%lg, xmax=%lg\n",xmin[i],xmax[i]);
        printf("ymin=%lg, ymax=%lg\n",ymin[i],ymax[i]);
        printf("zmin=%lg, zmax=%lg\n",zmin[i],zmax[i]);
        printf("rmin=%lg, rmax=%lg\n",rmin[i],rmax[i]);
        printf("d_horizontal=%lg, d_vertical=%lg\n",d_horizontal[i],d_vertical[i]);
        printf("--------------------------------------\n");
    }
}

/* Print refinement regions for Fibonacci sphere on screen. -----------------------*/

void Refinement_Regions::print_fibonacci()
{
    printf("Total number of refinement regions: %d\n",n);
    printf("--------------------------------------\n");
    
    for (int i=0; i<n; i++)
    {
        printf("n_points=%d\n",n_points[i]);
        printf("xmin=%lg, xmax=%lg\n",xmin[i],xmax[i]);
        printf("ymin=%lg, ymax=%lg\n",ymin[i],ymax[i]);
        printf("rmin=%lg, rmax=%lg\n",rmin[i],rmax[i]);
        printf("d_vertical=%lg\n",d_vertical[i]);
        printf("b=(%lg %lg %lg)\n",b[0][i],b[1][i],b[2][i]);
        printf("%lg\n",phi[i]);
        printf("--------------------------------------\n");
    }
}

/* Print refinement regions for regular spherical grid on screen. -----------------*/

void Refinement_Regions::print_regular()
{
    printf("Total number of refinement regions: %d\n",n);
    printf("--------------------------------------\n");
    
    for (int i=0; i<n; i++)
    {
        printf("xmin=%lg, xmax=%lg\n",xmin[i],xmax[i]);
        printf("ymin=%lg, ymax=%lg\n",ymin[i],ymax[i]);
        printf("zmin=%lg, zmax=%lg\n",zmin[i],zmax[i]);
        printf("d_horizontal=%lg\n",d_horizontal[i]);
        printf("d_vertical=%lg\n",d_vertical[i]);
        printf("b=(%lg %lg %lg)\n",b[0][i],b[1][i],b[2][i]);
        printf("%lg\n",phi[i]);
        printf("--------------------------------------\n");
    }
}


/*=================================================================================*/
/* Projection from the unit cube onto the sphere of radius R [km]. ----------------*/
/*=================================================================================*/

void project2sphere(Point* p, int n_points, double R)
{
    /* Local variables. */
    double r, phi, theta;
    
    /* March through all the points of the array. */
    for (int i=0; i<n_points; i++)
    {
        /* Compute colatitude and longitude of this point. */
        r=sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
        if (r!=0)
        {
            theta=acos(p[i].z/r);
            phi=atan2(p[i].y,p[i].x);
        }
        
        /* Project onto the sphere. */
        p[i].x=R*sin(theta)*cos(phi);
        p[i].y=R*sin(theta)*sin(phi);
        p[i].z=R*cos(theta);
    }
}

Point project2sphere(double x, double y, double z, double R)
{
    /* Local variables. */
    Point p;
    double r, phi, theta;

    /* Compute colatitude and longitude of this point. */
    r=sqrt(x*x+y*y+z*z);
    if (r!=0)
    {
        theta=acos(z/r);
        phi=atan2(y,x);
    }
        
    /* Project onto the sphere. */
    p.x=R*sin(theta)*cos(phi);
    p.y=R*sin(theta)*sin(phi);
    p.z=R*cos(theta);
    
    return p;
}

/*=================================================================================*/
/* Rotation matrix for points in Cartesian coordinates. ---------------------------*/
/*=================================================================================*/

void rotation_matrix(double R[][3], double b[], double phi)
{
    double A[3][3], B[3][3], C[3][3];
    
    A[0][0]=b[0]*b[0]; A[0][1]=b[0]*b[1]; A[0][2]=b[0]*b[2];
    A[1][0]=b[1]*b[0]; A[1][1]=b[1]*b[1]; A[1][2]=b[1]*b[2];
    A[2][0]=b[2]*b[0]; A[2][1]=b[2]*b[1]; A[2][2]=b[2]*b[2];
    
    B[0][0]=1.0; B[0][1]=0.0; B[0][2]=0.0;
    B[1][0]=0.0; B[1][1]=1.0; B[1][2]=0.0;
    B[2][0]=0.0; B[2][1]=0.0; B[2][2]=1.0;
    
    C[0][0]=0.0; C[0][1]=-b[2]; C[0][2]=b[1];
    C[1][0]=b[2]; C[1][1]=0.0; C[1][2]=-b[0];
    C[2][0]=-b[1]; C[2][1]=b[0]; C[2][2]=0.0;
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            R[i][j]=(1.0-cos(phi))*A[i][j]+cos(phi)*B[i][j]+sin(phi)*C[i][j];
        }
    }
}

/*=================================================================================*/
/* Sorting algorithms for 3D points. ----------------------------------------------*/
/*=================================================================================*/

/* Quicksort algorithm, recursive part. -------------------------------------------*/

void quicksort(Point* a, unsigned int start, unsigned int end)
{
    unsigned int idx;
    
    if (start<end)
    {
        idx=partition(a,start,end);
        if (idx>0) quicksort(a,start,idx-1);
        if (idx<end) quicksort(a,idx+1,end);
    }
}

/* Quicksort algorithm, partitioning part. ----------------------------------------*/

unsigned int partition(Point* a, unsigned int start, unsigned int end)
{
    Point dummy, pivot=a[end];
    unsigned int i=start;
    
    for (int j=start; j<end; j++)
    {
        if (a[j]<pivot)
        {
            dummy=a[i];
            a[i]=a[j];
            a[j]=dummy;
            i++;
        }
    }
    
    dummy=a[end];
    a[end]=a[i];
    a[i]=dummy;
    
    return i;
}


/*=================================================================================*/
/* Remove newline character from end of a string. ---------------------------------*/
/*=================================================================================*/

/* Remove newline character from a string. ----------------------------------------*/

void remove_newline(char *line)
{
    int new_line = strlen(line) -1;
    if (line[new_line] == '\n') line[new_line] = '\0';
}

/*=================================================================================*/
/* Compute number of grid points within a refinement region. ----------------------*/
/*=================================================================================*/
int number_of_grid_points(Refinement_Regions &r, int index)
{
    int nx, ny, nz;
    
    nx=floor((r.xmax[index]-r.xmin[index])/r.d_horizontal[index]);
    ny=floor((r.ymax[index]-r.ymin[index])/r.d_horizontal[index]);
    nz=floor((r.zmax[index]-r.zmin[index])/r.d_horizontal[index]);
    
    if (nx==0) return ny*nz;
    else if (ny==0) return nx*nz;
    else if (nz==0) return nx*ny;
    else return 0;
}

