// g++ test.c libeltopo_release.a -llapack -lblas -lstdc++ -lm -I../common -I../eltopo3d -o  test 
//C hello world example
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#define NSIZE 2000

#include <eltopo.h>
#include <subdivisionscheme.h>


struct ElTopoMesh input; 
struct ElTopoMesh output; 
struct ElTopoGeneralOptions general_otions;
struct ElTopoStaticOperationsOptions options; 
struct ElTopoDefragInformation defrag_info;

#include "mesh.c"

int printmesh(int i, int Nverticies, int Ntriangles, float * verticies, int * triangles) {

  char buf[256];
  snprintf(buf, sizeof buf, "%s%d%s", "output/", i, ".obj");

  FILE *fp;

  fp = fopen(buf, "w+");

  //fprintf(fp, "This is testing for fprintf...\n");

  for (int i=0;i<Nverticies;i++){
    fprintf(fp,"v %f %f %f \n",verticies[3*i],verticies[3*i+1],verticies[3*i+2]);
  }

  for (int i=0;i<Ntriangles;i++) {
    fprintf(fp,"f %d %d %d\n",triangles[3*i]+1,triangles[3*i+1]+1,triangles[3*i+2]+1);
  }
    

  //  fputs("This is testing for fputs...\n", fp);
  fclose(fp);
}


void enright_velocity(float t, float * pos, float * out )
{
    float x = pos[0]; 
    float y = pos[1]; 
    float z = pos[2];

    out[0] = 2.0 * sin(M_PI*x) * sin(M_PI*x) * sin(2.0*M_PI*y) * sin(2.0*M_PI*z);
    out[1] = -sin(2.0*M_PI*x) * sin(M_PI*y)*sin(M_PI*y) * sin(2.0*M_PI*z);
    out[2] = -sin(2.0*M_PI*x) * sin(2.0*M_PI*y) * sin(M_PI*z) * sin(M_PI*z);

    out[0] *= sin(M_PI * t * 2 / 3);
    out[1] *= sin(M_PI * t * 2 / 3);
    out[2] *= sin(M_PI * t * 2 / 3);
}

void set_predicted_vertex_positions(int Nverticies, float * current_positions,float * predicted_positions, float current_t, float dt )
{   

    float tmp[3];
    float v[3];
    float k1[3];
    float k2[3];
    float k3[3];
    float k4[4];

    float x,y,z;
    
    for(int i=0; i < Nverticies; ++i) 
    {

      x = current_positions[3*i];
      y = current_positions[3*i+1];
      z = current_positions[3*i+2];

      // RK4
      // k1 = dt * f( t, x );
      tmp[0] = x;
      tmp[1] = y;
      tmp[2] = z;
      enright_velocity( current_t, tmp, v);
      k1[0] = dt * v[0];
      k1[1] = dt * v[1];
      k1[2] = dt * v[2];
        
      // k2 = dt * f( t + 0.5*dt, x + 0.5*k1 );
      tmp[0] = x + 0.5 * k1[0];
      tmp[1] = y + 0.5 * k1[1];
      tmp[2] = z + 0.5 * k1[2];
      enright_velocity( current_t + 0.5*dt, tmp, v );
      k2[0] = dt * v[0];
      k2[1] = dt * v[1];
      k2[2] = dt * v[2];

      // k3 = dt * f( t + 0.5*dt, x + 0.5*k2 );
      tmp[0] = x + 0.5 * k2[0];
      tmp[1] = y + 0.5 * k2[1];
      tmp[2] = z + 0.5 * k2[2];
      enright_velocity( current_t + 0.5*dt, tmp, v );
      k3[0] = dt * v[0];
      k3[1] = dt * v[1];
      k3[2] = dt * v[2];

      // k4 = dt * f( t + dt, x + k3 );
      tmp[0] = x + k3[0];
      tmp[1] = y + k3[1];
      tmp[2] = z + k3[2];
      enright_velocity( current_t + dt, tmp, v );
      k4[0] = dt * v[0];
      k4[1] = dt * v[1];
      k4[2] = dt * v[2];


      predicted_positions[3*i] = x + 1./6. * ( k1[0] + k4[0] ) + 1./3. * ( k2[0] + k3[0] );
      predicted_positions[3*i+1] = y + 1./6. * ( k1[1] + k4[1] ) + 1./3. * ( k2[1] + k3[1] );
      predicted_positions[3*i+2] = z + 1./6. * ( k1[2] + k4[2] ) + 1./3. * ( k2[2] + k3[2] );
      
    }
    
}


int main()
{

/*   surftrack_parameters */
/* ( */
/*    use_fraction 1  !!!!!!!!!!!!!!
/*    min_edge_length 0.5 */
/*    max_edge_length 1.5 */
/*    max_volume_change 0.1 */
/*    min_curvature_multiplier 1.0 
/*    max_curvature_multiplier 1.0 */
/*    merge_proximity 0.001    !!!!!!
/*    repulsion_proximity 0.0001 !!!!!!!!
/*    friction_coefficient 0.0 */
/*    perform_improvement 1 */
/*    allow_topology_changes 0 */
/*    collision_safety 1 */
/*    subdivision_scheme butterfly 
/* ) */


  general_otions.m_verbose = 0;
  general_otions.m_collision_safety = 1;
//  general_otions.m_proximity_epsilon = 1e-3; //!!!!!!!!!

  options.m_perform_improvement = true;
  options.m_allow_topology_changes = false;
  options.m_max_volume_change = 0.1;
  options.m_min_edge_length = 0.5;
  options.m_max_edge_length = 1.5;
//  options.m_min_triangle_area = ;
//  options.m_min_triangle_angle = ;
//  options.m_max_triangle_angle = ;
  options.m_use_curvature_when_splitting = true;
  options.m_use_curvature_when_collapsing = true;
  options.m_min_curvature_multiplier = 1.0;
  options.m_max_curvature_multiplier = 1.0;
//  options.m_allow_vertex_movement = ;
//  options.m_edge_flip_min_length_change = ;
//  options.m_merge_proximity_epsilon = ;
//  options.m_subdivision_scheme = ; // buterfly
  options.m_subdivision_scheme = new ButterflyScheme();
  options.m_collision_safety = true;
  options.m_allow_non_manifold = false;

  
  printmesh(0,Nverticies,Ntriangles,verticies,triangles);

  float predicted[NSIZE];
  float dt = 0.1;
  float t = 0;

  for (int j=1;j<10;j++) {

    set_predicted_vertex_positions(Nverticies,verticies,predicted,t,dt);
    for(int i=0;i<Nverticies*3;i++) {
      verticies[i] = predicted[i];
    }

    printmesh(j,Nverticies,Ntriangles,verticies,triangles);
    t += dt;
  }
  
  printf("Hello world\n");
  return 0;
}
