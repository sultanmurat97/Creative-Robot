//
//  main.cpp
//  Assignment3
//
//  Created by Sultanmurat on 10/26/19.
//  Copyright © 2019 Sultanmurat. All rights reserved.
//
#include <GLUT/glut.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <array>
#include"/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/gnuplot_i.hpp"
using namespace std;
std::random_device rd; // obtain a random number from hardware
std::mt19937 eng(rd()); // seed the generator
std::uniform_real_distribution<> distr01(0.0, 1.0);

double overall_mass = 80;
int number_regions = 4;
double initial_mass = 0.1;
const double g = -9.81;           // gravity in z, [meter/ s^2]
const double dt = 0.0001;         // simulation time-step
const double mu = 1;           // coefficient of friction static
const double mu_k = 0.8;
const double k_ground = 100000;
const double damping = 0.999;
double T = 0;               // global time variable

bool dampeningOn = true;
bool breathingOn = true;
bool newtonianFriction = true;
bool shadedCube = true;



double w = 10;
const int dim = 5;
double errors[10000000] = {0};
int number_population = 10;
int number_generations = 100;

struct Force{
        double x;
        double y;
        double z;
};
struct Mass{
        double m; // mass, [kg]
        double px, py, pz;  // position, [meter]
        double vx, vy, vz;  // velocity, [meter/s]
        double ax, ay, az;  // acceleration, [meters/s^2]
};
struct Spring{
    double k;   // spring constant, [N/m]
    double L0;  // original rest length, [meters]
    int m1, m2; // indices of two masses that are connected by this spring
    double a, b, c;
};

struct Coef{
    double x;
    double y;
    double z;
    double r;
    double b;
    double c;
    double k;
};

class Individual2{
    public:
    //bool voxels[dim][dim][dim];
    vector <Coef> coefs;
};

double distanceCoefs(Coef c1, Coef c2);

bool isValidRobot(Individual2 ind);

vector<Mass> masses;
vector<Spring> springs;
vector<int> imMasses;

double normVector(double x, double y, double z);
void wait_for_key();
double distance(Mass m1, Mass m2);

void plotEnergy(void);
void plotLearningCurve(void);
void plotDot(void);
void plotDiversity(void);

ofstream kineticEnergyFile;
ofstream potentialEnergyFile;
ofstream totalEnergyFile;
bool theSame(Mass m1, Mass m2);

void animate(void){
    for (int u = 0; u < 100; u++){
        vector<Force> forces;
        double kinetic_energy = 0;
        double potential_energy = 0;
        double total_energy = 0;
        for (int i = 0; i < springs.size(); i++){
            // breathing
            if (breathingOn == true){
                    double a = springs[i].a;
                    double b = springs[i].b;
                    double c = springs[i].c;
                    springs[i].L0 = a+b*sin(T*w+c);
                }
        }
        for (int i = 0; i < masses.size(); i++)
        {
            Force current_force;
            current_force.x = 0; current_force.y = 0; current_force.z = 0;
            double F_c = 0;
            double F_h = 0;
            double F_v = 0;
            
            if (masses[i].pz < 0)
            {
                F_c = -k_ground * masses[i].pz;
                current_force.z = current_force.z + F_c;
                current_force.y = 0;
                current_force.x = 0;
                
                if (newtonianFriction == true){
                    F_h = sqrt( pow(current_force.x,2) + pow(current_force.y, 2) );
                    F_v = current_force.z;
                    if (F_h < mu*F_v){
                        masses[i].vx = 0;
                        masses[i].vy = 0;
                    } else {
                        current_force.x = current_force.x - F_v*mu_k;
                        current_force.y = current_force.y - F_v*mu_k;
                    }
                }
                
                potential_energy += abs(F_c*masses[i].pz);
            }
            
            for (int j = 0; j < springs.size(); j++)
            {
                if (springs[j].m1 == i || springs[j].m2 == i)
                {
                    Force spring_force;
                    spring_force.x = 0; spring_force.y = 0; spring_force.z = 0;
                    
                    double current_length = distance(masses[springs[j].m1], masses[springs[j].m2]);
                    double F = springs[j].k * (current_length - springs[j].L0);
                    
                    if (springs[j].m1 != i)
                    {
                        spring_force.x = F*(masses[springs[j].m1].px - masses[springs[j].m2].px)/current_length;
                        spring_force.y = F*(masses[springs[j].m1].py - masses[springs[j].m2].py)/current_length;
                        spring_force.z = F*(masses[springs[j].m1].pz - masses[springs[j].m2].pz)/current_length;
                    } else
                    {
                        spring_force.x = F*(masses[springs[j].m2].px - masses[springs[j].m1].px)/current_length;
                        spring_force.y = F*(masses[springs[j].m2].py - masses[springs[j].m1].py)/current_length;
                        spring_force.z = F*(masses[springs[j].m2].pz - masses[springs[j].m1].pz)/current_length;
                        
                    }
                    double delta_x = abs(current_length-springs[j].L0);
                    potential_energy += 0.5*springs[j].k*pow(delta_x, 2);
                    
                    current_force.x = current_force.x + spring_force.x;
                    current_force.y = current_force.y + spring_force.y;
                    current_force.z = current_force.z + spring_force.z;
                    
                    
                }
            }
            double gravity = masses[i].m*g;
            current_force.z = current_force.z + gravity;
        
            forces.push_back(current_force);
        }
        for (int i = 0; i < masses.size(); i++)
        {
            // velocity dampening
            if (dampeningOn == true){
                masses[i].vx *= damping;
                masses[i].vy *= damping;
                masses[i].vz *= damping;
            }
            
            // updating acceleration
            masses[i].ax = forces[i].x/masses[i].m;
            masses[i].ay = forces[i].y/masses[i].m;
            masses[i].az = forces[i].z/masses[i].m;
            
            // updating velocity
            masses[i].vx = masses[i].vx + masses[i].ax*dt;
            masses[i].vy = masses[i].vy + masses[i].ay*dt;
            masses[i].vz = masses[i].vz + masses[i].az*dt;
            
            // updating positions
            masses[i].px = masses[i].px + masses[i].vx*dt;
            masses[i].py = masses[i].py + masses[i].vy*dt;
            masses[i].pz = masses[i].pz + masses[i].vz*dt;
            
            //cout << masses[i].vx << " " << masses[i].vy <<  " " << masses[i].vz << endl;
            
            kinetic_energy = kinetic_energy + masses[i].m*pow(normVector(masses[i].vx, masses[i].vy, masses[i].vz),2)/2;
            potential_energy = potential_energy + abs(masses[i].m*g*masses[i].pz);
            
        }
        
        // visualize
        // draw cube
        glutPostRedisplay();
        
        // output energy values
        total_energy = kinetic_energy + potential_energy;
        
        kineticEnergyFile << T << " " << kinetic_energy << endl;
        potentialEnergyFile << T << " " <<potential_energy << endl;
        totalEnergyFile << T << " " << total_energy << endl;
        
        
        T = T + dt;
    }
}
void display(void){
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Draw a white grid "floor" for the tetrahedron to sit on.
     
    for (GLfloat i = -2; i <= 2; i += 0.1) {
         glColor3f(1.0, 1.0, 1.0);
     glBegin(GL_LINES);
         glVertex3f(i, 2, 0); glVertex3f(i, -2, 0);
               glVertex3f(2, i, 0); glVertex3f(-2, i, 0);
         glEnd();
        
         glColor3f(0.22, 0.78, 0.82);
         
         glBegin(GL_POLYGON);
       glVertex3f(i, 2, 0); glVertex3f(i, -2, 0);
       glVertex3f(2, i, 0); glVertex3f(-2, i, 0);
         glEnd();

     }
    
    
     //glEnd();
   
    // draw cube
   glPushMatrix();
      glRotatef(0, 1, 0, 0); //rotate alpha around the x axis
      glRotatef(0, 0, 1, 0); //rotate beta around the y axis
      glRotatef(0, 0, 0, 1); //rotate gamma around the z axis
    for (int i = 0; i < masses.size(); i++){
            
                glColor3f(1,1,0);
            GLUquadric *quad;
            quad = gluNewQuadric();
            glTranslatef(masses[i].px,masses[i].py,masses[i].pz);
            gluSphere(quad,0.01,100,20);
            glTranslatef(-masses[i].px, -masses[i].py, -masses[i].pz);
            
    }
     
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    for (int i = 0; i < springs.size(); i++)
    {
        
        glVertex3f(masses[springs[i].m1].px, masses[springs[i].m1].py, masses[springs[i].m1].pz);
        glVertex3f(masses[springs[i].m2].px, masses[springs[i].m2].py, masses[springs[i].m2].pz);
        
    }
    glEnd();
    
     
    if (shadedCube == true){
        glColor3f(1.0, 0.5, 0.0);
        glBegin(GL_TRIANGLES);
        //cout << masses.size() << endl;
        for (int i = 0; i < imMasses.size(); i=i+8){
            for (int k1 = i; k1 < i+6; k1++){
                for (int k2 = k1 +1; k2 < i+7; k2++){
                    for (int k3 = k2+1; k3 < i+8; k3++){
                        glVertex3f(masses[imMasses[k1]].px, masses[imMasses[k1]].py, masses[imMasses[k1]].pz);
                        glVertex3f(masses[imMasses[k2]].px, masses[imMasses[k2]].py, masses[imMasses[k2]].pz);
                        glVertex3f(masses[imMasses[k3]].px, masses[imMasses[k3]].py, masses[imMasses[k3]].pz);
                    }
                }
            }
            
        }
        glEnd();
    }
    
    
    glPopMatrix();
    
   
  
  glutSwapBuffers();
}
void init(void){
    glLineWidth(10);
    glPointSize(100.0);
 

  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);
        
    
  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 60.0,
    /* aspect ratio */ 1.0,
    /* Z near */ 0.1, /* Z far */ 4.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(2.5, 2.5, 1.5,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    -1, -1, 0.0);      /* up is in positive Y direction */

}
void displayRobot(void){
    
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowPosition(80, 80);
    glutInitWindowSize(800, 1000);
    glutCreateWindow("Cube");
    init();
    glutIdleFunc(animate);
    glutDisplayFunc(display);


    glutMainLoop();
}


void createCube(double x, double y, double z, vector<Coef> coefs){
    int n = masses.size();
    int n2 = springs.size();
    Coef tempCoef;
    tempCoef.x = x; tempCoef.y = y; tempCoef.z = z;
    vector<Coef> goodCoefs;
    
   
    for (int i = 0; i < coefs.size(); i++){
        if (distanceCoefs(tempCoef, coefs[i]) <= coefs[i].r){
            goodCoefs.push_back(coefs[i]);
        }
    }
    if (goodCoefs.size() == 0) return;
    else {
        double min_distance = distanceCoefs(tempCoef, goodCoefs[0]);
        Coef bestCoef = goodCoefs[0];
        for (int i = 0; i < goodCoefs.size(); i++){
            if (distanceCoefs(tempCoef, goodCoefs[i]) < min_distance){
                min_distance = distanceCoefs(tempCoef, goodCoefs[i]);
                bestCoef = goodCoefs[i];
            }
        }
    
        vector<Mass> tempMasses;
        Mass m1, m2, m3, m4, m5, m6, m7, m8;
        m1.m = initial_mass;
        m1.px = x; m1.py = y; m1.pz = z;
        
        m2.m = initial_mass;
        m2.px = x+0.1; m2.py = y+0; m2.pz = z;
        
        m3.m = initial_mass;
        m3.px = x+0.1; m3.py = y+0.1; m3.pz = z;
        
        m4.m = initial_mass;
        m4.px = x+0; m4.py = y+0.1; m4.pz = z;
        
        m5.m = initial_mass;
        m5.px = x+0; m5.py = y+0; m5.pz = z+0.1;
        
        m6.m = initial_mass;
        m6.px = x+0.1; m6.py = y+0; m6.pz = z+0.1;
        
        m7.m = initial_mass;
        m7.px = x+0.1; m7.py = y+0.1; m7.pz = z+0.1;
        
        m8.m = initial_mass;
        m8.px = x+0; m8.py = y+0.1; m8.pz = z+0.1;
        
        tempMasses.push_back(m1);
        tempMasses.push_back(m2);
        tempMasses.push_back(m3);
        tempMasses.push_back(m4);
        tempMasses.push_back(m5);
        tempMasses.push_back(m6);
        tempMasses.push_back(m7);
        tempMasses.push_back(m8);
        
        vector<int> coinciding;
        for (int j = 0; j < tempMasses.size(); j++){
            for (int i = 0; i < masses.size(); i++){
                if (theSame(tempMasses[j], masses[i]) == true ){
                    tempMasses[j].m = 0;
                    coinciding.push_back(i);
                    imMasses.push_back(i);
                }
            }
        }
        int count = 0;
        for (int i = 0; i < tempMasses.size(); i++){
            if (tempMasses[i].m != 0) {
                masses.push_back(tempMasses[i]);
                imMasses.push_back(n+count);
                count++;
            }
        }
        
        for (int i = n; i < masses.size()-1; i++)
        {
            for (int j = i+1; j <masses.size(); j++)
            {
                Spring s;
                s.k = bestCoef.k;
                s.m1 = j;
                s.m2 = i;
                s.L0 = distance(masses[i], masses[j]);
                springs.push_back(s);
            }
        }
        
        for (int i = 0; i < coinciding.size(); i++){
            for (int j = n; j < masses.size(); j++){
                Spring s;
                s.k = bestCoef.k;
                s.m1 = coinciding[i];
                s.m2 = j;
                s.L0 = distance(masses[s.m1], masses[s.m2]);
                springs.push_back(s);
            }
        }
        
        for (int i = n2; i < springs.size(); i++){
            springs[i].a = springs[i].L0;
            springs[i].b = bestCoef.b;
            springs[i].c = bestCoef.c;
        }
        for (int i = n; i < masses.size(); i++)
        {
            masses[i].vx = 0;
            masses[i].vy = 0;
            masses[i].vz = 0;
            masses[i].ax = 0;
            masses[i].ay = 0;
            masses[i].az = 0;
        }
    }
}
void deleteCube(void){
   
    int n1 = springs.size();
    for (int i = 0; i < n1; i++){
        springs.pop_back();
    }

    int n2 = masses.size();
    for (int i = 0; i < n2; i++){
        masses.pop_back();
    }
}

void createRobot(Individual2 ind2){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            for (int k = 0; k < dim; k++){
                createCube((double)(i)/10.0, (double)(j)/10.0, (double)(k)/10.0, ind2.coefs);
            }
        }
    }
    double eachMass = overall_mass/masses.size();
    for (int i = 0; i < masses.size(); i++){
        masses[i].m = eachMass;
    }
}

Individual2 mutateRobot(Individual2 ind){
    // mutation of parameters
    for (int i = 0; i < ind.coefs.size(); i++){
        ind.coefs[i].b = (0.95 + 0.1*distr01(eng) ) * ind.coefs[i].b;
        ind.coefs[i].c = (0.95 + 0.1*distr01(eng) ) * ind.coefs[i].c;
        
        ind.coefs[i].x = (0.95 + 0.1*distr01(eng) ) * ind.coefs[i].x;
        ind.coefs[i].y = (0.95 + 0.1*distr01(eng) ) * ind.coefs[i].y;
        ind.coefs[i].z = (0.95 + 0.1*distr01(eng) ) * ind.coefs[i].z;
        ind.coefs[i].r = (0.95 + 0.1*distr01(eng) ) * ind.coefs[i].r;
        if (ind.coefs[i].b > 0.1) ind.coefs[i].b = 0.1;
        if (ind.coefs[i].b < -0.1) ind.coefs[i].b = -0.1;
        if (ind.coefs[i].x > dim*0.1) ind.coefs[i].x = dim*0.1;
        if (ind.coefs[i].y > dim*0.1) ind.coefs[i].y = dim*0.1;
        if (ind.coefs[i].z > dim*0.1) ind.coefs[i].z = dim*0.1;
        if (ind.coefs[i].x < 0) ind.coefs[i].x = 0.1;
        if (ind.coefs[i].y < 0) ind.coefs[i].y = 0.1;
        if (ind.coefs[i].z < 0) ind.coefs[i].z = 0.1;
        
        if (ind.coefs[i].r < dim*0.1/3) ind.coefs[i].r = dim*0.1/3;
    }
    
    return ind;
}
vector<Individual2> crossoverRobot(Individual2 parent1, Individual2 parent2){
    int cutPoint1 = 1;
    int cutPoint2 = 3;
    vector<Individual2> children(2);
    if (parent1.coefs.size() != parent2.coefs.size()){
        cout << "Parents sizes are not equal!";
        return children;
    }
    for (int i = 0; i < cutPoint1; i++){
        children[0].coefs.push_back(parent1.coefs[i]);
        children[1].coefs.push_back(parent2.coefs[i]);
    }
    for (int i = cutPoint1; i < cutPoint2; i++){
        children[0].coefs.push_back(parent2.coefs[i]);
        children[1].coefs.push_back(parent1.coefs[i]);
    }
    for (int i = cutPoint2; i < parent1.coefs.size(); i++){
        children[0].coefs.push_back(parent1.coefs[i]);
        children[1].coefs.push_back(parent2.coefs[i]);
    }
    
//    int x1 = (int) (distr01(eng)*dim);
//    int x2;
//    if (x1 < dim/2) x2 = x1 + (int) (distr01(eng)*dim/2);
//    else x2 = x1 - (int) (distr01(eng)*dim/2);
//    if (x1 > x2) swap(x1,x2);
//
//    int y1 = (int) (distr01(eng)*dim);
//    int y2;
//    if (y1 < dim/2) y2 = y1 + (int) (distr01(eng)*dim/2);
//    else y2 = y1 - (int) (distr01(eng)*dim/2);
//    if (y1 > y2) swap(y1,y2);
//
//    int z1 = (int) (distr01(eng)*dim);
//    int z2;
//    if (z1 < dim/2) z2 = z1 + (int) (distr01(eng)*dim/2);
//    else z2 = z1 - (int) (distr01(eng)*dim/2);
//    if (z1 > z2) swap(z1,z2);
//
//    for (int i = 0; i < dim; i++){
//        for (int j = 0; j < dim; j++){
//            for (int k = 0; k < dim; k++){
//                if (i >= x1 && i <=x2 && j >= y1 && j <= y2 && k >= z1 && k <= z2  ){
//                    children[0].voxels[i][j][k] = parent2.voxels[i][j][k];
//                    children[1].voxels[i][j][k] = parent1.voxels[i][j][k];
//                } else {
//                    children[0].voxels[i][j][k] = parent1.voxels[i][j][k];
//                    children[1].voxels[i][j][k] = parent2.voxels[i][j][k];
//                }
//            }
//        }
//    }
    
    return children;
}

Individual2 readParametersRobot(void){
    Individual2 res;
    ifstream parametersFileInput;
    parametersFileInput.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/parametersFileRandom.txt");
    double x,y,z,r, b, c, k;
    int count = 0;
    while (count < number_regions){
        parametersFileInput >> x >> y >> z >> r >> b >> c >> k;
        Coef temp;
        temp.x = x;
        temp.y = y;
        temp.z = z;
        temp.r = r;
        temp.b = b;
        temp.c = c;
        temp.k = k;
        
        res.coefs.push_back(temp);
        count++;
    }
    parametersFileInput.close();
    return res;
}

double calculate_fitness_general(Individual2 ind2){
    createRobot(ind2);
    while (T < 10*(2*3.1415/w) ){
        vector<Force> forces;
        for (int i = 0; i < springs.size(); i++){
            // breathing
            if (breathingOn == true){
                    double a = springs[i].a;
                    double b = springs[i].b;
                    double c = springs[i].c;
                    springs[i].L0 = a+b*sin(T*w+c);
                }
        }
        
        for (int i = 0; i < masses.size(); i++)
        {
            Force current_force;
            current_force.x = 0; current_force.y = 0; current_force.z = 0;
            double F_c = 0;
            double F_h = 0;
            double F_v = 0;
            
            if (masses[i].pz <= 0)
            {
                F_c = -k_ground * masses[i].pz;
                current_force.z = current_force.z + F_c;
                current_force.y = 0;
                current_force.x = 0;
                
                if (newtonianFriction == true){
                    F_h = sqrt( pow(current_force.x,2) + pow(current_force.y, 2) );
                    F_v = current_force.z;
                    if (F_h < mu*F_v){
                        masses[i].vx = 0;
                        masses[i].vy = 0;
                    } else {
                        current_force.x = current_force.x - F_v*mu_k;
                        current_force.y = current_force.y - F_v*mu_k;
                    }
                }
                
            }
            
            for (int j = 0; j < springs.size(); j++)
            {
                if (springs[j].m1 == i || springs[j].m2 == i)
                {
                    Force spring_force;
                    spring_force.x = 0; spring_force.y = 0; spring_force.z = 0;
                    
                    double current_length = distance(masses[springs[j].m1], masses[springs[j].m2]);
                    
                        double F = springs[j].k * (current_length - springs[j].L0);
                    
                    if (springs[j].m1 != i)
                    {
                        spring_force.x = F*(masses[springs[j].m1].px - masses[springs[j].m2].px)/current_length;
                        spring_force.y = F*(masses[springs[j].m1].py - masses[springs[j].m2].py)/current_length;
                        spring_force.z = F*(masses[springs[j].m1].pz - masses[springs[j].m2].pz)/current_length;
                    } else
                    {
                        spring_force.x = F*(masses[springs[j].m2].px - masses[springs[j].m1].px)/current_length;
                        spring_force.y = F*(masses[springs[j].m2].py - masses[springs[j].m1].py)/current_length;
                        spring_force.z = F*(masses[springs[j].m2].pz - masses[springs[j].m1].pz)/current_length;
                        
                    }
                    
                    current_force.x = current_force.x + spring_force.x;
                    current_force.y = current_force.y + spring_force.y;
                    current_force.z = current_force.z + spring_force.z;
                    
                    
                }
            }
            double gravity = masses[i].m*g;
            current_force.z = current_force.z + gravity;
            
            forces.push_back(current_force);
        }
        
        for (int i = 0; i < masses.size(); i++)
        {
            // velocity dampening
            if (dampeningOn == true){
                masses[i].vx *= damping;
                masses[i].vy *= damping;
                masses[i].vz *= damping;
            }
            
            // updating acceleration
            masses[i].ax = forces[i].x/masses[i].m;
            masses[i].ay = forces[i].y/masses[i].m;
            masses[i].az = forces[i].z/masses[i].m;
            
            // updating velocity
            masses[i].vx = masses[i].vx + masses[i].ax*dt;
            masses[i].vy = masses[i].vy + masses[i].ay*dt;
            masses[i].vz = masses[i].vz + masses[i].az*dt;
            
            // updating positions
            masses[i].px = masses[i].px + masses[i].vx*dt;
            masses[i].py = masses[i].py + masses[i].vy*dt;
            masses[i].pz = masses[i].pz + masses[i].vz*dt;
        
        }
        T = T + dt;
    }
    T = 0;
    double sum = 0;
    for (int i = 0; i < masses.size(); i++){
        sum += masses[i].px;
    }
    double average_x = sum/masses.size();
    deleteCube();
    return average_x;
}

int main(int argc, char **argv) {

    //plotLearningCurve();
    //plotDot();
    //plotDiversity();
    //wait_for_key();
   
    Individual2 ind1 = readParametersRobot();
    
    cout << "Maximum fitness: " << calculate_fitness_general(ind1) << endl;
    createRobot(ind1);
    //if (isValidRobot(ind1) == true) cout << "Robot is valid!" << endl;
    //else cout << " Robot is invalid!" << endl;
    glutInit(&argc, argv);
    displayRobot();
    wait_for_key();


      // hill climber
//    ofstream parametersFileHC;
//    parametersFileHC.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/parametersFileHC.txt");
//    vector<Individual2> population;
//
//    //create initial population
//    vector<double> max_fitnesses;
//    for (int l = 0; l < number_population; l++){
//        Individual2 tempInd;
//        for (int i = 0; i < number_regions; i++){
//            Coef tempCoef;
//            tempCoef.x = distr01(eng)*dim*0.1;
//            tempCoef.y = distr01(eng)*dim*0.1;
//            tempCoef.z = distr01(eng)*dim*0.1;
//            tempCoef.r = distr01(eng)*(2*dim/3)*0.1;
//            tempCoef.b = -0.05 + distr01(eng)*0.1;
//            tempCoef.c = -0.05 + distr01(eng)*0.1;
//            tempCoef.k = 100000;
//            tempInd.coefs.push_back(tempCoef);
//        }
//
//        while ( isValidRobot(tempInd) == false){
//            for (int i = 0; i <number_regions; i++){
//                tempInd.coefs.pop_back();
//            }
//            for (int i = 0; i < number_regions; i++){
//                Coef tempCoef;
//                tempCoef.x = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.y = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.z = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.r = distr01(eng)*(dim*0.67)*0.1;
//                tempCoef.b = -0.05 + distr01(eng)*0.1;
//                tempCoef.c = -0.05 + distr01(eng)*0.1;
//                tempCoef.k = 100000;
//                tempInd.coefs.push_back(tempCoef);
//            }
//        }
//
//        population.push_back(tempInd);
//        double current_fitness = calculate_fitness_general(tempInd);
//        cout << current_fitness << endl;
//        max_fitnesses.push_back(current_fitness);
//    }
//    if (population.size() == number_population) cout << "Initial Population created!" << endl;
//
//    ofstream learningCurveHC;
//    learningCurveHC.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/learningCurveHC.txt");
//
//
//    for (int j = 0; j < number_generations; j++){
//        double total_fitness = 0;
//        double real_max = max_fitnesses[0];
//        int index_max = 0;
//        for (int i = 0; i < number_population; i++){
//            double fitness_max = max_fitnesses[i];
//            cout << fitness_max << endl;
//            Individual2 mutated = mutateRobot(population[i]);
//            double mutated_fitness = calculate_fitness_general(mutated);
//
//            if (mutated_fitness > fitness_max && isValidRobot(population[i])){
//                population[i] = mutated;
//                max_fitnesses[i] = mutated_fitness;
//            }
//
//            if (max_fitnesses[i] > real_max && isValidRobot(population[i])){
//                real_max = max_fitnesses[i];
//                index_max = i;
//            }
//
//            total_fitness += max_fitnesses[i];
//        }
//
//        // find error
//        double mean = total_fitness/number_population;
//        double error_deviation = 0;
//        for (int i = 0; i < number_population; i++){
//            error_deviation = error_deviation + sqrt((max_fitnesses[i]-mean)*(max_fitnesses[i]-mean));
//        }
//        error_deviation = error_deviation/sqrt(number_population);
//        error_deviation = error_deviation/number_population;
//        errors[j] = error_deviation;
//
//        cout << j << " " << real_max << " " << error_deviation << endl;
//        learningCurveHC << j << " " << real_max  << " " << error_deviation << endl;
//    }
//    double maxim = calculate_fitness_general(population[0]);
//    cout << maxim << endl;
//    int index_max = 0;
//    for (int i = 1; i < number_population; i++){
//        double current_fitness = calculate_fitness_general(population[i]);
//        cout << current_fitness << endl;
//        if (current_fitness > maxim){
//            maxim = current_fitness;
//            index_max = i;
//        }
//    }
//    cout << "Overall maximum: " << calculate_fitness_general(population[index_max]) << endl;
//
//     for (int i = 0 ; i < number_regions; i ++){
//            parametersFileHC  << setprecision(20) << population[index_max].coefs[i].x << " " << setprecision(20) <<   population[index_max].coefs[i].y << " " << setprecision(20) << population[index_max].coefs[i].z << " " << setprecision(20) <<   population[index_max].coefs[i].r << " " << setprecision(20) << population[index_max].coefs[i].b << " " << setprecision(20) << population[index_max].coefs[i].c << " " << setprecision(5) << population[index_max].coefs[i].k << endl;
//        }
//        parametersFileHC.close();
    
        
    
    
    // random search
//    ofstream parametersFileRandom;
//    parametersFileRandom.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/parametersFileRandom.txt");
//    ofstream randomCurve;
//    randomCurve.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/randomCurve.txt");
//    vector<Individual2> individuals;
//    int max_index = 0;
//    double fitness_max = -100;
//    for (int l = 0; l < number_population*number_generations; l++){
//        Individual2 tempInd;
//        for (int i = 0; i < number_regions; i++){
//            Coef tempCoef;
//            tempCoef.x =(dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//            tempCoef.y = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//            tempCoef.z = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//            tempCoef.r = distr01(eng)*(dim*0.67)*0.1;
//            tempCoef.b = -0.05 + distr01(eng)*0.1;
//            tempCoef.c = -0.05 + distr01(eng)*0.1;
//            tempCoef.k = 100000;
//            tempInd.coefs.push_back(tempCoef);
//        }
//        while ( isValidRobot(tempInd) == false){
//            for (int i = 0; i <number_regions; i++){
//                tempInd.coefs.pop_back();
//            }
//            for (int i = 0; i < number_regions; i++){
//                Coef tempCoef;
//                tempCoef.x = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.y = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.z = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.r = distr01(eng)*(dim*0.67)*0.1;
//                tempCoef.b = -0.05 + distr01(eng)*0.1;
//                tempCoef.c = -0.05 + distr01(eng)*0.1;
//                tempCoef.k = 100000;
//                tempInd.coefs.push_back(tempCoef);
//            }
//        }
//
//        individuals.push_back(tempInd);
//
//        double current_fitness = calculate_fitness_general(tempInd);
//
//        if (current_fitness > fitness_max && isValidRobot(tempInd)){
//            max_index = l;
//            fitness_max = current_fitness;
//        }
//        cout << (double) l/number_population << " " << fitness_max << endl;
//        randomCurve << (double) l/number_population << " " << fitness_max << endl;
//    }
//
//    for (int i = 0 ; i < number_regions; i ++){
//        //cout << individuals[max_index].coefs[i].b << " " <<  individuals[max_index].coefs[i].c << endl;
//        parametersFileRandom  << setprecision(20) << individuals[max_index].coefs[i].x << " " << setprecision(20) <<   individuals[max_index].coefs[i].y << " " << setprecision(20) << individuals[max_index].coefs[i].z << " " << setprecision(20) <<   individuals[max_index].coefs[i].r << " " <<    individuals[max_index].coefs[i].b << " " << setprecision(20) << individuals[max_index].coefs[i].c << " " << individuals[max_index].coefs[i].k <<  endl;
//    }
//    parametersFileRandom.close();
    
    // Genetic Algorithm
//    ofstream parametersFile;
//    parametersFile.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/parametersFileRobot.txt");
//
//    //create initial population
//    vector<Individual2> population;
//
//    for (int i = 0; i < number_population; i++){
//        Individual2 tempInd;
//        for (int i = 0; i < number_regions; i++){
//            Coef tempCoef;
//            tempCoef.x = distr01(eng)*(dim*0.1);
//            tempCoef.y = distr01(eng)*(dim*0.1);
//            tempCoef.z = distr01(eng)*(dim*0.1);
//            tempCoef.r = distr01(eng)*(dim*0.67)*0.1;
//            tempCoef.b = -0.05 + distr01(eng)*0.1;
//            tempCoef.c = -0.05 + distr01(eng)*0.1;
//            tempCoef.k = 100000;
//            tempInd.coefs.push_back(tempCoef);
//        }
//
//        while ( isValidRobot(tempInd) == false){
//            for (int i = 0; i <number_regions; i++){
//                tempInd.coefs.pop_back();
//            }
//            for (int i = 0; i < number_regions; i++){
//                Coef tempCoef;
//                tempCoef.x = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.y = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.z = (dim*0.1)/4 + distr01(eng)*(dim*0.1)/2;
//                tempCoef.r = distr01(eng)*(dim*0.67)*0.1;
//                tempCoef.b = -0.05 + distr01(eng)*0.1;
//                tempCoef.c = -0.05 + distr01(eng)*0.1;
//                tempCoef.k = 100000;
//                tempInd.coefs.push_back(tempCoef);
//            }
//        }
//        population.push_back(tempInd);
//    }
//    if (population.size() == number_population) cout << "Initial Population created!" << endl;
//
//    ofstream learningCurve;
//    ofstream diversity;
//    ofstream dot;
//    dot.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/dotRobot.txt");
//    diversity.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/diversityRobot.txt");
//    learningCurve.open("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/learningCurveRobot.txt");
//
//    for (int j = 0; j < number_generations; j++){
//
//        vector<Individual2> next_population;
//        array<double, 2> fitness[2*number_population];
//
//        double total_fitness = 0;
//        double min_fitness = 100;
//        double max_fitness = -100;
//        int max_index = -1;
//        int min_index = -1;
//
//        for (int i = 0; i < number_population; i++){
//            double current_fitness;
//            if (!isValidRobot(population[i]) ){
//                current_fitness = -1;
//            }
//            else {
//                current_fitness = calculate_fitness_general(population[i]);
//            }
//
//            cout << current_fitness << endl;
//            fitness[i][0] = current_fitness;
//            fitness[i][1] = i;
//            if (current_fitness > max_fitness){
//                max_fitness = current_fitness;
//                max_index = i;
//            }
//            if (current_fitness < min_fitness){
//                min_fitness = current_fitness;
//                min_index = i;
//            }
//
//            dot << j << " " << current_fitness << endl;
//
//            total_fitness += current_fitness;
//        }
//
//        // find error
//        double mean = total_fitness/number_population;
//        double error_deviation = 0;
//        for (int i = 0; i < number_population; i++){
//            error_deviation = error_deviation + sqrt((fitness[i][0]-mean)*(fitness[i][0]-mean));
//        }
//        error_deviation = error_deviation/sqrt(number_population);
//        error_deviation = error_deviation/number_population;
//        errors[j] = error_deviation;
//
//        next_population.push_back(population[max_index]);
//        // Producing new children
//        for (int m = 0; m < number_population/2; m++){
//            //Select using tournament selection
//            double k = 3; // set selection pressure
//            int r1 = -1;
//            int r2 = -1;
//            for (int i =0; i < k; i++){
//                int ind = (int) (distr01(eng)*number_population);
//                if (r1 == -1) {
//                    r1 = ind;
//                }else {
//                    double current = fitness[ind][0];
//                    double best = fitness[r1][0];
//                    //cout << ind << endl;
//                    if (current > best){
//                        r1 = ind;
//                    }
//                }
//            }
//            for (int i =0; i < k; i++){
//                int ind = (int) (distr01(eng)*number_population);
//                if (r2 == -1) {
//                    r2 = ind;
//                }else {
//                    double current = fitness[ind][0];
//                    double best = fitness[r2][0];
//                    //cout << ind << endl;
//                    if (current > best){
//                        r2 = ind;
//                    }
//                }
//            }
//
//            Individual2 parent1 = population[r1];
//            Individual2 parent2 = population[r2];
//
//            // pushing children values into new population
//            vector<Individual2> children;
//            children = crossoverRobot(parent1, parent2);
//            Individual2 child1;
//            Individual2 child2;
//            child1 = children[0];
//            child2 = children[1];
//            child1 = mutateRobot(child1);
//            child2 = mutateRobot(child2);
//            next_population.push_back(child1);
//            next_population.push_back(child2);
//        }
//        next_population.pop_back();
//        population = next_population;
//
//        double diversityNumber = abs(error_deviation/mean);
//
//        cout << "Gen: " << j*0.5 << "   Max: " << setprecision(4) << max_fitness << "   Diff: " << setprecision(5) << diversityNumber
//        <<endl;
//
//        // Plots
//        learningCurve << j*0.5 << " " << max_fitness  << " " << error_deviation << endl;
//        diversity << j*0.5 << " " << diversityNumber << endl;
//
//        if (diversityNumber < 0.0000000001 ) break;
//    }
//
//    cout << endl << "Finding the best individuals in the last population: " << endl;
//    double maxim = calculate_fitness_general(population[0]);
//    int index_max = 0;
//    cout << maxim << endl;
//    for (int i = 1; i < number_population; i++){
//        double current_fitness = calculate_fitness_general(population[i]);
//        cout << current_fitness << endl;
//        if (current_fitness > maxim){
//            maxim = current_fitness;
//            index_max = i;
//        }
//    }
//    cout << "Overall maximum: " << calculate_fitness_general(population[index_max]) << endl;
//
//    for (int i = 0 ; i < number_regions; i ++){
//        parametersFile  << setprecision(20) << population[index_max].coefs[i].x << " " << setprecision(20) <<   population[index_max].coefs[i].y << " " << setprecision(20) << population[index_max].coefs[i].z << " " << setprecision(20) <<   population[index_max].coefs[i].r << " " << setprecision(20) << population[index_max].coefs[i].b << " " << setprecision(20) << population[index_max].coefs[i].c << " " << setprecision(5) << population[index_max].coefs[i].k << endl;
//    }
//    parametersFile.close();


    return 0;
}

bool theSame(Mass m1, Mass m2){
    if (abs(m1.px-m2.px) < 0.000001 && abs(m1.py-m2.py)<0.000001 && abs(m1.pz-m2.pz)<0.000001 ){
        return true;
    }
    return false;
}

double normVector(double x, double y, double z){
    double res = 0;
    res = sqrt(x*x+y*y+z*z);
    return res;
}

bool isValidRobot(Individual2 ind){
    bool res = true;
    for (int i = 0; i < ind.coefs.size()-1; i++){
        for (int j = 0; j < ind.coefs.size(); j++){
            if (distanceCoefs(ind.coefs[i], ind.coefs[j]) >= (ind.coefs[i].r + ind.coefs[j].r) ) {
                res = false;
                break;
            }
        }
    }
   return res;
}

double distanceCoefs(Coef c1, Coef c2){
    double res = sqrt( (c1.x-c2.x)*(c1.x-c2.x) + (c1.y-c2.y)*(c1.y-c2.y) + (c1.z-c2.z)*(c1.z-c2.z) );
    return res;
}

double distance(Mass m1, Mass m2){
    double res = sqrt( (m1.px-m2.px)*(m1.px-m2.px) + (m1.py-m2.py)*(m1.py-m2.py) + (m1.pz-m2.pz)*(m1.pz-m2.pz) );
    return res;
}

void plotDiversity(void){
    Gnuplot gp("lines");
    gp.set_xlabel("Generations");
    gp.set_ylabel("Diversity");
    gp.set_title("Diversity chart");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/diversityRobot.txt", 1, 2, "Diversity Chart");
    gp.showonscreen();
    wait_for_key();
}
void plotDot(void){
    Gnuplot gp("points");
    gp.set_xlabel("Generations");
    gp.set_ylabel("Fitness");
    gp.set_title("Dot chart");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/dotRobot.txt", 1, 2, "Dot Chart");
    gp.showonscreen();
    wait_for_key();
  
}
void plotLearningCurve(void){
    Gnuplot gp("lines");
    gp.set_xlabel("Generations");
    gp.set_ylabel("Fitness");
    gp.set_title("Learning curve");
    gp.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/learningCurveRobot.txt", 1, 2, 3, "Learning Curve GA");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/randomCurve.txt", 1, 2, "Learning Curve Random");
    gp.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/learningCurveHC.txt", 1, 2, 3, "Learning Curve HC");
    gp.showonscreen();
    wait_for_key();
}
void plotEnergy(void){
    Gnuplot gp("lines");
    gp.set_xlabel("Time");
    gp.set_ylabel("Energy");
    gp.set_title("Potential, kinetic and total energies");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/kineticEnergy.txt", 1, 2, "Kinetic energy");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/potentialEnergy.txt", 1, 2, "Potential energy");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment3/Assignment3/totalEnergy.txt", 1, 2, "Total energy");
    gp.showonscreen();
    wait_for_key();
}
void wait_for_key (){
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;
    
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}
