#include <tuple>

#include "constants.h"
#include "vector3D.h"
#include "random64.h"

// ----- Class Declarations -----
class Body
{
private:
    // It has a position r, velocity v and a total force F acts on it.
    Vector3D r, v, F;
    // It has a mass and a total potential energy U.
    double m, U;
    // For this first implementation I'll ignore the geometry and assume it as a point mass.
    double q;

    bool oc;
        
public:
  void init(double x, double y, double z, double vx, double vy, double vz, double mass = 1.0, double charge = 1.0, bool ocup = false);
    void printState(double t, std::ofstream &File);

    // --Inline Functions --

    // Makes the total force equal to 0 for recalculation.
    inline void resetForce(void) { F.load(0, 0, 0); };

    // Adds a force to the total force on the body
    inline void addForce(Vector3D dF) { F += dF; };

    // Makes the total potential equal to 0 for recalculation.
    inline void resetU(void) { U = 0; };

    // Adds a potential to the total  Potential energy of the body
    inline void addU(double dU) { U += dU; };

    // update of velocity
    inline void moveV(double dt, double coef = 1.0) { v += (dt * coef / m) * F; };

    // update of position
    inline void moveR(double dt, double coef = 1.0) { r += (dt * coef) * v; };

    // Get position
    inline Vector3D getR(void) { return r; };

    // Get velocity
    inline Vector3D getV(void) { return v; };

    // Get charge
    inline double getQ(void) { return q; };

    // Get ocupation
    inline bool getoc(void) { return oc; };
        
    // Set position
    inline void setR(Vector3D newR) { r = newR; };

    // Set velocity
    inline void setV(Vector3D newV) { v = newV; };

    //set charge
    inline void setq(double newq) { q = newq; };

    //set ocupation
    inline void setoc(bool newoc) { oc = newoc; };    
        
    //Get kinetic energy
    inline double getK() { return 0.5 * m * norm2(v); };

    // print gnuplot animation
    inline void print() { std::cout << "," << r.x() << "+0.05*cos(t)," << r.y() << "+0.05*sin(t)"; };

    //-- Friends --

    friend Vector3D rMin(Body b1, Body b2);
};

// ----- Function implementations -----

/*
 * Initializes the molecule in a given initial position, with a given velocity and a given mass
 * Force and potential start being 0.
 * 
 * @param double x, y, z: Initial position of the body.
 * @param double vx, vy, vz: Initial velocity of the body.
 * @param double m: Mass of the body. Defaults to 1.
 */
void Body::init(double x, double y, double z, double vx, double vy, double vz, double mass, double charge, bool ocup)
{
    r.load(x, y, z);
    v.load(vx, vy, vz);
    m = mass;
    q = charge;
    oc = ocup;
    resetForce();
    resetU();
}

/*
 * Prints the current state of the body in csv format. For the moment it uses cout
 * a future implementation is set to make it write in a given file directly.
 * 
 * @param double t: The current timestep in the simulation
 */
void Body::printState(double t, std::ofstream &File)
{
    File << t << "," << r.x() << "," << r.y() << "," << r.z() << "," << v.x() << "," << v.y() << "," << v.z() << "," << U << '\n';
}

/*
 * Calculates the minimum image distance of 2 bodies. This is based in the image
 * strategy for implementing periodic boundaries.
 * 
 * @param Body b1, b2: the 2 bodies whose minimum distance we want to calculate.
 * 
 * @return Vector3D: Returns the minimum distance between the bodies.
 */
Vector3D rMin(Body b1, Body b2)
{
    Vector3D dr;
    dr = b1.r - b2.r;
    double x, y, z;
    x = dr.x() - Lx * std::round(dr.x() / Lx);
    y = dr.y() - Ly * std::round(dr.y() / Ly);
    Vector3D rmin(x, y, z);
    return rmin;
}

