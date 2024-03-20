#pragma once
#include <algorithm>
#include <limits>
#include <omp.h>
#include <execution>
#include "mesh.hpp"

// =================================================================================
// ============================== CUBIC SPLINE KERNEL ==============================
// =================================================================================

class CubicSpline
{
public:
    explicit CubicSpline(float const h = 1);
    void setSmoothingLen(const float h);
    float smoothingLen() const;
    float supportRadius() const;
    float f(const float l) const;
    float derivative_f(const float l) const;
    float w(const glm::vec3 & rij) const;
    glm::vec3 grad_w(const glm::vec3 & rij) const; 
    glm::vec3 grad_w(const glm::vec3 & rij, const float len) const;

private:
    float _h;
    float _sr;
    float _c;
    float _gc;
};

// ========================================================================
// ============================== SPH SOLVER ==============================
// ========================================================================

class SphSolver
{
public:
    explicit SphSolver(float nu = 0.08f, float h = 0.5f, float density = 1000.0f, glm::vec3 g = glm::vec3(0.0f, -9.8f, 0.0f), float eta = 0.01f, float gamma = 7.0f);
    void set_boundary(AABB iAABB);
    void set_grid_resolution(int res_x, int res_y, int res_z);
    void initScene();
    void update();
    int particleCount() const;
    const glm::vec3 & position(const int i) const;
    int resX() const;
    int resY() const;
    int resZ() const;
    float equationOfState(const float d, const float d0, const float k, const float gamma = 7.0f);
    glm::ivec3 get_particle_cell_coordinates(glm::vec3 const& pos, glm::vec3 const & cell_stride);
    void compute_particles_indices_in_grid_cell_omp();
    void compute_particles_neighbors_omp();
    void compute_particles_indices_in_grid_cell();
    void compute_particles_neighbors();
    void set_initial_data_for_position_SSBO();
    void reset_GPU_SSBOs();
    void buildNeighbor_sonic_boom();
    void sph_sonic_boom();
    void buildNeighbor_sonic();
    void buildNeighbor_eggman();
    void computeDensityPressure_omp();
    void computeDensityPressure();
    void applyForcesAndResolveCollisions_omp();
    void applyForcesAndResolveCollisions();
    int idx1d(int x, int y, int z);
    void set_neighbors_compute_method(NEIGHBORS_COMPUTE_METHOD method);
    void reset_spatial_lookup_arrays();

public:

    CubicSpline _kernel;

    // particle data
    std::vector<glm::vec3> _init_pos;   // initial position
    std::vector<glm::vec3> _pos;        // position
    std::vector<glm::vec3> _vel;        // velocity
    std::vector<glm::vec3> _acc;        // acceleration
    std::vector<float> _p;              // pressure
    std::vector<float> _d;              // density

    std::vector<std::vector<int>> _pidxNeighbor;    // neighbor particles
    std::vector<int> _pNumNeighbors;                // num neighbors per particle

    // simulation time step
    float _dt;

    // simulation grid
    int _resX;
    int _resY;
    int _resZ;
    glm::vec3 _domain_dimensions;
    glm::vec3 _cell_stride;
    int _num_cells;
    std::vector<int> _start_index;
    std::vector<glm::ivec2> _spatial_lookup;

    // performances
    NEIGHBORS_COMPUTE_METHOD m_neighbors_method;
    bool m_use_omp;
    double _prev_time;
    double _curr_time;
    float _milliseconds;

    // wall boundaries
    float _l;
    float _r;
    float _b;
    float _t;
    float _front;
    float _back;

    // SPH coefficients
    float _nu;                     // viscosity coefficient
    float _d0;                     // rest density
    float _h;                      // particle spacing (i.e., diameter)
    glm::vec3 _g;                  // gravity
    float _m0;                     // rest mass
    float _k;                      // EOS coefficient
    float _eta;
    float _c;                      // speed of sound
    float _gamma;                  // EOS power factor

    // compute shaders
    GLbitfield _writeMask;
    GLbitfield _readMask;
    Shader _cs_build_neighbors;
    Shader _cs_sph;
    GLuint _pos_SSBO;
    GLuint _vel_SSBO;
    GLuint _density_SSBO;
    GLuint _pressure_SSBO;
    GLuint _neighbors_SSBO;
    GLuint _numNeighbors_SSBO;
};

// ===========================================================================
// ============================== FLUID DISPLAY ==============================
// ===========================================================================

class Fluid
{
public:
	Fluid() = delete;
	Fluid(std::shared_ptr<Mesh> iFluidMesh, float iParticleRadius, std::shared_ptr<Mesh> iDomainMesh);
	void draw(std::shared_ptr<Shader> iShaderDomain, std::shared_ptr<Shader> iShaderFluid, bool iWireframe);
	void update_particles_radius(float iParticleRadius);
    void reset();
    void update();

public:
	std::shared_ptr<Mesh> m_fluid_mesh;
	float m_particle_radius;
    AABB m_domain_aabb;
	GLuint m_domain_vao;
	GLuint m_domain_vbo;
    GLuint m_particle_vao;
    GLuint m_particle_vbo;
    SphSolver m_solver;
    bool m_is_running;
};

class Ray
{
public:
	Ray(glm::vec3 const& iOrigin, glm::vec3 const& iDirection)
	{
		m_origin = iOrigin;
		m_direction = glm::normalize(iDirection);
	}

	glm::vec3 getOrigin() { return m_origin; }
	glm::vec3 const& getOrigin() const { return m_origin; }
	void setOrigin(glm::vec3 const& iOrigin) { m_origin = iOrigin; }
	glm::vec3 getDirection() { return m_direction; }
	glm::vec3 const& getDirection() const { return m_direction; }
	void setDirection(glm::vec3 const& iDirection) { m_direction = glm::normalize(iDirection); }

private:
	glm::vec3 m_origin;
	glm::vec3 m_direction;
};

bool rayTriangleIntersection(Ray const& iRay, glm::vec3 const& p0, glm::vec3 const& p1, glm::vec3 const& p2, glm::vec3& out);