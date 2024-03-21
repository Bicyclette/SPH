#version 460 core

// AMD WARP = 64 threads invocations per working group
layout(local_size_x = 64, local_size_y = 1, local_size_z = 1) in;

float PI = 3.14159f;

layout (location = 0) uniform int particleCount;
layout (location = 1) uniform float particleRadius;

float square(float a) { return a * a; }

float cube(float a) { return a * a * a; }

// =======================================================
// ==================== PARTICLE DATA ====================
// =======================================================
layout (binding = 0, std430) buffer position
{
    float _pos[];
};

layout (binding = 1, std430) buffer velocity
{
    float _vel[];
};

layout (binding = 2, std430) buffer density
{
    float _density[];
};

layout (binding = 3, std430) buffer pressure
{
    float _pressure[];
};

layout (binding = 4, std430) buffer neighbors
{
    int _neighbors[];
};

layout (binding = 5, std430) buffer numNeighbors
{
    int _numNeighbors[];
};

// ==============================================
// ==================== WALL ====================
// ==============================================

uniform float _l;
uniform float _r;
uniform float _b;
uniform float _t;
uniform float _front;
uniform float _back;

// =================================================================
// ==================== SPH cubic spline kernel ====================
// =================================================================

// SPH Kernel function: cubic spline
struct CubicSpline
{
    float _h;
    float _sr;
    float _c;
    float _gc;
};

CubicSpline _kernel;

void setSmoothingLen(float h)
{
    float h2 = square(h);
    float h3 = h2 * h;
    _kernel._h = h;
    _kernel._sr = 2.0f * h;
    _kernel._c = 1.0f / (PI * h3);
    _kernel._gc = _kernel._c / h;
}

float f(float l)
{
    float q = l / _kernel._h;
    if (q < 1.0f)
    {
        return _kernel._c * (1.0f - 1.5f * square(q) + 0.75f * cube(q));
    }
    else if (q < 2.0f)
    {
        return _kernel._c * (0.25f * cube(2.0f - q));
    }
    return 0.0f;
}

float derivative_f(float l)
{
    float q = l / _kernel._h;
    if (q <= 1.0f)
    {
        return _kernel._gc * (-3.0f * q + 2.25f * square(q));
    }
    else if (q < 2.0f)
    {
        return -_kernel._gc * 0.75f * square(2.0f - q);
    }
    return 0.0f;
}

float w(vec3 rij) { return f(length(rij)); }

vec3 grad_w(vec3 rij, float len)
{
    return derivative_f(len) * rij / len;
}

vec3 grad_w(vec3 rij) { return grad_w(rij, length(rij)); }

// ====================================================
// ==================== SPH solver ====================
// ====================================================

struct SphSolver
{
    float _nu;
    float _h;
    float _d0;
    vec3 _g;
    float _eta;
    float _gamma;
    float _dt;
    float _m0;
    float _c;
    float _k;
};

SphSolver _solver;

void init_solver()
{
    _solver._nu = 0.0075f;
    _solver._h = particleRadius * 2.0f;
    _solver._d0 = 1000.0f;
    _solver._g = vec3(0.0f, -9.81f, 0.0f);
    _solver._eta = 0.01f;
    _solver._gamma = 7.0f;
    _solver._dt = 0.00025f;
    _solver._m0 = _solver._d0 * cube(_solver._h);
    _solver._c = 9.81f / _solver._eta;
    _solver._k = _solver._d0 * square(_solver._c) / _solver._gamma;
}

// ========================================================
// ==================== SPH simulation ====================
// ========================================================

 float equationOfState(float d, float d0, float k)
{
    return k * (pow(d / d0, _solver._gamma) - 1.0f);
}

void computeDensity(int particle_index)
{
    _density[particle_index] = _solver._m0 * w(vec3(0.0f, 0.0f, 0.0f));
    vec3 pi = vec3(_pos[particle_index * 4], _pos[particle_index * 4 + 1], _pos[particle_index * 4 + 2]);
    
    // take all neighbors
    int neighbor_count = _numNeighbors[particle_index];
    
    // compute smoothing weight
    for (int j = 0; j < neighbor_count; ++j)
    {
        int neighbor_index = _neighbors[particle_index * 100 + j];
        vec3 pj = vec3(_pos[neighbor_index * 4], _pos[neighbor_index * 4 + 1], _pos[neighbor_index * 4 + 2]);
        vec3 rij = pi - pj;
        float weight = w(rij);
        _density[particle_index] += _solver._m0 * weight;
    }
}

void computePressure(int particle_index)
{
    float eos = equationOfState(_density[particle_index], _solver._d0, _solver._k);
    _pressure[particle_index] = (eos > 0.0f) ? eos : 0.0f;
}

vec3 applyForce(int particle_index)
{
    vec3 acc = vec3(0.0f, 0.0f, 0.0f);

    // gravity
    acc += _solver._g;

    vec3 pi = vec3(_pos[particle_index * 4], _pos[particle_index * 4 + 1], _pos[particle_index * 4 + 2]);
    int neighbor_count = _numNeighbors[particle_index];
    
    // pressure force and viscous force
    vec3 f_pressure = vec3(0.0f, 0.0f, 0.0f);
    vec3 f_viscous = vec3(0.0f, 0.0f, 0.0f);
    for (int j = 0; j < neighbor_count; ++j)
    {
        int neighbor_index = _neighbors[particle_index * 100 + j];
        vec3 pj = vec3(_pos[neighbor_index * 4], _pos[neighbor_index * 4 + 1], _pos[neighbor_index * 4 + 2]);
        vec3 rij = pi - pj;
        vec3 gw = grad_w(rij);
        vec3 uij = vec3(_vel[particle_index * 4], _vel[particle_index * 4 + 1], _vel[particle_index * 4 + 2]) - vec3(_vel[neighbor_index * 4], _vel[neighbor_index * 4 + 1], _vel[neighbor_index * 4 + 2]);

        float pi_o_di2 = _pressure[particle_index] / pow(_density[particle_index], 2.0f);
        float pj_o_dj2 = _pressure[neighbor_index] / pow(_density[neighbor_index], 2.0f);
        f_pressure -= _solver._m0 * (pi_o_di2 + pj_o_dj2) * gw;

        f_viscous += 2.0f * _solver._nu * _solver._m0 * (_solver._m0 / _density[neighbor_index]) * uij * (dot(rij, gw) / (dot(rij, rij) + 0.01f * pow(_solver._h, 2.0f)));
    }
    acc += f_pressure + f_viscous;

    return acc;
}

void updateVelocity(int particle_index, vec3 acc)
{
    _vel[particle_index * 4] = _vel[particle_index * 4] + _solver._dt * acc.x;
    _vel[particle_index * 4 + 1] = _vel[particle_index * 4 + 1] + _solver._dt * acc.y;
    _vel[particle_index * 4 + 2] = _vel[particle_index * 4 + 2] + _solver._dt * acc.z;
}

void updatePosition(int particle_index)
{
    _pos[particle_index * 4] = _pos[particle_index * 4] + _solver._dt * _vel[particle_index * 4];
    _pos[particle_index * 4 + 1] = _pos[particle_index * 4 + 1] + _solver._dt * _vel[particle_index * 4 + 1];
    _pos[particle_index * 4 + 2] = _pos[particle_index * 4 + 2] + _solver._dt * _vel[particle_index * 4 + 2];
}

void resolveCollision(int particle_index)
{
    float bias = _solver._h * 0.5f;
    float left = _l + bias;
    float right = _r - bias;
    float bottom = _b + bias;
    float top = _t - bias;
    float back = _back + bias;
    float front = _front - bias;

    vec3 particle_position = vec3(_pos[particle_index * 4], _pos[particle_index * 4 + 1], _pos[particle_index * 4 + 2]);
    
    if (particle_position.x < left || particle_position.x > right ||
        particle_position.y < bottom || particle_position.y > top ||
        particle_position.z < back || particle_position.z > front)
    {
        vec3 p = vec3(_pos[particle_index * 4], _pos[particle_index * 4 + 1], _pos[particle_index * 4 + 2]);
        _pos[particle_index * 4] = clamp(particle_position.x, left, right);
        _pos[particle_index * 4 + 1] = clamp(particle_position.y, bottom, top);
        _pos[particle_index * 4 + 2] = clamp(particle_position.z, back, front);
        _pos[particle_index * 4 + 3] = 0.0f;

        _vel[particle_index * 4] = _pos[particle_index * 4] - p.x;
        _vel[particle_index * 4 + 1] = _pos[particle_index * 4 + 1] - p.y;
        _vel[particle_index * 4 + 2] = _pos[particle_index * 4 + 2] - p.z;
        _vel[particle_index * 4 + 3] = 0.0f;
    }
}

void main()
{
    uint index = gl_GlobalInvocationID.x;
    if(index >= particleCount)
    {
        return;
    }
    int particle_index = int(index);

    setSmoothingLen(2.0f * particleRadius);
    init_solver();
    computeDensity(particle_index);
    computePressure(particle_index);
    vec3 acc = applyForce(particle_index);
    updateVelocity(particle_index, acc);
    updatePosition(particle_index);
    resolveCollision(particle_index);
}
