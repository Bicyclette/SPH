#version 460 core

// AMD WARP = 64 threads invocations per working group
layout(local_size_x = 64, local_size_y = 1, local_size_z = 1) in;

layout (location = 0) uniform int particleCount;

// =======================================================
// ==================== PARTICLE DATA ====================
// =======================================================
layout (binding = 0, std430) buffer position
{
    float _pos[];
};

layout (binding = 4, std430) buffer neighbors
{
    int _neighbors[];
};

layout (binding = 5, std430) buffer numNeighbors
{
    int _numNeighbors[];
};

// ==========================================================
// ==================== update particles ====================
// ==========================================================

float particle_radius;
float _h;
float _sr;

void buildNeighbors(int particle_index)
{
    float x = _pos[particle_index * 4];
    float y = _pos[particle_index * 4 + 1];
    float z = _pos[particle_index * 4 + 2];
    vec3 me = vec3(x, y, z);

    int numNeighbors = 0;
    for(int i = 0; i < particleCount; ++i)
    {
        if(i == particle_index) { continue; }
        x = _pos[i * 4];
        y = _pos[i * 4 + 1];
        z = _pos[i * 4 + 2];
        vec3 particle_position = vec3(x, y, z);

        float d = length(me - particle_position);
        if(d <= _sr)
        {
            _neighbors[particle_index * 100 + numNeighbors] = i;
            numNeighbors += 1;
        }
    }
    _numNeighbors[particle_index] = numNeighbors;
}

void main()
{
    uint index = gl_GlobalInvocationID.x;
    if(index >= particleCount)
    {
        return;
    }
    
    particle_radius = 0.25f;
    _h = 2.0f * particle_radius;
    _sr = 2.0f * _h;

    int particle_index = int(index);
    buildNeighbors(particle_index);
}
