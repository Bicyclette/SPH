#include "sph_fluid.hpp"

// =================================================================================
// ============================== CUBIC SPLINE KERNEL ==============================
// =================================================================================

CubicSpline::CubicSpline(float const h)
{
	setSmoothingLen(h);
}

void CubicSpline::setSmoothingLen(const float h)
{
	float const h2 = square(h);
	float const h3 = h2 * h;
	_h = h;
	_sr = 2.0f * h;
	_c = 1.0f / (M_PI * h3);
	_gc = _c / h;
}

float CubicSpline::smoothingLen() const { return _h; }

float CubicSpline::supportRadius() const { return _sr; }

float CubicSpline::f(const float l) const
{
	float const q = l / _h;
	if (q < 1.0f)
	{
		return _c * (1.0f - 1.5f * square(q) + 0.75f * cube(q));
	}
	else if (q < 2.0f)
	{
		return _c * (0.25f * cube(2.0f - q));
	}
	return 0.0f;
}

float CubicSpline::derivative_f(const float l) const
{
	float const q = l / _h;
	if (q <= 1.0f) return _gc * (-3.0f * q + 2.25f * square(q));
	else if (q < 2.0f) return -_gc * 0.75f * square(2.0f - q);
	return 0.0f;
}

float CubicSpline::w(const glm::vec3 & rij) const { return f(glm::length(rij)); }

glm::vec3 CubicSpline::grad_w(const glm::vec3 & rij) const { return grad_w(rij, glm::length(rij)); }

glm::vec3 CubicSpline::grad_w(const glm::vec3 & rij, const float len) const
{
	return derivative_f(len) * rij / len;
}

// ========================================================================
// ============================== SPH SOLVER ==============================
// ========================================================================

SphSolver::SphSolver(float nu, float h, float density, glm::vec3 g, float eta, float gamma) :
_kernel(h),
_nu(nu),
_h(h),
_d0(density),
_g(g),
_eta(eta),
_gamma(gamma),
_cs_build_neighbors("../shaders/compute/NV/build_neighbors.glsl"),
_cs_sph("../shaders/compute/NV/sph.glsl")
{
	_dt = 0.00025f;
	_m0 = _d0 * _h * _h * _h;
	_c = std::fabs(_g.y) / _eta;
	_k = _d0 * _c * _c / _gamma;

	m_neighbors_method = SONIC;
	m_use_omp = true;
	
	// create Shader Storage Buffer Objects
	_writeMask = GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT;
	_readMask = GL_MAP_READ_BIT;

	glGenBuffers(1, &_pos_SSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _pos_SSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, c_max_particle_count * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, _pos_SSBO);

	glGenBuffers(1, &_vel_SSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _vel_SSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, c_max_particle_count * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, _vel_SSBO);

	glGenBuffers(1, &_density_SSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _density_SSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, c_max_particle_count * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, _density_SSBO);

	glGenBuffers(1, &_pressure_SSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _pressure_SSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, c_max_particle_count * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, _pressure_SSBO);

	glGenBuffers(1, &_neighbors_SSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _neighbors_SSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, c_max_particle_count * c_max_particle_neighbors * sizeof(int), nullptr, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, _neighbors_SSBO);

	glGenBuffers(1, &_numNeighbors_SSBO);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _numNeighbors_SSBO);
	glBufferData(GL_SHADER_STORAGE_BUFFER, c_max_particle_count * sizeof(int), nullptr, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, _numNeighbors_SSBO);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void SphSolver::set_boundary(AABB iAABB)
{
	_l = iAABB.min_x;
	_r = iAABB.max_x;
	_b = iAABB.min_y;
	_t = iAABB.max_y;
	_back = iAABB.min_z;
	_front = iAABB.max_z;

	const float domain_width = fabs(_r - _l);
	const float domain_height = fabs(_t - _b);
	const float domain_depth = fabs(_front - _back);
	_domain_dimensions = glm::vec3(domain_width, domain_height, domain_depth);
	
	_cs_build_neighbors.use();
	_cs_build_neighbors.set("_l", _l);
	_cs_build_neighbors.set("_r", _r);
	_cs_build_neighbors.set("_b", _b);
	_cs_build_neighbors.set("_t", _t);
	_cs_build_neighbors.set("_front", _front);
	_cs_build_neighbors.set("_back", _back);

	_cs_sph.use();
	_cs_sph.set("_l", _l);
	_cs_sph.set("_r", _r);
	_cs_sph.set("_b", _b);
	_cs_sph.set("_t", _t);
	_cs_sph.set("_front", _front);
	_cs_sph.set("_back", _back);
}

void SphSolver::set_grid_resolution(int res_x, int res_y, int res_z)
{
	_resX = res_x;
	_resY = res_y;
	_resZ = res_z;

	const float stride_width = _domain_dimensions.x / static_cast<float>(resX());
	const float stride_height = _domain_dimensions.y / static_cast<float>(resY());
	const float stride_depth = _domain_dimensions.z / static_cast<float>(resZ());
	
	_cell_stride = glm::vec3(stride_width, stride_height, stride_depth);
	_num_cells = resX() * resY() * resZ();
	_start_index.clear();
	_start_index.resize(_num_cells, -1);
}

void SphSolver::initScene()
{
	_init_pos.clear();
	_pos.clear();

	// make sure for the other particle quantities
	_vel = std::vector<glm::vec3>(_pos.size(), glm::vec3(0.0f));
	_acc = std::vector<glm::vec3>(_pos.size(), glm::vec3(0.0f));
	_p = std::vector<float>(_pos.size(), 0);
	_d = std::vector<float>(_pos.size(), 0);
}

void SphSolver::update()
{
	_prev_time = glfwGetTime();

	if (m_neighbors_method == EGGMAN)
	{
		buildNeighbor_eggman();
	}
	else if (m_neighbors_method == SONIC)
	{
		buildNeighbor_sonic();
	}
	else if (m_neighbors_method == SONIC_BOOM)
	{
		buildNeighbor_sonic_boom();
		sph_sonic_boom();

		_curr_time = glfwGetTime();
		_milliseconds = static_cast<float>(_curr_time - _prev_time);

		return;
	}
	else
	{
		std::cerr << "Error : invalid method for neighbors computation." << std::endl;
		std::exit(-1);
	}

	if(m_use_omp)
	{
		computeDensityPressure_omp();
		_acc = std::vector<glm::vec3>(_pos.size(), glm::vec3(0.0f));
		applyForcesAndResolveCollisions_omp();
	}
	else
	{
		computeDensityPressure();
		_acc = std::vector<glm::vec3>(_pos.size(), glm::vec3(0.0f));
		applyForcesAndResolveCollisions();
	}

	_curr_time = glfwGetTime();
	_milliseconds = static_cast<float>(_curr_time - _prev_time);
}

int SphSolver::particleCount() const { return _pos.size(); }

const glm::vec3 & SphSolver::position(const int i) const { return _pos[i]; }

int SphSolver::resX() const { return _resX; }

int SphSolver::resY() const { return _resY; }

int SphSolver::resZ() const { return _resZ; }

float SphSolver::equationOfState(const float d, const float d0, const float k, const float gamma)
{
	return k * (pow(d / d0, gamma) - 1.0f);
}

glm::ivec3 SphSolver::get_particle_cell_coordinates(glm::vec3 const& pos, glm::vec3 const & cell_stride)
{
	glm::vec3 w_origin(_l, _b, _back);
	glm::vec3 w_new_pos = pos - w_origin;

	int x_water_tank_coord = static_cast<int>(floor(w_new_pos.x / cell_stride.x));
	int y_water_tank_coord = static_cast<int>(floor(w_new_pos.y / cell_stride.y));
	int z_water_tank_coord = static_cast<int>(floor(w_new_pos.z / cell_stride.z));

	glm::ivec3 coords(x_water_tank_coord, y_water_tank_coord, z_water_tank_coord);
	return coords;
}

void SphSolver::compute_particles_indices_in_grid_cell_omp()
{
	#pragma omp parallel for
	for (int p = 0; p < particleCount(); ++p)
	{
		glm::ivec3 cell_coordinates = get_particle_cell_coordinates(_pos[p], _cell_stride);

		// get index of cell
		int cell_index = idx1d(cell_coordinates.x, cell_coordinates.y, cell_coordinates.z);

		// fill spatial lookup array
		_spatial_lookup[p].x = p;
		_spatial_lookup[p].y = cell_index;
	}
}

void SphSolver::compute_particles_neighbors_omp()
{
	float sr = _kernel.supportRadius();
	int offset_x = static_cast<int>(ceil(sr / _cell_stride.x));
	int offset_y = static_cast<int>(ceil(sr / _cell_stride.y));
	int offset_z = static_cast<int>(ceil(sr / _cell_stride.z));

	#pragma omp parallel for
	for (int i = 0; i < particleCount(); ++i)
	{
		_pNumNeighbors[i] = 0;

		glm::vec3 pos = position(i);
		glm::ivec3 cell_coordinates = get_particle_cell_coordinates(pos, _cell_stride);
		for (int x = -offset_x; x <= offset_x; ++x)
		{
			for (int y = -offset_y; y <= offset_y; ++y)
			{
				for (int z = -offset_z; z <= offset_z; ++z)
				{
					glm::ivec3 cell = cell_coordinates + glm::ivec3(x, y, z);
					// check if we are inside the simulation domain
					if (cell.x < 0 || cell.y < 0 || cell.z < 0 || cell.x >= resX() || cell.y >= resY() || cell.z >= resZ())
					{
						continue;
					}

					// get cell index
					int cell_index = idx1d(cell.x, cell.y, cell.z);

					// get index from which to get data in spatial lookup array
					int start_index_lookup = _start_index[cell_index];
					if (start_index_lookup == -1) { continue; }

					// distance check
					while (start_index_lookup < particleCount() && _spatial_lookup[start_index_lookup].y == cell_index)
					{
						int particle_index = _spatial_lookup[start_index_lookup].x;
						if (particle_index == i)
						{
							start_index_lookup++;
							continue;
						}
						float d = glm::length(pos - _pos[particle_index]);
						if (d <= sr)
						{
							_pidxNeighbor[i][_pNumNeighbors[i]] = particle_index;
							_pNumNeighbors[i] += 1;
						}
						start_index_lookup++;
					}
				}
			}

		}
	}
}

void SphSolver::compute_particles_indices_in_grid_cell()
{
	for (int p = 0; p < particleCount(); ++p)
	{
		glm::ivec3 cell_coordinates = get_particle_cell_coordinates(_pos[p], _cell_stride);

		// get index of cell
		int cell_index = idx1d(cell_coordinates.x, cell_coordinates.y, cell_coordinates.z);

		// fill spatial lookup array
		_spatial_lookup[p].x = p;
		_spatial_lookup[p].y = cell_index;
	}
}

void SphSolver::compute_particles_neighbors()
{
	float sr = _kernel.supportRadius();
	int offset_x = static_cast<int>(ceil(sr / _cell_stride.x));
	int offset_y = static_cast<int>(ceil(sr / _cell_stride.y));
	int offset_z = static_cast<int>(ceil(sr / _cell_stride.z));

	for (int i = 0; i < particleCount(); ++i)
	{
		_pNumNeighbors[i] = 0;

		glm::vec3 pos = position(i);
		glm::ivec3 cell_coordinates = get_particle_cell_coordinates(pos, _cell_stride);
		for (int x = -offset_x; x <= offset_x; ++x)
		{
			for (int y = -offset_y; y <= offset_y; ++y)
			{
				for (int z = -offset_z; z <= offset_z; ++z)
				{
					glm::ivec3 cell = cell_coordinates + glm::ivec3(x, y, z);
					// check if we are inside the simulation domain
					if (cell.x < 0 || cell.y < 0 || cell.z < 0 || cell.x >= resX() || cell.y >= resY() || cell.z >= resZ())
					{
						continue;
					}

					// get cell index
					int cell_index = idx1d(cell.x, cell.y, cell.z);

					// get index from which to get data in spatial lookup array
					int start_index_lookup = _start_index[cell_index];
					if (start_index_lookup == -1) { continue; }

					// distance check
					while (start_index_lookup < particleCount() && _spatial_lookup[start_index_lookup].y == cell_index)
					{
						int particle_index = _spatial_lookup[start_index_lookup].x;
						if (particle_index == i)
						{
							start_index_lookup++;
							continue;
						}
						float d = glm::length(pos - _pos[particle_index]);
						if (d <= sr)
						{
							_pidxNeighbor[i][_pNumNeighbors[i]] = particle_index;
							_pNumNeighbors[i] += 1;
						}
						start_index_lookup++;
					}
				}
			}

		}
	}
}

void SphSolver::set_initial_data_for_position_SSBO()
{
	// send particle positions to position SSBO
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _pos_SSBO);
	float* positions_ptr = (float*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, c_max_particle_count * 4 * sizeof(float), _writeMask);
	for (size_t i = 0; i < particleCount(); ++i)
	{
		positions_ptr[i * 4] = _pos[i].x;
		positions_ptr[i * 4 + 1] = _pos[i].y;
		positions_ptr[i * 4 + 2] = _pos[i].z;
	}
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
}

void SphSolver::reset_GPU_SSBOs()
{
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _pos_SSBO);
	float* positions_ptr = (float*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, c_max_particle_count * 4 * sizeof(float), _writeMask);
	for (size_t i = 0; i < particleCount(); ++i)
	{
		positions_ptr[i * 4] = _init_pos[i].x;
		positions_ptr[i * 4 + 1] = _init_pos[i].y;
		positions_ptr[i * 4 + 2] = _init_pos[i].z;
	}
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, _vel_SSBO);
	float* velocity_ptr = (float*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, c_max_particle_count * 4 * sizeof(float), _writeMask);
	for (size_t i = 0; i < particleCount(); ++i)
	{
		velocity_ptr[i * 4] = 0.0f;
		velocity_ptr[i * 4 + 1] = 0.0f;
		velocity_ptr[i * 4 + 2] = 0.0f;
		velocity_ptr[i * 4 + 3] = 0.0f;
	}
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
}

void SphSolver::buildNeighbor_sonic_boom()
{
	_cs_build_neighbors.use();
	_cs_build_neighbors.set("particleCount", particleCount());

	glDispatchCompute((particleCount() / 32) + 32, 1, 1);
	glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void SphSolver::sph_sonic_boom()
{
	_cs_sph.use();
	_cs_sph.set("particleCount", particleCount());
	_cs_sph.set("particleRadius", _h * 0.5f);

	glDispatchCompute((particleCount() / 32) + 32, 1, 1);
	glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void SphSolver::buildNeighbor_sonic()
{
	reset_spatial_lookup_arrays();

	// first assign for each grid cell the index of particles that belongs into it
	if (m_use_omp)
	{
		compute_particles_indices_in_grid_cell_omp();
	}
	else
	{
		compute_particles_indices_in_grid_cell();
	}

	std::sort(std::execution::par_unseq, _spatial_lookup.begin(), _spatial_lookup.end(), [] (glm::ivec2 const & a, glm::ivec2 const & b) -> bool { return a.y < b.y; });

	// start index array
	int current_value = -1;
	for (size_t particle_index = 0; particle_index < particleCount(); ++particle_index)
	{
		int cell_index = _spatial_lookup[particle_index].y;
		if (cell_index != current_value)
		{
			_start_index[cell_index] = particle_index;
			current_value = cell_index;
		}
	}

	if (m_use_omp)
	{
		compute_particles_neighbors_omp();
	}
	else
	{
		compute_particles_neighbors();
	}
}

void SphSolver::buildNeighbor_eggman()
{
	float sr = _kernel.supportRadius();
	if(m_use_omp)
	{
		#pragma omp parallel for
		for (size_t i = 0; i < particleCount(); ++i)
		{
			_pNumNeighbors[i] = 0;
			glm::vec3 pos_i = position(i);
			for (size_t j = 0; j < particleCount(); ++j)
			{
				if (i == j) { continue; }
				glm::vec3 pos_j = position(j);
				float d = glm::length(pos_i - pos_j);
				if (d <= sr)
				{
					_pidxNeighbor[i][_pNumNeighbors[i]] = j;
					_pNumNeighbors[i] = 0;
				}
			}
		}
	}
	else
	{
		for (size_t i = 0; i < particleCount(); ++i)
		{
			_pNumNeighbors[i] = 0;
			glm::vec3 pos_i = position(i);
			for (size_t j = 0; j < particleCount(); ++j)
			{
				if (i == j) { continue; }
				glm::vec3 pos_j = position(j);
				float d = glm::length(pos_i - pos_j);
				if (d <= sr)
				{
					_pidxNeighbor[i][_pNumNeighbors[i]] = j;
					_pNumNeighbors[i] += 1;
				}
			}
		}
	}
}

void SphSolver::computeDensityPressure_omp()
{
	#pragma omp parallel for
	for (int i = 0; i < particleCount(); ++i)
	{
		_d[i] = _m0 * _kernel.w(glm::vec3(0));

		// compute smoothing weight
		for (int j = 0; j < _pNumNeighbors[i]; ++j)
		{
			int index_neighbor = _pidxNeighbor[i][j];
			glm::vec3 rij = position(i) - position(index_neighbor);
			float weight = _kernel.w(rij);
			_d[i] += _m0 * weight;
		}

		float eos = equationOfState(_d[i], _d0, _k);
		_p[i] = (eos > 0.0) ? eos : 0.0;
	}
}

void SphSolver::computeDensityPressure()
{
	for (int i = 0; i < particleCount(); ++i)
	{
		_d[i] = _m0 * _kernel.w(glm::vec3(0));

		// compute smoothing weight
		for (int j = 0; j < _pNumNeighbors[i]; ++j)
		{
			int index_neighbor = _pidxNeighbor[i][j];
			glm::vec3 rij = position(i) - position(index_neighbor);
			float weight = _kernel.w(rij);
			_d[i] += _m0 * weight;
		}

		float eos = equationOfState(_d[i], _d0, _k);
		_p[i] = (eos > 0.0) ? eos : 0.0;
	}
}

void SphSolver::applyForcesAndResolveCollisions_omp()
{
	const float bias = _h * 0.5f;
	float left = _l + bias;
	float right = _r - bias;
	float bottom = _b + bias;
	float top = _t - bias;
	float back = _back + bias;
	float front = _front - bias;

	#pragma omp parallel for
	for (int i = 0; i < particleCount(); ++i)
	{
		_acc[i] += _g;

		// pressure and viscous forces
		glm::vec3 f_pressure(0.0f, 0.0f, 0.0f);
		glm::vec3 f_viscosity(0.0f);
		for (int j = 0; j < _pNumNeighbors[i]; ++j)
		{
			int index_neighbor = _pidxNeighbor[i][j];
			glm::vec3 rij = position(i) - position(index_neighbor);
			glm::vec3 gw = _kernel.grad_w(rij);
			float pi_o_di2 = _p[i] / pow(_d[i], 2.0);
			float pj_o_dj2 = _p[index_neighbor] / pow(_d[index_neighbor], 2.0);
			f_pressure -= _m0 * (pi_o_di2 + pj_o_dj2) * gw;
			glm::vec3 uij = _vel[i] - _vel[index_neighbor];
			f_viscosity += 2.0f * _nu * _m0 * (_m0 / _d[index_neighbor]) * uij * (glm::dot(rij, gw) / (glm::dot(rij, rij) + 0.01f * pow(_h, 2.0f)));
		}
		_acc[i] += f_pressure + f_viscosity;

		_vel[i] = _vel[i] + _dt * _acc[i];
		_pos[i] = _pos[i] + _dt * _vel[i];
	
		if (_pos[i].x < left || _pos[i].x > right ||
			_pos[i].y < bottom || _pos[i].y > top ||
			_pos[i].z < back || _pos[i].z > front)
		{
			const glm::vec3 p = _pos[i];
			_pos[i].x = clamp(_pos[i].x, left, right);
			_pos[i].y = clamp(_pos[i].y, bottom, top);
			_pos[i].z = clamp(_pos[i].z, back, front);
			_vel[i] *= -0.7f;
		}
	}
}

void SphSolver::applyForcesAndResolveCollisions()
{
	const float bias = _h * 0.5f;
	float left = _l + bias;
	float right = _r - bias;
	float bottom = _b + bias;
	float top = _t - bias;
	float back = _back + bias;
	float front = _front - bias;

	for (int i = 0; i < particleCount(); ++i)
	{
		_acc[i] += _g;

		// pressure and viscous forces
		glm::vec3 f_pressure(0.0f, 0.0f, 0.0f);
		glm::vec3 f_viscosity(0.0f);
		for (int j = 0; j < _pNumNeighbors[i]; ++j)
		{
			int index_neighbor = _pidxNeighbor[i][j];
			glm::vec3 rij = position(i) - position(index_neighbor);
			glm::vec3 gw = _kernel.grad_w(rij);
			float pi_o_di2 = _p[i] / pow(_d[i], 2.0);
			float pj_o_dj2 = _p[index_neighbor] / pow(_d[index_neighbor], 2.0);
			f_pressure -= _m0 * (pi_o_di2 + pj_o_dj2) * gw;
			glm::vec3 uij = _vel[i] - _vel[index_neighbor];
			f_viscosity += 2.0f * _nu * _m0 * (_m0 / _d[index_neighbor]) * uij * (glm::dot(rij, gw) / (glm::dot(rij, rij) + 0.01f * pow(_h, 2.0f)));
		}
		_acc[i] += f_pressure + f_viscosity;

		_vel[i] = _vel[i] + _dt * _acc[i];
		_pos[i] = _pos[i] + _dt * _vel[i];
	
		if (_pos[i].x < left || _pos[i].x > right ||
			_pos[i].y < bottom || _pos[i].y > top ||
			_pos[i].z < back || _pos[i].z > front)
		{
			const glm::vec3 p = _pos[i];
			_pos[i].x = clamp(_pos[i].x, left, right);
			_pos[i].y = clamp(_pos[i].y, bottom, top);
			_pos[i].z = clamp(_pos[i].z, back, front);
			_vel[i] *= -0.7f;
		}
	}
}

int SphSolver::idx1d(int x, int y, int z)
{
	int slice = resX() * resY();
	return z * slice + x + y * resX();
}

void SphSolver::set_neighbors_compute_method(NEIGHBORS_COMPUTE_METHOD method)
{
	m_neighbors_method = method;
}

void SphSolver::reset_spatial_lookup_arrays()
{
	_spatial_lookup.assign(particleCount(), glm::ivec2(-1, -1));
	_start_index.assign(_num_cells, -1);
}

// ===========================================================================
// ============================== FLUID DISPLAY ==============================
// ===========================================================================

Fluid::Fluid(std::shared_ptr<Mesh> iFluidMesh, float iParticleRadius, std::shared_ptr<Mesh> iDomainMesh) :
m_solver(0.001f, iParticleRadius * 2.0f, 1000.0f, glm::vec3(0.0f, -9.8f, 0.0f), 0.01f, 7.0f)
{
	// fluid mesh (volume)
	m_fluid_mesh = iFluidMesh;

	m_particle_radius = iParticleRadius;

	// domain (container)
	m_domain_aabb = iDomainMesh->compute_axis_aligned_bounding_box();
	m_solver.set_boundary(m_domain_aabb);
	m_solver.set_grid_resolution(m_solver._domain_dimensions.x, m_solver._domain_dimensions.y, m_solver._domain_dimensions.z); // water tank dimensions

	float domain[72] =
	{
		// along x
		m_domain_aabb.min_x, m_domain_aabb.min_y, m_domain_aabb.min_z,
		m_domain_aabb.max_x, m_domain_aabb.min_y, m_domain_aabb.min_z,
		m_domain_aabb.min_x, m_domain_aabb.max_y, m_domain_aabb.min_z,
		m_domain_aabb.max_x, m_domain_aabb.max_y, m_domain_aabb.min_z,
		m_domain_aabb.min_x, m_domain_aabb.min_y, m_domain_aabb.max_z,
		m_domain_aabb.max_x, m_domain_aabb.min_y, m_domain_aabb.max_z,
		m_domain_aabb.min_x, m_domain_aabb.max_y, m_domain_aabb.max_z,
		m_domain_aabb.max_x, m_domain_aabb.max_y, m_domain_aabb.max_z,
		// along y
		m_domain_aabb.min_x, m_domain_aabb.min_y, m_domain_aabb.min_z,
		m_domain_aabb.min_x, m_domain_aabb.max_y, m_domain_aabb.min_z,
		m_domain_aabb.max_x, m_domain_aabb.min_y, m_domain_aabb.min_z,
		m_domain_aabb.max_x, m_domain_aabb.max_y, m_domain_aabb.min_z,
		m_domain_aabb.min_x, m_domain_aabb.min_y, m_domain_aabb.max_z,
		m_domain_aabb.min_x, m_domain_aabb.max_y, m_domain_aabb.max_z,
		m_domain_aabb.max_x, m_domain_aabb.min_y, m_domain_aabb.max_z,
		m_domain_aabb.max_x, m_domain_aabb.max_y, m_domain_aabb.max_z,
		// along z
		m_domain_aabb.min_x, m_domain_aabb.min_y, m_domain_aabb.min_z,
		m_domain_aabb.min_x, m_domain_aabb.min_y, m_domain_aabb.max_z,
		m_domain_aabb.min_x, m_domain_aabb.max_y, m_domain_aabb.min_z,
		m_domain_aabb.min_x, m_domain_aabb.max_y, m_domain_aabb.max_z,
		m_domain_aabb.max_x, m_domain_aabb.min_y, m_domain_aabb.min_z,
		m_domain_aabb.max_x, m_domain_aabb.min_y, m_domain_aabb.max_z,
		m_domain_aabb.max_x, m_domain_aabb.max_y, m_domain_aabb.min_z,
		m_domain_aabb.max_x, m_domain_aabb.max_y, m_domain_aabb.max_z
	};

	// DOMAIN's VAO
	glGenVertexArrays(1, &m_domain_vao);
	glBindVertexArray(m_domain_vao);

	// DOMAIN's VBO
	glGenBuffers(1, &m_domain_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, m_domain_vbo);
	glBufferData(GL_ARRAY_BUFFER, 72 * sizeof(float), domain, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

	glBindVertexArray(0);

	// load particle mesh
	m_particle_mesh = std::make_shared<Mesh>("../assets/sphere.obj");
	m_particle_mesh->m_material.m_albedo = glm::vec3(0.05f, 0.05f, 0.85f);

	// create particles
	update_particles_radius(iParticleRadius);
	m_solver.set_initial_data_for_position_SSBO();
	m_is_running = false;
}

void Fluid::draw(std::shared_ptr<Shader> iShader, bool iWireframe)
{	
	// draw domain
	glLineWidth(2.0f);
	glBindVertexArray(m_domain_vao);
	iShader->use();
	iShader->set("instanced_rendering", false);
	iShader->set("draw_domain", true);
	glDrawArrays(GL_LINES, 0, 72);
	glLineWidth(1.0f);

	// draw fluid
	iShader->set("draw_domain", false);
	if (iWireframe)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	if (m_solver.m_neighbors_method != SONIC_BOOM)
	{
		m_particle_mesh->set_instance_rendering(m_solver._pos, m_solver._h);
		m_particle_mesh->draw(iShader);
	}
	else
	{
		iShader->set("particleScale", m_solver._h);
		m_particle_mesh->m_instance_rendering = true;
		m_particle_mesh->m_instance_count = m_solver.particleCount();
		m_fluid_mesh->link_to_pos_SSBO(m_solver._pos_SSBO);
		m_particle_mesh->draw(iShader, true);
	}
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Fluid::update_particles_radius(float iParticleRadius)
{
	m_particle_radius = iParticleRadius;
	m_solver.initScene();
	m_solver._kernel.setSmoothingLen(iParticleRadius * 2.0f);

	// get mesh bounding box
	AABB aabb = m_fluid_mesh->compute_axis_aligned_bounding_box();

	// mesh data
	std::vector<uint32_t> triangles = m_fluid_mesh->m_cpu_geometry.m_triangleIndices;
	std::vector<Vertex> vertices = m_fluid_mesh->m_cpu_geometry.m_vertices;

	// discretize mesh by looking along negative x axis (front axis)
	float bias = 0.5f;
	float origin_x = aabb.max_x + bias;
	float step = (iParticleRadius * 2.0f);
	for (float z = aabb.min_z; z < aabb.max_z; z += step)
	{
		for (float y = aabb.min_y; y < aabb.max_y; y += step)
		{
			Ray ray(glm::vec3(origin_x, y, z), glm::vec3(-1.0f, 0.0f, 0.0f));

			// check against all triangles
			std::vector<glm::vec3> hit;
			for (size_t i = 0; i < triangles.size(); i += 3)
			{
				glm::vec3 v0 = vertices[triangles[i]].m_position;
				glm::vec3 v1 = vertices[triangles[i + 1]].m_position;
				glm::vec3 v2 = vertices[triangles[i + 2]].m_position;
				glm::vec3 intersection;
				if (rayTriangleIntersection(ray, v0, v1, v2, intersection))
				{
					if (std::find(hit.begin(), hit.end(), intersection) == hit.end())
					{
						hit.push_back(intersection);
					}
				}
			}

			if (hit.size() >= 2)
			{
				// sort hitted points from front to back (along x axis)
				std::sort(hit.begin(), hit.end(), [](glm::vec3 const& a, glm::vec3 const& b) -> bool {
					return a.x > b.x;
					});

				// add particles in between hit points
				for (size_t i = 0; i < hit.size(); i += 2)
				{
					if ((i + 1) < hit.size())
					{
						glm::vec3 from = hit[i];
						glm::vec3 to = hit[i + 1];

						for (float x = from.x; x >= to.x; x -= step)
						{
							m_solver._pos.emplace_back(x, from.y, from.z);
						}
					}
				}
			}
		}
	}

	m_solver._init_pos = m_solver._pos;
	m_solver._spatial_lookup.clear();
	m_solver._spatial_lookup.resize(m_solver.particleCount(), glm::ivec2(-1, -1));
	m_solver._pidxNeighbor.clear();
	m_solver._pidxNeighbor.resize(m_solver.particleCount());
	for(int i = 0; i < m_solver.particleCount(); ++i)
	{
		m_solver._pidxNeighbor[i].resize(c_max_particle_neighbors);
	}
	m_solver._pNumNeighbors.resize(m_solver.particleCount(), 0);

	// make sure for the other particle quantities
	m_solver._vel = std::vector<glm::vec3>(m_solver._pos.size(), glm::vec3(0.0f));
	m_solver._acc = std::vector<glm::vec3>(m_solver._pos.size(), glm::vec3(0.0f));
	m_solver._p = std::vector<float>(m_solver._pos.size(), 0);
	m_solver._d = std::vector<float>(m_solver._pos.size(), 0);
}

void Fluid::reset()
{
	if (m_solver.m_neighbors_method == SONIC_BOOM)
	{
		m_solver.reset_GPU_SSBOs();
	}
	else
	{
		m_solver._pos = m_solver._init_pos;
		// make sure for the other particle quantities
		m_solver._vel = std::vector<glm::vec3>(m_solver._pos.size(), glm::vec3(0.0f));
		m_solver._acc = std::vector<glm::vec3>(m_solver._pos.size(), glm::vec3(0.0f));
		m_solver._p = std::vector<float>(m_solver._pos.size(), 0);
		m_solver._d = std::vector<float>(m_solver._pos.size(), 0);
	}
}

void Fluid::update()
{
	if (m_is_running)
	{
		// solve 10 steps
		for (int i = 0; i < 10; ++i)
		{
			m_solver.update();
		}
	}
}

bool rayTriangleIntersection(Ray const& iRay, glm::vec3 const& p0, glm::vec3 const& p1, glm::vec3 const& p2, glm::vec3& out)
{
	std::vector<glm::vec3> hit;

	glm::vec3 o = iRay.getOrigin();
	glm::vec3 w = iRay.getDirection();
	glm::vec3 e0 = p1 - p0;
	glm::vec3 e1 = p2 - p0;
	glm::vec3 n = glm::normalize(glm::cross(e0, e1));
	glm::vec3 q = glm::cross(w, e1);
	float a = glm::dot(e0, q);
	const float epsilon = 0.001f;
	if (fabs(a) < epsilon)
	{
		return false;
	}
	glm::vec3 s = (o - p0) / a;
	glm::vec3 r = glm::cross(s, e0);
	glm::vec3 baryCoords(glm::dot(s, q), glm::dot(r, w), 1.0f - glm::dot(s, q) - glm::dot(r, w));
	if (baryCoords.x < 0.0f || baryCoords.y < 0.0f || baryCoords.z < 0.0f)
	{
		return false;
	}
	float t = glm::dot(e1, r);
	if (t >= 0.0f)
	{
		out = o + t * w;
		return true;
	}
	return false;
}