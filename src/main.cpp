#include "sph_fluid.hpp"
#include "camera.hpp"
#include "shader.hpp"
#include "light.hpp"
#include "grid_axis.hpp"

// ==================== FUNCTION HEADERS ====================
// ==========================================================
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void initScene();
void renderScene();
void draw_UI();
void ImGui_init(GLFWwindow* win);
void set_particle_shader();
void set_sonic_boom_shader();

// ==================== GLOBAL OBJECTS ====================
// ========================================================

// camera
std::shared_ptr<Camera> g_cam;

// mesh and fluid
std::shared_ptr<Mesh> g_mesh_water_tank;
std::shared_ptr<Fluid> g_fluid;

// grid axis
std::shared_ptr<GridAxis> g_grid_axis;

// shaders
std::shared_ptr<Shader> g_shader_particle;
std::shared_ptr<Shader> g_shader_sonic_boom;

// user interface
struct Mouse g_mouse;
struct UI g_ui;
double g_prev_time;
double g_curr_time;
float g_print_ms;

// ==================== FUNCTION DEFINITIONS ====================
// ==============================================================

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	g_cam->updateProjectionMatrix(width, height);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	g_cam->zoomView(static_cast<float>(yoffset));
}

void processInput(GLFWwindow* window)
{
	if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
		glfwSetWindowShouldClose(window, true);
	}

	if (glfwGetKey(window, GLFW_KEY_F12) == GLFW_PRESS) // recompile diffuse shader
	{
		g_shader_particle = std::make_shared<Shader>("../shaders/particle/vertex.glsl", "../shaders/particle/fragment.glsl");
	}

	// CAMERA RELATED
	double xpos;
	double ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	g_mouse.m_prevX = g_mouse.m_currX;
	g_mouse.m_prevY = g_mouse.m_currY;
	g_mouse.m_currX = xpos;
	g_mouse.m_currY = ypos;
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS && glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) != GLFW_PRESS)
	{
		g_cam->rotateView(g_mouse);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS && glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
	{
		g_cam->panView(g_mouse);
	}
}

void initScene()
{
	g_grid_axis = std::make_shared<GridAxis>(50);

	g_cam = std::make_shared<Camera>(glm::vec3(25.0f, 20.0f, 0.0f), glm::vec3(0.0f, 5.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), VIEWPORT_WIDTH, VIEWPORT_HEIGHT);
	
	std::shared_ptr<Mesh> mesh_water_volume = std::make_shared<Mesh>("../assets/water_volume.obj");
	g_mesh_water_tank = std::make_shared<Mesh>("../assets/water_tank.obj");
	g_fluid = std::make_shared<Fluid>(mesh_water_volume, g_ui.m_particle_radius, g_mesh_water_tank);

	g_shader_particle = std::make_shared<Shader>("../shaders/particle/vertex.glsl", "../shaders/particle/fragment.glsl");
	g_shader_sonic_boom = std::make_shared<Shader>("../shaders/sonic_boom/vertex.glsl", "../shaders/sonic_boom/fragment.glsl");
	
	g_fluid->m_solver.set_neighbors_compute_method(SONIC_BOOM);
}

void renderScene()
{
	// draw grid
	g_grid_axis->draw(g_cam->getViewMatrix(), g_cam->getProjectionMatrix());

	if (g_fluid->m_solver.m_neighbors_method == SONIC_BOOM)
	{
		set_sonic_boom_shader();
		g_fluid->draw(g_shader_sonic_boom, g_ui.m_draw_mode == 0);
	}
	else
	{
		set_particle_shader();
		g_fluid->draw(g_shader_particle, g_ui.m_draw_mode == 0);
	}

	// reset
	glBindVertexArray(0);
	glActiveTexture(GL_TEXTURE0);
}

void draw_UI()
{
	ImGui::SetNextWindowPos(ImVec2(0, 0));
	ImGui::SetNextWindowSize(ImVec2(300, 80));
	ImGui::Begin("Drawing mode");
	ImGui::RadioButton("wireframe", &g_ui.m_draw_mode, 0);
	ImGui::RadioButton("shaded", &g_ui.m_draw_mode, 1);
	ImGui::End();

	ImGui::SetNextWindowPos(ImVec2(0, 80));
	ImGui::SetNextWindowSize(ImVec2(300, 80));
	ImGui::Begin("Particle");
	ImGui::SliderFloat("radius", &g_ui.m_particle_radius, 0.250f, 0.330f);
	ImGui::End();

	ImGui::SetNextWindowPos(ImVec2(0, 160));
	ImGui::SetNextWindowSize(ImVec2(300, 240));
	ImGui::Begin("SPH simulation");
	if (ImGui::Button("run"))
	{
		g_fluid->m_is_running = true;
	}
	if (ImGui::Button("pause"))
	{
		g_fluid->m_is_running = false;
	}
	if (ImGui::Button("reset"))
	{
		g_fluid->reset();
	}
	ImGui::Checkbox("use OpenMP", &g_fluid->m_solver.m_use_omp);
	ImGui::RadioButton("Brute Force", reinterpret_cast<int*>(& g_fluid->m_solver.m_neighbors_method), EGGMAN);
	ImGui::RadioButton("CPU - 3D grid lookup", reinterpret_cast<int*>(&g_fluid->m_solver.m_neighbors_method), SONIC);
	ImGui::RadioButton("GPU", reinterpret_cast<int*>(&g_fluid->m_solver.m_neighbors_method), SONIC_BOOM);
	if (g_fluid->m_is_running)
	{
		std::string frame_ms = std::to_string(g_print_ms) + " ms";
		ImGui::Text(frame_ms.c_str());
	}
	else
	{
		ImGui::Text("0 ms");
	}
	std::string particle_count_txt = std::string("particle count = ") + std::to_string(g_fluid->m_solver.particleCount());
	ImGui::Text(particle_count_txt.c_str());
	ImGui::End();
}

void ImGui_init(GLFWwindow* win)
{
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	(void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(win, true);
	ImGui_ImplOpenGL3_Init("#version 460");
}

void set_particle_shader()
{
	g_shader_particle->use();
	g_shader_particle->set("view", g_cam->getViewMatrix());
	g_shader_particle->set("proj", g_cam->getProjectionMatrix());
}

void set_sonic_boom_shader()
{
	g_shader_sonic_boom->use();
	g_shader_sonic_boom->set("view", g_cam->getViewMatrix());
	g_shader_sonic_boom->set("proj", g_cam->getProjectionMatrix());
}

int main()
{
	// ==================== INIT GLFW ====================
	// ===================================================
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(VIEWPORT_WIDTH, VIEWPORT_HEIGHT, "SPH", NULL, NULL);
	if(window == NULL)
	{
		std::cerr << "Failed to create GLFW window !" << std::endl;
		glfwTerminate();
		std::exit(-1);
	}

	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetScrollCallback(window, scroll_callback);

	// ==================== INIT GLAD ====================
	// ===================================================
	if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cerr << "Failed to initialize GLAD !" << std::endl;
		std::exit(-1);
	}

	// ==================== OpenGL STATES ====================
	// =======================================================
	glViewport(0, 0, VIEWPORT_WIDTH, VIEWPORT_HEIGHT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// ==================== MAIN LOOP ====================
	// ===================================================
	ImGui_init(window);
	g_ui.m_draw_mode = 1;
	g_ui.m_particle_radius = 0.25f;

	initScene();

	g_prev_time = glfwGetTime();
	while(!glfwWindowShouldClose(window))
	{
		processInput(window);

		// ImGui new frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		// update particle radius
		if (g_ui.m_particle_radius != g_fluid->m_particle_radius)
		{
			g_fluid->update_particles_radius(g_ui.m_particle_radius);
		}

		g_fluid->update();

		renderScene();

		// draw UI
		draw_UI();
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
		glfwPollEvents();

		g_curr_time = glfwGetTime();
		if (g_curr_time - g_prev_time > 0.25)
		{
			g_print_ms = g_fluid->m_solver._milliseconds * 1000.0f;
			g_prev_time = g_curr_time;
		}
	}

	// ==================== CLEANUP ====================
	// =================================================
	glfwTerminate();
	return 0;
}