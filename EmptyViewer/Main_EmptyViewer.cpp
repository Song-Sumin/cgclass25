#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------

// class Ray
class Ray {
public:
	vec3 origin;
	vec3 direction;
};

// class Surface and virtual bool intersect
class Surface {
public:
	virtual bool intersect(const Ray& ray, float tMin, float tMax, float& t) const = 0;
};

// make Sphere using discriminant
class Sphere : public Surface {
public:
	vec3 center;
	float radius;

	Sphere(const vec3& c, float r) : center(c), radius(r) {}

	bool intersect(const Ray& ray, float tMin, float tMax, float& t) const override {
		vec3 oc = ray.origin - center;
		float a = dot(ray.direction, ray.direction);
		float b = 2.0f * dot(oc, ray.direction);
		float c = dot(oc, oc) - radius * radius;
		float discriminant = b * b - 4 * a * c;
		if (discriminant > 0) {
			float temp = (-b - std::sqrt(discriminant)) / (2.0f * a);
			if (temp < tMax && temp > tMin) {
				t = temp;
				return true;
			}
			temp = (-b + std::sqrt(discriminant)) / (2.0f * a);
			if (temp < tMax && temp > tMin) {
				t = temp;
				return true;
			}
		}
		return false;
	}
};

// class Plane
class Plane : public Surface {
public:
	vec3 point;
	vec3 normal;

	Plane(const vec3& p, const vec3& n) : point(p), normal(normalize(n)) {}

	bool intersect(const Ray& ray, float tMin, float tMax, float& t) const override {
		float denom = dot(normal, ray.direction);
		if (abs(denom) > 1e-6) {
			vec3 p0l0 = point - ray.origin;
			t = dot(p0l0, normal) / denom;
			if (t >= tMin && t <= tMax) {
				return true;
			}
		}
		return false;
	}
};

// set Camera and get Ray
class Camera {
public:
	vec3 eye;
	vec3 u, v, w;
	float l, r, b, t, d;

	Camera() {
		eye = vec3(0.0f, 0.0f, 0.0f);
		u = vec3(1.0f, 0.0f, 0.0f);
		v = vec3(0.0f, 1.0f, 0.0f);
		w = vec3(0.0f, 0.0f, 1.0f);
		l = -0.1f;
		r = 0.1f;
		b = -0.1f;
		t = 0.1f;
		d = 0.1f;
	}

	Ray getRay(float i, float j) const {
		float u_coord = l + (r - l) * (i + 0.5f) / Width;
		float v_coord = b + (t - b) * (j + 0.5f) / Height;
		vec3 direction = normalize(u_coord * u + v_coord * v - d * w);
		return Ray{ eye, direction };
	}
};

// class Scene if hit object return white else return black
class Scene {
public:
	std::vector<Surface*> surfaces;

	vec3 trace(const Ray& ray, float tMin, float tMax) const {
		float closest_t = tMax;
		const Surface* hit_surface = nullptr;
		for (const auto& surface : surfaces) {
			float t;
			if (surface->intersect(ray, tMin, closest_t, t)) {
				closest_t = t;
				hit_surface = surface;
			}
		}
		if (hit_surface) {
			return vec3(1.0f, 1.0f, 1.0f); // white
		}
		return vec3(0.0f, 0.0f, 0.0f); // black
	}
};

void render()
{
	//Create our image. We don't want to do this in 
	//the main loop since this may be too slow and we 
	//want a responsive display of our beautiful image.
	//Instead we draw to another buffer and copy this to the 
	//framebuffer using glDrawPixels(...) every refresh

	//Create camera and scene
	Camera camera;
	Scene scene;
	scene.surfaces.push_back(new Plane(vec3(0.0f, -2.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f))); // Plane P
	scene.surfaces.push_back(new Sphere(vec3(-4.0f, 0.0f, -7.0f), 1.0f)); // Sphere S1
	scene.surfaces.push_back(new Sphere(vec3(0.0f, 0.0f, -7.0f), 2.0f)); // Sphere S2
	scene.surfaces.push_back(new Sphere(vec3(4.0f, 0.0f, -7.0f), 1.0f)); // Sphere S3

	OutputImage.clear();
	for (int j = 0; j < Height; ++j)
	{
		for (int i = 0; i < Width; ++i)
		{
			// Create a ray from the camera and check plane and spheres with scene.trace
			Ray ray = camera.getRay(i, j);
			vec3 color = scene.trace(ray, 0.0f, std::numeric_limits<float>::max());

			OutputImage.push_back(color.x); // R
			OutputImage.push_back(color.y); // G
			OutputImage.push_back(color.z); // B
		}
	}
}

void resize_callback(GLFWwindow*, int nw, int nh)
{
	//This is called in response to the window resizing.
	//The new width and height are passed in so we make 
	//any necessary changes:
	Width = nw;
	Height = nh;
	//Tell the viewport to use all of our screen estate
	glViewport(0, 0, nw, nh);

	//This is not necessary, we're just working in 2d so
	//why not let our spaces reflect it?
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0.0, static_cast<double>(Width)
		, 0.0, static_cast<double>(Height)
		, 1.0, -1.0);

	//Reserve memory for our render so that we don't do 
	//excessive allocations and render the image
	OutputImage.reserve(Width * Height * 3);
	render();
}

int main(int argc, char* argv[])
{
	// -------------------------------------------------
	// Initialize Window
	// -------------------------------------------------

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//We have an opengl context now. Everything from here on out 
	//is just managing our window or opengl directly.

	//Tell the opengl state machine we don't want it to make 
	//any assumptions about how pixels are aligned in memory 
	//during transfers between host and device (like glDrawPixels(...) )
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//We call our resize function once to set everything up initially
	//after registering it as a callback with glfw
	glfwSetFramebufferSizeCallback(window, resize_callback);
	resize_callback(NULL, Width, Height);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// -------------------------------------------------------------
		//Rendering begins!
		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
		//and ends.
		// -------------------------------------------------------------

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		//Close when the user hits 'q' or escape
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
	}

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
