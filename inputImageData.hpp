//
//  inputImageData.hpp
//  Raytracer
//
//  Created by Maryam Imran on 30/03/2022.
//
#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <Eigen/Dense>


using std::endl;
using std::cout;
using std::vector;
using std::ios;
using std::ofstream;

struct color {
    int R;
    int G;
    int B;
};

struct coordinates {
    double x;
    double y;
    double z;
};

struct near_plane {
    double left;
    double right;
    double bottom;
    double top;
};

struct resolution {
    int width;
    int height;
};

class Ray {
private:
    Eigen::Vector3d pointA;
    Eigen::Vector3d pointB;
    Eigen::Vector3d join;
public:
    Ray(Eigen::Vector3d pointA, Eigen::Vector3d pointB) {
        this->pointA = pointA;
        this->pointB = pointB;
        this->join = this->pointB - this->pointA;
    }
    
    Eigen::Vector3d get_pointA(void) {
        return this->pointA;
    }
    
    Eigen::Vector3d get_pointB(void) {
        return this->pointB;
    }
    
    Eigen::Vector3d get_join(void) {
        return this->join;
    }
};

class Camera {
private:
    Eigen::Vector3d position;
    Eigen::Vector3d gaze;
    Eigen::Vector3d up;
    near_plane np;
    double near_distance;
    resolution res;
    
    //things to calculate
    Eigen::Vector3d m;
    Eigen::Vector3d u;
    Eigen::Vector3d q;
    
public:
    
    Camera() {}
    
    void calculate_camera_params(void) {
        this->m = this->position + (this->near_distance*this->gaze);
        this->u = this->gaze.cross(this->up);
        this->q = this->m + (this->np.left*this->u) + (this->np.top*this->up);
    }
    
    Ray compute_ray(int pixel_i, int pixel_j) {
        double su = calculate_su(pixel_i);
        double sv = calculate_sv(pixel_j);
        Eigen::Vector3d s = calculate_s(su, sv);
        
        return(Ray(this->position, s));
    }
    
    Eigen::Vector3d get_position(void) {
        return this->position;
    }
    
    Eigen::Vector3d get_gaze(void) {
        return this->gaze;
    }
    
    Eigen::Vector3d get_up(void) {
        return this->up;
    }
    
    near_plane get_np(void) {
        return this->np;
    }
    
    double get_near_distance(void) {
        return this->near_distance;
    }
    
    resolution get_res(void) {
        return this->res;
    }
    
    Eigen::Vector3d get_m(void) {
        return this->m;
    }
    
    Eigen::Vector3d get_u(void) {
        return this->u;
    }
    
    Eigen::Vector3d get_q(void) {
        return this->q;
    }
    
    void set_position(Eigen::Vector3d position) {
        this->position = position;
    }
    
    void set_gaze(Eigen::Vector3d gaze) {
        this->gaze = gaze;
    }
    
    void set_up(Eigen::Vector3d up) {
        this->up = up;
    }
    
    void set_np(near_plane np) {
        this->np = np;
    }
    
    void set_near_distance(double near_distance) {
        this->near_distance = near_distance;
    }
    
    void set_res(resolution res) {
        this->res = res;
    }
    
    void set_m(Eigen::Vector3d m) {
        this->m = m;
    }
    
    void set_u(Eigen::Vector3d u) {
        this->u = u;
    }
    
    void set_q(Eigen::Vector3d q) {
        this->q = q;
    }
    
    double calculate_su(int i) {
        double width = this->np.right - this->np.left;
        double su = (width*(i+0.5))/this->res.width;
        return su;
    }
    
    double calculate_sv(int j) {
        double height = this->np.top - this->np.bottom;
        double sv = (height*(j+0.5))/this->res.height;
        return sv;
    }
    
    Eigen::Vector3d calculate_s(double su, double sv) {
        return (this->q + su*this->u - sv*this->up);
    }
};

class Material {
private:
    int index;
    color ambient_reflectance;
    color diffuse_reflectance;
    color specular_reflectance;
    double phong_exponent;
    color mirror_reflectance;
public:
    Material() {}
    
    int get_index(void) {
        return this->index;
    }
    
    color get_ambient_reflectance(void) {
        return this->ambient_reflectance;
    }
    
    color get_diffuse_reflectance(void) {
        return this->diffuse_reflectance;
    }
    
    color get_specular_reflectance(void) {
        return this->specular_reflectance;
    }
    
    double get_phong_exponent(void) {
        return this->phong_exponent;
    }
    
    color get_mirror_reflectance(void) {
        return this->mirror_reflectance;
    }
    
    void set_index(int index) {
        this->index = index;
    }
    
    void set_ambient_reflectance(color ambient_reflectance) {
        this->ambient_reflectance = ambient_reflectance;
    }
    
    void set_diffuse_reflectance(color diffuse_reflectance) {
        this->diffuse_reflectance = diffuse_reflectance;
    }
    
    void set_specular_reflectance(color specular_reflectance) {
        this->specular_reflectance = specular_reflectance;
    }
    
    void set_phong_exponent(double phong_exponent) {
        this->phong_exponent = phong_exponent;
    }
    
    void set_mirror_reflectance(color mirror_reflectance) {
        this->mirror_reflectance = mirror_reflectance;
    }
};

class PointLight {
private:
    int index;
    coordinates position;
    color intensity;
public:
    PointLight() {}
    
    int get_index(void) {
        return this->index;
    }
    
    coordinates get_position(void) {
        return this->position;
    }
    
    color get_intensity(void) {
        return this->intensity;
    }
    
    void set_index(int index) {
        this->index = index;
    }
    
    void set_position(coordinates position) {
        this->position = position;
    }
    
    void set_intensity(color intensity) {
        this->intensity = intensity;
    }
};

class Sphere {
private:
    int index;
    int material_index;
    int center_vertex_id;
    double radius;
public:
    Sphere() {}
    
    bool test_intersection(Ray test_ray, Eigen::Vector3d center) {
        Eigen::Vector3d e = test_ray.get_pointA();
        Eigen::Vector3d s = test_ray.get_pointB();
        Eigen::Vector3d d = s - e;
        double a = d.dot(d);
        double b = (2.0 * d).dot(e - center);
        double c = (e - center).dot((e - center)) - pow(this->radius, 2);
        double D = pow(b,2) - 4.0 * a * c;
        //cout << "D = " << D << endl;
        if(D > 0.0)
            return true;
        else
            return false;
    }
    
    int get_index(void) {
        return this->index;
    }
    
    int get_material_index(void) {
        return this->material_index;
    }
    
    int get_center_vertex_id(void) {
        return this->center_vertex_id;
    }
    
    double get_radius(void) {
        return this->radius;
    }
    
    void set_index(int index) {
        this->index = index;
    }
    
    void set_material_index(int material_index) {
        this->material_index = material_index;
    }
    
    void set_center_vertex_id(int center_vertex_id) {
        this->center_vertex_id = center_vertex_id;
    }
    
    void set_radius(double radius) {
        this->radius = radius;
    }
};

class Triangle {
private:
    int index;
    int material_index;
    int *vertexes;
public:
    Triangle() {}
    
    int get_index(void) {
        return this->index;
    }
    
    int get_material_index(void) {
        return this->material_index;
    }
    
    int* get_vertexes(void) {
        return this->vertexes;
    }
    
    void set_index(int index) {
        this->index = index;
    }
    
    void set_material_index(int material_index) {
        this->material_index = material_index;
    }
    
    void set_vertexes(int *vertexes) {
        int i;
        this->vertexes = new int[3];
        for(i=0; i<3; i++)
            this->vertexes[i] = vertexes[i];
    }
    bool test_intersection(Ray test_ray, Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2)
    {
        // compute plane's normal
        Eigen::Vector3d v0v1 = v1 - v0;
        Eigen::Vector3d v0v2 = v2 - v0;
        // no need to normalize
        Eigen::Vector3d N = v0v1.cross(v0v2); // N

        float t;
     
        // Step 1: finding P
     
        // check if ray and plane are parallel ?
        double NdotRayDirection = N.dot(test_ray.get_join());
        if (fabs(NdotRayDirection) < pow(2,-52)) // almost 0
            return false; // they are parallel so they don't intersect !
     
        // compute d parameter using equation 2
        float d = -N.dot(v0);
     
        // compute t (equation 3)
        t = -(N.dot(test_ray.get_pointA()) + d) / NdotRayDirection;
     
        // check if the triangle is in behind the ray
        if (t < 0) return false; // the triangle is behind
     
        // compute the intersection point using equation 1
        Eigen::Vector3d P = test_ray.get_pointA() + t * test_ray.get_join();
     
        // Step 2: inside-outside test
        Eigen::Vector3d C; // vector perpendicular to triangle's plane
     
        // edge 0
        Eigen::Vector3d edge0 = v1 - v0;
        Eigen::Vector3d vp0 = P - v0;
        C = edge0.cross(vp0);
        if (N.dot(C) < 0) return false; // P is on the right side
     
        // edge 1
        Eigen::Vector3d edge1 = v2 - v1;
        Eigen::Vector3d vp1 = P - v1;
        C = edge1.cross(vp1);
        if (N.dot(C) < 0)  return false; // P is on the right side
     
        // edge 2
        Eigen::Vector3d edge2 = v0 - v2;
        Eigen::Vector3d vp2 = P - v2;
        C = edge2.cross(vp2);
        if (N.dot(C) < 0) return false; // P is on the right side;
     
        return true; // this ray hits the triangle
    }
};

class Mesh {
private:
    int index;
    int material_index;
    int total_triangles;
    vector<vector<int>> triangle_vertex_indices;
public:
    Mesh() {}
    
    int get_index(void) {
        return this->index;
    }
    
    int get_material_index(void) {
        return this->material_index;
    }
    
    vector<vector<int>> get_triangle_vertex_indices(void) {
        return this->triangle_vertex_indices;
    }
    
    int get_total_triangles(void) {
        return this->total_triangles;
    }
    
    void set_index(int index) {
        this->index = index;
    }
    
    void set_material_index(int material_index) {
        this->material_index = material_index;
    }
    
    void set_triangle_vertex_indices(int total_triangles, vector<vector<int>> triangle_vertex_indices) {
        int i;
        for(i=0; i<total_triangles; i++) {
            this->triangle_vertex_indices.push_back(triangle_vertex_indices[i]);
        }
        this->total_triangles = total_triangles;
    }
};

class Image {
private:
    color backgroundColor;
    int maxResursionDepth;
    double shadowRayEpsilon;
    Camera cam;
    vector<Material> materials;
    color ambient_light;
    vector<PointLight> pointLight;
    vector<Eigen::Vector3d> vertexList;
    vector<Sphere> spheres;
    vector<Triangle> triangles;
    vector<Mesh> meshes;
public:
    Image(color backgroundColor, int maxResursionDepth, double shadowRayEpsilon, Camera cam, int total_materials, vector<Material> materials, color ambient_light, int total_pointlights, vector<PointLight> pointLight, int total_vertexes, vector<Eigen::Vector3d> vertexList, int total_spheres, vector<Sphere> spheres, int total_triangles, vector<Triangle> triangles, int total_meshes, vector<Mesh> meshes) {
        
        int i;
        
        this->backgroundColor = backgroundColor;
        this->maxResursionDepth = maxResursionDepth;
        this->shadowRayEpsilon = shadowRayEpsilon;
        
        this->cam.set_gaze(cam.get_gaze());
        this->cam.set_near_distance(cam.get_near_distance());
        this->cam.set_np(cam.get_np());
        this->cam.set_up(cam.get_up());
        this->cam.set_res(cam.get_res());
        this->cam.set_position(cam.get_position());
        this->cam.set_m(cam.get_m());
        this->cam.set_u(cam.get_u());
        this->cam.set_q(cam.get_q());
        
        for(i=0; i<total_materials; i++) {
            this->materials.push_back(materials[i]);
        }
        
        this->ambient_light = ambient_light;
        
        for(i=0; i<total_pointlights; i++) {
            this->pointLight.push_back(pointLight[i]);
        }
        
        for(i=0; i<total_vertexes; i++) {
            this->vertexList.push_back(vertexList[i]);
        }
        
        for(i=0; i<total_spheres; i++) {
            this->spheres.push_back(spheres[i]);
        }
        
        for (i=0; i<total_triangles; i++) {
            this->triangles.push_back(triangles[i]);
        }
        
        for(i=0; i<total_meshes; i++) {
            this->meshes.push_back(meshes[i]);
        }
        
    }
    
    color get_backgroundColor(void) {
        return this->backgroundColor;
    }
    
    int get_maxResursionDepth(void) {
        return this->maxResursionDepth;
    }
    
    double get_shadowRayEpsilon(void) {
        return this->shadowRayEpsilon;
    }
    
    Camera get_cam(void) {
        return this->cam;
    }
    
    vector<Material> get_materials(void) {
        return this->materials;
    }
    
    color get_ambient_light(void) {
        return this->ambient_light;
    }
    
    vector<PointLight> get_pointLight(void) {
        return this->pointLight;
    }
    
    vector<Eigen::Vector3d> get_vertexList(void) {
        return this->vertexList;
    }
    
    vector<Sphere> get_spheres(void) {
        return this->spheres;
    }
    
    vector<Triangle> get_triangles(void) {
        return this->triangles;
    }
    
    vector<Mesh> get_meshes(void) {
        return this->meshes;
    }
};

class outputImage {
private:
    vector<vector<color>> imagePixels;
    int x_res;
    int y_res;
public:
    outputImage() {
        this->x_res = 0;
        this->y_res = 0;
    }
    
    void initialize(int x_res, int y_res) {
        imagePixels.resize(x_res, vector<color>(y_res));
        this->x_res = x_res;
        this->y_res = y_res;
    }
    
    void setPixel(int x, int y, double r, double g, double b) {
        imagePixels.at(x).at(y).R = r;
        imagePixels.at(x).at(y).G = g;
        imagePixels.at(x).at(y).B = b;
    }
    
    vector<vector<color>> get_image(void) {
        return imagePixels;
    }
    
    int get_x_res(void) {
        return this->x_res;
    }
    
    int get_y_res(void) {
        return this->y_res;
    }
    
    void get_output_file() {
        
        ofstream outputImageFile("/Users/maryamimran/Desktop/output.ppm", ios::binary);
        if(outputImageFile.is_open()) {
            outputImageFile << "P3" << endl;
            outputImageFile << this->x_res << endl;
            outputImageFile << this->y_res << endl;
            outputImageFile << 255 << endl;
            
            int i, j;
            
            for(i=0; i<this->x_res; i++) {
                //cout << "i = " << i << endl;
                for(j=0; j<this->y_res; j++) {
                    outputImageFile << (int) imagePixels[i][j].R << ' ';
                    outputImageFile << (int) imagePixels[i][j].G << ' ';
                    outputImageFile << (int) imagePixels[i][j].B << '\n';
                }
            }
            outputImageFile.close();
        }
    }
};

