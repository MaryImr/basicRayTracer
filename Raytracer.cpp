//
//  main.cpp
//  Raytracer
//
//  Created by Maryam Imran on 30/03/2022.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "inputImageData.hpp"
#include <Eigen/Dense>


using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::ifstream;
using std::stoi;
using std::stod;
using std::ios_base;

Image readFile(string filename) {
    
    color backgroundColor;
    int maxResursionDepth;
    double shadowRayEpsilon;
    Camera cam;
    vector<Material> materials;
    Material tempMat;
    color ambient_light;
    vector<PointLight> pointLight;
    PointLight tempPL;
    vector<Eigen::Vector3d> vertexList;
    vector<Sphere> spheres;
    Sphere tempSp;
    vector<Triangle> triangles;
    Triangle tempTriangle;
    vector<Mesh> meshes;
    Mesh tempMesh;
    
    ifstream imageFile (filename);
    string readString;
    if(!imageFile.is_open()) {
        cout << "File not opened!" << endl;
        exit(1);
    }
    
    while(!imageFile.eof()) {
        imageFile >> readString;
        //cout << readString << endl;
        if (readString == "#BackgroundColor") {
            //cout << readString << endl;
            imageFile >> readString;
            backgroundColor.R = stoi(readString);
            imageFile >> readString;
            backgroundColor.G = stoi(readString);
            imageFile >> readString;
            backgroundColor.B = stoi(readString);
        }
        else if(readString == "#MaxRecursionDepth") {
            imageFile >> readString;
            maxResursionDepth = stoi(readString);
        }
        else if(readString == "#ShadowRayEpsilon") {
            imageFile >> readString;
            shadowRayEpsilon = stod(readString);
        }
        else if(readString == "#Camera") {
            
            //setting camera position
            Eigen::Vector3d temp;
            imageFile >> readString;
            temp(0) = stod(readString);
            imageFile >> readString;
            temp(1) = stod(readString);
            imageFile >> readString;
            temp(2) = stod(readString);
            cam.set_position(temp);
            
            //setting camera gaze
            imageFile >> readString;
            temp(0) = stod(readString);
            imageFile >> readString;
            temp(1) = stod(readString);
            imageFile >> readString;
            temp(2) = stod(readString);
            cam.set_gaze(temp);
            
            //setting camera up
            imageFile >> readString;
            temp(0) = stod(readString);
            imageFile >> readString;
            temp(1) = stod(readString);
            imageFile >> readString;
            temp(2) = stod(readString);
            cam.set_up(temp);
            
            //setting camera near plane
            near_plane np;
            imageFile >> readString;
            np.left = stod(readString);
            imageFile >> readString;
            np.right = stod(readString);
            imageFile >> readString;
            np.bottom = stod(readString);
            imageFile >> readString;
            np.top = stod(readString);
            cam.set_np(np);
            
            imageFile >> readString;
            cam.set_near_distance(stod(readString));
            
            resolution res;
            imageFile >> readString;
            res.width = stoi(readString);
            imageFile >> readString;
            res.height = stoi(readString);
            cam.set_res(res);
            
            cam.calculate_camera_params();
        }
        else if(readString == "#Material") {
            
            color tempCol;
            
            imageFile >> readString;
            tempMat.set_index(stoi(readString));
            
            //setting ambient reflectence
            imageFile >> readString;
            tempCol.R = stod(readString);
            imageFile >> readString;
            tempCol.G = stod(readString);
            imageFile >> readString;
            tempCol.B = stod(readString);
            tempMat.set_ambient_reflectance(tempCol);
            
            //setting diffuse reflectance
            imageFile >> readString;
            tempCol.R = stod(readString);
            imageFile >> readString;
            tempCol.G = stod(readString);
            imageFile >> readString;
            tempCol.B = stod(readString);
            tempMat.set_diffuse_reflectance(tempCol);
            
            //setting specular reflectance
            imageFile >> readString;
            tempCol.R = stod(readString);
            imageFile >> readString;
            tempCol.G = stod(readString);
            imageFile >> readString;
            tempCol.B = stod(readString);
            tempMat.set_specular_reflectance(tempCol);
            
            imageFile >> readString;
            tempMat.set_phong_exponent(stod(readString));
            
            //setting mirror reflectance
            imageFile >> readString;
            tempCol.R = stod(readString);
            imageFile >> readString;
            tempCol.G = stod(readString);
            imageFile >> readString;
            tempCol.B = stod(readString);
            tempMat.set_mirror_reflectance(tempCol);
            
            materials.push_back(tempMat);
        }
        else if(readString == "#AmbientLight") {
            
            imageFile >> readString;
            ambient_light.R = stod(readString);
            imageFile >> readString;
            ambient_light.G = stod(readString);
            imageFile >> readString;
            ambient_light.B = stod(readString);
        }
        else if (readString == "#PointLight") {
            
            coordinates tempPos;
            color tempCol;
            
            imageFile >> readString;
            tempPL.set_index(stoi(readString));
            
            //setting position
            imageFile >> readString;
            tempPos.x = stod(readString);
            imageFile >> readString;
            tempPos.y = stod(readString);
            imageFile >> readString;
            tempPos.z = stod(readString);
            tempPL.set_position(tempPos);
            
            //setting intensity
            imageFile >> readString;
            tempCol.R = stod(readString);
            imageFile >> readString;
            tempCol.G = stod(readString);
            imageFile >> readString;
            tempCol.B = stod(readString);
            tempPL.set_intensity(tempCol);
            
            pointLight.push_back(tempPL);
        }
        else if (readString == "#VertexList") {
            
            Eigen::Vector3d tempVert;
            imageFile >> readString;
            do{
                tempVert(0) = stod(readString);
                imageFile >> readString;
                tempVert(1) = stod(readString);
                imageFile >> readString;
                tempVert(2) = stod(readString);
                vertexList.push_back(tempVert);
                imageFile >> readString;
            }
            while(readString[0] != '#');
                
            imageFile.seekg(-readString.length(), ios_base::cur);
        }
        else if (readString == "#Sphere") {
            
            imageFile >> readString;
            tempSp.set_index(stoi(readString));
            imageFile >> readString;
            tempSp.set_material_index(stoi(readString));
            imageFile >> readString;
            tempSp.set_center_vertex_id(stoi(readString));
            imageFile >> readString;
            tempSp.set_radius(stod(readString));
            
            spheres.push_back(tempSp);
        }
        else if (readString == "#Triangle") {
            
            int tempVert[3];
            
            imageFile >> readString;
            tempTriangle.set_index(stoi(readString));
            imageFile >> readString;
            tempTriangle.set_material_index(stoi(readString));
            imageFile >> readString;
            tempVert[0] = stoi(readString);
            imageFile >> readString;
            tempVert[1] = stoi(readString);
            imageFile >> readString;
            tempVert[2] = stoi(readString);
            tempTriangle.set_vertexes(tempVert);
            
            triangles.push_back(tempTriangle);
        }
        else if (readString == "#Mesh") {
            
            vector<vector<int>> triangle_vertex_indices;
            vector<int> vertex_ind;
            
            imageFile >> readString;
            tempMesh.set_index(stoi(readString));
            imageFile >> readString;
            tempMesh.set_material_index(stoi(readString));
            
            imageFile.seekg(1, ios_base::cur);
            vertex_ind.resize(3);
            while(imageFile.peek() != '\n') {
                imageFile >> readString;
                vertex_ind.at(0) = (stoi(readString));
                imageFile >> readString;
                vertex_ind.at(1) = (stoi(readString));
                imageFile >> readString;
                vertex_ind.at(2) = (stoi(readString));
                triangle_vertex_indices.push_back(vertex_ind);
                imageFile.seekg(1, ios_base::cur);
            }
            
            tempMesh.set_triangle_vertex_indices(triangle_vertex_indices.size(), triangle_vertex_indices);
            meshes.push_back(tempMesh);
        }
    }
    
    Image inputImage(backgroundColor, maxResursionDepth, shadowRayEpsilon, cam, materials.size(), materials, ambient_light, pointLight.size(), pointLight, vertexList.size(), vertexList, spheres.size(), spheres, triangles.size(), triangles, meshes.size(), meshes);
    
    return inputImage;
}

int main(int argc, const char * argv[]) {
    
    Image inputImage = readFile("/Users/maryamimran/Documents/Undergraduate/3rd Year/6th Semester/CNG 477/HWs/homework1/hw1/input3.txt");
    outputImage myImage;
    myImage.initialize(inputImage.get_cam().get_res().width, inputImage.get_cam().get_res().height);
    int i, j, shapeCount, meshTriCount;
    bool valid;
    
    for(i=0; i<inputImage.get_cam().get_res().width; i++) {
        for(j=0; j<inputImage.get_cam().get_res().height; j++) {
            Ray camR = inputImage.get_cam().compute_ray(i, j);
            //intersecting with all spheres
            for(shapeCount = 0; shapeCount < inputImage.get_spheres().size(); shapeCount++) {
                valid = inputImage.get_spheres()[shapeCount].test_intersection(camR, inputImage.get_vertexList()[inputImage.get_spheres()[shapeCount].get_center_vertex_id()-1]);
                if(valid) {
                    myImage.setPixel(i, j, 255.0, 0.0, 0.0);
                }
                else {
                    myImage.setPixel(i, j, inputImage.get_backgroundColor().R, inputImage.get_backgroundColor().G, inputImage.get_backgroundColor().B);
                }
            }
            //intersecting with all triangles
            for(shapeCount = 0; shapeCount < inputImage.get_triangles().size(); shapeCount++) {
                valid = inputImage.get_triangles()[shapeCount].test_intersection(camR, inputImage.get_vertexList()[inputImage.get_triangles()[shapeCount].get_vertexes()[0]-1], inputImage.get_vertexList()[inputImage.get_triangles()[shapeCount].get_vertexes()[1]-1], inputImage.get_vertexList()[inputImage.get_triangles()[shapeCount].get_vertexes()[2]-1]);
                if(valid) {
                    myImage.setPixel(i, j, 255.0, 0.0, 0.0);
                }
            }
            //intersecting with all meshes
            for(shapeCount = 0; shapeCount < inputImage.get_meshes().size(); shapeCount++) {
                valid = false;
                for(meshTriCount = 0; meshTriCount < inputImage.get_meshes()[shapeCount].get_triangle_vertex_indices().size(); meshTriCount++) {
                    if(!valid) {
                        valid = inputImage.get_triangles()[meshTriCount].test_intersection(camR, inputImage.get_vertexList()[inputImage.get_meshes()[shapeCount].get_triangle_vertex_indices()[meshTriCount][0]-1], inputImage.get_vertexList()[inputImage.get_meshes()[shapeCount].get_triangle_vertex_indices()[meshTriCount][1]-1], inputImage.get_vertexList()[inputImage.get_meshes()[shapeCount].get_triangle_vertex_indices()[meshTriCount][2]-1]);
                        if(valid) {
                            myImage.setPixel(i, j, 255.0, 0.0, 0.0);
                        }
                    }
                }
            }
        }
    }
    
    myImage.get_output_file();
    
    cout << "Successfull!!" << endl;
    
    return 0;
}
