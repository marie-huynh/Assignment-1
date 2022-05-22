#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <tuple>

#include <iostream>
#include <random>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <list>

static std::default_random_engine engine(10) ; // random seed=10
static std::uniform_real_distribution<double> uniform (0,1);


class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector times(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}


Vector random_cos(const Vector &N){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2*M_PI*r1) * sqrt(1-r2);
    double y = sin(2*M_PI*r1) * sqrt(1-r2);
    double z = sqrt(r2);
    Vector T1;
    if ((abs(N[0]) < abs(N[1])) && (abs(N[0]) < abs(N[2]))){
        T1[0] = 0;
        T1[1] = -N[2]; 
        T1[2] = N[1];
        
    }
    else if ((abs(N[1]) < abs(N[0])) && (abs(N[1]) < abs(N[2]))){
        T1[0] = -N[2];
        T1[1] = 0; 
        T1[2] = N[0];
    }
    else {
        T1[0] = -N[1];
        T1[1] = N[0]; 
        T1[2] = 0;
    }
    T1.normalize();
    Vector T2 = cross(T1, N);
    Vector V = x*T1 + y*T2 + z*N;
    return V;
}

Vector boxMuller (double stdev) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = sqrt(-2 * log(r1)) *cos(2 * M_PI*r2) * stdev;
    double y = sqrt(-2 * log(r1)) *sin(2 * M_PI*r2) * stdev;
    return Vector(x, y, 0.);
}

//#############################################################################
class Ray {
public:
    //origin vector O and unit direction u
    Vector O;
    Vector u;
    //We define the constructor
    Ray(const Vector& origin, const Vector& unit){
        O = origin;
        u = unit;
    }
}; 


//#############################################################################

class Geometry {
    public :
        Vector albedo;
        bool mirror;
        virtual bool intersect(const Ray& ray, Vector &P, Vector &N, double &t, int &index) = 0;
};

//#############################################################################

class Sphere : public Geometry{
public:
    // center vector C and a double radius R
    Vector C;
    double R;
    
    
    Sphere(const Vector& center, const double radius, const Vector& a, bool m){
        C = center;
        R = radius;
        albedo = a;
        mirror = m;
    }
    // we define the function intersect that computes the point of intersection
    // between a ray and a sphere if any and returns a bool
    bool intersect(const Ray& ray, Vector &P, Vector &N, double &t, int &index)  {
        //We want to solve t^2 + 2 t (u, O − C) + ‖O − C‖^2 − R^2 = 0
        //We compute the discriminant delta.
        //double a = dot(ray.u, ray.O - C);
        double delta = pow(dot(ray.u, ray.O - C),2) - ((ray.O - C).norm2() - R*R);
        //If delta < 0, we have no intersection.
        if (delta < 0){
            return false;
        }
       // double sqrt_delta = sqrt(delta);
        double t1, t2;
        t1 = -dot(ray.u, ray.O - C) - sqrt(delta);
        t2 = -dot(ray.u, ray.O - C) + sqrt(delta);
        if (t2 < 0) {
            return false;
        }
        if (t1 >= 0){
            t = t1;
        }
        else{
            t = t2;
        }
        P = ray.O + t*ray.u;

        //We retrieve unit normal N at P
        N = P - C;
        N.normalize();
        return true;
    }

};


//#############################################################################

class Scene {
public : 
    std::vector<Geometry*> objects;
    Vector L;
    double I;

    Scene(){};

    void add_an_object(Geometry& object){
        objects.push_back(&object);
    }

    bool intersect(const Ray& ray, Vector &P, Vector &N, double &t, int &index){
        t = 1e9;
        bool intersect = false;
        for (int i = 0; i < objects.size(); i++){
            Vector P_sphere, N_sphere; 
            double t_sphere;
            int useless;
            int useless2;
            bool flag = objects[i]->intersect(ray, P_sphere, N_sphere, t_sphere, useless);
            if (flag){
                intersect = true;
                if ((t_sphere < t) and (0 < t_sphere)){
                    t = t_sphere;
                    P = P_sphere;
                    N = N_sphere;
                    index = i;
                }
            } 
        }
        return  intersect;
    }

    Vector get_color_intersection(const Ray &ray, int ray_depth){
        if (ray_depth == 0){
            return Vector(0,0,0);
        } //to terminate the recursion

        Vector P, N;
        Vector L0(0,0,0);
        int index;
        double t;
        if (intersect(ray, P, N, t, index)){
            //reflections
            if (objects[index]->mirror){
                Vector u = ray.u - 2 * dot(ray.u, N) * N; //std::max(0., dot(ray.u, N)
                Ray reflected_ray = Ray(P+0.001*N, u);
                return get_color_intersection(reflected_ray, ray_depth - 1);
            }
            // handle diffuse surfaces
            Vector R = L-P;
            int d2 = (R).norm2();
            R.normalize();
            Ray light_ray = Ray(P+0.01*N, R);
            Vector P_light_ray;
            Vector N_light_ray;
            int index_light_ray;
            double t_light_ray;
            int visibility = 1;

            if (intersect(light_ray, P_light_ray, N_light_ray, t_light_ray, index_light_ray)){
                //Shading and shadows computation 
                //an intersection exists but further than the light source
                if (t_light_ray*t_light_ray <= d2){
                    visibility = 0;
                }
            }

            double attenuation_with_distance = I/(4 * M_PI * (L - P).norm2());
            Vector material = objects[index]->albedo/M_PI;
            double solid_angle = std::max(0., dot(N, R));
            L0 = attenuation_with_distance * material * visibility * solid_angle;

            //We add indirect lighting
            Vector V = random_cos(N);
            Ray randomRay = Ray(P+0.001*N, V);
            L0 = L0 + times(objects[index]->albedo, get_color_intersection(randomRay, ray_depth-1));
        }
        return L0;
    }
};



//#############################################################################

class BoundingBox {
    public :
    Vector B_min, B_max;  

    BoundingBox(){}

    BoundingBox(Vector B_min, Vector B_max){
        this->B_max = B_max;
        this->B_min = B_min;
    }

    bool ray_bbox_intersect(const Ray&ray, double &t){
        //We follow the procedure explained in the video on moodle
        double t1_x, t1_y, t1_z;
        double t0_x, t0_y, t0_z;

        t0_x = (B_min[0] - ray.O[0])/ray.u[0];
        t0_y = (B_min[1] - ray.O[1])/ray.u[1];
        t0_z = (B_min[2] - ray.O[2])/ray.u[2];

        t1_x = (B_max[0] - ray.O[0])/ray.u[0];
        t1_y = (B_max[1] - ray.O[1])/ray.u[1];
        t1_z = (B_max[2] - ray.O[2])/ray.u[2];

        double a = std::min(t0_x, t1_x);
        double b = std::min(t0_y, t1_y);
        double c = std::min(t0_z, t1_z);

        double d = std::max(t0_x, t1_x);
        double e = std::max(t0_y, t1_y);
        double f = std::max(t0_z, t1_z);

        double minimum = min_of_three(d, e, f); 
        double maximum = max_of_three(a, b, c);

        if ((maximum < minimum) and (0 < minimum)){
            t = maximum;
            return true;
        }
        return false;
    }

    double max_of_three(double x, double y, double z){
        double max1 = std::max(x, y);
        return std::max(max1, z);
    }

    double min_of_three(double x, double y, double z){
        double min1 = std::min(x, y);
        return std::min(min1, z);
    }
};

class Node{
    public : 
    Node* left;
    Node* right;
    BoundingBox bbox;
    int starting_triangle;
    int ending_triangle;

    Node(){
        left = NULL;
        right = NULL;
    }
};


//#############################################################################

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};


class TriangleMesh : public Geometry{
public:
    ~TriangleMesh() {}

	TriangleMesh(const Vector& a, bool m) {
        albedo = a;
        mirror = m;
        root = new Node;
    };
	

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);
        scale_and_translate();
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    BoundingBox bbox;
    Node* root;
	
    bool ray_triangle_intersect(const Ray &ray, const Vector &A, const Vector &B, const Vector &C, double &t, Vector &P, Vector &N){
        Vector e1 = B - A;
        Vector e2 = C - A;
        N = cross(e1, e2);
        double det = dot(ray.u, N);

        Vector v = cross((A-ray.O), ray.u);
        double beta = dot(e2, v)/det;
        double gamma = -dot(e1, v)/det;
        double alpha = 1 - beta - gamma; 

        if ((beta < 0) ||  (beta > 1)){
            return false;
        }

        if ((gamma < 0) || (gamma > 1)){
            return false;
        }

        if (alpha < 0){
            return false;
        }

        t = dot(A - ray.O, N)/det;
        P = ray.O + ray.u*t;//A + beta * e1 + gamma * e2;
        N.normalize();
        return true;
    }

	//We build a ray mesh intersection function
	//bool intersect(const Ray& ray, Vector &P, Vector &N, double &t, int &index){
		// t = 1e9;
		// bool intersect = false;
		// for (int i=0; i < indices.size(); i++){
		// 	TriangleIndices triangle_indices = indices[i];
		// 	Vector A(vertices[triangle_indices.vtxi]);
		// 	Vector B(vertices[triangle_indices.vtxj]);
		// 	Vector C(vertices[triangle_indices.vtxk]);
		// 	double t_triangle;
        //     Vector P_triangle, N_triangle;
		// 	bool flag = ray_triangle_intersect(ray, A, B, C, t_triangle, P_triangle, N_triangle);
		// 	if (flag){
        //         intersect = true;
        //         if ((t_triangle < t) and (0 < t_triangle)){
        //             t = t_triangle;
        //             P = P_triangle;
        //             N = N_triangle;
        //             index = i;
        //         }
		// }
	//}
    //return intersect;
    //}

    void scale_and_translate() {
        for (int v=0; v<vertices.size(); v++) {
            vertices[v] = Vector(0,-10,0) + vertices[v] * 0.6;
        }
    }


    void compute_bbox(BoundingBox& bbox, int starting_triangle, int ending_triangle){
        double min_x, min_y, min_z;
        double max_x, max_y, max_z;
        min_x = MAXFLOAT;
        min_y = MAXFLOAT;
        min_z = MAXFLOAT;
        
        max_x = -MAXFLOAT;
        max_y = -MAXFLOAT;
        max_z = -MAXFLOAT;

        for (int i = starting_triangle; i < ending_triangle; i++) {
            //We are looking for the minimum coordinates and the maximum coordinates
            // of the x y and z of all the vertices to define the bounding box
            std::vector<Vector> vertices;
            Vector vertex1 = this->vertices[this->indices[i].vtxi];
            vertices.push_back(vertex1);
            Vector vertex2 = this->vertices[this->indices[i].vtxj];
            vertices.push_back(vertex2);
            Vector vertex3 = this->vertices[this->indices[i].vtxk];
            vertices.push_back(vertex3);
            for (int i=0; i < vertices.size(); i++){
                if (vertices[i][0] < min_x){
                    min_x = vertices[i][0];
                }
                if (vertices[i][0] > max_x){
                    max_x = vertices[i][0];
                }
                if (vertices[i][1] < min_y){
                    min_y = vertices[i][1];
                }
                if (vertices[i][1] > max_y){
                    max_y = vertices[i][1];
                }
                if (vertices[i][2] < min_z){
                    min_z = vertices[i][2];
                }
                if (vertices[i][2] > max_z){
                    max_z = vertices[i][2];
                }
            }
        }
            bbox.B_min = Vector(min_x, min_y, min_z);
            bbox.B_max = Vector(max_x, max_y, max_z);
    }

    int get_longest(Vector diag){
        int longest = 0;
        double max = -MAXFLOAT;
        for (int i = 0; i < 3; i++) {
            if (abs(diag[i]) > max) {
            max = abs(diag[i]);
            longest = i;
            }
        }
        return longest;
    }

    void build_tree(Node* node, int starting_triangle, int ending_triangle){
        //We follow the pseudo-code given in the lecture notes
        compute_bbox(node->bbox, starting_triangle, ending_triangle) ; //BBox from 
        //starting-triangle included to endingtriangle excluded 
        node->starting_triangle = starting_triangle;
        node->ending_triangle = ending_triangle;

        Vector diag = node->bbox.B_max - node->bbox.B_min;
        Vector middle_diag = node->bbox.B_min + diag * 0.5;

        int longest_axis = get_longest(diag);

        int pivot_index = starting_triangle;

        for ( int i=starting_triangle ; i<ending_triangle ; i++) {
            Vector vertex1 = this->vertices[this->indices[i].vtxi];
            Vector vertex2 = this->vertices[this->indices[i].vtxj]; 
            Vector vertex3 = this->vertices[this->indices[i].vtxk]; 

            Vector barycenter = 1/3 * (vertex1 + vertex2 + vertex3);
            
            // the swap below guarantees triangles whose barycenter are smaller than
            //middle_diag are before ”pivot index”
            if (barycenter[longest_axis] < middle_diag[longest_axis]){
                std::swap(indices[i], indices[pivot_index]);
                pivot_index++;
            }
        }

        //stopping criterion
        if ((pivot_index<=starting_triangle) || (pivot_index>=ending_triangle-1) || (ending_triangle - starting_triangle< 5)){
            return;
        } 

        Node* left = new Node();
        Node* right =  new Node();
        node->left = left;
        node->right = right;
        build_tree(left, starting_triangle, pivot_index);
        build_tree(right, pivot_index, ending_triangle);
        }


        bool intersect(const Ray& ray, Vector &P, Vector &N, double &t, int &index){
            bool intersect = false;
            if (!this->root->bbox.ray_bbox_intersect(ray, t)){
                return false;
            } 
            std::list<Node*> nodes_to_visit;
            nodes_to_visit.push_front(this->root);
            t = std::numeric_limits<double>:: max();
            while (!nodes_to_visit.empty()){
                Node* curNode = nodes_to_visit.back();
                nodes_to_visit.pop_back() ;
                // if there is one child, then it is not a leaf, so test the bounding box
                if (curNode->left) {
                    double inter_distance;
                    if (curNode->left->bbox.ray_bbox_intersect(ray, inter_distance)){
                            nodes_to_visit.push_back(curNode->left ) ;
                    }
                    if (curNode->right->bbox.ray_bbox_intersect(ray, inter_distance)){
                            nodes_to_visit.push_back(curNode->right);
                    }
                }
                else {
                // test all triangles between curNode−>starting triangle
                // and curNode−>ending triangle as before .
                // if an intersection is found, update best inter distance if needed
                    for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++){
                        TriangleIndices triangle_indices = indices[i];
                        Vector A(vertices[triangle_indices.vtxi]);
                        Vector B(vertices[triangle_indices.vtxj]);
                        Vector C(vertices[triangle_indices.vtxk]);
                        double t_triangle;
                        Vector P_triangle, N_triangle;
                        bool flag = ray_triangle_intersect(ray, A, B, C, t_triangle, P_triangle, N_triangle);
                        if (flag){
                            intersect = true;
                            if ((t_triangle < t) and (0 < t_triangle)){
                                t = t_triangle;
                                P = P_triangle;
                                N = N_triangle;
                                index = i;
                            }
                        }
                    }
                }
            }
            return intersect;
        }
};


//#############################################################################

int main() {
    int W = 512;
    int H = 512;
    double alpha = 60;

    //Converting alpha to radians;
    double alpha_radians = alpha*M_PI/180;
    //Defining z
    double z = -W/(2*tan(alpha_radians/2));
    //Defining the camera center as a vector Q
    Vector Q(0, 0, 55); 
    //We define the center spheres of radius 10 and albedo (0.7, 0.3, 0.1)
    Sphere S(Vector(-10,0,0), 10, Vector(0.7, 0.3, 0.1), false);
    //Sphere S2(Vector(10,0,0), 10, Vector(0.7, 0.3, 0.1), true);

    //We define the other spheres : walls, a ground and a ceiling
    Sphere S_ceiling(Vector(0,1000,0), 940, Vector(0.8, 0.9, 0.4), false);
    Sphere S_back(Vector(0,0,1000), 940, Vector(0.2, 0.3, 0.7), false);
    Sphere S_floor(Vector(0,-1000,0), 990, Vector(0.4, 0.7, 0.7), false);
    Sphere S_front(Vector(0,0,-1000), 940, Vector(0.5, 0.7, 0.1), false);
    Sphere S_right(Vector(1000,0,0), 940, Vector(0.2, 0.8, 0.6), false);
    Sphere S_left(Vector(-1000,0,0), 940, Vector(0.7, 0.2, 0.4), false);

    //We define a scene of position L, intensity I, and we add the spheres
    Scene scene;


    TriangleMesh mesh(Vector(0.7, 0.3, 0.1), false); 
    mesh.readOBJ("cat.obj"); 
    //mesh.readOBJ("test.obj");
    //we scale and translate in readOBJ
    mesh.build_tree(mesh.root, 0, mesh.indices.size());
    scene.add_an_object(mesh);

    //scene.add_an_object(S);
    //scene.add_an_object(S2);
    scene.add_an_object(S_ceiling);
    scene.add_an_object(S_back);
    scene.add_an_object(S_floor);
    scene.add_an_object(S_front);
    scene.add_an_object(S_right);
    scene.add_an_object(S_left);

    scene.L = Vector(-10, 20, 40);
    scene.I = 1e10;

    //define the number of paths and the max_length_path
    int nb_paths = 10;
    int max_length_path = 3;

    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) { //
        for (int j = 0; j < W; j++) {
            double x = j-W/2+0.5;
            double y = H/2-i-0.5;
            Vector u(x, y, z);
            //We normalize u 
            //u.normalize();
            //We define a ray of camera center Q and unit direction u. 
            Ray ray(Q, u);
            Vector color(0., 0., 0.);

            //We add antialiasing
            for (int k = 0; k<nb_paths; k++){
                //We define a random direction vector 
                double x2;
                double y2;
                Vector random_dir =  u + boxMuller(0.5);
                random_dir.normalize();
                Ray ray_2(Q, random_dir); // cast a ray from the camera center Q with rand dir
                //We set the color vector
                color = color + scene.get_color_intersection(ray_2, max_length_path);
                //color = scene.get_color_intersection(ray, max_length_path);

            }

            //We add the Gamma Correction
            image[(i * W + j) * 3 + 0] = std::min(255., pow(color[0]/nb_paths, 1./2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., pow(color[1]/nb_paths, 1./2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., pow(color[2]/nb_paths, 1./2.2));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}