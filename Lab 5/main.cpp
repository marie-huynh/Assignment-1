
#include "vector.cpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <algorithm>

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);
static std::normal_distribution<double> normal(0,1);

Vector random_direction() {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    return Vector(x, y, z);
}


void color_matching(double* image_double, double* image_target_double, int W, int H, double* image_result_double){
	memcpy(image_result_double, image_double, W*H*3*sizeof(double));
	//We use 100 iterations to be safe
	for (int i = 0; i < 100; i++){
		//We create a random unit vector 
		Vector v = random_direction();

		std::vector< std::pair<double, int> > PI(W*H);
		std::vector< std::pair<double, int> > PM(W*H);

		//We compute the projections and 
		//store the dot product and pixel index as a pair of values
		for (int j=0; j < W*H; j++) {
			double projection_1 = image_double[j*3+0]*v[0] + image_double[j*3+1]*v[1] + image_double[j*3+2]*v[2];
			double projection_2 = image_target_double[j*3+0]*v[0] + image_target_double[j*3+1]*v[1] + image_target_double[j*3+2]*v[2];
			PI[j] = std::pair<double, int>(projection_1, j);
			PM[j] = std::pair<double, int>(projection_2, j);
		}
		//We can sort the vectors by the projection values
		std::sort(PI.begin(), PI.end());
		std::sort(PM.begin(), PM.end());

		// Advect initial point cloud 
		for (int j=0; j < W*H; j++) {
			image_result_double[PI[j].second * 3 + 0] += (PM[j].first - PI[j].first)*v[0];
			image_result_double[PI[j].second * 3 + 1] += (PM[j].first - PI[j].first)*v[1];
			image_result_double[PI[j].second * 3 + 2] += (PM[j].first - PI[j].first)*v[2];
		}
}

int main() {

	int W, H, C;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	
	unsigned char *target_image = stbi_load("4852775794_c671f133d0_b.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);

	std::vector<double> image_double(W*H*3);
	std::vector<double> target_image_double(W*H*3);
	std::vector<double> image_result_double(W*H*3);
	
	color_matching(&image_double[0], &image_target_double[0], W, H, &image_result_double[0]);

	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			image_result[(i*W + j) * 3 + 0] = std::min(255., std::max(0., image_result_double[(i*W+j)*3+0]));
			image_result[(i*W + j) * 3 + 1] = std::min(255., std::max(0.,image_result_double[(i*W+j)*3+1]));
			image_result[(i*W + j) * 3 + 2] = std::min(255., std::max(0.,image_result_double[(i*W+j)*3+2]));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}