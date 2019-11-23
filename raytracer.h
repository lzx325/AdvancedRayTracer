#ifndef RAYTRACER_H
#define RAYTRACER_H
#include <iostream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>
#include <omp.h>
#include <ctime>
#include "surface.h"
arma::fvec3 raycolor(const Ray& ray,const Group& scene, const std::vector<std::shared_ptr<LightBase>>& lights,
                     const arma::fvec3& Ia, const arma::fvec3& bg, float t0, float t1, const DRParams& drparams);
std::vector<arma::fvec3> render(const Group & scene, const std::vector<std::shared_ptr<LightBase> > &lights,float tan_vfov);
bool write_ppm(std::vector<arma::fvec3> img, size_t width, size_t height);
extern size_t img_width;
extern size_t img_height;
#endif // RAYTRACER_H
