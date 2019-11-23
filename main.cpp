#include <armadillo>
#include "surface.h"
#include "raytracer.h"
std::vector<arma::fvec3> get_image(){
    Material grey_material(0.1, 0.3,{0.5,0.5,0.5}, 0.0, 100.f);
    Material red_material(0.2, 0.6,{1,0,0}, 0.3, 25.f);
    Material glass_material(1,1,1.5,{0,0,0});
    Material green_glass_material(1,1,1.5,{0.2,0,0.2});
    Material blue_glass_material(1,1,1.5,{0.2, 0.2, 0.});
    Material blue_material(0.2,0.6,{0,0,1},0.3,25.f);
    Material green_material(0.2,0.6,{0,1,0},0.3,25.f);
    Group scene;
    scene.add(std::shared_ptr<Surface>(new Sphere({-1.f, -1.f, -6},1.5f,red_material,"Red Sphere")));
    scene.add(std::shared_ptr<Surface>(new Sphere({0.5f, -1.f, -4},1.f,blue_material,"Blue Sphere")));
//    scene.add(std::shared_ptr<Surface>(new Sphere({3.3f, -1.5f, -4},1.f,green_material,{0.5,0,0},"Moving Green Sphere")));
    scene.add(std::shared_ptr<Surface>(new Sphere({3.3f, -1.5f, -4},1.f,green_material,"Moving Green Sphere")));
    scene.add(std::shared_ptr<Surface>(new Ellipsoid({1.3f, -1.5f, -3.f},1.f,1.f,0.3f,blue_glass_material,"Blue Glass Sphere")));
    scene.add(std::shared_ptr<Surface>(new Sphere({-0.25f, -1.f, -2.f},0.7f,green_glass_material,"Green Glass Sphere")));
    scene.add(std::shared_ptr<Surface>(new Plane({0,1,0},2.5f,grey_material,"Grey Plane")));
    std::vector<std::shared_ptr<LightBase>> lights;
    lights.push_back(std::shared_ptr<LightBase>(new RectangleLight{arma::fvec3{0.5f,5,-4}, 3.0*arma::fvec3{1,1,1},{0,0,-0.5},{0.5,0,0}}));
    lights.push_back(std::shared_ptr<LightBase>(new RectangleLight{arma::fvec3{0.f,5,-4}, 3.0*arma::fvec3{1,1,1},{0,0,-0.5},{0.5,0,0}}));
//    lights.push_back(std::shared_ptr<LightBase>(new PointLight{arma::fvec3{0.5f,5,-4}, 3.0*arma::fvec3{1,1,1}}));
//    lights.push_back(std::shared_ptr<LightBase>(new PointLight{arma::fvec3{0.f,5,-4}, 3.0*arma::fvec3{1,1,1}}));
    float vfov=1.5f;
    auto img=render(scene,lights,vfov);
    return img;
}
int main(int argc, char *argv[])
{
    auto start1 = std::chrono::high_resolution_clock::now();
    auto img=get_image();
    write_ppm(img,img_width,img_height);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff1 = end1-start1;
    std::cout<<"time: "<<diff1.count()<<std::endl;


}
