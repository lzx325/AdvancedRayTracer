#include "utils.h"
#include <cassert>
#include <vector>
arma::fmat44 translation_matrix(const arma::fvec3& x){
    arma::fmat44 m(arma::fill::eye);
    m(0,3)=x[0];
    m(1,3)=x[1];
    m(2,3)=x[2];
    return m;
}

arma::fmat44 scale_matrix(float alpha, const arma::fvec3& center){
    assert(alpha>0);
    arma::fmat44 m1=translation_matrix(-center);
    arma::fmat44 m2=arma::diagmat(arma::fvec4{alpha,alpha,alpha,1});
    arma::fmat44 m3=translation_matrix(center);
    return m3*m2*m1;
}

arma::fmat44 rotation_matrix_by_axis(float theta, char axis, const arma::fvec3& center){
    arma::fmat44 m1=translation_matrix(-center);
    arma::fmat44 m2(arma::fill::zeros);
    if(axis=='z'){
        m2={{cosf(theta),-sinf(theta),0,0},
           {sinf(theta),cosf(theta),0,0},
           {0,0,1,0},
           {0,0,0,1}};
    } else if(axis=='x'){
        m2={{1,0,0,0},
           {0,cosf(theta),-sinf(theta),0},
           {0,sinf(theta),cosf(theta),0},
           {0,0,0,1}};

    } else if(axis=='y'){
        m2={{cosf(theta),0,sinf(theta),0},
           {0,1,0,0},
           {-sinf(theta),0,cosf(theta),0},
           {0,0,0,1}};
    } else{
        assert(false);
    }
    arma::fmat44 m3=translation_matrix(center);
    return m3*m2*m1;
}
arma::fmat33 cross_product_matrix(const arma::fvec3& vec){
    return {{0,-vec[2],vec[1]},
            {vec[2],0,-vec[0]},
            {-vec[1],vec[0],0}};
}
arma::fmat44 rotation_matrix_by_direction(float theta, const arma::fvec3& dir, const arma::fvec3& center){
    assert(fabsf(arma::norm(dir)-1)<1e-6f);
    arma::fvec3 u={dir[0],dir[1],dir[2]};
    arma::fmat44 m1=translation_matrix(-center);
    arma::fmat44 m2(arma::fill::zeros);
    m2(arma::span(0,2),arma::span(0,2))=\
            cosf(theta)*arma::fmat(3,3,arma::fill::eye)+\
            sinf(theta)*cross_product_matrix(u)+\
            (1-cosf(theta))*u*u.t();
    m2(3,3)=1;
    arma::fmat44 m3=translation_matrix(center);
    return m3*m2*m1;
}

arma::fvec3 reflect(const arma::fvec3& I, const arma::fvec3 &N)
{
    assert(fabsf(arma::norm(I)-1)<1e-6);

    if(! (fabsf(arma::norm(N)-1)<1e-6)){
        std::cout<<fabsf(arma::norm(N)-1)<<std::endl;
    }
    assert(fabsf(arma::norm(N)-1)<1e-6);
    if(! (arma::dot(I,N)<=0)){
        std::cout<<N<<std::endl;
        std::cout<<arma::dot(I,N)<<std::endl;
    }
    assert(arma::dot(I,N)<=1e-4);
    arma::fvec3 R=I - N * 2.f * arma::dot(I,N);
    return R/arma::norm(R);
}

bool refract(const arma::fvec3 &I, const arma::fvec3 &N, const float eta_t,arma::fvec3& t)
{ // Snell's law
    assert(fabsf(arma::norm(I)-1)<1e-6);
    assert(fabsf(arma::norm(N)-1)<1e-6);
    if(! (arma::dot(I,N)<=0)){
        std::cout<<N<<std::endl;
        std::cout<<arma::dot(I,N)<<std::endl;
    }
    assert(arma::dot(I,N)<=1e-4);
    float cosi = -std::max(-1.f, arma::dot(I,N));
    float eta = 1 / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    if(k<0){
        // total internal reflection
        return false;
    } else {
        t=I * eta + N * (eta * cosi - sqrtf(k));
        t/=arma::norm(t);
        return true;
    }
}

bool write_ppm(std::vector<arma::fvec3> img, size_t width, size_t height){
    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm", std::ios::binary);
    if(ofs){
        ofs << "P6\n"
            << width << " " << height << "\n255\n";
        for (size_t i = 0; i < height * width; ++i)
        {
            arma::fvec3 &c = img[i];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max > 1)
                c = c * (1.f / max);
            for (size_t j = 0; j < 3; j++)
            {
                ofs << char(255 * std::max(0.f, std::min(1.f, img[i][j])));
            }
        }
        ofs.close();
        return true;
    } else{
        return false;
    }

}


