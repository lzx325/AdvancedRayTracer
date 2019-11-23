#ifndef UTILS_H
#define UTILS_H
#include <armadillo>
#include <cassert>
arma::fmat44 translation_matrix(const arma::fvec3& x);
arma::fmat44 scale_matrix(float alpha, const arma::fvec3& center={0,0,0});
arma::fmat44 rotation_matrix_by_axis(float theta, char axis='z', const arma::fvec3& center={0,0,0});
arma::fmat33 cross_product_matrix(const arma::fvec3& vec);
arma::fmat44 rotation_matrix_by_direction(float theta, const arma::fvec3& dir, const arma::fvec3& center={0,0,0});
inline arma::fvec4 to_arma_fvec4(const arma::fvec3& v1,float pad=1){
    arma::fvec4 v2={v1[0],v1[1],v1[2],pad};
    return v2;
}
inline arma::fvec3 to_arma_fvec3(const arma::fvec4& v1,bool homo=false){

    arma::fvec3 v2;
    if(homo){
        if(!(fabsf(v1[3])>1e-10f)){
            std::cout<<fabsf(v1[3])<<std::endl;
        }
        assert(fabsf(v1[3])>1e-10f);
        v2={v1[0]/v1[3],v1[1]/v1[3],v1[2]/v1[3]};
    } else{
        v2={v1[0],v1[1],v1[2]};
    }
    return v2;
}
arma::fvec3 reflect(const arma::fvec3& I, const arma::fvec3 &N);
bool refract(const arma::fvec3 &I, const arma::fvec3 &N, const float eta_t,arma::fvec3& t);
#endif // UTILS_H
