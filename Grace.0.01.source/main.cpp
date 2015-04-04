#include <cmath>
#include "amp_fft.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "amp.h"
#include "amp_short_vectors.h"
#include "amp_math.h"
#include "script.h"
#include "Configuration.h"
#include "display.h"
#include "util.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

using namespace Concurrency;
using namespace Concurrency::fast_math;

extern ExternH* myExtH;

#define Index(x,y,z) ((x)+ (y) * nz + (z) * ny * nz)
#define pIndex(x,y,z) ((x)+(y) * nz/2 + (z) * ny * nz/4)
#define INDEX(X,Y,Z,dimX,dimY,dimZ) (Z)*(dimX)*(dimY)+(Y)*(dimX)+X 
#define CPU
#define GPU
#define CONSOLE 1
#define FFTW 1
#define CPPAMP 1
#define TEST 0
#define TESTEXCH 0
#define NAKATANI 1
#define WRITEOUT 1
#define DEMAG 1
#define EXCH 1
#define TWODIM 0 // two-dimensional exchange field; ignore contributions from z- and z+ cells
#define ANISOTROPY 1
#define EXTERN 1
#define FINITEDAMPING 1
#define READSCRIPT 1
#define POSTDISPLAY 0
#define READK 0

#define gamma (0.221f) // single precision precision mandated for accelerator compatiability
//#define alpha (.02f)
//#define Ms  (800.f)
//#define PI (3.14159265f)
#define mu_0 (1.256637f) // = 4 * PI / 1

#define WRITE_FAILURE -1

void default_properties() {
    accelerator default_acc;
	std::wcout << default_acc.description << "\n"; // why wcout here?
    std::wcout << default_acc.device_path << "\n";
    std::wcout << default_acc.dedicated_memory << "\n";
   // std::wcout << (default_acc.supports_cpu_shared_memory ? "CPU shared memory: true" : "CPU shared memory: false") << "\n";
    std::wcout << (default_acc.supports_double_precision ? 
        "double precision: true" : "double precision: false") << "\n";
    std::wcout << (default_acc.supports_limited_double_precision ? 
        "limited double precision: true" : "limited double precision: false") << "\n";
}

void dipoleK(float *Kxx, float *Kxy, float *Kxz, float *Kyy, float *Kyz, float *Kzz, int sizex, int sizey, int sizez) {
	int nx = sizex/2; // nx, ny, nz is the physical size
	int ny = sizey/2; 
	int nz = sizez/2; 
	float r;
	float dx = 1.f;
	float dy = 1.f;
	float dz = 1.f; 
	float prefactor = 1.f / 4 / PI;

	for (int K = - nz + 1; K < nz; K++) {
	for (int J = - ny + 1; J < ny; J++) {
	for (int I = - nx + 1; I < nx; I++) {   // calculate demag tensor
		if (I == 0 && J == 0 && K == 0) continue;
		int N = I + nx - 1; //shift the indices of demag tensors to 0, 0, 0
		int M = J + ny - 1;
		int L = K + nz - 1;
		r = sqrt (I*I*dz*dz + J*J*dy*dy + K*K*dx*dx);
		float r2 = r*r;
		float r3 = r*r*r;
		float r5 = r2*r3;
		Kxx[L * sizex * sizey + M * sizex + N] = -1.f / r3 * (1.f - 3.f * I * I * dx * dx / r2) * prefactor;
		Kyy[L * sizex * sizey + M * sizex + N] = -1.f / r3 * (1.f - 3.f * J * J * dy * dy / r2) * prefactor;
		Kzz[L * sizex * sizey + M * sizex + N] = -1.f / r3 * (1.f - 3.f * K * K * dz * dz / r2) * prefactor;
		Kxy[L * sizex * sizey + M * sizex + N] = +3.f * I * dx * J * dy / r5 * prefactor;
		Kyz[L * sizex * sizey + M * sizex + N] = +3.f * J * dy * K * dz / r5 * prefactor;
		Kxz[L * sizex * sizey + M * sizex + N] = +3.f * I * dx * K * dz / r5 * prefactor;
	}
	}
	}
	return;
}

int powi (int x, int p) {
  int i = 1;
  for (int j = 1; j <= p; j++)  i *= x;
  return i;
}

void NakataniK (float *Kxx, float *Kxy, float *Kxz, float *Kyy, float *Kyz, float *Kzz, int nx, int ny, int nz) {
	float r;
	float dx = 1.f;
	float dy = 1.f;
	float dz = 1.f;

	int sizex = nx/2;
	int sizey = ny/2;
	int sizez = nz/2;

	float prefactor = 1.f / 4 / PI;

	for (int K = - sizez + 1; K < sizez; K++) {
	for (int J = - sizey + 1; J < sizey; J++) {
	for (int I = - sizex + 1; I < sizex; I++) {
		if (I == 0 && J == 0 && K == 0) continue;
		int L = K + sizez - 1; //shift the indices of demag tensors to 0, 0, 0
		int M = J + sizey - 1;
		int N = I + sizex - 1;
		for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			r = sqrt( (I + i - 0.5f) * dx * (I + i - 0.5f) * dx
					+ (J + j - 0.5f) * dy * (J + j - 0.5f) * dy
					+ (K + k - 0.5f) * dz * (K + k - 0.5f) * dz
					);

			Kxx[L * nx * ny + M * nx + N] += float (powi(-1,i+j+k)) * atan (
				  
					(K + k - 0.5f) * (J + j - 0.5f) * dz * dy 
					/ (r * (I + i - 0.5f) * dx)
				 
				);
			Kxy[L * nx * ny + M * nx + N] += float (powi(-1,i+j+k)) * log (
				
					(K + k - 0.5f) * dz + r
				
				);
			Kxz[L * nx * ny + M * nx + N] += float (powi(-1,i+j+k)) * log ( // Kxz = Kzx
				
					(J + j - 0.5f) * dy + r
				
				);
			Kyy[L * nx * ny + M * nx + N] += float (powi(-1,i+j+k)) * atan (
				//abs( 
					(I + i - 0.5f) * (K + k - 0.5f) * dz * dx 
					/ (r * (J + j - 0.5f) * dy)
				//)
				);
			Kyz[L * nx * ny + M * nx + N] += float (powi(-1,i+j+k)) * log (
				abs(
					(I + i - 0.5f) * dx + r
				)
				);
			Kzz[L * nx * ny + M * nx + N] += float (powi(-1,i+j+k)) * atan (
				//abs( 
					(J + j - 0.5f) * (I + i - 0.5f) * dy * dx 
					/ (r * (K + k - 0.5f) * dz)
				//)
				);

			} // end for k
			} // end for j
			} // end for i
		    Kxx[L * nx * ny + M * nx + N] *= prefactor;
			Kyy[L * nx * ny + M * nx + N] *= prefactor;
			Kzz[L * nx * ny + M * nx + N] *= prefactor;
			Kxy[L * nx * ny + M * nx + N] *= -1.f * prefactor;
			Kxz[L * nx * ny + M * nx + N] *= -1.f * prefactor;
			Kyz[L * nx * ny + M * nx + N] *= -1.f * prefactor;
#if 0
			std::cout << N << " "
				<< M << " "
				<< L << " " 
				<< Kxx [L * nx * ny + M * nx + N] << " "
				<< Kxy [L * nx * ny + M * nx + N] << " "
				<< Kxz [L * nx * ny + M * nx + N] << " "
				<< Kyy [L * nx * ny + M * nx + N] << " "
				<< Kyz [L * nx * ny + M * nx + N] << " "
				<< Kzz [L * nx * ny + M * nx + N] << std::endl;
#endif
				
	} // end for K
	} // end for J
	} // end for I
					
	return;
}

void setM (float *M, float val, int dimX, int dimY, int dimZ) { // set uniform magnetization
	std::cout << "dimX = " << dimX << " dimY = " << dimY << " dimZ = " << dimZ << std::endl;
	for (int i = 0; i < dimZ/2; i++) {
	for (int j = 0; j < dimY/2; j++) {
	for (int k = 0; k < dimX/2; k++) {
		//if (i==1 && j==1 && k==1) // for test purpose only
		M[INDEX(k,j,i,dimX,dimY,dimZ)] = val;
		//std::cout << "INDEX = " << INDEX(k,j,i,dimX,dimY,dimZ) << std::endl;
	}
	}
	}
}

#if DEMAG

template <int dims>
void calcDemag(fft<float,dims> &transform,
			   array<float, dims> &Mx_GPU,
			   array<float, dims> &My_GPU,
			   array<float, dims> &Mz_GPU,
			   array<Concurrency::graphics::float_2, dims> &Mx_fft,
			   array<Concurrency::graphics::float_2, dims> &My_fft,
			   array<Concurrency::graphics::float_2, dims> &Mz_fft,
			   array<float, dims> &Kxx_GPU,
			   array<float, dims> &Kxy_GPU,
			   array<float, dims> &Kxz_GPU,
			   array<float, dims> &Kyy_GPU,
			   array<float, dims> &Kyz_GPU,
			   array<float, dims> &Kzz_GPU,
			   array<Concurrency::graphics::float_2, dims> &Kxx_fft,
			   array<Concurrency::graphics::float_2, dims> &Kxy_fft,
			   array<Concurrency::graphics::float_2, dims> &Kxz_fft,
			   array<Concurrency::graphics::float_2, dims> &Kyy_fft,
			   array<Concurrency::graphics::float_2, dims> &Kyz_fft,
			   array<Concurrency::graphics::float_2, dims> &Kzz_fft,
			   array<Concurrency::graphics::float_2, dims> &Hx_fft,
			   array<Concurrency::graphics::float_2, dims> &Hy_fft,
			   array<Concurrency::graphics::float_2, dims> &Hz_fft,
			   array<float, dims> &Hx,
			   array<float, dims> &Hy,
			   array<float, dims> &Hz,
			   array<float, dims> &Hx_raw,
			   array<float, dims> &Hy_raw,
			   array<float, dims> &Hz_raw,
			   int nx,
			   int ny,
			   int nz
			   ) {

		transform.forward_transform(Mx_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Mx_fft));
		transform.forward_transform(My_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (My_fft));
		transform.forward_transform(Mz_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Mz_fft));
#if TEST
		std::vector<Concurrency::graphics::float_2> Mx_fft_CPU = Mx_fft;
//		My_fft_CPU = My_fft;
//		Mz_fft_CPU = Mz_fft;
		for (int z = 0; z < nx; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nz; x++) {
			std::cout << "(" << Mx_fft_CPU[Index(x,y,z)].x << " + i" << Mx_fft_CPU[Index(x,y,z)].y << " )";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif

#if 0
		std::vector <float> Mx = Mx_GPU;
		std::vector <float> My = My_GPU;
		std::vector <float> Mz = Mz_GPU;
		std::cout << "In calcDemag, M = " << std::endl;
		for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			std::cout  << " ("
				<< Mx[z * nx * ny + y * nx + x] << " "
				<< My[z * nx * ny + y * nx + x] << " "
				<< Mz[z * nx * ny + y * nx + x] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif 

#if 0
		std::vector<Concurrency::graphics::float_2> Mx_fft_CPU = Mx_fft;
		std::vector<Concurrency::graphics::float_2> My_fft_CPU = My_fft;
		std::vector<Concurrency::graphics::float_2> Mz_fft_CPU = Mz_fft;
//		My_fft_CPU = My_fft;
//		Mz_fft_CPU = Mz_fft;
		std::cout << " Mx_fft = " << std::endl;
		for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			//std::cout << "(" << Mx_fft_CPU[Index(x,y,z)].x << " + i " << Mx_fft_CPU[Index(x,y,z)].y << " )";
			std::cout << "(" << Mx_fft_CPU[z * nx * ny + y * nx + x].x << " + i " 
				<< Mx_fft_CPU[z * nx * ny + y * nx + x].y << " )";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}

		std::cout << " My_fft = " << std::endl;
		for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			//std::cout << "(" << Mx_fft_CPU[Index(x,y,z)].x << " + i " << Mx_fft_CPU[Index(x,y,z)].y << " )";
			std::cout << "(" << My_fft_CPU[z * nx * ny + y * nx + x].x << " + i " 
				<< My_fft_CPU[z * nx * ny + y * nx + x].y << " )";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
		
		std::cout << " Mz_fft = " << std::endl;
		for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			//std::cout << "(" << Mx_fft_CPU[Index(x,y,z)].x << " + i " << Mx_fft_CPU[Index(x,y,z)].y << " )";
			std::cout << "(" << Mz_fft_CPU[z * nx * ny + y * nx + x].x << " + i " 
				<< Mz_fft_CPU[z * nx * ny + y * nx + x].y << " )";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif

		parallel_for_each (Hx_fft.extent, [&] (index<3> idx) restrict(amp) { // element wise product: Hx~ = Kxx~ Mx~ + Kxy~ My~ + Kxz~ Mz~ 
			/* (a + bi) * (c + di) = (a * c - b * d) + (a * d + c * b)i */
			Hx_fft[idx].x = Mx_fft[idx].x * Kxx_fft[idx].x - Mx_fft[idx].y * Kxx_fft[idx].y
						  + My_fft[idx].x * Kxy_fft[idx].x - My_fft[idx].y * Kxy_fft[idx].y
						  + Mz_fft[idx].x * Kxz_fft[idx].x - Mz_fft[idx].y * Kxz_fft[idx].y;
			Hx_fft[idx].y = Mx_fft[idx].x * Kxx_fft[idx].y + Mx_fft[idx].y * Kxx_fft[idx].x
						  + My_fft[idx].x * Kxy_fft[idx].y + My_fft[idx].y * Kxy_fft[idx].x
						  + Mz_fft[idx].x * Kxz_fft[idx].y + Mz_fft[idx].y * Kxz_fft[idx].x;
		}); 

		parallel_for_each (Hy_fft.extent, [&] (index<3> idx) restrict(amp) { 
			/* (a + bi) * (c + di) = (a * c - b * d) + (a * d + c * b)i */
			Hy_fft[idx].x = Mx_fft[idx].x * Kxy_fft[idx].x - Mx_fft[idx].y * Kxy_fft[idx].y
						  + My_fft[idx].x * Kyy_fft[idx].x - My_fft[idx].y * Kyy_fft[idx].y
						  + Mz_fft[idx].x * Kyz_fft[idx].x - Mz_fft[idx].y * Kyz_fft[idx].y;
			Hy_fft[idx].y = Mx_fft[idx].x * Kxy_fft[idx].y + Mx_fft[idx].y * Kxy_fft[idx].x
						  + My_fft[idx].x * Kyy_fft[idx].y + My_fft[idx].y * Kyy_fft[idx].x
						  + Mz_fft[idx].x * Kyz_fft[idx].y + Mz_fft[idx].y * Kyz_fft[idx].x;
		}); 

		parallel_for_each (Hz_fft.extent, [&] (index<3> idx) restrict(amp) { 
			/* (a + bi) * (c + di) = (a * c - b * d) + (a * d + c * b)i */
			Hz_fft[idx].x = Mx_fft[idx].x * Kxz_fft[idx].x - Mx_fft[idx].y * Kxz_fft[idx].y
						  + My_fft[idx].x * Kyz_fft[idx].x - My_fft[idx].y * Kyz_fft[idx].y
						  + Mz_fft[idx].x * Kzz_fft[idx].x - Mz_fft[idx].y * Kzz_fft[idx].y;
			Hz_fft[idx].y = Mx_fft[idx].x * Kxz_fft[idx].y + Mx_fft[idx].y * Kxz_fft[idx].x
						  + My_fft[idx].x * Kyz_fft[idx].y + My_fft[idx].y * Kyz_fft[idx].x
						  + Mz_fft[idx].x * Kzz_fft[idx].y + Mz_fft[idx].y * Kzz_fft[idx].x;
		}); 

		transform.inverse_transform(reinterpret_cast<array<std::complex<float>, dims>&> (Hx_fft), Hx_raw);
		transform.inverse_transform(reinterpret_cast<array<std::complex<float>, dims>&> (Hy_fft), Hy_raw);
		transform.inverse_transform(reinterpret_cast<array<std::complex<float>, dims>&> (Hz_fft), Hz_raw);

		try {

			parallel_for_each (Hx.extent, [=, &Hx, &Hx_raw] (index<3> idx) restrict(amp) { 

				/* assign Hx_raw to Hx */
				Hx [idx] = Hx_raw [idx[0]+nz/2-1][idx[1]+ny/2-1][idx[2]+nx/2-1] ;
			}); 
		} catch (std::exception &ex) {
			std::cout << ex.what() << std::endl;
		}

			parallel_for_each (Hy.extent, [=, &Hy, &Hy_raw] (index<3> idx) restrict(amp) { 
			/* assign Hx_raw to Hx */
			Hy [idx] = Hy_raw [idx[0]+nz/2-1][idx[1]+ny/2-1][idx[2]+nx/2-1] ;
			}); 

			parallel_for_each (Hz.extent, [=, &Hz, &Hz_raw] (index<3> idx) restrict(amp) { 
			/* assign Hx_raw to Hx */
			Hz [idx] = Hz_raw [idx[0]+nz/2-1][idx[1]+ny/2-1][idx[2]+nx/2-1] ;
			}); 


//#if TEST
#if 0
		std::vector<float> Hx_raw_CPU = Hx_raw;
		std::vector<float> Hy_raw_CPU = Hy_raw;
		std::vector<float> Hz_raw_CPU = Hz_raw;

		std::cout << "H_raw = " << std::endl;
		
		for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			std::cout << "(" << Hx_raw_CPU[z * nx * ny + y * nx + x] 
			<< " " << Hy_raw_CPU[z * nx * ny + y * nx + x] 
			<< " " << Hz_raw_CPU[z * nx * ny + y * nx + x] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}

		std::cout << " after added H_demag, H = " << std::endl;
		std::vector<float> Hx_vec = Hx; std::vector<float> Hy_vec  = Hy; std::vector<float> Hz_vec = Hz;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			//std::cout << Hx_raw_CPU[Index(x,y,z)] << " ";

			std::cout << "(" << Hx_vec[z * nx * ny /4 + y * nx/2 + x] << " " 
				<< Hy_vec[z * nx * ny /4 + y * nx/2 + x] << " " 
				<< Hz_vec[z * nx * ny /4 + y * nx/2 + x] << ") ";
			/* std::cout << Hx_vec[z * nx * ny + y * nx + x] << " "
				<< Hy_vec[z * nx * ny + y * nx + x] << " "
				<< Hz_vec[z * nx * ny + y * nx + x] << std::endl;
				*/
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
		std::cout << "end outputing H." << std::endl;
#endif
}

#endif // if DEMAG

#if ANISOTROPY
template <int dims>
void calcAnis(array<float, dims> &Hx, array<float, dims> &Hy, array<float, dims> &Hz,
		 array <float, dims> &Mx_GPU, array <float, dims> &My_GPU, array <float, dims> &Mz_GPU,
		 float H_anisX, float H_anisY, float H_anisZ,
		 float Ms,
		 int nx, int ny, int nz) {

        parallel_for_each (Hx.extent, [=, &Hx, &Hy, &Hz, &Mx_GPU, &My_GPU, &Mz_GPU] (index<3> idx) restrict (amp) {
            Hx[idx] += H_anisX * Mx_GPU[idx] / Ms;	Hy[idx] += H_anisY * My_GPU[idx] / Ms;	Hz[idx] += H_anisZ * Mz_GPU[idx] / Ms;
        });
        /*
		parallel_for_each (Hx.extent, [=, &Hx, &Hy, &Hz] (index<3> idx) restrict (amp) {
		Hx[idx] += H_anisX;
		Hy[idx] += H_anisY;
		Hz[idx] += H_anisZ;
		});
        */

#if TEST
		std::cout << "After No. "<< i << " iteration and added H_anis, H = " << std::endl;
		Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
		for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			std::cout << "(" << Hx_vec[pIndex(x,y,z)] << " " << Hy_vec[pIndex(x,y,z)] << " " << Hz_vec[pIndex(x,y,z)] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif

}
#endif // ANISOTROPY

#if EXTERN

template <int dims>
void addHextern (array<float, dims> &Hx, array<float, dims> &Hy, array<float, dims> &Hz,
				 float Hx_extern, float Hy_extern, float Hz_extern,
				 std::vector<float> &Hx_vec, std::vector<float> &Hy_vec, std::vector<float> &Hz_vec,  
				 int nx, int ny, int nz,
				 float Ms
				 ) {
#if READSCRIPT
#else
		float Hx_extern = Ms; float Hy_extern = 0.f; float Hz_extern = 0.f;
#endif

		parallel_for_each (Hx.extent, [=, &Hx, &Hy, &Hz] (index<3> idx) restrict (amp) {
			Hx[idx] += Hx_extern; Hy[idx] += Hy_extern; Hz[idx] += Hz_extern;
		});

#if TEST
		std::cout << " after added H_extern, H = " << std::endl;
		Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" << Hx_vec[pIndex(x,y,z)] << " " << Hy_vec[pIndex(x,y,z)] << " " << Hz_vec[pIndex(x,y,z)] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif
}

#endif // EXTERN

#if EXCH
template <int dims>
void calcExch (array <float, dims> &Hexchx, array <float, dims> &Hexchy, array <float, dims> &Hexchz, 
			   array <float, dims> &Mx_GPU, array <float, dims> &My_GPU, array <float, dims> &Mz_GPU, 
			   std::vector<float> &Hx_vec, std::vector<float> &Hy_vec, std::vector<float> &Hz_vec,
			   array<float, dims> &Hx, array<float, dims> &Hy, array<float, dims> &Hz,
			   float Ms,
			   int nx, int ny, int nz,
			   float exchConstant) {
#if TESTEXCH
			Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
			std::cout << "After No. "<< i << " iteration, H =:" << std::endl;
			for (int z = 0; z < nz/2; z++) {
			for (int y = 0; y < ny/2; y++) {
			for (int x = 0; x < nx/2; x++) {
				std::cout << "(" << Hx_vec[pIndex(x,y,z)]  << " " << Hy_vec[pIndex(x,y,z)]  << " " << Hz_vec[pIndex(x,y,z)]  << ") ";
			}
			std::cout << std::endl; 
			}
			std::cout << std::endl;
			}
#endif

			//NULL;

			//std::cout << "in exch...";

#if EXCH
		
		parallel_for_each (Hexchx.extent, [=, &Hexchx, &Hexchy, &Hexchz] (index<3> idx) restrict (amp) {
			Hexchx[idx] = 0.f;
			Hexchy[idx] = 0.f;
			Hexchz[idx] = 0.f;
		});

		parallel_for_each (Hexchx.extent, [=, &Hexchx, &Hexchy, &Hexchz, &Mx_GPU, &My_GPU, &Mz_GPU] (index<3> idx) restrict (amp) { 

			if (nz/2 != 1) {
				if (idx[0] != 0 && idx[0] != nz/2 -1) {
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx[0]-1][idx[1]][idx[2]] + Mx_GPU[idx[0]+1][idx[1]][idx[2]] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx[0]-1][idx[1]][idx[2]] + My_GPU[idx[0]+1][idx[1]][idx[2]] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx[0]-1][idx[1]][idx[2]] + Mz_GPU[idx[0]+1][idx[1]][idx[2]] - 2.f * Mz_GPU[idx]
					);
				} else if (idx[0] == 0) {
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx] + Mx_GPU[idx[0]+1][idx[1]][idx[2]] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx] + My_GPU[idx[0]+1][idx[1]][idx[2]] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx] + Mz_GPU[idx[0]+1][idx[1]][idx[2]] - 2.f * Mz_GPU[idx]
					);
				} else { // idx[0] == nz/2 -1
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx[0]-1][idx[1]][idx[2]] + Mx_GPU[idx] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx[0]-1][idx[1]][idx[2]] + My_GPU[idx] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx[0]-1][idx[1]][idx[2]] + Mz_GPU[idx] - 2.f * Mz_GPU[idx]
					);
				}
			}

			if (ny/2 != 1) {
				if (idx[1] != 0 && idx[1] != ny/2 -1) {
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx[0]][idx[1]-1][idx[2]] + Mx_GPU[idx[0]][idx[1]+1][idx[2]] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx[0]][idx[1]-1][idx[2]] + My_GPU[idx[0]][idx[1]+1][idx[2]] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx[0]][idx[1]-1][idx[2]] + Mz_GPU[idx[0]][idx[1]+1][idx[2]] - 2.f * Mz_GPU[idx]
					);
				} else if (idx[1] == 0) {
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx] + Mx_GPU[idx[0]][idx[1]+1][idx[2]] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx] + My_GPU[idx[0]][idx[1]+1][idx[2]] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx] + Mz_GPU[idx[0]][idx[1]+1][idx[2]] - 2.f * Mz_GPU[idx]
					);
				} else { // idx[1] == ny/2 -1
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx[0]][idx[1]-1][idx[2]] + Mx_GPU[idx] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx[0]][idx[1]-1][idx[2]] + My_GPU[idx] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx[0]][idx[1]-1][idx[2]] + Mz_GPU[idx] - 2.f * Mz_GPU[idx]
					);
				}
			}
			
			if (nx/2 != 1) {
				if (idx[2] != 0 && idx[2] != nx/2 -1) {
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx[0]][idx[1]][idx[2]-1] + Mx_GPU[idx[0]][idx[1]][idx[2]+1] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx[0]][idx[1]][idx[2]-1] + My_GPU[idx[0]][idx[1]][idx[2]+1] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx[0]][idx[1]][idx[2]-1] + Mz_GPU[idx[0]][idx[1]][idx[2]+1] - 2.f * Mz_GPU[idx]
					);
				} else if (idx[2] == 0) {
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx] + Mx_GPU[idx[0]][idx[1]][idx[2]+1] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx] + My_GPU[idx[0]][idx[1]][idx[2]+1] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx] + Mz_GPU[idx[0]][idx[1]][idx[2]+1] - 2.f * Mz_GPU[idx]
					);
				} else { // idx[2] == nx/2 -1
					Hexchx [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mx_GPU[idx[0]][idx[1]][idx[2]-1] + Mx_GPU[idx] - 2.f * Mx_GPU[idx]
					);
					Hexchy [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						My_GPU[idx[0]][idx[1]][idx[2]-1] + My_GPU[idx] - 2.f * My_GPU[idx]
					);
					Hexchz [idx] += 2 * exchConstant / mu_0 / Ms / Ms * (
						Mz_GPU[idx[0]][idx[1]][idx[2]-1] + Mz_GPU[idx] - 2.f * Mz_GPU[idx]
					);
				}
			}
			

        }); // end parallel for each Hexchx.extent

		parallel_for_each (Hx.extent, [&] (index<3> idx) restrict (amp) {
			Hx[idx] += Hexchx[idx];
			Hy[idx] += Hexchy[idx];
			Hz[idx] += Hexchz[idx];
		});

#if TESTEXCH
		Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
		std::cout << "After one iteration and added H_exch, H =:" << std::endl;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" << Hx_vec[pIndex(x,y,z)] << " " << Hy_vec[pIndex(x,y,z)] << " " << Hz_vec[pIndex(x,y,z)] << ") ";
			// std::cout << Hx_vec[pIndex(x,y,z)]  << " ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif

#endif // EXCH
		}
#endif

template <int dims>
void updateSystem (array<Concurrency::graphics::float_3, dims> &MxH,
				   array<float, dims> &Mx_GPU, array<float, dims>&My_GPU, array<float, dims> &Mz_GPU,
				   array<float, dims> &Hx, array<float, dims> &Hy, array<float, dims> &Hz,
				   int nx, int ny, int nz,
				   float Ms,
				   array<float, dims> &deltaMx, array<float, dims> &deltaMy, array<float, dims> &deltaMz,
				   float prefactor1, float prefactor2) {

#if 0
std::cout << "After one iteration and normalization, M =:" << std::endl;
		std::vector <float>Mx_vec = Mx_GPU; 
		std::vector <float>My_vec = My_GPU; 
		std::vector <float>Mz_vec = Mz_GPU;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" 
				<< Mx_vec[z * nx * ny / 4 + y * nx / 2 + x] << " "
				<< My_vec[z * nx * ny / 4 + y * nx / 2 + x] << " "
				<< Mz_vec[z * nx * ny / 4 + y * nx / 2 + x] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
		std::cout << "output of M done." << std::endl;
#endif

parallel_for_each (MxH.extent, [&] (index<3> idx) restrict(amp) { // # of writable array in one parallel_for_each cannot exceed 8, strange rule
		/* calculate dM = M ^ H */
		/*MxH[idx].x = My_GPU[idx] * Hx_raw[idx] - Mz_GPU[idx] * Hx_raw[idx];
		MxH[idx].y = Mz_GPU[idx] * Hx_raw[idx] - Mx_GPU[idx] * Hx_raw[idx];
		MxH[idx].z = Mx_GPU[idx] * Hx_raw[idx] - My_GPU[idx] * Hx_raw[idx];*/
		MxH[idx].x = My_GPU[idx] * Hz[idx] - Mz_GPU[idx] * Hy[idx]; // is dims.e compatiable with dims.pe?
		MxH[idx].y = Mz_GPU[idx] * Hx[idx] - Mx_GPU[idx] * Hz[idx];
		MxH[idx].z = Mx_GPU[idx] * Hy[idx] - My_GPU[idx] * Hx[idx];
		});  

#if TEST
		std::vector<Concurrency::graphics::float_3> MxH_CPU = MxH;
		std::cout << " After No. " << i << " iteration, MxH = " << std::endl;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" << MxH_CPU[pIndex(x,y,z)].x << " " << MxH_CPU[pIndex(x,y,z)].y << " " << MxH_CPU[pIndex(x,y,z)].z << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif
		//parallel_for_each (MxH.extent, [&] (index<3> idx) restrict(amp) { 
		//	/* calculate dM = M ^ H */
		//	deltaMx[idx] = MxH[idx].x;
		//	deltaMy[idx] = MxH[idx].y;
		//	deltaMz[idx] = MxH[idx].z;
		//}); 

		parallel_for_each (MxH.extent, [=, &deltaMx, &deltaMy, &deltaMz, &MxH, &Mx_GPU, &My_GPU, &Mz_GPU] (index<3> idx) restrict(amp) { 
		/* calculate dM = M ^ H */
#if FINITEDAMPING
		deltaMx[idx] = prefactor1 * MxH[idx].x  + prefactor2 * (My_GPU[idx] * MxH[idx].z - Mz_GPU[idx] * MxH[idx].y);
		deltaMy[idx] = prefactor1 * MxH[idx].y  + prefactor2 * (Mz_GPU[idx] * MxH[idx].x - Mx_GPU[idx] * MxH[idx].z);
		deltaMz[idx] = prefactor1 * MxH[idx].z  + prefactor2 * (Mx_GPU[idx] * MxH[idx].y - My_GPU[idx] * MxH[idx].x);
#else
        deltaMx[idx] = prefactor2 * (My_GPU[idx] * MxH[idx].z - Mz_GPU[idx] * MxH[idx].y);
		deltaMy[idx] = prefactor2 * (Mz_GPU[idx] * MxH[idx].x - Mx_GPU[idx] * MxH[idx].z);
		deltaMz[idx] = prefactor2 * (Mx_GPU[idx] * MxH[idx].y - My_GPU[idx] * MxH[idx].x);
#endif
		}); 
#if TESTEXCH
		std::cout << "*** Before iteration, M =:" << std::endl;
		std::vector<float> Mx_vec = Mx_GPU; std::vector<float> My_vec = My_GPU; std::vector<float> Mz_vec = Mz_GPU;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" << Mx_vec[Index(x,y,z)] << " " << My_vec[Index(x,y,z)] << " " << Mz_vec[Index(x,y,z)] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif
		
		NULL;
#if TEST
		std::cout << "*** After *** No. "<< i << " iteration, deltaM =:" << std::endl;
		std::vector<float> deltaMx_vec = deltaMx; std::vector<float> deltaMy_vec = deltaMy; std::vector<float> deltaMz_vec = deltaMz;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
		//std::cout << C[INDEX(x,y,z)] << " ";
		//std::cout << Mx_f[Index(x,y,z)] << " ";
			std::cout << "(" << deltaMx_vec[pIndex(x,y,z)] << " " << deltaMy_vec[pIndex(x,y,z)] << " " << deltaMz_vec[pIndex(x,y,z)] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif

		parallel_for_each (MxH.extent, [&] (index<3> idx) restrict(amp) {
			Mx_GPU[idx] += deltaMx[idx];
			My_GPU[idx] += deltaMy[idx];
			Mz_GPU[idx] += deltaMz[idx];
		});
#if (TESTEXCH || TEST)
		std::cout << "After one iteration, added deltaM, M =:" << std::endl;
		Mx_vec = Mx_GPU; My_vec = My_GPU; Mz_vec = Mz_GPU;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" << Mx_vec[Index(x,y,z)] << " " << My_vec[Index(x,y,z)] << " " << Mz_vec[Index(x,y,z)] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
#endif
		parallel_for_each (MxH.extent, [=, &Mx_GPU, &My_GPU, &Mz_GPU] (index<3> idx) restrict(amp) {
		float prefactor = Ms/concurrency::fast_math::sqrt(Mx_GPU[idx] * Mx_GPU[idx] + My_GPU[idx] * My_GPU[idx] + Mz_GPU[idx] * Mz_GPU[idx]);
		Mx_GPU[idx] *= prefactor;
		My_GPU[idx] *= prefactor;
		Mz_GPU[idx] *= prefactor;
		});

#if 0
		std::cout << "After one iteration and normalization, M =:" << std::endl;
		Mx_vec = Mx_GPU; 
		My_vec = My_GPU; 
		Mz_vec = Mz_GPU;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" 
				<< Mx_vec[z * nx * ny / 4 + y * nx / 2 + x] << " "
				<< My_vec[z * nx * ny / 4 + y * nx / 2 + x] << " "
				<< Mz_vec[z * nx * ny / 4 + y * nx / 2 + x] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}
		std::cout << "output of M done." << std::endl;
#endif

#if (TESTEXCH || TEST)
		std::cout << "After one iteration and normalization, M =:" << std::endl;
		Mx_vec = Mx_GPU; My_vec = My_GPU; Mz_vec = Mz_GPU;
		for (int z = 0; z < nz/2; z++) {
		for (int y = 0; y < ny/2; y++) {
		for (int x = 0; x < nx/2; x++) {
			std::cout << "(" << Mx_vec[Index(x,y,z)] << " " << My_vec[Index(x,y,z)] << " " << Mz_vec[Index(x,y,z)] << ") ";
		}
		std::cout << std::endl; 
		}
		std::cout << std::endl;
		}

		std::cout << " ------------ one iteration done -------------- " << std::endl;

#endif
		
}

#if WRITEOUT
template <int dims>
void writeOut (int &writeCounter, int &writeInterval,
			   std::vector <float> &Mx_vec, std::vector <float> &My_vec, std::vector <float> &Mz_vec,
			   array<float, dims> &Mx_GPU, array<float, dims> &My_GPU, array<float, dims> &Mz_GPU,
			   std::vector <float> &Hx_vec, std::vector <float> &Hy_vec, std::vector <float> &Hz_vec,
			   array<float, dims> &Hx, array<float, dims> &Hy, array<float, dims> &Hz,
			   std::ofstream &datafile, std::ofstream &datafile2,
			   int &currentFrame,
			   int nx, int ny, int nz,
			   float Ms,
			   int sampleX, int sampleY, int sampleZ, bool alwaysWrite) {
		writeCounter++;
		if (writeCounter == writeInterval || alwaysWrite) {
			Mx_vec = Mx_GPU; My_vec = My_GPU; Mz_vec = Mz_GPU;
			Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
			datafile << "t " << currentFrame << std::endl;
			for (int z = 0; z < nz/2; z+=sampleZ) {
			for (int y = 0; y < ny/2; y+=sampleY) {
			for (int x = 0; x < nx/2; x+=sampleX) {
				try {
				//datafile << x << " " << y << " " << z << " "
				datafile << x << " " << y << " " << z << " "
						 //<< Mx_vec[Index(x,y,z)] / Ms << " " << My_vec[Index(x,y,z)] / Ms << " " << Mz_vec[Index(x,y,z)] / Ms << " "
						 //<< Hx_vec[pIndex(x,y,z)] << " " << Hy_vec[pIndex(x,y,z)] << " " << Hz_vec[pIndex(x,y,z)] 
						 << Mx_vec[z * nx * ny + y * nx + x] / Ms << " " 
						 << My_vec[z * nx * ny + y * nx + x] / Ms << " " 
						 << Mz_vec[z * nx * ny + y * nx + x] / Ms << " "
						 << Hx_vec[z * nx * ny / 4 + y * nx / 2 + x] << " " 
						 << Hy_vec[z * nx * ny / 4 + y * nx / 2 + x] << " " 
						 << Hz_vec[z * nx * ny / 4 + y * nx / 2 + x] 
						 << std::endl;
				} catch (std::exception &ex) {
					std::cout << ex.what() << std::endl;
				}
			}
			}
			}
			datafile << std::endl;

			float Mx = 0.;
			float My = 0.;
			float Mz = 0.;
			datafile2 << currentFrame+1 << " "; // so that the first frame is frame 1, not 0, for compatiability
			for (int z = 0; z < nz/2; z++) {
			for (int y = 0; y < ny/2; y++) {
			for (int x = 0; x < nx/2; x++) {
				Mx += Mx_vec[z * nx * ny + y * nx + x];
				My += My_vec[z * nx * ny + y * nx + x];
				Mz += Mz_vec[z * nx * ny + y * nx + x];
			}
			}
			}
			Mx /= nx * ny * nz / 8 * Ms;
			My /= nx * ny * nz / 8 * Ms;
			Mz /= nx * ny * nz / 8 * Ms;
			datafile2 << Mx << " " << My << " " << Mz << std::endl;

			currentFrame++;
			writeCounter = 0; // reset counter
		}
}
		
#endif

void readState(ifstream& stateFile, std::vector<float> &Mx, std::vector<float> &My,
			   std::vector<float> &Mz, int nx, int ny, int nz, float Ms) {
	float junk;
	char c;
	stateFile >> c >> junk;
	assert(c == 't');
	for (int z = 0; z < nz/2; z++) {
	for (int y = 0; y < ny/2; y++) {
	for (int x = 0; x < nx/2; x++) {
		try {
		//datafile << x << " " << y << " " << z << " "
		stateFile >> x >> y >> z 
				 >> Mx[z * nx * ny + y * nx + x] 
				 >> My[z * nx * ny + y * nx + x]
				 >> Mz[z * nx * ny + y * nx + x]
				 >> junk
				 >> junk
				 >> junk;
				 
		} catch (std::exception &ex) {
			std::cout << ex.what() << std::endl;
		}
		Mx[z * nx * ny + y * nx + x] *= Ms;
		My[z * nx * ny + y * nx + x] *= Ms;
		Mz[z * nx * ny + y * nx + x] *= Ms;
		/*
		std::cout << x << " "
			<< y << " "
			<< z << " "
			<< Mx[z * nx * ny + y * nx + x] << " "
			<< My[z * nx * ny + y * nx + x] << " "
			<< Mz[z * nx * ny + y * nx + x] << std::endl;
			*/

	}
	}
	}

	//std::cout << "last number from .fc.txt:" << junk << std::endl;
	
}

template <int dims>
void idle(fft<float,dims> & transform, 
		  array<float, dims> &Mx_GPU, array<float, dims> &My_GPU, array<float, dims> &Mz_GPU, 
			array<Concurrency::graphics::float_2, dims> &Mx_fft, 
			array<Concurrency::graphics::float_2, dims> &My_fft, 
			array<Concurrency::graphics::float_2, dims> &Mz_fft,
			array<float, dims> &Kxx_GPU, 
			array<float, dims> &Kxy_GPU, 
			array<float, dims> &Kxz_GPU, 
			array<float, dims> &Kyy_GPU, 
			array<float, dims> &Kyz_GPU, 
			array<float, dims> &Kzz_GPU,
			array<Concurrency::graphics::float_2, dims> &Kxx_fft, 
			array<Concurrency::graphics::float_2, dims> &Kxy_fft, 
			array<Concurrency::graphics::float_2, dims> &Kxz_fft, 
			array<Concurrency::graphics::float_2, dims> &Kyy_fft, 
			array<Concurrency::graphics::float_2, dims> &Kyz_fft, 
			array<Concurrency::graphics::float_2, dims> &Kzz_fft,
			array<Concurrency::graphics::float_2, dims> &Hx_fft, 
			array<Concurrency::graphics::float_2, dims> &Hy_fft, 
			array<Concurrency::graphics::float_2, dims> &Hz_fft,
			array<float, dims> &Hx, 
			array<float, dims> &Hy, 
			array<float, dims> &Hz,
			array<float, dims> &Hx_raw, array<float, dims> &Hy_raw, array<float, dims> &Hz_raw,
			std::vector<float> &Hx_vec, std::vector<float> &Hy_vec, std::vector<float> &Hz_vec,
			array <float, dims> &Hexchx, array <float, dims> &Hexchy, array <float, dims> &Hexchz, 
			float H_anisX, float H_anisY, float H_anisZ,
			float Hx_extern, float Hy_extern, float Hz_extern,
			array<float, dims> &deltaMx, array<float, dims> &deltaMy, array<float, dims> &deltaMz, 
			float prefactor1, float prefactor2,
			array<Concurrency::graphics::float_3, dims> &MxH,
			float Ms,
			int nx, int ny, int nz,
			float exchConstant,
			int writeCounter, int writeInterval, 
			std::vector<float> &Mx_vec, std::vector<float> &My_vec, std::vector<float> &Mz_vec,
			std::ofstream &datafile, int currentFrame,
			int sampleX, int sampleY, int sampleZ
		) {
#if DEMAG
	calcDemag<3> (transform, Mx_GPU, My_GPU, Mz_GPU, 
		Mx_fft, My_fft, Mz_fft,
		Kxx_GPU, Kxy_GPU, Kxz_GPU, Kyy_GPU, Kyz_GPU, Kzz_GPU,
		Kxx_fft, Kxy_fft, Kxz_fft, Kyy_fft, Kyz_fft, Kzz_fft,
		Hx_fft, Hy_fft, Hz_fft,
		Hx, Hy, Hz,
		Hx_raw, Hy_raw, Hz_raw,
		nx, ny, nz);
#endif 

#if EXCH
	calcExch<3> (Hexchx, Hexchy, Hexchz, 
		Mx_GPU, My_GPU, Mz_GPU,
		Hx_vec, Hy_vec, Hz_vec,
		Hx, Hy, Hz,
		Ms,
		nx, ny, nz,
		exchConstant);
#endif

#if ANISOTROPY
		calcAnis(Hx, Hy, Hz, Mx_GPU, My_GPU, Mz_GPU, H_anisX, H_anisY, H_anisZ, Ms, nx, ny, nz);
#endif

#if EXTERN
		

		addHextern(Hx, Hy, Hz, Hx_extern, Hy_extern, Hz_extern, Hx_vec, Hy_vec, Hz_vec, nx, ny, nz, Ms);
#endif

		updateSystem(MxH, Mx_GPU, My_GPU, Mz_GPU, Hx, Hy, Hz, nx, ny, nz, Ms, deltaMx, deltaMy, deltaMz, prefactor1, prefactor2);

#if WRITEOUT
		writeOut(writeCounter, writeInterval, Mx_vec, My_vec, Mz_vec, Mx_GPU, My_GPU, Mz_GPU, Hx_vec, Hy_vec, Hz_vec, 
			Hx, Hy, Hz, datafile, currentFrame, nx, ny, nz, Ms, sampleX, sampleY, sampleZ);
#endif
}

template <int dims>
void test_AMP(int argc, char** argv) {
	int nx, ny, nz; // dimensions used in fft, larger than physical dimensions by factors of 2
	int timesteps;

#if READSCRIPT
    Script myScript;
    myScript.readScript(argc, argv);

#if TEST
    myScript.print();
#endif

    nx = myScript.nx;
    ny = myScript.ny;
    nz = myScript.nz;
    timesteps = myScript.timesteps;
	std::cout << "total timesteps = " << timesteps << std::endl;
#else
	nx = 480;
	ny = 128;
	nz = 4;
	timesteps = 100;
#endif // if READSCRIPT

	assert (dims == 3);
	Concurrency::extent<dims> e;
	/*e[0] = nx; // e[0] is equivalent to nx? TODO: check this
	e[1] = ny;
	e[2] = nz;
	*/

	e[0] = nz;
	e[1] = ny;
	e[2] = nx;

	Concurrency::extent<dims> pe; // physical dim
/*	pe[0] = nx/2;
	pe[1] = ny/2;
	pe[2] = nz/2;
	*/
	pe[0] = nz/2;
	pe[1] = ny/2;
	pe[2] = nx/2;

	fft<float,dims> transform(e); // fft object, size given by e
	float* Mx = new float[nx*ny*nz] (); // CPU Magnetizations
	float* My = new float[nx*ny*nz] (); 
	/*float* Mz = new float[2]();
	try {*/
	float* Mz = new float[nx*ny*nz] ();
	/*}
	catch (std::exception &ex) {
		std::cout << ex.what() << std::endl;
	}*/
	//std::vector<float> Kxx(n*n*n, 0.);
	float* Kxx = new float[nx*ny*nz] (); // CPU demagnetization tensors
	float* Kxy = new float[nx*ny*nz] ();
	float* Kxz = new float[nx*ny*nz] ();
	float* Kyy = new float[nx*ny*nz] ();
	float* Kyz = new float[nx*ny*nz] ();
	float* Kzz = new float[nx*ny*nz] ();

    float Ms; // saturation magnetization
    Ms = sqrtf (myScript.initMx * myScript.initMx + myScript.initMy * myScript.initMy + myScript.initMz * myScript.initMz);

    array<float, dims> Mx_GPU(e); // GPU magnetizations
    array<float, dims> My_GPU(e);
    array<float, dims> Mz_GPU(e);
	parallel_for_each(e, [=, &Mx_GPU, &My_GPU, &Mz_GPU] (index<3> idx) restrict(amp) { // initialize magnetizations
            Mx_GPU[idx] = 0.f;
            My_GPU[idx] = 0.f;
            Mz_GPU[idx] = 0.f;
	});

#if READSCRIPT
    float initMx = myScript.initMx; // initial magnetizations
    float initMy = myScript.initMy;
    float initMz = myScript.initMz;

        try {
		parallel_for_each(pe, [=, &Mx_GPU, &My_GPU, &Mz_GPU] (index<3> idx) restrict(amp) { // initialize magnetizations
            Mx_GPU[idx] = initMx;
            My_GPU[idx] = initMy;
            Mz_GPU[idx] = initMz;
		});
			
		} catch (std::exception &ex) {
            std::cout << ex.what(); 
		}

		//std::cout << "INDEX = " << INDEX(k,j,i,dimX,dimY,dimZ) << std::endl;
#else
	for (int i = 0; i < nx/2; i++) {
	for (int j = 0; j < ny/2; j++) {
	for (int k = 0; k < nz/2; k++) {
		//if (i==1 && j==1 && k==1) // for test purpose only
        Mx[INDEX(k,j,i,nx,ny,nz)] = Ms;
	}
	}
	}
#endif

	std::cout << "Calculating demag tensor... ";
#if NAKATANI
	NakataniK(Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, nx, ny, nz); // initialize demag tensors, according to NAKATANI 1989
#else 
	dipoleK(Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, nx, ny, nz); // dipole approximation
#endif
	std::cout << "done." << std::endl;


//#if 0
#if READK
	ifstream ifs;
	ifs.open("kernels.txt");
	if (ifs.good()) cout << "read kernel from kernels.txt successful." << endl;
	while (!ifs.eof()) {
		int i, j, k;
		float xx, xy, xz, yy, yz, zz;
		ifs >> i >> j >> k >> xx >> xy >> xz >> yy >> yz >> zz;
		/*
		cout << "read in: " 
			<< i << " "
			<< j << " "
			<< k << " "
			<< xx << " "
			<< xy << " "
			<< xz << " "
			<< yy << " "
			<< yz << " "
			<< zz << " "
			<< endl;
			*/
		if (i < nx/2
			&& i >= -nx/2 +1
			&& j < ny/2
			&& j >= -ny/2 +1
			&& k < nz/2
			&& k >= -nz/2 +1
			) {
				int N = i + nx/2 -1;
				int M = j + ny/2 -1;
				int L = k + nz/2 -1;
		Kxx[Index(L, M, N)] = xx;
		Kxy[Index(L, M, N)] = xy;
		Kxz[Index(L, M, N)] = xz;
		Kyy[Index(L, M, N)] = yy;
		Kyz[Index(L, M, N)] = yz;
		Kzz[Index(L, M, N)] = zz;

		cout << "in: " 
			<< i << " "
			<< j << " "
			<< k << " ["
			<< L * ny * nx + M * nx + N
			<< "] = "
			<< xx << " "
			<< xy << " "
			<< xz <<" "
			<< yy <<" "
			<< yz <<" "
			<< zz << endl;
		}
	}
	/*
	for (int N = 0; N < nx; N++)
	for (int M = 0; M < ny; M++)
	for (int L = 0; L < nz; L++)
	{
		int i = N - nx/2 + 1;
		int j = M - ny/2 + 1;
		int k = L - nz/2 + 1;
		cout 
			<< i << " "
			<< j << " "
			<< k << " "
			<< Kxx[Index(L, M, N)] << " "
			<< Kxy[Index(L, M, N)] << " "
			<< Kxz[Index(L, M, N)] << " "
			<< Kyy[Index(L, M, N)] << " "
			<< Kyz[Index(L, M, N)] << " "
			<< Kzz[Index(L, M, N)] << endl;
	}
	*/

	ifs.close();
#endif
    
	

#if WRITEOUT
	int writeCounter = 0;

#if READSCRIPT
    int writeInterval = myScript.writeInterval; // set interval to write output data
#else
	int writeInterval = 1;
#endif 

	std::cout << "writeInterval = " << writeInterval << std::endl;

	int currentFrame = 0;
#if 1
	int sampleX = nx < 32 ? 1 : nx / 100;
	int sampleY = ny < 32 ? 1 : ny / 50;
	int sampleZ = nz < 32 ? 1 : nz / 20; // why 8 for z?
#else
	int sampleX = 1;
	int sampleY = 1;
	int sampleZ = 1;

#endif
	std::cout << "Sample X, Y, Z = " << sampleX << " "
		<< sampleY << " "
		<< sampleZ << std::endl;
	std::ofstream datafile;
	std::ofstream datafile2;
    std::string datafilename;
//    char filename[128];
 //   sprintf(filename, "%f data.txt", myScript.runID);
    datafilename = myScript.runID + ".data.txt"; // runID is specified by default.txt file
	datafile.open (datafilename);
	if (!datafile.is_open()) { std::cout << "cannot open output file. " << std::endl; exit(WRITE_FAILURE);}
#if READSCRIPT
    myScript.print(datafile);
#endif
	datafile << "dims " << e[0]/2 << " " << e[1]/2 << " " << e[2]/2 << std::endl;
	datafile << "sampling" << " " << sampleX << " " << sampleY << " " << sampleZ << std::endl;
	datafile << "frames " << timesteps/writeInterval << std::endl;

	std::string datafile2name = myScript.runID + ".mvst.txt"; // Magnetization vs. time file
	datafile2.open(datafile2name);
	if (!datafile2.is_open()) { std::cout << "cannot open output file2. " << std::endl; exit (-4);}

	
	
#endif
	//using namespace concurrency::graphics;

	
    /*
	parallel_for_each (Mx_GPU.extent, [&] (index<3> idx) restrict (amp) {
		Mx_GPU[idx] = Ms;
	});
    */
	std::vector<float> Mx_vec = Mx_GPU; std::vector<float> My_vec = My_GPU; std::vector<float> Mz_vec = Mz_GPU; // copy magnetizations from CPU to GPU

	if (myScript.readInitState) {
		ifstream statefile;
		statefile.open(myScript.runID + ".fc.txt");
		if (!statefile.good()) { 
			std::cout << "cannot open init state file." << std::endl;
			exit(-3);
		}

		std::cout << "reading initial state file... ";
		readState(statefile, Mx_vec, My_vec, Mz_vec, nx, ny, nz, Ms);
		copy(Mx_vec.begin(), Mx_vec.end(), Mx_GPU);
		copy(My_vec.begin(), My_vec.end(), My_GPU);
		copy(Mz_vec.begin(), Mz_vec.end(), Mz_GPU);
		std::cout << "done. " << std::endl;
		/*
		Mx_GPU = Mx_vec;
		My_GPU = My_vec;
		Mz_GPU = Mz_vec;
		*/
	}

#if TESTEXCH
	std::cout << "M =:" << std::endl;

	for (int z = 0; z < nx/2; z++) {
	for (int y = 0; y < ny/2; y++) {
	for (int x = 0; x < nz/2; x++) {
		std::cout << "(" << Mx_vec[Index(x,y,z)] << " " << My_vec[Index(x,y,z)] << " " << Mz_vec[Index(x,y,z)] << ") ";
	}
	std::cout << std::endl; 
	}
	std::cout << std::endl;
	}
#endif

	array<float, dims> deltaMx(pe); // MxH is the cross product of M and H
	array<float, dims> deltaMy(pe);
	array<float, dims> deltaMz(pe);
	array<Concurrency::graphics::float_3, dims> MxH (pe);
	array<float, dims> Hx_raw(e), Hy_raw(e), Hz_raw(e);
	array<float, dims> Hx (pe), Hy (pe), Hz (pe);
	std::vector<float> Hx_vec; std::vector<float> Hy_vec; std::vector<float> Hz_vec;
	std::vector <Concurrency::graphics::float_2> Kxx_fft_CPU;

#if DEMAG
	array<float, dims> Kxx_GPU(e, Kxx), Kxy_GPU(e, Kxy), Kxz_GPU(e, Kxz), Kyy_GPU(e, Kyy), Kyz_GPU(e, Kyz), Kzz_GPU(e, Kzz); // demag tensors in GPU
	array<Concurrency::graphics::float_2, dims> Mx_fft(e), My_fft(e), Mz_fft(e); // magnetizations in GPU after FFT
	array<Concurrency::graphics::float_2, dims> Kxx_fft(e), Kxy_fft(e), Kxz_fft(e), Kyy_fft(e), Kyz_fft(e), Kzz_fft(e); // demag tensors in GPU after FFT
	array<Concurrency::graphics::float_2, dims> Hx_fft(e), Hy_fft(e), Hz_fft(e); // demag fields in GPU
	
	transform.forward_transform(Kxx_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Kxx_fft));
	transform.forward_transform(Kxy_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Kxy_fft));
	transform.forward_transform(Kxz_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Kxz_fft));
	transform.forward_transform(Kyy_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Kyy_fft));
	transform.forward_transform(Kyz_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Kyz_fft));
	transform.forward_transform(Kzz_GPU, reinterpret_cast<array<std::complex<float>, dims>&> (Kzz_fft));
	
#endif

#if TEST
	Kxx_fft_CPU = Kxx_fft;
	for (int z = 0; z < nz; z++) {
	for (int y = 0; y < ny; y++) {
	for (int x = 0; x < nx; x++) {
		std::cout << "(" << Kxx_fft_CPU[Index(x,y,z)].x << " + " << Kxx_fft_CPU[Index(x,y,z)].y  << "i ) ";
	}
	std::cout << std::endl; 
	}
	std::cout << std::endl;
	}

#endif
	
#if EXCH

	array <float, dims> Hexchx(e), Hexchy(e), Hexchz(e);  
	float exchConstant;
	std::cout << "size = ";
	std::cout << nx << " x " << ny << " x " << nz << " = " << Hexchx.extent.size() << std::endl;

#if READSCRIPT
    exchConstant = myScript.exchConstant; 
    exchConstant *= 1E18f; // A in J/m converted to J'/nm; 1 J' = 1E-27 J
#else
	exchConstant = 1.3E-11f * 1E18f; // A in J/m converted to J'/nm; 1 J' = 1E-27 J
#endif

#if TEST
	std::vector <float> Hexchx_CPU (nx*ny*nz, 0.f), Hexchy_CPU (nx*ny*nz, 0.f), Hexchz_CPU (nx*ny*nz, 0.f);
#endif
	std::cout << "L_exch2 = " << sqrtf(2 * exchConstant / mu_0 / Ms / Ms ) << " nm." << std::endl;
#endif

#if ANISOTROPY
	float H_anisX;
	float H_anisY;
	float H_anisZ;
	H_anisX = myScript.Hkx;
    H_anisY = myScript.Hky;
    H_anisZ = myScript.Hkz;
    float H_k = sqrtf(H_anisX * H_anisX + H_anisY * H_anisY + H_anisZ * H_anisZ);

#if EXCH
	std::cout << "L_exch1 = " << sqrtf(2 * exchConstant / mu_0 / H_k / Ms) << " nm." << std::endl;
#endif

#if TEST
	std::cout << "H_anis = " << std::endl;
	std::cout << "(" << H_anisX << " " << H_anisY << " " << H_anisZ << ") " << std::endl;
	/*std::cout << " after added H_anis, H = " << std::endl;
	Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
	for (int z = 0; z < nz/2; z++) {
	for (int y = 0; y < ny/2; y++) {
	for (int x = 0; x < nx/2; x++) {
		std::cout << "(" << Hx_vec[pIndex(x,y,z)] << " " << Hy_vec[pIndex(x,y,z)] << " " << Hz_vec[pIndex(x,y,z)] << ") ";
	}
	std::cout << std::endl; 
	}
	std::cout << std::endl;
	}*/
#endif

#else
	float EAx = 1.f; float EAy = 0.f; float EAz = 0.f;
	float H_k = 10.f;

	float H_anisX = EAx * H_k; float H_anisY = EAy * H_k; float H_anisZ = EAz * H_k;
#endif

#if EXTERN
	float Hx_extern;
	float Hy_extern;
	float Hz_extern;

#if READSCRIPT
    Hx_extern = myScript.externFieldx;
    Hy_extern = myScript.externFieldy;
    Hz_extern = myScript.externFieldz;
#else
	Hx_extern = 2000.f; Hy_extern = 2000.f; Hz_extern = 2000.f;
#endif	

#if TEST
	std::cout << " H_extern = " << std::endl;
	std::cout << "(" << Hx_extern << " " << Hy_extern << " " << Hz_extern << ")" << std::endl;
	/*Hx_vec = Hx; Hy_vec = Hy; Hz_vec = Hz;
	for (int z = 0; z < nx/2; z++) {
	for (int y = 0; y < ny/2; y++) {
	for (int x = 0; x < nz/2; x++) {
		std::cout << "(" << Hx_vec[pIndex(x,y,z)] << " " << Hy_vec[pIndex(x,y,z)] << " " << Hz_vec[pIndex(x,y,z)] << ") ";
	}
	std::cout << std::endl; 
	}
	std::cout << std::endl;
	}*/
#endif

#endif


#if READSCRIPT
        float dt = myScript.dt;
        float alpha = myScript.alpha;
        float prefactor1 = (-0.221f) * dt;
        float prefactor2 = prefactor1 * alpha / Ms;
#else

        float alpha = 1.f;
        float dt = .001f;
        float prefactor1 = (-0.221f) * dt;
        float prefactor2 = prefactor1 * alpha / Ms; // = prefactor1 * alpha * Ms
#endif

#if CONSOLE
#else
	GLInit(argc, argv); // glut init. functions
	glutKeyboardFunc(keyboardFunc); // press esc or q to quit, r to redraw display
	glutDisplayFunc(displayFunc); // set vector positions and magnitudes
	glutIdleFunc(idle<3>); // update system parameters
	glutReshapeFunc(reShape); // handles display resizing
	glutMouseFunc(mouseButton); // handles mouse button clicking
	glutMotionFunc(mouseMove); // handles mouse moving

	// OpenGL init
	glEnable(GL_DEPTH_TEST);
#endif
	std::cout << "Calculating... " << std::endl;
	/*for (int iters = 100; iters < 100001; iters *= 2) {
	clock_t start = clock();
	for (int i = 0; i < iters; i++ ) {*/
    int updateSteps = 0;
    int updateProgressCounter = 0;
    if (timesteps >= 100) updateSteps = timesteps / 100;
	clock_t start = clock();
#if CONSOLE
	
	for (int i = 0; i < timesteps; i++ ) {

        if (updateProgressCounter == updateSteps) {
            printf("\b\b\b\b\b\b %3.0f%%",(float)i / timesteps * 100 + 1);
			//std::cout << "\b\b\b\b" << (float)i / timesteps * 100 << "%";
            updateProgressCounter = 0;
		} else
			updateProgressCounter ++;

#if DEMAG
	calcDemag<3> (transform, Mx_GPU, My_GPU, Mz_GPU, 
		Mx_fft, My_fft, Mz_fft,
		Kxx_GPU, Kxy_GPU, Kxz_GPU, Kyy_GPU, Kyz_GPU, Kzz_GPU,
		Kxx_fft, Kxy_fft, Kxz_fft, Kyy_fft, Kyz_fft, Kzz_fft,
		Hx_fft, Hy_fft, Hz_fft,
		Hx, Hy, Hz,
		Hx_raw, Hy_raw, Hz_raw,
		nx, ny, nz);
#endif 

#if EXCH
	calcExch<3> (Hexchx, Hexchy, Hexchz, 
		Mx_GPU, My_GPU, Mz_GPU,
		Hx_vec, Hy_vec, Hz_vec,
		Hx, Hy, Hz,
		Ms,
		nx, ny, nz,
		exchConstant);
#endif

#if ANISOTROPY
		calcAnis(Hx, Hy, Hz, Mx_GPU, My_GPU, Mz_GPU, H_anisX, H_anisY, H_anisZ, Ms, nx, ny, nz);
#endif

#if EXTERN
		Hx_extern = 0.;
		Hy_extern = 0.;
		Hz_extern = 0.;

		for (int F = 0; F < myExtH->numOfFields; F++) { // add all ext. fields
					  if (i >= myExtH[F].startTime && i < myExtH[F].endTime) {
						  if (i < myExtH[F].decayTime) {
							  Hx_extern += myExtH[F].externHx;
							  Hy_extern += myExtH[F].externHy;
							  Hz_extern += myExtH[F].externHz;
						  } else {
							  Hx_extern += myExtH[F].externHx * 
								  (myExtH[F].endTime - i) / (myExtH[F].endTime - myExtH[F].decayTime);
							  Hy_extern += myExtH[F].externHy * 
								  (myExtH[F].endTime - i) / (myExtH[F].endTime - myExtH[F].decayTime);
							  Hz_extern += myExtH[F].externHz * 
								  (myExtH[F].endTime - i) / (myExtH[F].endTime - myExtH[F].decayTime);
						  }
					  } else {
						  Hx_extern += 0.;
						  Hy_extern += 0.;
						  Hz_extern += 0.;
					  }
					  
					  /*if (i % 100 == 0) {
						  cout << "# ext H = (" 
							  << Hx_extern << " "
							  << Hy_extern << " "
							  << Hz_extern << ") " << endl;
					  }*/
					  
		}
		addHextern(Hx, Hy, Hz, Hx_extern, Hy_extern, Hz_extern, Hx_vec, Hy_vec, Hz_vec, nx, ny, nz, Ms);
#endif

		updateSystem(MxH, Mx_GPU, My_GPU, Mz_GPU, Hx, Hy, Hz, nx, ny, nz, Ms, deltaMx, deltaMy, deltaMz, prefactor1, prefactor2);

#if WRITEOUT
		writeOut(writeCounter, writeInterval, Mx_vec, My_vec, Mz_vec, Mx_GPU, My_GPU, Mz_GPU, Hx_vec, Hy_vec, Hz_vec, 
			Hx, Hy, Hz, datafile, datafile2, currentFrame, nx, ny, nz, Ms, sampleX, sampleY, sampleZ, false);
#endif
	} // end for int i

#else  // if CONSOLE
	glutMainLoop();

#endif
	std::cout << " \a... done. " << std::endl;
	std::cout << "time spent = " << (clock() - start ) / float (CLOCKS_PER_SEC) * 1000 << " ms " << std::endl;
	//std::cout << iters << "\t";
	//std::cout << "time spent per iter = " << (clock() - start ) / float (CLOCKS_PER_SEC) / iters *1000 << " ms " << std::endl;
	//
	//} // end while

    datafile.close();

	if (myScript.writeFinalState) {
		ofstream finalCondition;
		finalCondition.open(myScript.runID + ".fc.txt");
		if (!finalCondition.good()) std::cout << "can't open final condition file." << std::endl; 
		//writeCounter = 0; // so that the write can be successful
		//writeInterval = 1; 
		writeOut(writeCounter, writeInterval, Mx_vec, My_vec, Mz_vec, Mx_GPU, My_GPU, Mz_GPU, Hx_vec, Hy_vec, Hz_vec, 
				Hx, Hy, Hz, finalCondition, datafile2, currentFrame, nx, ny, nz, Ms, 1, 1, 1, true);
		finalCondition.close();
		std::cout << "final condition file written to " << myScript.runID << ".fc.txt" << std::endl;
	}
#if POSTDISPLAY
    std::cout << "display results? Y/N: ";
    char answer;
    std::cin >> answer;
    switch (answer) {
        case 'Y':
        case 'y':
            extern Configuration * myConfig;
            myConfig = new Configuration("example.txt");
            Display myDisplay;

            myDisplay.Run(argc, argv);
            break;

        case 'N':
        case 'n':
        default:
            break;
	} // end switch
#endif

	ofstream plotFile;
	plotFile.open (myScript.runID + ".plt");
	if (!plotFile.good()) { cout << " cannot open gnuplot script file..." << endl; }
	else {
		plotFile << "set xrange [-1:" << nx/2 + 1 << "]" << endl; // physical size is a half of fft size
		plotFile << "set yrange [-1:" << ny/2 + 1 << "]" << endl; 
		plotFile << "set zrange [-1:" << nz/2 + 1 << "]" << endl; 
		plotFile << "set xlabel\"x\"" << endl;
		plotFile << "set xlabel\"y\"" << endl;
		plotFile << "set xlabel\"z\"" << endl;
		plotFile << "set terminal wxt size 600, 600" << endl; // for best resolution compatiability
		plotFile << "set view equal xyz" << endl; // squal scale for x y z
		plotFile << "set view 0, 0, 1.2" << endl; // zoom in a bit
		plotFile << "set ticslevel 0" << endl; // place system in screen center
		plotFile << "do for [ii=0:" << timesteps / writeInterval - 1 << "] {" << endl
			<< "set label 1 sprintf('frame = %d', ii) at -1, -2 right front" << endl
			<< "splot '" << myScript.runID + ".data.txt" << "' every :::ii::ii using 1:2:3:($4*"<< sampleY * 0.5 
			<< "):($5*" << sampleY * 0.5
			<< "):($6*" << sampleY * 0.5
			<< ") with vectors title 'magnetizations'" << endl
			<< "pause 0.1" << endl
			<< "}" << endl;
	};

	plotFile.close();

    //Sleep(1000);
	std::cin.get();

}

int main (int argc, char** argv) {
	//test_fft_3d<3>();
	//test_FFTW();

	//default_properties();
	test_AMP<3>(argc, argv);
	return 0;
}