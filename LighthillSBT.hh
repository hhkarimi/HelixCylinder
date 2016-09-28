/**
 * \file LighthillSBT.hh
 *
 * \author Khalid Jawed
 * \date 09/25/2013
 *
 */

#ifndef LIGHTHILLSBT_HH
#define LIGHTHILLSBT_HH
#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"
#include <fstream> // file stream
#include <string> // string manipulations
#include <boost/lexical_cast.hpp>
#include <iostream>     // std::cout
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort
#include <vector>       // std::vector

namespace BASim
{
    // Implements Lighthill Slender Body Theory
    class LighthillSBT : public RodExternalForce
    {
    public:
        
        /** Constructor         
         \param[in] mu medium viscosity coefficient
         */
        LighthillSBT(Scalar & mu, Scalar & m_time, 
			Scalar & m_forceP, Scalar & m_torqueP, Scalar & m_percentResidue, Scalar & contactFlag,
			Scalar& eta_per, Scalar & eta_par, int SBT_method, 
			int SBTpolyN, int SBTpreconditioner, bool boundaryEffect, Scalar topBoundary, Scalar rCyl) :
        m_eta_per(eta_per),
        m_eta_par(eta_par)
        {
            m_mu = &mu,
            ttime = &m_time;            
            forceP = &m_forceP; //total propulsive force            
            torqueP = &m_torqueP; //total torque
            percentResidue = &m_percentResidue;
            m_contactFlag = &contactFlag;
            m_method = SBT_method;            
            polyN = SBTpolyN;
            preconditioner = SBTpreconditioner;
            beffect = boundaryEffect; //whether to include boundary at z=0
            m_top = topBoundary;
            m_rCyl = rCyl;
            if (m_top <= 0.0)
                chosenAxis = 0;
            else
                chosenAxis = 2;
            
            // load solution flow field data due to a stokeslet in a cylinder
                // store file names
            std::string filex = "../assets/FlowFieldData/ux.txt";
            std::string filey = "../assets/FlowFieldData/uy.txt";
            std::string filez = "../assets/FlowFieldData/uz.txt";
                // load files
            std::ifstream fileX (filex.c_str());
            std::ifstream fileY (filey.c_str());
            std::ifstream fileZ (filez.c_str());
            if ( (!fileX) || (!fileY) || (!fileZ) )
            {
                std::cerr << "Could not open data files";
            }
                // size of data set and mesh
            Nr = 11; Nf = 11; Nz = 17; Nb = 11; sd = 3;
            Ntot = Nr*Nf*Nz*Nb*sd;
//            rArray[0] = 1; phiArray[0] = 1; zArray[0] = 1; bArray[0] = 1;
            rArray[0] = 1.0e-04; rArray[1] = 1.0009e-01; rArray[2] = 2.0008e-01; rArray[3] = 3.0007e-01; rArray[4] = 4.0006e-01; rArray[5] = 5.0005e-01; rArray[6] = 6.0004e-01; rArray[7] = 7.00030000e-01; rArray[8] = 8.0002e-01; rArray[9] = 9.0001e-01; rArray[10] = 1.0e+00;
            phiArray[0] = 0.0; phiArray[1] = 0.62831853; phiArray[2] = 1.25663706; phiArray[3] = 1.88495559; phiArray[4] = 2.51327412; phiArray[5] = 3.14159265; phiArray[6] = 3.76991118; phiArray[7] = 4.39822972; phiArray[8] = 5.02654825; phiArray[9] = 5.65486678; phiArray[10] = 6.28318531;
            zArray[0] = -1.0; zArray[1] = -0.875; zArray[2] = -0.75; zArray[3] = -0.625; zArray[4] = -0.5; zArray[5] = -0.375; zArray[6] = -0.25; zArray[7] = -0.125; zArray[8] = 0.0; zArray[9] =  0.125; zArray[10] = 0.25; zArray[11] = 0.375; zArray[12] = 0.5; zArray[13] = 0.625; zArray[14] = 0.75; zArray[15] =  0.875; zArray[16] = 1.0;
            bArray[0] = 0.0; bArray[1] = 0.1; bArray[2] = 0.2; bArray[3] = 0.3; bArray[4] = 0.4; bArray[5] = 0.5; bArray[6] = 0.6; bArray[7] = 0.7; bArray[8] = 0.8; bArray[9] = 0.9; bArray[10] = 1.0;

                // go through files, line by line, and assign string value to double element in MatXd
            std::string line;
            while (getline(fileX,line))
            {
                std::stringstream ss(line);
                std::string ssString;
                ss >> ssString;
                double ssDouble = boost::lexical_cast<double>(ssString); //std::stod(ssString)
                uxBar.push_back(ssDouble);
            }
            while (getline(fileY,line))
            {
                std::stringstream ss(line);
                std::string ssString;
                ss >> ssString;
                double ssDouble = boost::lexical_cast<double>(ssString); //std::stod(ssString)
                uyBar.push_back(ssDouble);
            }
            while (getline(fileZ,line))
            {
                std::stringstream ss(line);
                std::string ssString;
                ss >> ssString;
                double ssDouble = boost::lexical_cast<double>(ssString); //std::stod(ssString)
                uzBar.push_back(ssDouble);
            }
        }
        
        void computeForce(ElasticRod& rod, VecXd& force)
        {
            Vec3d f;
            for (int i = 0; i < N; i++)
            {
                f = m_force.segment<3>( i * 3);
                force.segment<3>(i * 4) += f;
            }
        }
        
        void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}
        
        void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}
        
        VecXd updateForce(const ElasticRod & rod)
        {
            if (*ttime <= 0.0)
            {
                N = rod.nv();
                if (preconditioner == 2)
                    makeMatrixC_Average();
                else if(preconditioner == 3)
                    makeMatrixC_reducedAverage(); // best performance
                else
                    makeMatrixC();
            }
            
            VecXd u = VecXd::Zero(N * 3);
            VecXd positions = VecXd::Zero(N * 3);
            VecXd tangents = VecXd::Zero(N * 3);
            Vec3d T;
            
            for (int i = 0; i < N; i++)
            {
                u.segment<3>(i * 3) = rod.getVelocity(i);
                positions.segment<3>(i * 3) = rod.getVertex(i);
                if (i==0)
					T = rod.getVertex(1) - rod.getVertex(0);
                else if (i == N-1)
					T = rod.getVertex(N-1) - rod.getVertex(N-2);
				else
					T = rod.getVertex(i+1) - rod.getVertex(i-1);
                T.normalize();
                tangents.segment<3>(i * 3) = T;
			}
            // u.segment<3>(0 * 3) = rod.getVelocity(1);
            // u.segment<3>( (N-1) * 3) = rod.getVelocity(N-2);
            
            /* Populate the matrix A that connect the velocity vector U with
			 * the force vector F through the relation U = A F. */
			MatXd A(N * 3, N * 3);
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Vec3d x_u, x_f;
					/* Get r */
					x_u = positions.segment<3>(i * 3);
					x_f = positions.segment<3>(j * 3);
					Vec3d r = x_u - x_f;
            
					if ( i!=j ) // r.norm() >= 2.0 * deltaCutoff
					{
//						A.block<3, 3>(i * 3, j * 3) = Bij(x_u, x_f, beffect);
						A.block<3, 3>(i * 3, j * 3) = BijCyl(x_u, x_f, beffect);
					}
					else //if (i==j) //(abs(i - j) <= Ndelta/2.0)
					{
						Vec3d T = tangents.segment<3>(i * 3);
						Mat3d Stokeslet1 = (Mat3d::Identity() - T * T.transpose());
						Scalar len = ( j == N-1 ) ? rod.getEdgeLength(j-1) : rod.getEdgeLength(j);
						A.block<3, 3>(i * 3, j * 3) = 2.0 * Stokeslet1 / len; // original
//						A.block<3, 3>(i * 3, j * 3) = 2.0 * Stokeslet1 / len * (8.0 * M_PI * (*m_mu)); // with scaling
//						A.block<3, 3>(i * 3, j * 3) = -2.0 * Stokeslet1 / len * (8.0 * M_PI * (*m_mu)); // scaling and sign
					}
                }

            }
                        
            VecXd f;
            /*
             Straight forward way. Prone to ill-conditioned matrix.
             */
            if(m_method==1)
            {
                A.svd().solve(u, &f);
            }
            /*
             Polynomial approximation
             */
            else if(m_method==2)
            {
                MatXd AC = A * C;

                VecXd P = VecXd::Zero(polyN * 3);
                AC.svd().solve(u, &P);
                f = C * P;
            }
            /*
             Tikhonov regularization
             */
            // Tikhonov regularization: A^-1 = (A^T A + h^2 I)^{-1} A^T with ad hoc choice of h
            // reference: https://www.mat.univie.ac.at/~neum/ms/regtutorial.pdf
            else if(m_method==3)
            {
                MatXd ATA = A.transpose() * A;
                double h = 10;
                h = 1e-4;
                do {
                    h *= 10;
                    for (int i = 0; i < N; i++)
                        ATA(i, i) += h * h;
                        //f = ATA.partialPivLu().solve(A.transpose() * u);
                        ATA.lu().solve(A.transpose() * u, &f);
                    }
                while (f.norm() > 1.0e10 || (ATA * f - A.transpose() * u).norm() > 1e-6);
                // The while condition should be updated
            }
            else
            {
                std::cout << "Error: cannot read SBT method type\n";
                exit(1);
            }
//            m_force = -f * (8.0 * M_PI * (*m_mu));
            m_force = -f * ((*m_mu) * m_rCyl);
//            m_force = f * ((*m_mu) * m_rCyl);
            
            // Check residue
            *percentResidue = (u - A * f).norm() / u.norm() * 100;
            // std::cout << "Residue " << *percentResidue << std::endl;
            
            
			// Calculate the propulsive force and torque
            *forceP = 0.0;
            *torqueP = 0.0;
            Vec3d localTorque;
            for (int i = 0; i < N; i++)
            {
                Vec3d ft = m_force.segment<3>(i * 3);
                *forceP += ft(chosenAxis);
                localTorque = (rod.getVertex(i)).cross(ft);
                *torqueP += localTorque(0);
            }
            return m_force;
        }
        
        Mat3d Bij(Vec3d x_u, Vec3d x_f, int be)
        {
            Vec3d r = x_u - x_f;
            Vec3d x_fi = x_f;
            x_fi.z() = - x_f.z();
            Vec3d R = x_u - x_fi;
            Scalar h = x_f.z();
            
            Scalar rnorm = r.norm();
            Scalar Rnorm = R.norm();
            Vec3d rhat = r/rnorm;
            Vec3d Rhat = R/Rnorm;
                        
            Mat3d S1 = (Mat3d::Identity() + rhat * rhat.transpose()) / rnorm;
            if (be==0) return S1;
            Mat3d S2 = (Mat3d::Identity() + Rhat * Rhat.transpose()) / Rnorm;
            Mat3d B = S1 - S2;
            Mat3d Bimg;
            for (int  i = 1; i <= 3; i++)
                for (int  j = 1; j <= 3; j++)
                {
                    Scalar term1 = Delta(j,1) * (del(i,1)* h/cube(Rnorm) -
                                                 3.0 *h*Rhat(i-1)*Rhat(0)/cube(Rnorm) +
                                                 del(i,3) * Rhat(0)/ square(Rnorm) -
                                                 del(i,1) * Rhat(2) / square(Rnorm) +
                                                 3.0 * Rhat(i-1) * Rhat(0) * Rhat(2) / square(Rnorm));
                    Scalar term2 = Delta(j,2) * (del(i,2)* h/cube(Rnorm) -
                                                 3.0 *h*Rhat(i-1)*Rhat(1)/cube(Rnorm) +
                                                 del(i,3) * Rhat(1)/ square(Rnorm) -
                                                 del(i,2) * Rhat(2) / square(Rnorm) +
                                                 3.0 * Rhat(i-1) * Rhat(1) * Rhat(2) / square(Rnorm));
                    Scalar term3 = Delta(j,3) * (del(i,3)* h/cube(Rnorm) -
                                                 3.0 *h*Rhat(i-1)*Rhat(2)/cube(Rnorm) +
                                                 del(i,3) * Rhat(2)/ square(Rnorm) -
                                                 Rhat(i-1) / square(Rnorm) -
                                                 del(i,3) * Rhat(2) / square(Rnorm) +
                                                 3.0 * Rhat(i-1) * Rhat(2) * Rhat(2) / square(Rnorm));
                    Bimg(i-1, j-1) = 2.0* h * (term1 + term2 + term3);
                }
            B = B + Bimg;
            return B;
        }
        
        Mat3d BijCyl(Vec3d x_uBA, Vec3d x_fBA, int be)
        {
            // x_u is the position of the vertex where velocity is to be evaluated
            // x_f is the position of the stokeslet/source vertex
         
            // coordinate transformation from BASim to numerics from Liron and Shahar
            Mat3d rotM;
            rotM << 0,  0, 1,
                 0, -1, 0,
                 1,  0, 0;
            Vec3d x_u = rotM * x_uBA;
            Vec3d x_f = rotM * x_fBA;
            
            // set geometry
            Vec2d r( x_u.x()/m_rCyl , x_u.y()/m_rCyl);
            Scalar z = ( x_u.z() - x_f.z() )/m_rCyl;
            Vec2d b( x_f.x()/m_rCyl , x_f.y()/m_rCyl);
            Scalar rnorm = r.norm();
            if (rnorm < rArray[0])
            {
                Vec2d r(rArray[0],0.0);
                Scalar rnorm = r.norm();
            }
            Scalar bnorm = b.norm();
            Scalar phi = acos( b.dot(r) / (rnorm*bnorm) );
            if (phi < 0.0)
            {
                phi = phi + 2*M_PI;
            }
            if ( std::isnan(phi) ) // i.e. r is parallel to b
                {
                phi = 0.001;
                }
            
//            std::cout << "rnorm, phi, z, bnorm = " << rnorm << ", " << phi << ", " << z << ", " << bnorm << "\n";
            
            // populate induced velocity matrix if z < zArray.end()
            Mat3d B = Mat3d::Zero();
            if ( std::abs(z) < zArray[Nz-1] )
            {
                // interpolate indices
                int nr = lowestindex( rArray, Nr, rnorm);
                Scalar Dr = rArray[nr+1] - rArray[nr];
                Scalar dr = rnorm - rArray[nr];
                
                int nf = lowestindex( phiArray, Nf, phi);
                Scalar Df = phiArray[nf+1] - phiArray[nf];
                Scalar df = phi - phiArray[nf];
                
                int nz = lowestindex( zArray, Nz, z);
                Scalar Dz = zArray[nz+1] - zArray[nz];
                Scalar dz = z - zArray[nz];
                
                int nb = lowestindex( bArray, Nb, bnorm);
                Scalar Db = bArray[nb+1] - bArray[nb];
                Scalar db = bnorm - bArray[nb];
                
                // iterate over the 3 stokeslet directions for ux,uy,uz
                for (int sd = 0; sd < 3; sd++)
                    {
                    // interpolate velocities
                    int n = nr + nf*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz;
                    int npr = (nr+1) + nf*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz;
                    int npf = nr + (nf+1)*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz;
                    int npz = nr + nf*Nr + (nz+1)*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz;
                    int npb = nr + nf*Nr + nz*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz;
                    Scalar ux = uxBar[n] + (uxBar[npr]-uxBar[n])/Dr * dr + (uxBar[npf]-uxBar[n])/Df * df + (uxBar[npz]-uxBar[n])/Dz * dz + (uxBar[npb]-uxBar[n])/Db * db;
                    Scalar uy = uyBar[n] + (uyBar[npr]-uyBar[n])/Dr * dr + (uyBar[npf]-uyBar[n])/Df * df + (uyBar[npz]-uyBar[n])/Dz * dz + (uyBar[npb]-uyBar[n])/Db * db;
                    Scalar uz = uzBar[n] + (uzBar[npr]-uzBar[n])/Dr * dr + (uzBar[npf]-uzBar[n])/Df * df + (uzBar[npz]-uzBar[n])/Dz * dz + (uzBar[npb]-uzBar[n])/Db * db;
                    
                    // populate B matix: u = B*f (velocity = B * force)
                    B(0,sd) = ux;
                    B(1,sd) = uy;
                    B(2,sd) = uz;
                    
//                    std::cout << "rnorm,phi,z,bnorm,sd = " << rnorm << ", " << phi << ", " << z << ", " << bnorm << ", " sd << "\n";
//                    std::cout << "ux,uy,uz = " << ux << ", " << uy << ", " << uz << "\n";
//                    /****************** pause code ********************/
//                    char pauseWord; // do not release terminal when pausing
//                    std::cout << "enter any letter to continue. \n";
//                    std::cin >> pauseWord;
                    /**************************************************/
                    }
            }
            
            // undo rotation matrix:
            B = rotM.inverse() * B * rotM;
            return B;
        }
        
        
        int lowestindex(double array[], int arrayLength, double target)
        {
            int N = arrayLength;
            int index = -1;
            for (int i = 0; i < N; i++)
            {
                if ( target < array[i] )
                {
                    index = i-1;
                    break;
                }
            }
            if (index == -1)
            {
                std::cout << "Error: target value is not within array.\n";
                std::cout << "target = " << target << ", array[end] = " << array[N] << "\n";
                exit(1);
            }
            return index;
        }
        
        Scalar Delta(int j, int k)
        {
            if (j==k)
            {
                if (j==1 || j == 2) return 1;
                else return -1;
            }
            return 0;
        }
        Scalar del(int i, int j)
        {
            if (i==j) return 1;
            return 0;
        }
        
        void makeMatrixC()
        {
            // Make the C matrix in U = S * C * p
            //number of polynomial terms in functional form of f
            MatXd Csmall = MatXd::Zero(N, polyN);
            for (int i=0; i < N; i++)
            {
                for (int j=0; j < polyN; j++)
                {
                    Csmall(i,j) = pow( (double)i / (double)N, (double)j );
                }
            }
            C = MatXd::Zero(N*3, polyN*3);
            for (int i=0; i < N; i++)
            {
                for (int j=0; j < polyN; j++)
                {
                    C.block<3,3>(i*3, j*3) = Mat3d::Identity() * Csmall(i,j);
                }
            }
        }

        void makeMatrixC_Average()
        {
            polyN = N;
            int bandWidth = 10;
            // Make the C matrix in U = S * C * p
            //number of polynomial terms in functional form of f
            MatXd Csmall = MatXd::Zero(N, N);
            for (int i=0; i < N; i++)
            {
                for (int j=i-bandWidth; j <= i+bandWidth ; j++)
                {
                    if (j >= 0 && j < N)
                        Csmall(i,j) = (double)(bandWidth - abs(i-j)) / (double)(bandWidth);
                }
            }
            C = MatXd::Zero(N*3, N*3);
            for (int i=0; i < N; i++)
            {
                for (int j=0; j < N; j++)
                {
                    C.block<3,3>(i*3, j*3) = Mat3d::Identity() * Csmall(i,j);
                }
            }
        }
        
        void makeMatrixC_reducedAverage()
        {
            Scalar bandWidth = (double)N / (double)polyN;
            // Make the C matrix in U = S * C * p
            //number of polynomial terms in functional form of f
            MatXd Csmall = MatXd::Zero(N, polyN);
            for (int i=0; i < N; i++)
            {
                for (int j=0; j < polyN; j++)
                {
                    Scalar prefactor = (double)( N - abs(bandWidth*j - i) ) / (double) N;
                    Csmall(i,j) = prefactor;
                }
            }
            C = MatXd::Zero(N*3, polyN*3);
            for (int i=0; i < N; i++)
            {
                for (int j=0; j < polyN; j++)
                {
                    C.block<3,3>(i*3, j*3) = Mat3d::Identity() * Csmall(i,j);
                }
            }
        }
        
    protected:
        Scalar *m_mu;
        Scalar *ttime;
        Scalar *forceP;
        Scalar *torqueP;
        Scalar *m_contactFlag;
        Scalar *percentResidue;
        int m_method;
        int polyN;
        int preconditioner;
        bool beffect;
        int N;
        MatXd C;
        
        VecXd m_force;
        VecXd u;
        Scalar m_eta_per;
        Scalar m_eta_par;
        Scalar m_top;
        Scalar m_rCyl;
        int chosenAxis;
            // declare non-dimensional velocity data set
        std::vector<double> uxBar;
        std::vector<double> uyBar;
        std::vector<double> uzBar;
        int Ntot; int Nr; int Nf; int Nz; int Nb; int sd;
        double rArray[11]; double phiArray[11]; double zArray[17]; double bArray[11];
    };
    
} // namespace BASim

#endif // LIGHTHILLSBT_HH

