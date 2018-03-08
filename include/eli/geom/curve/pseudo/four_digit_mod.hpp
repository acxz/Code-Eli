/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_curve_pseudo_four_digit_mod_hpp
#define eli_geom_curve_pseudo_four_digit_mod_hpp

#include <string>    // std::string
#include <sstream>   // std::ostringstream, std::istringstream
#include <iomanip>   // std::setw
#include <algorithm> // std::transform

#include "eli/code_eli.hpp"

namespace eli
{
  namespace geom
  {
    namespace curve
    {
      namespace pseudo
      {
        template<typename data__>
        class four_digit_mod
        {
          public:
            typedef data__ data_type;
            typedef Eigen::Matrix<data_type, 1, 2> point_type;
            typedef Eigen::Matrix<data_type, 4, 1> coefficient_type;
            typedef typename point_type::Index index_type;

          public:
            four_digit_mod() : thickness(00), camber(0), camber_loc(0), thickness_loc(5), sharp_te(false)
            {
              recalc_params();
              recalc_coefficients();
            }

            four_digit_mod(const four_digit_mod<data_type> &fs)
             : thickness(fs.thickness), camber(fs.camber), camber_loc(fs.camber_loc), thickness_loc(fs.thickness_loc),
               sharp_te(fs.sharp_te)
            {
              recalc_params();
              recalc_coefficients();
            }

            coefficient_type get_thickness_coefficients_front() const {return a;}

            data_type get_t0() const {return static_cast<data_type>(-1);}
            data_type get_tmax() const {return static_cast<data_type>(1);}

            void set_sharp_trailing_edge(bool fl)
            {
              sharp_te=fl;
              recalc_coefficients();
            }
            bool sharp_trailing_edge() const {return sharp_te;}

            // Valid values of thickness are greater than 0 and less than 100
            bool set_thickness(const data_type &t)
            {
              if ((t>=0) && (t<=100))
              {
                thickness=t;
                recalc_params();
                return true;
              }
              return false;
            }
            data_type get_thickness() const {return thickness;}

            // valid camber values are greater than or equal to zero and less or equal to 9
            // valid camber location values are zero or greater than or equal to 1 and less than or equal to 9
            //
            // note that if either is zero then they both should be zero
            bool set_camber(const data_type &cam, const data_type &cam_loc)
            {
              if ((cam == 0) || (cam_loc == 0))
              {
                camber = 0;
                camber_loc = 0;
                recalc_params();
                return true;
              }

              if ((cam<=0) || (cam>=9))
              {
                return false;
              }
              if ((cam_loc<1) || (cam_loc>9))
              {
                return false;
              }

              camber=cam;
              camber_loc=cam_loc;
              recalc_params();
              return true;
            }
            data_type get_maximum_camber() const {return camber;}
            data_type get_maximum_camber_location() const {return camber_loc;}

            bool set_name(const std::string &name)
            {
              std::istringstream istr(name);
              std::string buffer;

              // get first character and eat preceding white space
              istr >> std::setw(4) >> buffer;
              if ((buffer == "NACA") || (buffer == "naca") || (buffer=="Naca"))
              {
                unsigned short mm, pp;
                char buf[3] = "\n\n";

                istr >> std::setw(4) >> buffer;
                buf[0]=buffer[0]; mm=std::atoi(buf);
                buf[0]=buffer[1]; pp=std::atoi(buf);

                if (set_camber(mm, pp))
                {
                  buf[0]=buffer[2];
                  buf[1]=buffer[3];
                  if (set_thickness(std::atoi(buf)))
                    return true;
                }
              }

              return false;
            }

            std::string get_name() const
            {
              std::ostringstream str;

              str.fill('0');
              str << "NACA " << std::setw(1) << camber << std::setw(1) << camber_loc << std::setw(2) << thickness;
              return str.str();
            }

            point_type f(const data_type &xi) const
            {
              point_type x, xp, xpp;

              evaluate(x, xp, xpp, xi);

              return x;
            }

            point_type fp(const data_type &xi) const
            {
              point_type x, xp, xpp;

              evaluate(x, xp, xpp, xi);

              return xp;
            }

            point_type fpp(const data_type &xi) const
            {
              point_type x, xp, xpp;

              evaluate(x, xp, xpp, xi);

              return xpp;
            }

            void evaluate(point_type &x, point_type &xp, point_type &xpp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=-1) && (xi<=1));

              data_type xc, yc, ycp, ycpp, ycppp, delta, deltap, deltapp;
              const data_type one(1), two(2), three(3);
              index_type surf_sign;

              // calculate the lower surface
              if (xi<0)
              {
                xc=-xi;
                surf_sign=-1;
              }
              // calculate the upper surface
              else
              {
                xc=xi;
                surf_sign=1;
              }

              // calculate the supporting quantities needed
              calc_camber(yc, ycp, ycpp, ycppp, xc);
              calc_thickness(delta, deltap, deltapp, xc);

              data_type tmp1, tmp13, cos_theta, sin_theta, curv, curvp;

              tmp1=std::sqrt(one+ycp*ycp);
              tmp13=tmp1*tmp1*tmp1;
              cos_theta=one/tmp1;
              sin_theta=ycp/tmp1;
              curv=ycpp/tmp13;
              curvp=ycppp/tmp13-three*ycp*tmp1*curv*curv;

              // calculate the info
              x(0)=surf_sign*(xi-delta*sin_theta);
              x(1)=yc+surf_sign*delta*cos_theta;
              xp(0)=surf_sign-deltap*sin_theta-delta*curv;
              xp(1)=ycp*(surf_sign-delta*curv)+deltap*cos_theta;
              xpp(0)=-surf_sign*(deltapp*sin_theta+two*deltap*curv+delta*curvp);
              xpp(1)=ycpp*(one-surf_sign*delta*curv)+surf_sign*(deltapp*cos_theta-ycp*(two*deltap*curv+delta*curvp));
            }

            point_type tangent(const data_type &xi) const
            {
              point_type tgt(fp(xi));

              tgt.normalize();
              return tgt;
            }

          protected:
            void recalc_params()
            {
              data_type ten(10), one_hundred(100);

              p=camber/one_hundred;
              m=camber_loc/ten;
              mt=thickness_loc/ten;
              t=thickness/one_hundred;
            }

            void recalc_coefficients()
            {
              typedef Eigen::Matrix<data_type, 4, 4> coefficient_matrix_type;

              coefficient_matrix_type coef_mat_A, coef_mat_B;
              coefficient_type rhs;
              coefficient_type orig_a, orig_d;

              // set the specified coefficients
              orig_a << static_cast<data_type>(0.989665),
                        static_cast<data_type>(-0.239250),
                        static_cast<data_type>(-0.041000),
                        static_cast<data_type>(-0.559400);

              // set the specified coefficients
              orig_d << static_cast<data_type>(0.010000),
                      static_cast<data_type>(2.325000),
                      static_cast<data_type>(-3.420000),
                      static_cast<data_type>(1.460000);


              // if blunt trailing edge, then use specified coefficients
              if (!sharp_trailing_edge())
              {
                a=orig_a;
                return;
              }

              // if want sharp trailing edge then find "actual" constraints that
              // are placed on thickness distribution

              // calculate the constraint coefficients
              coef_mat_A(0,0) = static_cast<data_type>(1);                                             // (1) Set LE radius
              calc_four_digit_mod_front_args(coef_mat_A, 1, static_cast<data_type>(mt));      // (2) thickness at x/c=mt
              calc_four_digit_mod_front_der_args(coef_mat_A, 2, static_cast<data_type>(mt));  // (3) slope at x/c=mt
              calc_four_digit_mod_front_der2_args(coef_mat_A, 3, static_cast<data_type>(mt)); // (4) radius at x/c=mt
              // calculate the corresponding constraints for the blunt trailing edge
              rhs=coef_mat_A*orig_a;

              printf("\nRHS\n");
              for(int i = 0; i < rhs.size(); i++ )
              {
                printf( "%d %g\n", i, rhs[i] );
              }

              printf("Rle %g\n", rhs[0] * rhs[0] / 2.0 );


              data_type Rref = (1.0-mt)*(1.0-mt)/(2.0*0.1*orig_d(1)*(1.0-mt)-0.588);

              printf("Rmt %g Rref %g\n", 1.0 / rhs[3], Rref );


//              0 0.989665
//              1 0.451493
//              2 0.488548
//              3 -2.59465

//              RHS
//              i Actual              Theory
//              0 0.0021              0.002   TE thickness
//              1 -0.23385           -0.234   TE slope
//              2 0.0780462           0.078   Nose shape
//              3 0.100029            0.1     Max thickness
//              4 -0.000129621        0.0     Max thickness slope

              // correct the trailing edge thickness constraint to zero while leaving the rest un changed
//              rhs(0)=static_cast<data_type>(0);
//              a=coef_mat_A.lu().solve(rhs);
            }

            template<typename Derived1>
            static void calc_four_digit_mod_front_args(Eigen::MatrixBase<Derived1> &A, const typename Derived1::Index &i, const data_type &xi)
            {
              data_type xi2(xi*xi), xi3(xi2*xi);

              A(i,0)=std::sqrt(xi);
              A(i,1)=xi;
              A(i,2)=xi2;
              A(i,3)=xi3;
            }

            template<typename Derived1>
            static void calc_four_digit_mod_rear_args(Eigen::MatrixBase<Derived1> &B, const typename Derived1::Index &i, const data_type &xi)
            {
              const data_type omx( static_cast<data_type>(1) - xi );
              const data_type omx2(omx*omx), omx3(omx*omx2);

              B(i,0)=static_cast<data_type>(1);
              B(i,1)=omx;
              B(i,2)=omx2;
              B(i,3)=omx3;
            }

            template<typename Derived1>
            static void calc_four_digit_mod_front_der_args(Eigen::MatrixBase<Derived1> &A, const typename Derived1::Index &i, const data_type &xi)
            {
              data_type xi2(xi*xi);

              A(i,0)=static_cast<data_type>(0.5)/std::sqrt(xi);
              A(i,1)=static_cast<data_type>(1);
              A(i,2)=static_cast<data_type>(2)*xi;
              A(i,3)=static_cast<data_type>(3)*xi2;
            }

            template<typename Derived1>
            static void calc_four_digit_mod_rear_der_args(Eigen::MatrixBase<Derived1> &B, const typename Derived1::Index &i, const data_type &xi)
            {
              const data_type omx( static_cast<data_type>(1) - xi );
              const data_type omx2(omx*omx), omx3(omx*omx2);

              B(i,0)=static_cast<data_type>(0);
              B(i,1)=static_cast<data_type>(-1);
              B(i,2)=static_cast<data_type>(-2)*omx;
              B(i,3)=static_cast<data_type>(-3)*omx2;
            }

            template<typename Derived1>
            static void calc_four_digit_mod_front_der2_args(Eigen::MatrixBase<Derived1> &A, const typename Derived1::Index &i, const data_type &xi)
            {
              A(i,0)=static_cast<data_type>(-1) / ( std::sqrt(xi) * xi * static_cast<data_type>(4) );
              A(i,1)=static_cast<data_type>(0);
              A(i,2)=static_cast<data_type>(2);
              A(i,3)=static_cast<data_type>(6)*xi;
            }


            void calc_camber(data_type &y, data_type &yp, data_type &ypp, data_type &yppp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              data_type zero(0), one(1), two(2);

              // short circuit if no camber
              if (camber==0)
              {
                y=zero;
                yp=zero;
                ypp=zero;
                yppp=zero;

                return;
              }

              if (xi<=m)
              {
                data_type pm2(p/(m*m));

                y=pm2*(xi*(two*m-xi));
                yp=two*pm2*(m-xi);
                ypp=-two*pm2;
                yppp=zero;
              }
              else
              {
                data_type p1m2(p/(one+m*(m-two)));

                y=p1m2*(one-two*m+xi*(two*m-xi));
                yp=two*p1m2*(m-xi);
                ypp=-two*p1m2;
                yppp=zero;
              }
            }

            void calc_thickness(data_type &y, data_type &yp, data_type &ypp, const data_type &xi) const
            {
              // check to make sure given valid parametric value
              assert((xi>=0) && (xi<=1));

              const data_type zero(0), one(1), two(2), three(3), four(4), six(6), half(one/two), quarter(one/four);
              const data_type trat(t/static_cast<data_type>(0.20));

              // short circuit for no thickness
              if (thickness==0)
              {
                y=zero;
                yp=zero;
                ypp=zero;
                return;
              }

              if (xi==0)
              {
                y=zero;
                yp=one/std::numeric_limits<data_type>::epsilon();
                ypp=yp;
                return;
              }
              else if ((xi==1) && sharp_trailing_edge())
              {
                y=zero;
                yp=trat*(a.sum()-half*a(0));
                ypp=trat*(-quarter*a(0)+two*a(2)+six*a(3));
                return;
              }

              if (xi<mt)
              {
                const data_type xi2(xi*xi), xi3(xi*xi2), sqrtxi(std::sqrt(xi));

                y=trat*(a(0)*sqrtxi+a(1)*xi+a(2)*xi2+a(3)*xi3);
                yp=trat*(half*a(0)/sqrtxi+a(1)+two*a(2)*xi+three*a(3)*xi2);
                ypp=trat*(-quarter*a(0)/sqrtxi/xi+two*a(2)+six*a(3)*xi);
              }
              else
              {
                const data_type omx( one - xi );
                const data_type omx2(omx*omx), omx3(omx*omx2);

                y=trat*( d(0)+d(1)*omx+d(2)*omx2+d(3)*omx3 );
                yp=trat*( -d(1) - d(2)*two*omx - d(3)*three*omx2 );
                ypp=trat*( d(2)*two + d(3)*six*omx );
              }

            }

          private:
            data_type thickness;    // thickness index (integer [00,99] with practical limit of [00, 30].
                                    // Index is interpreted as 100 times the maximum thickness.
            data_type camber;       // max camber index (integer [0,9]). Index will be
                                    // interpreted as 100 times the maximum camber.
            data_type camber_loc;   // chord-wise location of maximum camber index (integer [0,9]).
                                    // Index is interpreted as 10 times the maximum camber location.
            data_type thickness_loc; // chord-wise location of maximum thickness index (integer [0,9]).
                                    // Index is interpreted as 10 times the maximum camber location.

            data_type m;  // location of maximum camber
            data_type mt;  // location of maximum thickness
            data_type p;  // maximum camber
            data_type t;  // maximum thickness
            bool sharp_te; // flag to indicate if the trailing edge should be sharp
            coefficient_type a; // coefficients for thickness distribution
            coefficient_type d; // coefficients for thickness distribution
        };
      }
    }
  }
}
#endif
