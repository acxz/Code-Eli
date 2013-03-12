/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef eli_geom_surface_piecewise_hpp
#define eli_geom_surface_piecewise_hpp

#include <vector>

#include "eli/util/tolerance.hpp"

#include "eli/geom/general/continuity.hpp"

namespace eli
{
  namespace geom
  {
    namespace surface
    {
      template<template<typename, unsigned short, typename> class surface__, typename data__, unsigned short dim__, typename tol__=eli::util::tolerance<data__> >
      class piecewise
      {
        public:
          typedef surface__<data__, dim__, tol__> surface_type;
          typedef typename surface_type::index_type index_type;
          typedef typename surface_type::point_type point_type;
          typedef typename surface_type::control_point_type control_point_type;
          typedef data__ data_type;
          typedef unsigned short dimension_type;
          typedef tol__ tolerance_type;
          enum error_code
          {
            NO_ERROR=0,
            INVALID_INDEX=1,
            INDEX_NOT_FOUND=2,
            INVALID_PARAM=50,
            INVALID_PARAM_DIFFERENCE=51,
            PATCH_NOT_CONNECTED=100,
            UNKNOWN_ERROR=999
          };

        public:
          piecewise() : u0(0), v0(0), nu(0), nv(0) {}
          piecewise(const piecewise<surface__, data_type, dim__, tol__> &p)
            : patches(p.patches), u0(p.u0), v0(p.v0), nu(p.nu), nv(p.nv) {}
          ~piecewise() {}

          bool operator==(const piecewise<surface__, data_type, dim__> &p) const
          {
            if (this==&p)
              return true;
            if (u0!=p.u0)
              return false;
            if (v0!=p.v0)
              return false;
            if (nu!=p.nu)
              return false;
            if (nv!=p.nv)
              return false;
            if (number_u_patches()!=p.number_u_patches())
              return false;
            if (number_v_patches()!=p.number_v_patches())
              return false;
            typename patch_collection_type::const_iterator scit, it;
            for (scit=patches.begin(), it=p.patches.begin(); scit!=patches.end(); ++scit, ++it)
            {
              if ((*it)!=(*scit))
                return false;
            }

            return true;
          }

          bool operator!=(const piecewise<surface__, data_type, dim__> &p) const
          {
            return !operator==(p);
          }

          static dimension_type dimension() {return dim__;}

          const data_type & get_u0() const {return u0;}
          void set_u0(const data_type &u0_in) {u0=u0_in;}

          const data_type & get_v0() const {return v0;}
          void set_v0(const data_type &v0_in) {v0=v0_in;}

          index_type number_u_patches() const {return nu;}
          index_type number_v_patches() const {return nv;}

          void get_parameter_min(data_type &umin, data_type &vmin) const
          {
            umin=u0;
            vmin=v0;
          }

          void get_parameter_max(data_type &umax, data_type &vmax) const
          {
            index_type i, j;
            typename patch_collection_type::const_iterator pcit;

            umax=u0;
            vmax=v0;

            for (i=0, pcit=patches.begin(); i<nu; ++i, ++pcit)
            {
              umax+=pcit->delta_u;
            }
            for (j=0, pcit=patches.begin(); j<nv; ++j, pcit+=nu)
            {
              vmax+=pcit->delta_v;
            }
          }

          void resize(const index_type &nu_in, const index_type &nv_in)
          {
            if ( (nu_in<=0) || (nv_in<=0) )
            {
              return;
            }

            patches.resize(nu_in*nv_in);
            nu=nu_in;
            nv=nv_in;
          }

          bool open() const
          {
#if 0
            return check_continuity(patches.rbegin()->s, patches.begin()->s, eli::geom::general::C0);
#else
            return false;
#endif
          }
          bool closed() const
          {
            return !open();
          }

          void get_bounding_box(point_type &pmin, point_type &pmax) const
          {
            typename patch_collection_type::const_iterator it=patches.begin();
            point_type pmintmp, pmaxtmp;

            // cycle through all patches to get each bounding box to compare
            it->s.get_bounding_box(pmin, pmax);
            for (++it; it!=patches.end(); ++it)
            {
              it->s.get_bounding_box(pmintmp, pmaxtmp);
              for (index_type i=0; i<dim__; ++i)
              {
                if (pmintmp(i)<pmin(i))
                {
                  pmin(i)=pmintmp(i);
                }
                if (pmaxtmp(i)>pmax(i))
                {
                  pmax(i)=pmaxtmp(i);
                }
              }
            }
          }

          void reverse()
          {
#if 0
            // reverse order of patches
            patches.reverse();

            // reverse each patch
            for (typename patch_collection_type::iterator it=patches.begin(); it!=patches.end(); ++it)
            {
              it->s.reverse();
            }

            // check if still connected
            assert(check_continuity(eli::geom::general::C0));
#else
            assert(false);
#endif
          }

          void swap_uv()
          {
            assert(false);
          }

          void clear() {patches.clear();}

          error_code get(surface_type &surf, const index_type &ui, const index_type &vi) const
          {
            data_type du, dv;
            return get(surf, du, dv, ui, vi);
          }

          error_code get(surface_type &surf, data_type &du, data_type &dv, const index_type &ui, const index_type &vi) const
          {
            if ( (ui>=number_u_patches()) || (vi>=number_v_patches()) )
              return INVALID_INDEX;

            // advance to desired index
            index_type i, j;
            typename patch_collection_type::const_iterator pcit;
            for (i=0, pcit=patches.begin(); i<ui; ++i, ++pcit) {}
            for (j=0; j<vi; ++j, pcit+=nu) {}

            surf=pcit->s;
            du=pcit->delta_u;
            dv=pcit->delta_v;

            return NO_ERROR;
          }

          error_code replace(const surface_type &surf, const index_type &ui, const index_type &vi)
          {
            if ( (ui>=number_u_patches()) || (vi>=number_v_patches()) )
              return INVALID_INDEX;

            // advance to desired index
            index_type i, j;
            typename patch_collection_type::iterator pcit, pcito;
            for (i=0, pcit=patches.begin(); i<ui; ++i, ++pcit) {}
            for (j=0; j<vi; ++j, pcit+=nu) {}

            // check the connectivity on adjacent nodes (if available)
            if (ui>0)
            {
              pcito=pcit;
              --pcito;
              if (!check_u_continuity(pcito->s, surf, eli::geom::general::C0))
              {
                return PATCH_NOT_CONNECTED;
              }
            }
            if ((ui+1)<number_u_patches())
            {
              pcito=pcit;
              ++pcito;
              if (!check_u_continuity(surf, pcito->s, eli::geom::general::C0))
              {
                return PATCH_NOT_CONNECTED;
              }
            }
            if (vi>0)
            {
              pcito=pcit;
              pcito-=nu;
              if (!check_v_continuity(pcito->s, surf, eli::geom::general::C0))
              {
                return PATCH_NOT_CONNECTED;
              }
            }
            if ((vi+1)<number_v_patches())
            {
              pcito=pcit;
              pcito+=nu;
              if (!check_u_continuity(surf, pcito->s, eli::geom::general::C0))
              {
                return PATCH_NOT_CONNECTED;
              }
            }

            // set the new surf
            pcit->s=surf;

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
          }

#if 0
          error_code split_u(const data_type &u0)
          {
// This needs to split the entire column of surfaces that correspond to the u0 location
#if 0
            // find patch that corresponds to given t
            typename patch_collection_type::iterator it;
            data_type uu(0), vv(0), u, v;
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
              return INVALID_PARAM;

            // split the patch and replace
            surface_type cl, cr;
            patch_info stl, str;
            it->s.split(cl, cr, tt);
            stl.c=cl;
            stl.delta_t=it->delta_t*tt;
            str.c=cr;
            str.delta_t=it->delta_t*(1-tt);
            (*it)=str;
            patches.insert(it, stl);

            assert(check_continuity(eli::geom::general::C0));

            return NO_ERROR;
#else
            return UNKNOWN_ERROR;
#endif
          }

          error_code split_v(const data_type &v0)
          {
// This needs to split the entire row of surfaces that correspond to the v0 location
            return UNKNOWN_ERROR;
          }
#endif

          point_type f(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f(uu, vv);
          }

          point_type f_u(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_u(uu, vv)/it->delta_u;
          }

          point_type f_v(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_v(uu, vv)/it->delta_v;
          }

          point_type f_uu(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_uu(uu, vv)/(it->delta_u*it->delta_u);
          }

          point_type f_uv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_uv(uu, vv)/(it->delta_u*it->delta_v);
          }

          point_type f_vv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_vv(uu, vv)/(it->delta_v*it->delta_v);
          }

          point_type f_uuu(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_uuu(uu, vv)/(it->delta_u*it->delta_u*it->delta_u);
          }

          point_type f_uuv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_uuv(uu, vv)/(it->delta_u*it->delta_u*it->delta_v);
          }

          point_type f_uvv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_uvv(uu, vv)/(it->delta_u*it->delta_v*it->delta_v);
          }

          point_type f_vvv(const data_type &u, const data_type &v) const
          {
            // find patch that corresponds to given t
            typename patch_collection_type::const_iterator it;
            data_type uu(0), vv(0);
            find_patch(it, uu, vv, u, v);

            if (it==patches.end())
            {
              assert(false);
              --it;
            }

            return it->s.f_vvv(uu, vv)/(it->delta_v*it->delta_v*it->delta_v);
          }

          point_type normal(const data_type &u, const data_type &v) const
          {
            point_type n=f_u(u, v).cross(f_v(u, v));
            n.normalize();
            return n;
          }

          // TODO: NEED TO IMPLEMENT
          //       * fit
          //       * interpolate

        private:
//           template<template<typename, unsigned short, typename> class surf1__,
//                    typename data1__, unsigned short dim1__, typename tol1__>
//           friend void area(typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &len,
//                            const piecewise<surf1__, data1__, dim1__, tol1__> &pc,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &tol);
//           template<template<typename, unsigned short, typename> class surf1__,
//                             typename data1__, unsigned short dim1__, typename tol1__>
//           friend void area(typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &len,
//                            const piecewise<surf1__, data1__, dim1__, tol1__> &pc,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &t0,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &t1,
//                            const typename piecewise<surf1__, data1__, dim1__, tol1__>::data_type &tol);

          struct patch_info
          {
            surface_type s;
            data_type delta_u, delta_v;

            patch_info() : delta_u(1), delta_v(1) {}
            patch_info(const patch_info &si) : s(si.s), delta_u(si.delta_u), delta_v(si.delta_v) {}
            ~patch_info() {}

            bool operator==(const patch_info &si) const
            {
              if (this==&si)
                return true;
              if (delta_u!=si.delta_u)
                return false;
              if (delta_v!=si.delta_v)
                return false;
              if (s!=si.s)
                return false;

              return true;
            }

            bool operator!=(const patch_info &si) const
            {
              return !operator==(si);
            }
          };
          typedef std::vector<patch_info> patch_collection_type;

          patch_collection_type patches;
          data_type u0, v0;
          index_type nu, nv;

        protected:
          bool check_continuity(const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

          bool check_u_continuity(const surface_type &/*s1*/, const surface_type &/*s2*/, const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

          bool check_v_continuity(const surface_type &/*s1*/, const surface_type &/*s2*/, const eli::geom::general::continuity &/*cont*/) const
          {
            // TODO: Need to implement this
            return true;
          }

        private:
          void find_patch(typename patch_collection_type::const_iterator &it,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in) const
          {
            data_type u(u0), v(v0);
            index_type i, j;

            // check to see if have invalid u_in or v_in
            if ((u_in<u0) || (v_in<v0))
            {
              it=patches.end();
              return;
            }

            // cycle through the u-coordinates to find match
            for (i=0, it=patches.begin(); i<nu; ++i, ++it)
            {
              if (u_in<=u+it->delta_u)
              {
                uu=(u_in-u)/it->delta_u;
                break;
              }
              u+=it->delta_u;
            }
            if (i==nu)
            {
              assert(false);
              it=patches.end();
              return;
            }

            // cycle through the v-coordinates to find match
            for (j=0; j<nv; ++j, it+=nu)
            {
              if (v_in<=v+it->delta_v)
              {
                vv=(v_in-v)/it->delta_v;
//               std::cout << "searching for (u,v)=(" << u_in << "," << v_in << ")" << std::endl;
//               std::cout << "  found (i,j)=(" << i << "," << j << ")" << std::endl;
//               std::cout << "  local (u,v)=(" << uu << "," << vv << ")" << std::endl;
//               std::cout << "  distance=" << std::distance(patches.begin(), it) << std::endl;
                return;
              }
              v+=it->delta_v;
            }
            if (j==nv)
            {
              assert(false);
              it=patches.end();
              return;
            }
          }

          void find_patch(typename patch_collection_type::iterator &it,
                          data_type &uu, data_type &vv,
                          const data_type &u_in, const data_type &v_in)
          {
            data_type u(u0), v(v0);
            index_type i, j;

            // check to see if have invalid u_in or v_in
            if ((u_in<u0) || (v_in<v0))
            {
              it=patches.end();
              return;
            }

            // cycle through the u-coordinates to find match
            for (i=0, it=patches.begin(); i<nu; ++i, ++it)
            {
              if (u_in<=u+it->delta_u)
              {
                uu=(u_in-u)/it->delta_u;
                break;
              }
            }
            if (i==nu)
            {
              it=patches.end();
              return;
            }

            // cycle through the v-coordinates to find match
            for (j=0; j<nv; ++j, it+=nu)
            {
              if (v_in<=v+it->delta_v)
              {
                vv=(v_in-v)/it->delta_v;
                return;
              }
            }
            if (j==nv)
            {
              it=patches.end();
              return;
            }
          }
      };
    }
  }
}
#endif