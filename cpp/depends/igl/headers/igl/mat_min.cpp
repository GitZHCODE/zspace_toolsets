// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mat_min.h"

template <typename DerivedX, typename DerivedY, typename DerivedI>
IGL_INLINE void igl::mat_min(
  const Eigen::DenseBase<DerivedX> & X,
  const int dim,
  Eigen::PlainObjectBase<DerivedY> & Y,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  assert(dim==1||dim==2);

  // output size
  int n = (dim==1?X.cols():X.rows());
  // resize output
  Y.resize(n,1);
  I.resize(n,1);

  // loop over dimension opposite of dim
  for(int j = 0;j<n;j++)
  {
    typename DerivedX::Index PHONY,i;
    typename DerivedX::Scalar  m;
    if(dim==1)
    {
      m = X.col(j).minCoeff(&i,&PHONY);
    }else
    {
      m = X.row(j).minCoeff(&PHONY,&i);
    }
    Y(j) = m;
    I(j) = i;
  }
}

//template <typename T>
//IGL_INLINE Eigen::Matrix<T,Eigen::Dynamic,1> igl::mat_min(
//  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
//  const int dim)
//{
//  Eigen::Matrix<T,Eigen::Dynamic,1> Y;
//  Eigen::Matrix<int,Eigen::Dynamic,1> I;
//  mat_min(X,dim,Y,I);
//  return Y;
//}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::mat_min<Eigen::Array<bool, -1, 3, 0, -1, 3>, Eigen::Array<bool, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::DenseBase<Eigen::Array<bool, -1, 3, 0, -1, 3> > const&, int, Eigen::PlainObjectBase<Eigen::Array<bool, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
// generated by autoexplicit.sh
template void igl::mat_min<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::DenseBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif