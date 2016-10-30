/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of the ETH Zurich nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __ARRAY_H__
#define __ARRAY_H__

template<typename T>
class Array2D {
protected:
  size_t nx, ny;
  T* data;
public:
  Array2D(size_t x, size_t y) :
    nx(x), ny(y) 
  {
    data = (T*) malloc(x*y*sizeof(T));
  }
  virtual ~Array2D() { free(data); }
  inline T* pt() { return data; }
  inline T* pt(size_t x, size_t y) { return data + (x*ny + y); }
  inline const T* col_pt(size_t r) const { return data + r*ny; }
  inline T* operator()(size_t r) { return data + r*ny; }
  inline const T* operator()(size_t r) const { return data + r*ny; }
  inline T& operator()(size_t x, size_t y) { return data[x*ny + y]; }
  inline T& operator[](int i) { return data[i]; }
  inline void set_all(const T& val) { for (size_t i(0); i < nx*ny; ++i) data[i] = val; }
  inline size_t n_x() const { return nx; }
  inline size_t n_y() const { return ny; }
  inline size_t size() const { return nx*ny; }
};

template<typename T>
class Array3D {
protected:
  size_t nx, ny, nz;
  T* data;
public:
  Array3D(size_t x, size_t y, size_t z) :
    nx(x), ny(y), nz(z) 
  {
    data = (T*) malloc(x*y*z*sizeof(T));
  }
  virtual ~Array3D() { free(data); }
  inline T* pt() { return data; }
  inline T* pt(size_t x, size_t y, size_t z = 0) { return data + (x*ny*nz + y*nz + z); }
  inline T& operator()(size_t x, size_t y, size_t z) { return data[x*ny*nz + y*nz + z]; }
  inline void set_all(const T& val) { for (size_t i(0); i < nx*ny*nz; ++i) data[i] = val; }
  inline size_t n_x() const { return nx; }
  inline size_t n_y() const { return ny; }
  inline size_t n_z() const { return nz; }
};

// VIEWS =====================================================================

template<typename T>
class ArrayView {
protected:
  size_t n;
  T* data;
public:
  ArrayView(size_t n, double* x) : n(n), data(x) {}
  ArrayView(const ArrayView& a) : n(a.n), data(a.data) {}
  virtual ~ArrayView() {}
  inline T& operator[](int i) { return data[i]; }
};

template<typename T>
class Array2DView {
protected:
  size_t nx, ny;
  T* data;
public:
  Array2DView(size_t x, size_t y, T* a) : nx(x), ny(y), data(a) {}
  Array2DView(const Array2DView& a) : nx(a.nx), ny(a.ny), data(a.data) {}
  virtual ~Array2DView() {}
  inline T& operator()(size_t x, size_t y) { return data[x*ny + y]; }
  inline const T* col_pt(size_t r) const { return data + r*ny; }
  inline T* col_pt(size_t r) { return data + r*ny; }
};

#endif // __ARRAY_H__
