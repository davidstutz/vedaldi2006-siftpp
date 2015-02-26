// file:        sift.hpp
// author:      Andrea Vedaldi
// description: Sift declaration

// AUTORIGHTS
// Copyright (c) 2006 The Regents of the University of California
// All Rights Reserved.
// 
// Created by Andrea Vedaldi (UCLA VisionLab)
// 
// Permission to use, copy, modify, and distribute this software and its
// documentation for educational, research and non-profit purposes,
// without fee, and without a written agreement is hereby granted,
// provided that the above copyright notice, this paragraph and the
// following three paragraphs appear in all copies.
// 
// This software program and documentation are copyrighted by The Regents
// of the University of California. The software program and
// documentation are supplied "as is", without any accompanying services
// from The Regents. The Regents does not warrant that the operation of
// the program will be uninterrupted or error-free. The end-user
// understands that the program was developed for research purposes and
// is advised not to rely exclusively on the program for any reason.
// 
// This software embodies a method for which the following patent has
// been issued: "Method and apparatus for identifying scale invariant
// features in an image and use of same for locating an object in an
// image," David G. Lowe, US Patent 6,711,293 (March 23,
// 2004). Provisional application filed March 8, 1999. Asignee: The
// University of British Columbia.
// 
// IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
// FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
// INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
// ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
// CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
// BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
// MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

#ifndef VL_SIFT_HPP
#define VL_SIFT_HPP

#include<valarray>
#include<vector>
#include<ostream>
#include<assert.h>

namespace VL {

/** @brief Generic exception */
struct Exception
{
  Exception(std::string _msg) : msg(_msg) { }
  std::string msg ;
} ;

#define VL_THROW(x)                             \
  {                                             \
    std::ostringstream oss ;                    \
    oss << x ;                                  \
    throw VL::Exception(oss.str()) ;            \
  }

/** @brief PGM buffer descriptor
 **
 ** The structure describes a buffer for a gray scale image
 ** which is used by PGM input/output functions.
 **/
struct PgmBuffer
{
  int width ;    ///< Image width
  int height ;   ///< Image hegith
  float* data ;  ///< Image data
} ;

std::ostream& insertPgm(std::ostream&, float const* im, int width, int height) ;
std::istream& extractPgm(std::istream&, PgmBuffer& buffer) ;

/** @brief SIFT Gaussian Scale Space Filter
 **
 ** This class implements a filter computing the SIFT Gaussian scale
 ** space of an image.
 **/
class Sift
{

public:

  struct Keypoint
  {
    int o ;    ///< Keypoint octave index

    int ix ;   ///< Keypoint integer X coordinate (unnormalized)
    int iy ;   ///< Keypoint integer Y coordinate (unnormalized)
    int is ;   ///< Keypoint integer 

    float x  ;  ///< Keypoint fractional X coordinate
    float y  ;  ///< Keypoint fractional Y coordinate
    float s ;   ///< Keypoint fractional scale index

    float sigma ;  ///< Keypoint scale
  } ; 
  
  typedef std::vector<Keypoint>     Keypoints ;
  typedef Keypoints::iterator       KeypointsIter ;
  typedef Keypoints::const_iterator KeypointsConstIter ;

  Sift(const float* _im_pt, int _width, int _height,
                 float _sigman,
                 float _sigma0,
                 int _O, int _S,
                 int _omin, int _smin, int _smax) ;
  ~Sift() ;

  void process(const float* _im_pt, int _width, int _height) ;

  /** @brief Querying the Gaussian scale space */
  /*@{*/
  float* getOctave(int o) ;
  float* getLevel(int o, int s) ;
  int    getWidth() const ;
  int    getHeight() const ;
  int    getOctaveWidth(int o) const ;
  int    getOctaveHeight(int o) const ;
  float  getOctaveSamplingPeriod(int o) const ;
  float  getScaleFromIndex(float o, float s) const ;
  /*@{*/

  /** @brief Detector */
  void  detectKeypoints() ;
  int computeKeypointOrientations(float angles [4], Keypoint keypoint) ; 
  void  computeKeypointDescriptor(float* descr_pt, Keypoint keypoint, float angle) ;

  KeypointsIter keypointsBegin() { return keypoints.begin() ; } 
  KeypointsIter keypointsEnd()   { return keypoints.end() ; }
    
private:
  void prepareBuffers() ;
  void freeBuffers() ;
  void smooth(float* dst, float* temp, 
              float const* src, int width, int height, 
              float s) ;

  void prepareGrad(int o) ;
  
  // scale space parameters
  float sigman ;
  float sigma0 ;
  float sigmak ;

  int O ;
  int S ; 
  int omin ;
  int smin ; 
  int smax ;

  int width ;
  int height ;

  // buffers
  float*  temp ;
  int     tempReserved ;
  bool    tempIsGrad  ;
  int     tempOctave ;
  float** octaves ;
  
  float* filter ;
  int    filterReserved ;

  // keypoints
  float threshold ;
  float r ;

  Keypoints keypoints ;
  
} ;


}

#include "sift.ipp"

#endif
