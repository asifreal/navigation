/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, 2013, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Eitan Marder-Eppstein
 *         David V. Lu!!
 *********************************************************************/
#ifndef _PRM_H
#define _PRM_H

#include <global_planner/planner_core.h>
#include <global_planner/expander.h>
#include <global_planner/astar.h>
#include <global_planner/kdtree.h>
#include <vector>
#include <algorithm>
#include <set>

namespace global_planner {

class Point2d {
    public:
    Point2d(int x, int y, int idx) {
        x_ = x;
        y_ = y;
        idx_ = idx;
    }
    int x_, y_, idx_;
};

class PRMExpansion : public Expander {
    public:
        PRMExpansion(PotentialCalculator* p_calc, int nx, int ny);
        virtual ~PRMExpansion() {}
        bool calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y, int cycles,
                                float* potential);
    private:
        bool generateCandidatePoint(unsigned char* costs, int start_x, int start_y, int end_x, int end_y, 
                                                                      std::vector<Point2d >& points, std::map<int, int> &index);
        bool isGoodPoint(unsigned char* costs, int x, int y);
        void putPointInArray(int x, int y, int idx, std::vector<Point2d>& arr) {
            Point2d point(x, y, idx);
            arr.push_back(point);
        }
        bool isInSight(unsigned char* costs, std::vector<Point2d >& points, std::vector<Point2d >& path,int cur, int nxt);
        void add(unsigned char* costs, float* potential, std::vector<Point2d >& points, int cur, int nxt);
        std::vector<Index> queue_; // for A*
        float ratio_;
        int knn_;
};

} //end namespace global_planner
#endif

