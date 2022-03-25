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
#include<global_planner/prm.h>
#include<costmap_2d/cost_values.h>
#include <map>
#include <algorithm>    // std::sort, std::stable_sort
#include <stdlib.h>
#include <time.h>

namespace global_planner {

PRMExpansion::PRMExpansion(PotentialCalculator* p_calc, int xs, int ys) :
        Expander(p_calc, xs, ys), ratio_(0.01), knn_(8) {     // 0.01 * 500 * 500 = 2500
}

/**
 * @brief 随机生成点，自动过滤掉不合适的点，上限为xs * ys * ratio = 2500。同时生成坐标到数组的mapping关系
 * 
 * @param costs 
 * @param start_x 
 * @param start_y 
 * @param end_x 
 * @param end_y 
 * @param points 
 * @param index 
 * @return true 
 * @return false 
 */
bool PRMExpansion::generateCandidatePoint(unsigned char* costs, int start_x, int start_y, int end_x, int end_y, 
                                                                      std::vector<Point2d >& points, std::map<int, int> &index) {
    srand((unsigned)time(NULL));
    points.clear();
    index.clear();
    int idx = toIndex(end_x, end_y);
    putPointInArray(end_x, end_y, idx, points);
    index[idx] = 0;
    idx = toIndex(start_x, start_y);
    putPointInArray(start_x, start_y, idx, points);
    index[idx] = 1;
    int num = ns_ * ratio_;
    while(num-- > 0) {
        //generate
        int x = rand() % nx_;
        int y = rand() % ny_;
        idx = toIndex(x, y);
        //check valid
        if(index.find(idx) !=index.end()) continue;  //duplicate
        if(!isGoodPoint(costs, x, y)) continue;
        //put in array
        index[idx] = points.size();
        putPointInArray(x, y, idx, points);
    }
    return true;
}

/**
 * @brief 该点与周围8个点是否不行
 * 
 * @param costs 
 * @param x 
 * @param y 
 * @return true 
 * @return false 
 */
bool PRMExpansion::isGoodPoint(unsigned char* costs, int x, int y) {
    for(int i = -1; i < 2; i++){
        for (int j = -1; j < 2; j++){
            int idx = toIndex((x+i), (y+j));
            if(costs[idx]>=lethal_cost_ && !(unknown_ && costs[idx]==costmap_2d::NO_INFORMATION))
                return false;
        }
    }
    return true;
}


bool PRMExpansion::calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y,
                                        int cycles, float* potential) {
    ROS_ERROR("Start PRM calculate!");
    clock_t t0 = clock();
    std::vector<Point2d > rand_points;
    std::map<int, int> mapping;
    generateCandidatePoint(costs,start_x, start_y, end_x, end_y, rand_points, mapping);
    Kdtree::KdNodeVector tree_nodes = Kdtree::KdNodeVector();
    for(int i=0; i<rand_points.size(); i++){
        Kdtree::CoordPoint pos = Kdtree::CoordPoint();
        pos.push_back(rand_points[i].x_);
        pos.push_back(rand_points[i].y_);
        tree_nodes.push_back(Kdtree::KdNode(pos, i));
    }
    Kdtree::KdTree tree_index = Kdtree::KdTree(&tree_nodes);

    queue_.clear();
    int start_i = toIndex(start_x, start_y);
    queue_.push_back(Index(start_i, 0));

    std::fill(potential, potential + ns_, POT_HIGH);
    potential[start_i] = 0;

    int goal_i = toIndex(end_x, end_y);
    int cycle = 0;
    std::set<int> visited = std::set<int>();
    std::vector<Point2d> path = std::vector<Point2d>();
    visited.insert(start_i);
    Kdtree::KdNodeVector neighbors = Kdtree::KdNodeVector();
    while (queue_.size() > 0 && cycle < cycles) { // 最小堆
        Index top = queue_[0];
        std::pop_heap(queue_.begin(), queue_.end(), greater1());
        queue_.pop_back();

        int i = top.i;
        if (i == goal_i)
            return true;
        path.clear();
        if (isInSight(costs, rand_points, path, mapping[i], 0)) {
            for(int i=1; i<path.size();i++){
                Point2d &prev = path[i-1];
                Point2d &p = path[i];
                potential[p.idx_] = std::min(std::max(0, 0+costs[p.idx_]) + neutral_cost_ + potential[prev.idx_], potential[p.idx_]);
            }
            return true;
        }
            

        neighbors.clear();
        Kdtree::CoordPoint pos = Kdtree::CoordPoint();
        pos.push_back(rand_points[mapping[i]].x_);
        pos.push_back(rand_points[mapping[i]].y_);
        tree_index.k_nearest_neighbors(pos, knn_, &neighbors, NULL);
        for (int j=0; j<neighbors.size(); j++) {
            int nearby = neighbors[j]._id;
            if(visited.find(nearby) == visited.end()) {
                add(costs, potential, rand_points, mapping[i], nearby);
                visited.insert(nearby);
            } 
        }
        cycle++;
    }
    clock_t t1 = clock();
    ROS_ERROR("End PRM calculate, time=%.4f!", (double)(t1 - t0) / CLOCKS_PER_SEC);
    return false;
}

/**
 * @brief 检查是否可以直线到达终点，没有碰撞
 * 
 * @param costs 
 * @param points 
 * @param path 计算出的路径
 * @param cur from
 * @param nxt to
 * @return true 
 * @return false 
 */
bool PRMExpansion::isInSight(unsigned char* costs, std::vector<Point2d >& points, std::vector<Point2d >& path, int cur, int nxt) {
    Point2d &center = points[cur];
    Point2d &target = points[nxt];
    int dx = target.x_ - center.x_;
    int dy = target.y_ - center.y_;
    int ax = abs(dx);
    int ay = abs(dy);
    int signx = dx > 0 ? 1 : -1;
    int signy = dy > 0 ? 1 : -1;
    Point2d p = Point2d(center.x_, center.y_, center.idx_);
    path.push_back(p);
    for(int ix = 0, iy=0; ix < ax || iy < ay; ){
        if ((0.5 + ix) / ax < (0.5 + iy) / ay) { // horizontal
            p.x_ += signx;
            ix++;
        } else { // vertical
            p.y_ += signy;
            iy++;
        }
        p.idx_ = toIndex(p.x_, p.y_);
        
        if(costs[p.idx_]>=lethal_cost_ && !(unknown_ && costs[p.idx_]==costmap_2d::NO_INFORMATION)) // 检测碰撞
            return false;

        //potential[p.idx_] = p_calc_->calculatePotential(potential, costs[p.idx_] + neutral_cost_, p.idx_, potential[prev.idx_]);
        //potential[p.idx_] = std::min(std::max(0, 0+costs[p.idx_]) + neutral_cost_ + potential[prev.idx_], potential[p.idx_]);
        path.push_back(Point2d(p.x_, p.y_, p.idx_));  // for debug
    }
    return true;
}

/**
 * @brief 探索路径, 要检测碰撞
 * 
 * @param costs 
 * @param potential 
 * @param points 
 * @param cur 是在points中的下标
 * @param nxt 是在points中的下标
 */
void PRMExpansion::add(unsigned char* costs, float* potential, std::vector<Point2d >& points, int cur, int nxt) {
    Point2d &target = points[nxt];
    Point2d &endpos = points[0];
    
    std::vector<Point2d> path;
    if(!isInSight(costs, points, path, cur, nxt)) 
        return;

    for(int i=1; i<path.size();i++){
        Point2d &prev = path[i-1];
        Point2d &p = path[i];
        potential[p.idx_] = std::min(std::max(0, 0+costs[p.idx_]) + neutral_cost_ + potential[prev.idx_], potential[p.idx_]);
        //ROS_WARN("PATH %d: (%d, %d) = %f", i, path[i].x_, path[i].y_, potential[path[i].idx_]);
    }

    int ddx = (endpos.x_ - target.x_);
    int ddy = (endpos.y_ - target.y_);
    float distance = sqrt(ddx * ddx + ddy * ddy);
    queue_.push_back(Index(target.idx_, potential[target.idx_] + distance * neutral_cost_));
    std::push_heap(queue_.begin(), queue_.end(), greater1());
}


} //end namespace global_planner
