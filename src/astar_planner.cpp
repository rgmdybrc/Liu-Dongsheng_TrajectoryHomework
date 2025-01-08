#include <ros/ros.h>
#include <utility>
#include <vector>
#include <queue>
#include <cmath>
#include <Eigen/Dense>
#include "visualization_msgs/MarkerArray.h"
#include <geometry_msgs/Point.h>
#include "nav_msgs/Path.h"
#include <Eigen/Dense>
#include <vector>
#include <algorithm>

struct Node
{
    int x, y;                     // 节点所在的网格坐标
    double g_cost;                // 从起点到当前节点的代价
    double h_cost;                // 从当前节点到终点的估计代价
    std::shared_ptr<Node> parent; // 父节点，用于回溯路径

    Node(int x, int y, double g_cost, double h_cost, std::shared_ptr<Node> parent = nullptr)
        : x(x), y(y), g_cost(g_cost), h_cost(h_cost), parent(std::move(parent)) {}

    double f() const { return g_cost + h_cost; } // 总代价值
};
// 比较器，用于优先队列
struct cmp
{
    bool operator()(std::shared_ptr<Node> a, std::shared_ptr<Node> b)
    {
        return a->f() > b->f();
    }
};
struct GridMap
{
    int width;
    int height;
    double map_max;
    double map_min;
    double grid_resolution;
    std::vector<std::vector<int>> grid; // 0: 空闲, 1: 占用

    GridMap(int w, int h, double map_min_, double map_max_, double res) : width(w), height(h), map_min(map_min_), map_max(map_max_), grid_resolution(res), grid(w, std::vector<int>(h, 0)) {}

    void markObstacle(double cx, double cy, double radius)
    {
        int grid_cx = std::round((cx - map_min) / grid_resolution);
        int grid_cy = std::round((cy - map_min) / grid_resolution);
        int grid_radius = std::round(radius / grid_resolution);
        // Step 1: 将圆形区域标记为占用
        for (int dx = -grid_radius; dx <= grid_radius; ++dx)
        {
            for (int dy = -grid_radius; dy <= grid_radius; ++dy)
            {
                if (dx * dx + dy * dy <= grid_radius * grid_radius)
                {
                    int x = grid_cx + dx;
                    int y = grid_cy + dy;
                    if (x >= 0 && x < width && y >= 0 && y < height)
                    {
                        grid[x][y] = 1;
                    }
                }
            }
        }
        // your code
        // finish
    }
};
class AStarPlanner
{
public:
    AStarPlanner(int width, int height, double m_min, double m_max, double res, double safety_distance, double line_weight)
        : width_(width), height_(height), map_min_(m_min), map_max_(m_max),
          grid_resolution_(res), grid_map_(width, height, map_min_, map_max_, grid_resolution_),
          num_of_obs_(0), safety_distance_(safety_distance), line_weight_(line_weight) {}

    void setObstacle(double cx, double cy, double radius)
    {
        num_of_obs_++;
        double inflated_radius = radius + safety_distance_;
        grid_map_.markObstacle(cx, cy, inflated_radius);
    }

    void printGridMap()
    {
        for (int i = 0; i < width_; i++)
        {
            for (int j = 0; j < height_; j++)
            {
                std::cout << grid_map_.grid[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "num of obstacles: " << num_of_obs_ << std::endl;
    }

    // 修改findPath函数，使用新的评估函数
    std::vector<Eigen::Vector2d> findPath(Eigen::Vector2d start, Eigen::Vector2d goal)
    {
        auto gridStart = worldToGrid(start);
        auto gridGoal = worldToGrid(goal);

        std::vector<std::vector<bool>> closed_list(width_, std::vector<bool>(height_, false));
        auto current = std::make_shared<Node>(Node(gridStart.first, gridStart.second, 0.0, heuristic(gridStart, gridGoal)));

        while (current->x != gridGoal.first || current->y != gridGoal.second)
        {
            closed_list[current->x][current->y] = true;

            std::vector<std::shared_ptr<Node>> open_list;

            for (const auto &neighbor : getNeighbors(*current))
            {
                if (grid_map_.grid[neighbor.x][neighbor.y] == 1 || closed_list[neighbor.x][neighbor.y])
                {
                    continue;
                }

                double tentative_g_cost = current->g_cost + distance(*current, neighbor);

                auto neighbor_node = std::make_shared<Node>(
                    neighbor.x, neighbor.y,
                    tentative_g_cost,
                    heuristic({neighbor.x, neighbor.y}, gridGoal),
                    current);

                open_list.push_back(neighbor_node);
            }

            if (open_list.empty())
            {
                return {}; // No path found
            }

            // 使用新的评估函数选择最佳节点
            current = *std::min_element(open_list.begin(), open_list.end(),
                                        [this, gridStart, gridGoal](const std::shared_ptr<Node> &a, const std::shared_ptr<Node> &b)
                                        {
                                            return this->evaluateNode(*a, gridStart, gridGoal) < this->evaluateNode(*b, gridStart, gridGoal);
                                        });

            // 检查是否到达目标
            if (current->x == gridGoal.first && current->y == gridGoal.second)
            {
                std::vector<Eigen::Vector2d> path = reconstructPath(current);
                return path;
                // return smoothPathBSpline(path, 3, 0.01);
            }
        }

        return {}; // No path found
    }

    std::vector<Eigen::Vector2d> smoothPathBSpline(const std::vector<Eigen::Vector2d> &path, int degree = 3, double resolution = 0.01)
    {
        if (path.size() < degree + 1)
            return path;

        int n = path.size() - 1;
        int m = n + degree + 1;

        // 生成均匀的节点向量
        std::vector<double> knots(m + 1);
        for (int i = 0; i <= m; ++i)
        {
            if (i < degree + 1)
                knots[i] = 0;
            else if (i > m - degree - 1)
                knots[i] = 1;
            else
                knots[i] = (i - degree) / double(n - degree + 1);
        }

        std::vector<Eigen::Vector2d> smoothed_path;
        for (double t = 0; t <= 1; t += resolution)
        {
            Eigen::Vector2d point = bsplinePoint(t, degree, path, knots);

            // 检查平滑后的点是否满足安全距离要求
            bool is_safe = true;
            for (int i = 0; i < width_; i++)
            {
                for (int j = 0; j < height_; j++)
                {
                    if (grid_map_.grid[i][j] == 1)
                    {
                        Eigen::Vector2d obstacle = gridToWorld(i, j);
                        if ((point - obstacle).norm() < safety_distance_)
                        {
                            is_safe = false;
                            break;
                        }
                    }
                }
                if (!is_safe)
                    break;
            }

            if (is_safe)
            {
                smoothed_path.push_back(point);
            }
            else
            {
                // 如果不安全，可以尝试将点推离障碍物
                // 这里使用一个简单的方法，可能需要更复杂的逻辑来处理
                Eigen::Vector2d safe_point = point;
                for (const auto &original_point : path)
                {
                    if ((original_point - point).norm() < safety_distance_)
                    {
                        safe_point = original_point;
                        break;
                    }
                }
                smoothed_path.push_back(safe_point);
            }
        }

        return smoothed_path;
    }
    double safety_distance_;
    double line_weight_;

    void reset()
    {
        num_of_obs_ = 0;
        grid_map_.grid = std::vector<std::vector<int>>(width_, std::vector<int>(height_, 0));
    }

private:
    // 计算启发式代价（使用欧几里得距离）
    double heuristic(const std::pair<int, int> &from, const std::pair<int, int> &to)
    {
        return std::sqrt(std::pow(from.first - to.first, 2) + std::pow(from.second - to.second, 2));
    }

    // 计算两节点之间的距离（用于邻居代价计算）
    double distance(const Node &a, const Node &b)
    {
        return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
    }

    // 计算点到直线的距离
    double pointToLineDistance(const Eigen::Vector2d &point, const Eigen::Vector2d &lineStart, const Eigen::Vector2d &lineEnd)
    {
        Eigen::Vector2d line = lineEnd - lineStart;
        Eigen::Vector2d pointVector = point - lineStart;
        double lineLengthSquared = line.squaredNorm();

        if (lineLengthSquared == 0.0)
            return (point - lineStart).norm();

        double t = pointVector.dot(line) / lineLengthSquared;
        t = std::max(0.0, std::min(1.0, t));

        Eigen::Vector2d projection = lineStart + t * line;
        return (point - projection).norm();
    }

    // 从世界坐标转换到栅格坐标
    std::pair<int, int> worldToGrid(const Eigen::Vector2d &position)
    {
        int x = std::round((position.x() - map_min_) / grid_resolution_);
        int y = std::round((position.y() - map_min_) / grid_resolution_);
        return {x, y};
    }

    // 从栅格坐标转换到世界坐标（主要用于路径结果显示）
    Eigen::Vector2d gridToWorld(int x, int y)
    {
        double wx = x * grid_resolution_ + map_min_;
        double wy = y * grid_resolution_ + map_min_;
        return Eigen::Vector2d(wx, wy);
    }

    // 获取当前节点的所有邻居节点
    std::vector<Node> getNeighbors(const Node &current)
    {
        std::vector<Node> neighbors;

        // 32方向搜索
        for (int i = 0; i < 32; ++i)
        {
            double angle = i * (2 * M_PI / 32);
            int dx = std::round(std::cos(angle) * 4); // 将间距设为0.25的整数倍
            int dy = std::round(std::sin(angle) * 4);

            int new_x = current.x + dx;
            int new_y = current.y + dy;

            if (new_x >= 0 && new_x < width_ && new_y >= 0 && new_y < height_)
            {
                neighbors.emplace_back(new_x, new_y, 0, 0);
            }
        }

        return neighbors;
    }

    Eigen::Vector2d bsplinePoint(double t, int degree, const std::vector<Eigen::Vector2d> &controlPoints, const std::vector<double> &knots)
    {
        int n = controlPoints.size() - 1;
        int span = findSpan(n, degree, t, knots);
        std::vector<double> basisFuncs = basisFunctions(span, t, degree, knots);

        Eigen::Vector2d point(0, 0);
        for (int i = 0; i <= degree; ++i)
        {
            point += basisFuncs[i] * controlPoints[span - degree + i];
        }

        return point;
    }

    int findSpan(int n, int degree, double t, const std::vector<double> &knots)
    {
        if (t == knots[n + 1])
            return n;

        int low = degree;
        int high = n + 1;
        int mid = (low + high) / 2;

        while (t < knots[mid] || t >= knots[mid + 1])
        {
            if (t < knots[mid])
                high = mid;
            else
                low = mid;
            mid = (low + high) / 2;
        }

        return mid;
    }

    std::vector<double> basisFunctions(int span, double t, int degree, const std::vector<double> &knots)
    {
        std::vector<double> left(degree + 1);
        std::vector<double> right(degree + 1);
        std::vector<double> N(degree + 1);

        N[0] = 1.0;

        for (int j = 1; j <= degree; ++j)
        {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            double saved = 0.0;

            for (int r = 0; r < j; ++r)
            {
                double temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }

            N[j] = saved;
        }

        return N;
    }

    // 修改评估函数，考虑安全距离
    double evaluateNode(const Node &node, const std::pair<int, int> &start, const std::pair<int, int> &goal)
    {
        double g_cost = node.g_cost;
        double h_cost = heuristic({node.x, node.y}, goal);

        // 增加与障碍物距离的评估
        double obstacle_cost = 0;
        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                int nx = node.x + dx;
                int ny = node.y + dy;
                if (nx >= 0 && nx < width_ && ny >= 0 && ny < height_)
                {
                    if (grid_map_.grid[nx][ny] == 1)
                    {
                        obstacle_cost += safety_distance_ / (std::abs(dx) + std::abs(dy) + 0.1);
                    }
                }
            }
        }

        // 计算到直线路径的距离
        Eigen::Vector2d nodePos = gridToWorld(node.x, node.y);
        Eigen::Vector2d startPos = gridToWorld(start.first, start.second);
        Eigen::Vector2d goalPos = gridToWorld(goal.first, goal.second);
        double lineDistance = pointToLineDistance(nodePos, startPos, goalPos);

        // 直线引导成本
        double line_cost = line_weight_ * lineDistance;

        return g_cost + h_cost + obstacle_cost + line_cost;
    }

    // 回溯路径
    std::vector<Eigen::Vector2d> reconstructPath(std::shared_ptr<Node> node)
    {
        std::vector<Eigen::Vector2d> path;
        while (node)
        {
            path.push_back(gridToWorld(node->x, node->y));
            node = node->parent;
        }
        std::reverse(path.begin(), path.end());
        reset();
        return path;
    }

    double cross2D(const Eigen::Vector2d &a, const Eigen::Vector2d &b)
    {
        return a.x() * b.y() - a.y() * b.x();
    }

    // 地图数据
    int width_, height_;
    double map_min_, map_max_, grid_resolution_;
    GridMap grid_map_; // 栅格地图，0: 空闲，1: 障碍物
    int num_of_obs_;
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "astar_planner");
    ros::NodeHandle nh;
    double map_min_, map_max_, grid_resolution_;
    double start_x_, start_y_, goal_x_, goal_y_;
    double safety_distance, line_weight;

    nh.param("astar_planner/map_min", map_min_, -5.0);
    nh.param("astar_planner/map_max", map_max_, 5.0);
    nh.param("astar_planner/grid_resolution", grid_resolution_, 0.05);
    nh.param("astar_planner/start_x", start_x_, -4.5);
    nh.param("astar_planner/start_y", start_y_, -4.5);
    nh.param("astar_planner/goal_x", goal_x_, 4.5);
    nh.param("astar_planner/goal_y", goal_y_, 4.5);
    nh.param("astar_planner/safety_distance", safety_distance, 0.2);
    nh.param("astar_planner/line_weight", line_weight, 0.1);

    // 计算网格的宽度和高度
    int grid_width = std::round((map_max_ - map_min_) / grid_resolution_);
    int grid_height = grid_width;  // 假设地图是正方形的

    AStarPlanner planner(grid_width, grid_height, map_min_, map_max_, grid_resolution_, safety_distance, line_weight);


    // 障碍物订阅
    ros::Subscriber obstacle_sub = nh.subscribe<visualization_msgs::MarkerArray>("obstacles", 1,
                                                                                 [&planner, &grid_resolution_, &map_min_](const visualization_msgs::MarkerArray::ConstPtr &msg)
                                                                                 {
                                                                                     for (const auto &marker : msg->markers)
                                                                                     {
                                                                                         planner.setObstacle(marker.pose.position.x, marker.pose.position.y, marker.scale.x / 2.0);
                                                                                     }
                                                                                 });

    // 发布路径
    ros::Rate rate(10);
    ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("path", 1);
    // 起点和终点参数
    Eigen::Vector2d start(start_x_, start_y_);
    Eigen::Vector2d goal(goal_x_, goal_y_);
    while (ros::ok())
    {
        planner.reset();
        //        // 等待障碍物加载
        //        ros::Duration(1.0).sleep();
        ros::spinOnce();
        // 执行路径搜索
        std::vector<Eigen::Vector2d> path = planner.findPath(start, goal);

        // 路径可视化
        if (path.empty())
        {
            continue;
        }
        nav_msgs::Path path_msg;
        path_msg.header.frame_id = "map";
        path_msg.header.stamp = ros::Time::now();
        for (const auto &point : path)
        {
            geometry_msgs::PoseStamped pose;
            pose.pose.position.x = point.x();
            pose.pose.position.y = point.y();
            pose.pose.position.z = 0.0; // 平面路径，z 设置为 0
            path_msg.poses.push_back(pose);
        }
        path_pub.publish(path_msg);
        rate.sleep();
    }

    return 0;
}
