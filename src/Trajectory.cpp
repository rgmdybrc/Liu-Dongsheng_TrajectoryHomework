#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>

typedef std::vector<Eigen::Vector2d> Path;

class TrajectoryGenerator
{
public:
    TrajectoryGenerator() = default;

    std::vector<Eigen::Vector2d> generateBSpline(const Path &path, int degree = 3, double tension = 0.5)
    {
        Path interpolatedPath = interpolatePath(path);
        int n = interpolatedPath.size() - 1;
        int m = n + degree + 1;

        std::vector<double> knots = generateKnots(n, degree, tension);

        std::vector<Eigen::Vector2d> trajectory;
        for (double t = 0.0; t <= 1.0; t += 0.0005)  // 增加采样密度
        {
            Eigen::Vector2d point = deBoor(degree, n, knots, interpolatedPath, t);
            trajectory.push_back(point);
        }

        return smoothTrajectory(trajectory, 7);  // 增加平滑窗口大小
    }

private:
    Path interpolatePath(const Path& originalPath)
    {
        Path interpolatedPath;
        for (size_t i = 0; i < originalPath.size() - 1; ++i)
        {
            interpolatedPath.push_back(originalPath[i]);
            Eigen::Vector2d midpoint1 = originalPath[i] + (originalPath[i+1] - originalPath[i]) * 0.33;
            Eigen::Vector2d midpoint2 = originalPath[i] + (originalPath[i+1] - originalPath[i]) * 0.66;
            interpolatedPath.push_back(midpoint1);
            interpolatedPath.push_back(midpoint2);
        }
        interpolatedPath.push_back(originalPath.back());
        return interpolatedPath;
    }

    std::vector<double> generateKnots(int n, int degree, double tension)
    {
        int m = n + degree + 1;
        std::vector<double> knots(m + 1);
        for (int i = 0; i <= m; ++i)
        {
            if (i < degree + 1)
            {
                knots[i] = 0.0;
            }
            else if (i > m - degree - 1)
            {
                knots[i] = 1.0;
            }
            else
            {
                double t = static_cast<double>(i - degree) / (m - 2 * degree);
                knots[i] = std::pow(t, tension);
            }
        }
        return knots;
    }

    Eigen::Vector2d deBoor(int k, int n, const std::vector<double> &knots, const Path &controlPoints, double t)
    {
        int s = findSpan(n, k, t, knots);
        std::vector<Eigen::Vector2d> d(k + 1);

        for (int j = 0; j <= k; ++j)
        {
            d[j] = controlPoints[s - k + j];
        }

        for (int r = 1; r <= k; ++r)
        {
            for (int j = k; j >= r; --j)
            {
                double alpha = (t - knots[s - k + j]) / (knots[s + 1 + j - r] - knots[s - k + j]);
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }

        return d[k];
    }

    int findSpan(int n, int k, double t, const std::vector<double> &knots)
    {
        if (t >= knots[n + 1]) return n;
        if (t <= knots[k]) return k;

        int low = k;
        int high = n + 1;
        int mid = (low + high) / 2;

        while (t < knots[mid] || t >= knots[mid + 1])
        {
            if (t < knots[mid]) high = mid;
            else low = mid;
            mid = (low + high) / 2;
        }

        return mid;
    }

    std::vector<Eigen::Vector2d> smoothTrajectory(const std::vector<Eigen::Vector2d>& trajectory, int windowSize)
    {
        std::vector<Eigen::Vector2d> smoothedTrajectory;
        for (size_t i = 0; i < trajectory.size(); ++i)
        {
            Eigen::Vector2d sum(0, 0);
            int count = 0;
            for (int j = -windowSize/2; j <= windowSize/2; ++j)
            {
                int index = i + j;
                if (index >= 0 && index < static_cast<int>(trajectory.size()))
                {
                    sum += trajectory[index];
                    count++;
                }
            }
            smoothedTrajectory.push_back(sum / count);
        }
        return smoothedTrajectory;
    }
};

int main()
{
    TrajectoryGenerator generator;
    Path path = {{0, 0}, {1, 2}, {3, 3}, {4, 1}, {5, 0}};

    std::cout << "开始生成轨迹。" << std::endl;
    auto trajectory = generator.generateBSpline(path, 3, 0.5);  // 修改度数和张力

    for (const auto &point : trajectory)
    {
        std::cout << "Point: (" << point.x() << ", " << point.y() << ")\n";
    }

    return 0;
}
