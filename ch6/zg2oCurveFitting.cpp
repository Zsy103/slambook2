//$y = e^{ap^2 + bp + c}$
#include <iostream> // C++标准输入输出库
#include <g2o/core/g2o_core_api.h> // g2o核心库
#include <g2o/core/base_vertex.h> // g2o顶点基类
#include <g2o/core/base_unary_edge.h> // g2o一元边基类
#include <g2o/core/block_solver.h> // g2o块求解器
#include <g2o/core/optimization_algorithm_levenberg.h> // g2o Levenberg求解器
#include <g2o/core/optimization_algorithm_gauss_newton.h> // g2o Gauss-Newton求解器
#include <g2o/core/optimization_algorithm_dogleg.h> // g2o Dogleg求解器
#include <g2o/solvers/dense/linear_solver_dense.h> // g2o稠密线性求解器
#include <Eigen/Core> // Eigen核心库
#include <opencv2/core/core.hpp> // OpenCV核心库
#include <cmath> // 数学库
#include <chrono> // 时间库

using namespace std; // 标准命名空间
// preStep1: 定义图的顶点，模板参数：优化变量维度和数据类型
class CurveFittingVertex : public g2o::BaseVertex<3, Eigen::Vector3d> // 3: 优化变量维度(在流形空间manifold的最小表示), Eigen::Vector3d: 优化变量类型
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // 宏定义，确保new操作符对齐内存
    CurveFittingVertex() {} // 构造函数
    // 重置
    virtual void setToOriginImpl() override // 重置函数，设置优化变量的初始值
    {
        _estimate << 0, 0, 0; // 优化变量初始化
    }

    // 更新
    virtual void oplusImpl(const double *update) override // 更新函数，用于将增量update($\Delta x$)加到优化变量($x$)上，即$x \leftarrow x + \Delta x$
    {
        _estimate += Eigen::Vector3d(update); // 优化变量更新
    }

    // 存盘和读盘：留空
    virtual bool read(istream &in) {} // 存盘操作
    virtual bool write(ostream &out) const {} // 读盘操作
};
// preStep2: 定义图的边（误差模型），模板参数：观测值维度，类型，连接顶点类型
class CurveFittingEdge : public g2o::BaseUnaryEdge<1, double, CurveFittingVertex> // 1: 观测值维度, double: 观测值类型, CurveFittingVertex: 连接顶点类型
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // 宏定义，确保new操作符对齐内存
    CurveFittingEdge(double p) : BaseUnaryEdge(), _p(p) {} // 构造函数
    // 计算曲线模型误差
    virtual void computeError() override
    {   // 计算误差函数，即计算$e = y - \hat{y} = y - e^{ap^2 + bp + c}$，其中$\hat{y} = e^{ap^2 + bp + c}$
        const CurveFittingVertex *v = static_cast<const CurveFittingVertex *> (_vertices[0]); // 顶点指针：获取 _vertices 数组中的第一个元素，并将其转换为 CurveFittingVertex* 类型的指针，然后将这个指针赋值给 v
        const Eigen::Vector3d abc = v->estimate(); // 顶点的优化变量
        _error(0, 0) = _measurement - exp(abc(0, 0) * _p * _p + abc(1, 0) * _p + abc(2, 0)); // 计算误差 $e = y - \hat{y} = y - e^{ap^`2 + bp + c}$
    }

    // 计算雅可比矩阵
    virtual void linearizeOplus() override   
    {   // 计算雅可比矩阵函数，即计算$\frac{\partial e}{\partial x}$
        const CurveFittingVertex *v = static_cast<const CurveFittingVertex *> (_vertices[0]); // 顶点指针
        const Eigen::Vector3d abc = v->estimate(); // 顶点的优化变量
        double y = exp(abc(0, 0) * _p * _p + abc(1, 0) * _p + abc(2, 0)); // 计算$\hat{y} = e^{ap^2 + bp + c}$
        _jacobianOplusXi[0] = -_p * _p * y; // 计算$\frac{\partial e}{\partial a} = -p^2\hat{y}$
        _jacobianOplusXi[1] = -_p * y; // 计算$\frac{\partial e}{\partial b} = -p\hat{y}$
        _jacobianOplusXi[2] = -y; // 计算$\frac{\partial e}{\partial c} = -\hat{y}$
    }
    virtual bool read(istream &in) {} // 存盘操作
    virtual bool write(ostream &out) const {} // 读盘操作
public:
    double _p; // p 值，即$p$, y值为_measurement
};
int main()
{
    // 生成数据
    double a_gt = 1.0, b_gt = 2.0, c_gt = 1.0; // 真实参数值
    double a_est = 2.0, b_est = -1.0, c_est = 5.0; // 估计参数值
    int N = 100; // 数据点个数
    double w_sigma = 1.0; // 噪声Sigma值
    cv::RNG rng; // OpenCV随机数生成器
    vector<double> p_data, y_data; // 数据容器
    for (int i = 0; i < N; i++) // 循环生成数据
    {
        double p = i / 100.0; // p值
        p_data.push_back(p); // p值存入数据容器
        y_data.push_back(exp(a_gt * p * p + b_gt * p + c_gt) + rng.gaussian(w_sigma * w_sigma)); // y值存入数据容器
    }
    //g2o编程流程
    //step1: 设置BlockSolver的优化变量维度和误差值维度
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>> Block; // 3: 优化变量维度, 1: 误差值维度
    
    //step2: 设置线性求解器LinearSolver类型
    typedef g2o::LinearSolverDense<Block::PoseMatrixType> LinearSolver; // 线性求解器类型
    
    //step3: 创建总求解器Solver,梯度下降从GN,LM,Dogleg中选择一种
    //Gauss-Newton方法
    auto solver = new g2o::OptimizationAlgorithmGaussNewton(
        g2o::make_unique<Block>(g2o::make_unique<LinearSolver>())); // GN求解器
    //Levenberg方法
    // auto solver = new g2o::OptimizationAlgorithmLevenberg(
    //     g2o::make_unique<Block>(g2o::make_unique<LinearSolver>())); // LM求解器
    //Dogleg方法
    // auto solver = new g2o::OptimizationAlgorithmDogleg(
    //     g2o::make_unique<Block>(g2o::make_unique<LinearSolver>())); // Dogleg求解器
    
    //step4: 创建稀疏优化器SparseOptimizer
    g2o::SparseOptimizer optimizer; // 创建稀疏优化器
    optimizer.setAlgorithm(solver); // 设置求解器
    optimizer.setVerbose(true); // 打开调试输出

    //step5: 添加顶点
    CurveFittingVertex *v = new CurveFittingVertex(); // 创建顶点
    v->setEstimate(Eigen::Vector3d(a_est, b_est, c_est)); // 设置顶点的初始值
    v->setId(0); // 设置顶点的id
    optimizer.addVertex(v); // 添加顶点
    //step6: 添加边
    for (int i = 0; i < N; i++) // 循环添加边
    {
        CurveFittingEdge *edge = new CurveFittingEdge(p_data[i]); // 创建边
        edge->setId(i); // 设置边的id
        edge->setVertex(0, v); // 设置边的连接顶点
        edge->setMeasurement(y_data[i]); // 设置边的观测值
        edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity() * 1 / (w_sigma * w_sigma)); // 设置边的信息矩阵
        optimizer.addEdge(edge); // 添加边
    }
    //step7: 设置优化参数，开始执行优化
    cout << "start optimization" << endl; // 输出开始优化
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now(); // 计时开始
    optimizer.initializeOptimization(); // 初始化优化器
    optimizer.optimize(100); // 执行优化100次
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now(); // 计时结束
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t2 - t1); // 计算时间差
    cout << "solve time cost = " << time_used.count() << " seconds. " << endl; // 输出时间差
    //step8: 输出优化值
    Eigen::Vector3d abc_estimate = v->estimate(); // 优化值
    cout << "estimated model: " << abc_estimate.transpose() << endl; // 输出优化值
    return 0; // 程序结束
}
