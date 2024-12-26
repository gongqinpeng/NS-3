#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>

// 定义个体
class DOA_Unit
{
public:
    std::vector<double> position; // 节点位置
    double value; // 适应度值
    double survival_rate;
    std::vector<std::vector<bool>> connections; // 节点连接矩阵
};

class DOA_Base
{
public:
    DOA_Base(int num_nodes, int size, int iter_max, const std::vector<double>& range_min_list, const std::vector<double>& range_max_list)
            : m_num_nodes(num_nodes), m_size(size), m_iter_max(iter_max), m_range_min_list(range_min_list), m_range_max_list(range_max_list)
    {
        m_name = "DOA";
        m_P = 0.5;
        m_Q = 0.7;
        m_na_ini = 2;
        m_na_end = std::floor(m_size / 2.0);
    }

    void init();
    void run();
    void update_survival_rate();
    void update_position();
    std::vector<double> get_sumatory(int na, int id);
    int get_best_id();
    void export_topology(const std::string& filename);

protected:
    std::string m_name;
    double m_P, m_Q;
    int m_na_ini, m_na_end;
    int m_num_nodes, m_size, m_iter_max;
    std::vector<double> m_range_min_list, m_range_max_list;
    std::vector<DOA_Unit> m_unit_list;
    std::vector<double> m_position_best;

private:
    double cal_fitfunction(const std::vector<double>& position, const std::vector<std::vector<bool>>& connections);
    std::vector<double> get_out_bound_value(const std::vector<double>& position);
    std::vector<double> levy(int d);
};

void DOA_Base::init()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < m_size; ++i)
    {
        DOA_Unit unit;
        unit.position.resize(m_num_nodes * 2); // 二维位置
        unit.connections.resize(m_num_nodes, std::vector<bool>(m_num_nodes, false));

        // 随机生成节点位置
        for (int j = 0; j < m_num_nodes * 2; ++j)
        {
            std::uniform_real_distribution<> dis(m_range_min_list[j % 2], m_range_max_list[j % 2]);
            unit.position[j] = dis(gen);
        }

        // 随机生成连接
        for (int j = 0; j < m_num_nodes; ++j)
        {
            for (int k = j + 1; k < m_num_nodes; ++k)
            {
                if (std::uniform_real_distribution<>(0, 1)(gen) < 0.5)
                {
                    unit.connections[j][k] = unit.connections[k][j] = true;
                }
            }
        }

        unit.value = cal_fitfunction(unit.position, unit.connections);
        unit.survival_rate = 0.0; // 初始化为0
        m_unit_list.push_back(unit);
    }
    m_position_best = m_unit_list[get_best_id()].position;
}

void DOA_Base::run()
{
    for (int iter = 0; iter < m_iter_max; ++iter)
    {
        update_survival_rate();
        update_position();
        std::cout << "Iteration " << iter << ": Best Value = " << m_unit_list[get_best_id()].value << std::endl;
    }
}

void DOA_Base::update_survival_rate()
{
    std::vector<double> values;
    for (const auto& unit : m_unit_list)
    {
        values.push_back(unit.value);
    }
    std::sort(values.begin(), values.end(), std::greater<double>());

    for (auto& unit : m_unit_list)
    {
        double survival_rate = std::abs(values[0] - unit.value) / std::abs(values[0] - values.back()) + std::numeric_limits<double>::min();
        unit.survival_rate = survival_rate;
    }
}

void DOA_Base::update_position()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis_na(m_na_ini, m_na_end);
    std::uniform_real_distribution<> dis(-2, 2);
    std::uniform_real_distribution<> dis1(-1, 1);
    std::uniform_real_distribution<> dis_scavenging(-1, 1);
    std::uniform_real_distribution<> dis_binary(0, 1);

    for (int i = 0; i < m_size; ++i)
    {
        int na = dis_na(gen);
        if (std::uniform_real_distribution<>(0, 1)(gen) < m_P)
        {
            if (std::uniform_real_distribution<>(0, 1)(gen) < m_Q)
            {
                auto pos_sumatory = get_sumatory(na, i);
                std::vector<double> new_pos;
                for (int j = 0; j < m_num_nodes * 2; ++j)
                {
                    new_pos.push_back(pos_sumatory[j] * dis(gen) - m_position_best[j]);
                }
                new_pos = get_out_bound_value(new_pos);
                std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
                double new_value = cal_fitfunction(new_pos, new_connections);
                if (new_value > m_unit_list[i].value)
                {
                    m_unit_list[i].position = new_pos;
                    m_unit_list[i].value = new_value;
                }
            }
            else
            {
                // 这里简化了逻辑，实际上应包括连接的更新
                int r_id = std::uniform_int_distribution<>(0, m_size - 1)(gen);
                double beta1 = dis(gen);
                double beta2 = dis1(gen);
                std::vector<double> new_pos;
                for (int j = 0; j < m_num_nodes * 2; ++j)
                {
                    new_pos.push_back(m_position_best[j] + beta1 * std::exp(beta2) * (m_unit_list[r_id].position[j] - m_unit_list[i].position[j]));
                }
                new_pos = get_out_bound_value(new_pos);
                std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
                double new_value = cal_fitfunction(new_pos, new_connections);
                if (new_value > m_unit_list[i].value)
                {
                    m_unit_list[i].position = new_pos;
                    m_unit_list[i].value = new_value;
                }
            }
        }
        else
        {
            // 食腐行为
            int r_id = std::uniform_int_distribution<>(0, m_size - 1)(gen);
            double m = dis_scavenging(gen);
            double n = dis_binary(gen) <= 0.5 ? 0 : 1;
            std::vector<double> new_pos;
            for (int j = 0; j < m_num_nodes * 2; ++j)
            {
                new_pos.push_back(0.5 * (std::exp(m) * m_unit_list[r_id].position[j] - std::pow(-1, n) * m_unit_list[i].position[j]));
            }
            new_pos = get_out_bound_value(new_pos);
            std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
            double new_value = cal_fitfunction(new_pos, new_connections);
            if (new_value > m_unit_list[i].value)
            {
                m_unit_list[i].position = new_pos;
                m_unit_list[i].value = new_value;
            }
        }

        // 对于存活率低于0.3的个体，使用特殊策略更新位置
        if (m_unit_list[i].survival_rate < 0.3)
        {
            int r_id1 = std::uniform_int_distribution<>(0, m_size - 1)(gen);
            int r_id2 = std::uniform_int_distribution<>(0, m_size - 1)(gen);
            double binary = (dis_binary(gen) < 0.5) ? -1.0 : 1.0;
            std::vector<double> new_pos;
            for (int j = 0; j < m_num_nodes * 2; ++j)
            {
                new_pos.push_back(m_position_best[j] + (m_unit_list[r_id1].position[j] + binary * m_unit_list[r_id2].position[j]) / 2);
            }
            new_pos = get_out_bound_value(new_pos);
            std::vector<std::vector<bool>> new_connections = m_unit_list[i].connections; // 保留原连接
            double new_value = cal_fitfunction(new_pos, new_connections);
            if (new_value > m_unit_list[i].value)
            {
                m_unit_list[i].position = new_pos;
                m_unit_list[i].value = new_value;
            }
        }
    }
}

std::vector<double> DOA_Base::get_sumatory(int na, int id)
{
    std::vector<int> A(m_size);
    std::iota(A.begin(), A.end(), 0);
    std::random_shuffle(A.begin(), A.end());
    std::vector<int> r_ids;
    for (int i = 0; i < na + 1; ++i)
    {
        if (A[i] != id)
            r_ids.push_back(A[i]);
    }

    std::vector<double> pos_sumatory(m_num_nodes * 2, 0.0);
    for (auto rid : r_ids)
    {
        for (int i = 0; i < m_num_nodes * 2; ++i)
        {
            pos_sumatory[i] += m_unit_list[id].position[i] - m_unit_list[rid].position[i];
        }
    }
    for (auto& p : pos_sumatory)
    {
        p /= r_ids.size();
    }
    return pos_sumatory;
}
\
int DOA_Base::get_best_id()
{
    return std::max_element(m_unit_list.begin(), m_unit_list.end(),
                            [](const DOA_Unit& a, const DOA_Unit& b) { return a.value < b.value; }) - m_unit_list.begin();
}

// 适应度函数
double DOA_Base::cal_fitfunction(const std::vector<double>& position, const std::vector<std::vector<bool>>& connections)
{
    // 这里简化了适应度计算，实际应用中需要根据具体的网络拓扑需求来实现
    double fitness = 0.0;
    for (int i = 0; i < m_num_nodes; ++i)
    {
        for (int j = i + 1; j < m_num_nodes; ++j)
        {
            if (connections[i][j])
            {
                double distance = std::sqrt(std::pow(position[i * 2] - position[j * 2], 2) + std::pow(position[i * 2 + 1] - position[j * 2 + 1], 2));
                fitness += 1.0 / distance; // 假设更短的距离带来更高的适应度
            }
        }
    }
    return fitness;
}

std::vector<double> DOA_Base::get_out_bound_value(const std::vector<double>& position)
{
    std::vector<double> new_pos = position;
    for (int i = 0; i < m_num_nodes * 2; ++i)
    {
        if (new_pos[i] < m_range_min_list[i % 2] || new_pos[i] > m_range_max_list[i % 2])
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(m_range_min_list[i % 2], m_range_max_list[i % 2]);
            new_pos[i] = dis(gen);
        }
    }
    return new_pos;
}

// Levy函数的实现
std::vector<double> DOA_Base::levy(int d)
{
    double beta = 1.5;
    double sigma = std::pow((std::tgamma(1 + beta) * std::sin(M_PI * beta / 2) / (std::tgamma((1 + beta) / 2) * beta * std::pow(2, (beta - 1) / 2))), 1 / beta);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> normal(0, sigma);

    std::vector<double> step(d);
    for (int i = 0; i < d; ++i)
    {
        double u = normal(gen);
        double v = std::normal_distribution<>(0, 1)(gen);
        step[i] = u / std::pow(std::abs(v), 1 / beta);
    }
    return step;
}

// 导出拓扑结构到文件
void DOA_Base::export_topology(const std::string& filename)
{
    int best_id = get_best_id();
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    outFile<<"["<<std::endl;
    for (int i = 0; i < m_num_nodes; ++i)
    {
        outFile<<"   [";
        for (int j = 0; j < m_num_nodes; ++j)
        {
            outFile << (m_unit_list[best_id].connections[i][j] ? "1" : "0");
            if (j < m_num_nodes - 1)
                outFile << ",";
        }
        outFile<<"]";
        if(i < m_num_nodes -1) outFile<<",";
        outFile << std::endl;
    }
    outFile<<"]"<<std::endl;
    outFile.close();
}

int main()
{
    int num_nodes = 50;
    int size = 30;
    int iter_max = 100;
    std::vector<double> range_min_list = {-10000, -10000}; // 节点位置范围
    std::vector<double> range_max_list = {10000, 10000};

    DOA_Base doa(num_nodes, size, iter_max, range_min_list, range_max_list);
    doa.init();
    doa.run();
    doa.export_topology("C:\\Users\\Administrator\\CLionProjects\\aliyun1\\topology_output50.json");

    return 0;
}
